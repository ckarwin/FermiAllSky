"""Repository regression tests (no Fermipy/Fermitools or LAT data needed).

Run from the repository root:  python -m pytest -q tests
Requires NumPy, SciPy, Astropy, Matplotlib, PyYAML and pytest.
"""

import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pytest
import yaml
from astropy.io import fits

import fit_isotropic as driver
import loopI
from build_templates import build
from galprop_padding import extend_arrays
from native_spectra import check_seed
from reuse_products import (
    assemble_maps,
    map_files,
    source_index,
    source_key,
    sparse_problems,
    stage_products,
    install_source_reuse,
)
from solar_array import merge_gti, chunk_specs
from solar_diagnostics import car_solid_angles
from spectral_core import ProfileSpectra, parse_models, ratio

ROOT = Path(__file__).resolve().parents[1]


def config_file(tmp_path, **updates):
    """Write a copy of Isotropic/config.yaml with top-level overrides; return its path."""
    c = yaml.safe_load((ROOT / "Isotropic/config.yaml").read_text())
    c.update(updates)
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(c))
    return path


def test_default_directory_and_independent_norms(tmp_path, monkeypatch):
    """Shipped config: outdir resolves next to the YAML; Loop I shells come from LoopI/."""
    monkeypatch.chdir(tmp_path)
    c = driver.read_config(ROOT / "Isotropic/config.yaml")
    assert Path(c["outdir"]) == ROOT / "Isotropic/output_with_solar_and_loopI_separate"
    assert set(loopI.source_specs(c)) == {"LoopI_S1", "LoopI_S2"}
    for item in (c["additional_diffuse"][s] for s in loopI.source_specs(c)):
        assert Path(item["file"]).parent == ROOT / "LoopI"
        assert item["free"] and item["normalization"] == 0.14


def test_directory_override_and_shell_controls(tmp_path, monkeypatch):
    """loopI.directory (with $VARS) and per-shell overrides are honoured."""
    monkeypatch.setenv("LOOPI_TEST_DIR", "assets")
    path = config_file(
        tmp_path,
        loopI=dict(
            enabled=True,
            directory="$LOOPI_TEST_DIR",
            mode="separate",
            spectrum="empirical",
            free_norm=True,
            normalization=0.14,
            shells={"S1": {"free_norm": False, "normalization": 0.2}},
        ),
    )
    c = driver.read_config(path)
    assert loopI.asset_directory(c["loopI"]) == tmp_path / "assets"
    assert c["additional_diffuse"]["LoopI_S1"]["normalization"] == 0.2
    assert not c["additional_diffuse"]["LoopI_S1"]["free"]
    assert c["additional_diffuse"]["LoopI_S2"]["free"]


@pytest.mark.parametrize(
    "text", ["outdir: a\noutdir: b\n", "loopI:\n  mode: combined\n  mode: separate\n"]
)
def test_duplicate_yaml_rejected(tmp_path, text):
    """Duplicate keys (top level or nested) raise instead of being silently dropped."""
    p = tmp_path / "duplicate.yaml"
    p.write_text(text)
    with pytest.raises(yaml.constructor.ConstructorError, match="duplicate key"):
        driver.read_config(p)


def test_yaml_merge_override_is_explicit():
    """YAML merge keys (<<) remain allowed."""
    result = yaml.load(
        "base: &base {free: false}\nitem: {<<: *base, free: true}\n", Loader=driver.UniqueKeyLoader
    )
    assert result["item"]["free"] is True


@pytest.mark.parametrize("mode", ["combined", "separate"])
def test_missing_assets_explains_build(tmp_path, mode):
    """Missing Loop I files give an error that names build_templates.py."""
    c = driver.read_config(
        config_file(
            tmp_path, loopI=dict(enabled=True, mode=mode, spectrum="empirical", directory="missing")
        )
    )
    with pytest.raises(FileNotFoundError, match="build_templates.py"):
        loopI.ensure_cubes(c)


def test_disabled_loopi_needs_no_assets(tmp_path):
    """With Loop I disabled, no assets are needed."""
    c = driver.read_config(config_file(tmp_path, loopI=dict(enabled=False, directory="missing")))
    loopI.ensure_cubes(c)
    assert not loopI.source_specs(c)


@pytest.fixture(scope="module")
def assets(tmp_path_factory):
    """Run the real Loop I builder once into a temporary folder."""
    folder = tmp_path_factory.mktemp("loopI")
    # Exercise the actual full-resolution builder and configured output path.
    cfg = config_file(folder, loopI={"enabled": True, "directory": "assets"})
    build(cfg)
    return folder / "assets"


@pytest.mark.parametrize("spectrum", ["empirical", "powerlaw"])
def test_generated_shells_sum_and_seed_validation(assets, spectrum):
    """S1 + S2 cubes equal the combined cube; seed cubes pass check_seed; no needless rebuild."""
    c = {
        "galprop": {},
        "loopI": dict(enabled=True, mode="separate", spectrum=spectrum, directory=str(assets)),
    }
    loopI.configure(c)
    loopI.ensure_cubes(c)
    suffix = "empirical" if spectrum == "empirical" else "powerlaw_seed"
    with (
        fits.open(assets / f"loopI_{suffix}_mapcube.fits") as combined,
        fits.open(c["additional_diffuse"]["LoopI_S1"]["file"]) as s1,
        fits.open(c["additional_diffuse"]["LoopI_S2"]["file"]) as s2,
    ):
        for k in range(len(combined[0].data)):
            np.testing.assert_allclose(
                s1[0].data[k] + s2[0].data[k], combined[0].data[k], rtol=3e-6, atol=1e-44
            )
        omega = car_solid_angles(combined[0].header, combined[0].data.shape[1:])
        assert np.isclose(omega.sum(), 4 * np.pi)
    if spectrum != "empirical":
        for name, item in c["additional_diffuse"].items():
            check_seed(item["file"], loopI.analytic_model(c["loopI"], name))
    for item in c["additional_diffuse"].values():
        driver.validate_extra_cube(item["file"], 100, 300000)
    before = [Path(q["file"]).stat().st_mtime_ns for q in c["additional_diffuse"].values()]
    loopI.ensure_cubes(c)
    assert before == [Path(q["file"]).stat().st_mtime_ns for q in c["additional_diffuse"].values()]


def test_cache_invalidation_and_fit_only_changes(assets, tmp_path):
    """Fit-only edits keep the preparation signature; touching Loop I inputs changes it."""
    c = driver.read_config(
        config_file(
            tmp_path,
            loopI=dict(
                enabled=True,
                mode="separate",
                spectrum="empirical",
                directory=str(assets),
                free_norm=True,
            ),
        )
    )
    original = driver.template_signature(c)
    edited = copy.deepcopy(c)
    edited["loopI"]["normalization"] = 0.2
    edited["additional_diffuse"]["LoopI_S1"]["normalization"] = 0.2
    edited["mask"]["plane_abs_b_min"] = 20
    edited["plots"]["map_nlon"] = 360
    assert driver.template_signature(edited) == original
    assert driver.result_dir(edited) != driver.result_dir(c)
    for name in ("loopI_empirical_mapcube.fits", "model_metadata.json", "empirical_intensity.txt"):
        p = assets / name
        stat = p.stat()
        try:
            os.utime(p, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1000000))
            assert driver.template_signature(c) != original
        finally:
            os.utime(p, ns=(stat.st_atime_ns, stat.st_mtime_ns))


def test_export_uses_custom_directory_and_seed_index(tmp_path):
    """Analytic Loop I export uses the configured directory and the fitted index."""
    settings = {
        "enabled": True,
        "mode": "combined",
        "spectrum": "powerlaw",
        "directory": str(tmp_path),
        "parameters": {"index": {"value": 2.1, "free": False}},
    }
    (tmp_path / "model_metadata.json").write_text(
        json.dumps(
            dict(seed_index=2.8, intensity_at_1000_MeV=2e-9, allsky_map_integral_sr=4 * np.pi)
        )
    )
    model = loopI.analytic_model(settings)
    energy = np.array([500.0, 1000.0, 5000.0])
    loopI.export_spectrum(
        {"loopI": settings},
        energy,
        [{"name": "LoopI", "normalization": 0.4}],
        {"models": {"LoopI": model}},
        tmp_path,
    )
    np.testing.assert_allclose(
        np.loadtxt(tmp_path / "loopI_spectrum.txt")[:, 1], 0.4 * 2e-9 * (energy / 1000) ** -2.1
    )


def test_changed_file_content_and_relocation_fingerprints(tmp_path):
    """Source fingerprints follow file content, not file location."""
    p = tmp_path / "model.fits"
    p.write_bytes(b"one")
    q = tmp_path / "moved.fits"
    q.write_bytes(b"one")

    def source(file):
        return ET.fromstring(f'<source name="IC"><spatialModel file="{file}"/></source>')

    xml = tmp_path / "model.xml"
    key = source_key(source(p), xml)
    assert key == source_key(source(q), xml)
    p.write_bytes(b"two")
    assert key != source_key(source(p), xml)


def sparse_hdu(name, pixels=(-2147483648, 1, -2147483648, 2)):
    """Small sparse source-map HDU in the Fermitools layout (PIX nulls mimic the writer bug)."""
    h = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="PIX", format="J", null=-2147483648, array=pixels),
            fits.Column(name="CHANNEL", format="I", array=[0, 0, 1, 1]),
            fits.Column(name="VALUE", format="E", array=[1.0, 2.0, 3.0, 4.0]),
        ],
        name=name,
    )
    h.header.update(INDXSCHM="SPARSE", HPX_CONV="FGST_SRCMAP_SPARSE", NSIDE=1, N_ENERGY=2)
    return h


def test_chunk_reuse_and_pixel_zero_repair_without_mutating_old(tmp_path):
    """Maps in numbered chunks are found and repaired in the copy; originals untouched."""
    base = tmp_path / "srcmap_00.fits"
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(name="EBOUNDS")]).writeto(base)
    chunk = tmp_path / "srcmap_00.fits_02.fits"
    fits.HDUList([fits.PrimaryHDU(), sparse_hdu("IC")]).writeto(chunk)
    original = chunk.read_bytes()
    (tmp_path / "srcmap_00.fits_backup.fits").write_bytes(b"not a chunk")
    assert map_files(base) == [base, chunk]
    index = source_index(base, ["IC"])
    assert not sparse_problems(index)
    assert sparse_problems(index, allow_pixel_zero_repair=False)
    target = tmp_path / "new.fits"
    audit = assemble_maps(target, base, list(index.values()), ["IC"])
    assert audit[0]["rows"] == 2
    assert chunk.read_bytes() == original
    with fits.open(target) as h:
        np.testing.assert_array_equal(h["IC"].data["PIX"], [0, 1, 0, 2])
        np.testing.assert_array_equal(h["IC"].data["VALUE"], [1, 2, 3, 4])


def test_invalid_sparse_map_is_not_reused(tmp_path):
    """Sparse maps with invalid pixel indices are flagged for regeneration."""
    p = tmp_path / "srcmap_00.fits"
    fits.HDUList([fits.PrimaryHDU(), sparse_hdu("IC", [0, -1, 0, 2])]).writeto(p)
    assert "ic" in sparse_problems(source_index(p, ["IC"]))


def test_selective_reuse_omits_changed_diffuse(tmp_path):
    """Unchanged sources are reused; a changed diffuse source is regenerated."""
    old = tmp_path / "old" / "fermipy"
    old.mkdir(parents=True)
    new = tmp_path / "new"
    new.mkdir()
    same = '<source name="catalog"><spectrum type="PowerLaw"/></source>'
    before = '<source name="IC"><spectrum value="1"/></source>'
    after = '<source name="IC"><spectrum value="2"/></source>'
    (old / "srcmdl_00.xml").write_text("<source_library>" + same + before + "</source_library>")
    xml = new / "srcmdl_00.xml"
    xml.write_text("<source_library>" + same + after + "</source_library>")
    fits.HDUList([fits.PrimaryHDU(), sparse_hdu("catalog"), sparse_hdu("IC")]).writeto(
        old / "srcmap_00.fits"
    )
    ccube = new / "ccube.fits"
    fits.PrimaryHDU().writeto(ccube)

    class Component:
        files = {"srcmdl": str(xml), "srcmap": str(new / "srcmap.fits"), "ccube": str(ccube)}

        def _create_srcmaps(self, **kwargs):
            self.loaded = source_index(Path(self.files["srcmap"]), ["catalog", "IC"])

    component = Component()
    report = {}
    install_source_reuse(component, old, report, new / "report.json")
    component._create_srcmaps()
    assert set(component.loaded) == {"catalog"}
    assert report["regenerated_sources"] == ["IC"]
    assert report["regeneration_reasons"]["IC"] == "source_definition_changed"


def test_reuse_rejects_changed_response_before_copy(tmp_path):
    """A changed selection is rejected before anything is copied."""
    old = tmp_path / "old"
    old.mkdir()
    (old / "fermipy_config.yaml").write_text(
        yaml.safe_dump({"data": {}, "selection": {"emin": 300}})
    )
    with pytest.raises(ValueError, match="selection differs"):
        stage_products(old, tmp_path / "new" / "fermipy", {"data": {}, "selection": {"emin": 20}})
    assert not (tmp_path / "new").exists()


def test_padding_preserves_original_planes_and_zeros():
    """GALPROP padding keeps original planes and zero pixels exactly."""
    energy = np.array([100.0, 200.0, 400.0])
    data = (energy / 100.0)[:, None, None] ** -2 * np.array([[[1.0, 0.0], [2.0, 3.0]]])
    e, a, _ = extend_arrays(energy, data, 50, 800, {"low": True})
    for k, original in enumerate(energy):
        np.testing.assert_array_equal(a[np.flatnonzero(e == original)[0]], data[k])
    np.testing.assert_allclose(a, (e / 100.0)[:, None, None] ** -2 * data[0])


def test_solar_chunks_partition_observation():
    """Sun/Moon chunks cover the GTIs exactly once."""
    gti = merge_gti([(20, 30), (0, 10), (9, 15), (40, 40)])
    tasks = chunk_specs(gti, sun_days=10 / 86400, moon_days=20 / 86400)
    for body in ("SUN", "MOON"):
        intervals = [pair for task in tasks if task["body"] == body for pair in task["gti"]]
        np.testing.assert_allclose(merge_gti(intervals), gti)
        assert sum(hi - lo for lo, hi in intervals) == 25


def test_tied_spectral_shape_with_independent_shell_amplitudes():
    """Toy shared-shape fit recovers the index and both shell amplitudes."""
    energy = np.geomspace(300, 30000, 5)
    base = np.array(
        [
            np.full((5, 4), 20.0),
            (energy / 1000)[:, None] ** -2.4 * np.array([[10.0, 3.0, 1.0, 2.0]]),
            (energy / 1000)[:, None] ** -2.4 * np.array([[1.0, 2.0, 4.0, 12.0]]),
        ]
    )
    models = parse_models(
        {
            name: {"parameters": {"index": {"value": 2.4, "min": 1.0, "max": 4.0}}}
            for name in ("S1", "S2")
        }
    )

    def predict(name, s):
        return base[["iso", "S1", "S2"].index(name)] * ratio(energy, s)[:, None]

    truth = copy.deepcopy(models)
    for model in truth.values():
        model["parameters"]["index"]["value"] = 2.05
    fixed = np.full((5, 4), 2.0)
    counts = (
        fixed + 1.3 * base[0] + 0.7 * predict("S1", truth["S1"]) + 1.8 * predict("S2", truth["S2"])
    )
    problem = ProfileSpectra(
        counts,
        fixed,
        base,
        ["S1", "S2"],
        [True, True],
        [1.0, 1.0],
        models,
        predict,
        driver.JointPoisson,
        tied_sources=["S1", "S2"],
    )
    shape, (_, amplitudes, _, _) = problem.solve()
    np.testing.assert_allclose(shape, [2.05], atol=2e-4)
    np.testing.assert_allclose(amplitudes, [1.3] * 5 + [0.7, 1.8], rtol=2e-4)


def test_legacy_solar_entrypoint_reads_current_config(tmp_path):
    """Sun_and_Moon/fit_isotropic.py link reads the isotropic config correctly."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import fit_isotropic as f; c=f.read_config('../Isotropic/config.yaml'); "
            "print(c['additional_diffuse']['LoopI_S1']['file'])",
        ],
        cwd=ROOT / "Sun_and_Moon",
        text=True,
        capture_output=True,
        check=True,
    )
    assert str(ROOT / "LoopI" / "loopI_S1_empirical_mapcube.fits") in result.stdout


@pytest.mark.parametrize(
    "include_dm,include_bubbles", [(False, False), (False, True), (True, False), (True, True)]
)
def test_main_exports_only_enabled_components(include_dm, include_bubbles):
    """Main writes model maps only for the enabled optional components."""
    # Exercise the actual export block without importing the site-specific 3ML stack.
    import ast

    tree = ast.parse((ROOT / "Main/run_3ML_master.py").read_text())
    start = next(
        i
        for i, n in enumerate(tree.body)
        if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "optional_diffuse_names" for t in n.targets)
    )
    end = next(
        i
        for i, n in enumerate(tree.body)
        if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "ps_map" for t in n.targets)
    )

    class Maps:
        def __init__(self):
            self.written = []

        def write_model_map(self, path, name):
            self.written.append(name)

        def model_counts_map(self, exclude):
            self.excluded = exclude

    maps = Maps()
    fixed = ["Sun", "Moon", "Loop_I_S1", "Loop_I_S2"]
    namespace = dict(
        gta=maps,
        fixed_template_names=fixed,
        include_dm=include_dm,
        include_bubbles=include_bubbles,
        savename="test",
    )
    exec(
        compile(
            ast.Module(body=tree.body[start : end + 1], type_ignores=[]), "<export block>", "exec"
        ),
        namespace,
    )
    expected = (
        (["GCE_NFW"] if include_dm else [])
        + (["bubbles_low", "bubbles_high"] if include_bubbles else [])
        + fixed
    )
    assert maps.written == expected
    assert maps.excluded == ["pion_decay", "ic", "brem", "isodiff"] + expected
