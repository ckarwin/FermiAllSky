#!/usr/bin/env python3
"""Make full-sky Sun and Moon template mapcubes with the Fermitools Solar System Tools.

Outputs ``sun_template.fits`` (solar disk + solar inverse Compton) and
``moon_template.fits``: full-sky Galactic CAR differential-intensity cubes
(ph cm^-2 s^-1 sr^-1 MeV^-1) with an ENERGIES extension. They can be used
directly as ``additional_diffuse`` components in Isotropic/config.yaml.

The time-averaged appearance of a moving body depends on where it was while
the LAT observed, so the templates must be made for the same event selection
(GTIs, zenith cut, event type, IRFs) as the analysis. The steps are:

1. ordinary livetime cube (``gtltcube``) and exposure (``gtexpcube2``);
2. livetime in the moving body's frame (``gtltcubesun``) for time chunks of
   ``--sun-chunk-days`` / ``--moon-chunk-days``, summed with ``gtltsumsun``;
3. body-frame exposure (``gtexphpsun``);
4. template = exposure-weighted emission profile (``gtsuntemp``), using the
   FSSC baseline profiles (downloaded automatically) unless others are given;
5. diagnostic plots, ``validation.json`` and an example XML model.

Run modes:

* no stage flag      : do everything in one job (step 2 can take very long);
* ``--prepare-array``: write ``array_plan.json`` / ``array_range.txt`` for a
  Slurm array in which each task runs one ``gtltcubesun`` chunk
  (see ``submit_solar_chunks.sbatch`` and ``solar_array.py``);
* ``--finalize-array``: after all array tasks succeeded, sum the chunks and
  make the templates (refuses to run if any chunk is missing).

The config may be the isotropic wrapper YAML (Isotropic/config.yaml) or a
standard Fermipy YAML. Every step writes a receipt, so rerunning with the same
settings skips finished work; changed settings require a new ``--outdir``.
"""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("config", type=Path, help="Isotropic wrapper YAML or Fermipy YAML")
    parser.add_argument("--events", type=Path, required=True, help="Final GTI-filtered FT1 FITS")
    parser.add_argument("--spacecraft", help="FT2 FITS or @list; default: config scfile")
    parser.add_argument("--ltcube", type=Path, help="Matching ordinary livetime cube")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("solar_lunar_templates"),
        help="Output folder (one per selection/settings)",
    )
    parser.add_argument("--component", help="Standard Fermipy component index or name")
    parser.add_argument("--sun-profile", help="Alternative SST intensity profile FITS")
    parser.add_argument("--moon-profile", help="Alternative SST intensity profile FITS")
    parser.add_argument(
        "--map-binsz", type=float, default=0.5, help="Template pixel size in deg (default 0.5)"
    )
    parser.add_argument(
        "--sun-chunk-days", type=float, default=30.0, help="Sun livetime chunk length (days)"
    )
    parser.add_argument(
        "--moon-chunk-days", type=float, default=180.0, help="Moon livetime chunk length (days)"
    )
    parser.add_argument(
        "--padding-bins",
        type=int,
        default=None,
        help="Default abs(edisp_bins)+2; includes true-energy margins",
    )
    stage = parser.add_mutually_exclusive_group()
    stage.add_argument(
        "--prepare-array",
        action="store_true",
        help="Prepare shared products and write array_plan.json; do not run chunks.",
    )
    stage.add_argument(
        "--finalize-array",
        action="store_true",
        help="Require all completed chunks, then sum, make templates and plots.",
    )
    args = parser.parse_args()
    # ------------------------------------------------------------------
    # Settings. Command-line values are copied into these UPPER-CASE names;
    # the remaining constants are expert settings that rarely need changing.
    # ------------------------------------------------------------------
    CONFIG_YAML = args.config.resolve()
    EVENT_FILE = args.events.resolve()
    SCFILE = args.spacecraft
    LIVETIME_CUBE = args.ltcube.resolve() if args.ltcube else None
    YAML_BASE_DIR = None
    COMPONENT = (
        int(args.component) if args.component and args.component.isdigit() else args.component
    )
    OUTDIR = args.outdir.resolve()

    # None -> read from YAML. IRFS must be explicit, e.g. P8R3_SOURCE_V3.
    IRFS_OVERRIDE = None
    EVTYPE_OVERRIDE = None
    LIVETIME_ZMAX_OVERRIDE = None  # Match the ordinary LAT exposure's zenith prescription
    EMIN_OVERRIDE = None  # MeV, measured analysis band
    EMAX_OVERRIDE = None

    EXPOSURE_BINS_PER_DECADE = 10
    ENERGY_PADDING_BINS = None  # Auto: |edisp_bins|+2 if edisp=True, otherwise 0
    LIVETIME_BINSZ = 0.5  # Keep 0.5 deg unless you have validated another choice
    DCOSTHETA = 0.025
    SUN_CHUNK_DAYS = 30
    MOON_CHUNK_DAYS = 180
    MAP_BINSZ = 0.5  # Full-sky Galactic CAR grid (720 x 360 by default)
    PREVIEW_ENERGY_MEV = 1000.0

    # Optional alternative *intensity profile* files; NOT precomputed sky templates.
    SUN_PROFILE = None
    MOON_PROFILE = None
    # Otherwise download the FSSC tutorial's baseline profiles automatically.

    MAP_BINSZ = args.map_binsz
    SUN_CHUNK_DAYS = args.sun_chunk_days
    MOON_CHUNK_DAYS = args.moon_chunk_days
    SUN_PROFILE = args.sun_profile
    MOON_PROFILE = args.moon_profile
    ENERGY_PADDING_BINS = args.padding_bins

    # ------------------------------------------------------------------
    # Imports, signal handling and Fermitools availability
    # ------------------------------------------------------------------
    import os, sys, json, copy, hashlib, shutil, subprocess, time, math
    import urllib.request, importlib.metadata
    import threading, signal, fcntl

    stop_workers = threading.Event()

    def stop_on_signal(signum, frame):
        stop_workers.set()
        raise KeyboardInterrupt(f"Received signal {signum}; stopping child tools.")

    signal.signal(signal.SIGTERM, stop_on_signal)
    import numpy as np
    import yaml
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    TOOLS = ["gtltcube", "gtltcubesun", "gtltsumsun", "gtexpcube2", "gtexphpsun", "gtsuntemp"]
    missing = [t for t in TOOLS if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            "Activate your Fermitools environment first. Missing tools: " + ", ".join(missing)
        )

    def merge_dict(base, extra):
        result = copy.deepcopy(base)
        for key, value in extra.items():
            if isinstance(value, dict) and isinstance(result.get(key), dict):
                result[key] = merge_dict(result[key], value)
            else:
                result[key] = copy.deepcopy(value)
        return result

    # ------------------------------------------------------------------
    # Read the selection from the YAML (wrapper configs are translated to the
    # equivalent Fermipy sections first)
    # ------------------------------------------------------------------
    CONFIG_YAML = CONFIG_YAML.expanduser().resolve()
    base_dir = Path(YAML_BASE_DIR).expanduser().resolve() if YAML_BASE_DIR else CONFIG_YAML.parent
    raw_cfg = yaml.safe_load(CONFIG_YAML.read_text())
    if isinstance(raw_cfg, dict) and "galprop" in raw_cfg and "selection" not in raw_cfg:
        from fit_isotropic import read_config

        wrapper = read_config(CONFIG_YAML)
        raw_cfg = {
            "selection": {
                k: wrapper[k] for k in ["emin", "emax", "evtype", "zmax", "tmin", "tmax"]
            },
            "gtlike": {"irfs": wrapper["irfs"], "edisp": True, "edisp_bins": wrapper["edisp_bins"]},
            "binning": {"binsperdec": wrapper["binsperdec"]},
            "data": {"scfile": wrapper["scfile"]},
        }

    if not isinstance(raw_cfg, dict):
        raise ValueError("Expected a Fermipy YAML mapping.")
    base_cfg = {k: v for k, v in raw_cfg.items() if k != "components"}
    components = raw_cfg.get("components")
    if components:
        if COMPONENT is None:
            raise ValueError("Set COMPONENT to the component index (list) or name (mapping).")
        selected = components[COMPONENT]
        if not isinstance(selected, dict):
            raise ValueError("Inline the selected component configuration in the YAML first.")
        cfg = merge_dict(base_cfg, selected)
    else:
        cfg = base_cfg

    def resolve_path(p):
        p = Path(os.path.expandvars(str(p))).expanduser()
        return (base_dir / p).resolve() if not p.is_absolute() else p.resolve()

    EVENT_FILE = resolve_path(EVENT_FILE)
    OUTDIR = OUTDIR.expanduser().resolve()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    # Hold this lock for main()'s lifetime; only one driver may own an output.
    run_lock = open(OUTDIR / ".run.lock", "a")
    try:
        fcntl.flock(run_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("Another generator is using this outdir. Choose a new outdir.")
    for name in ["profiles", "chunks", "logs", "pfiles"]:
        (OUTDIR / name).mkdir(exist_ok=True)
    selection = cfg.get("selection", {})
    gtlike = cfg.get("gtlike", {})

    def choose(override, mapping, name):
        value = override if override is not None else mapping.get(name)
        if value is None:
            raise ValueError(f"Missing {name}: supply it in the YAML or configuration override.")
        return value

    IRFS = str(choose(IRFS_OVERRIDE, gtlike, "irfs"))
    if IRFS.upper() == "CALDB":
        raise ValueError("Set the explicit IRF name, rather than CALDB, for reproducibility.")
    EVTYPE = int(choose(EVTYPE_OVERRIDE, selection, "evtype"))
    ZMAX = float(choose(LIVETIME_ZMAX_OVERRIDE, selection, "zmax"))
    EMIN = float(choose(EMIN_OVERRIDE, selection, "emin"))
    EMAX = float(choose(EMAX_OVERRIDE, selection, "emax"))
    if not (0 < EMIN < EMAX and 0 < ZMAX <= 180 and EVTYPE > 0):
        raise ValueError("Invalid energy, zenith, or event-type selection.")
    if SUN_CHUNK_DAYS <= 0 or MOON_CHUNK_DAYS <= 0 or EXPOSURE_BINS_PER_DECADE <= 0:
        raise ValueError("Chunk durations and energy sampling must be positive.")

    # Convert a spacecraft list into an explicit, absolute-path @list.
    sc_input = SCFILE if SCFILE is not None else cfg.get("data", {}).get("scfile")
    if sc_input is None:
        raise ValueError("Set SCFILE or data.scfile in the YAML.")
    if isinstance(sc_input, (list, tuple)):
        sc_paths = [resolve_path(p) for p in sc_input]
    elif str(sc_input).startswith("@"):
        listpath = resolve_path(str(sc_input)[1:])
        # Relative entries in an @list are resolved relative to the list file itself.
        sc_paths = []
        for line in listpath.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                p = Path(os.path.expandvars(line)).expanduser()
                sc_paths.append((listpath.parent / p).resolve() if not p.is_absolute() else p)
    else:
        sc_paths = [resolve_path(sc_input)]
    for p in [EVENT_FILE, *sc_paths]:
        if not p.is_file():
            raise FileNotFoundError(p)
    if not sc_paths:
        raise ValueError("Spacecraft list is empty.")
    if len(sc_paths) == 1:
        scarg = str(sc_paths[0])
    else:
        sclist = OUTDIR / "spacecraft_files.txt"
        sclist.write_text("".join(str(p) + "\n" for p in sc_paths))
        scarg = "@" + str(sclist)

    # Each Fermitools call gets private parameter files and single-threaded math
    # libraries, so parallel jobs cannot interfere with each other.
    # Use a private writable parameter directory to avoid changing other jobs' pfiles.
    pfiles = os.environ.get("PFILES", "")
    if ";" not in pfiles:
        raise RuntimeError(
            "PFILES is not initialized. Run this script from the activated Fermitools environment."
        )
    system_pfiles = pfiles.split(";", 1)[1]
    TOOL_ENV = os.environ.copy()
    TOOL_ENV["PFILES"] = str(OUTDIR / "pfiles") + ";" + system_pfiles
    for variable in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        TOOL_ENV[variable] = "1"
    versions = {}
    for pkg in ["fermitools", "fermipy", "astropy", "numpy"]:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = "not reported by Python metadata"
    print(
        json.dumps(
            {
                "IRFS": IRFS,
                "evtype": EVTYPE,
                "livetime_zmax": ZMAX,
                "energy_band_MeV": [EMIN, EMAX],
                "event_file": str(EVENT_FILE),
                "spacecraft_files": [str(p) for p in sc_paths],
                "versions": versions,
            },
            indent=2,
        )
    )

    # ------------------------------------------------------------------
    # Good-time intervals (GTIs) and spacecraft coverage. The event file's GTIs
    # define the observing time; every chunk and sum is checked against them.
    # ------------------------------------------------------------------
    def merge_intervals(intervals):
        a = np.asarray(intervals, dtype=float).reshape(-1, 2)
        if not np.isfinite(a).all() or np.any(a[:, 1] <= a[:, 0]):
            raise ValueError("GTIs must be finite intervals with STOP > START.")
        result = []
        for start, stop in a[np.argsort(a[:, 0])]:
            if result and start <= result[-1][1] + 1e-6:
                result[-1][1] = max(result[-1][1], stop)
            else:
                result.append([float(start), float(stop)])
        return np.asarray(result, dtype=float).reshape(-1, 2)

    def read_gti(path):
        with fits.open(path, memmap=True) as f:
            g = f["GTI"].data
            return merge_intervals(np.column_stack([g["START"], g["STOP"]]))

    def clipped_gti(gti, start, stop):
        a = np.maximum(gti[:, 0], start)
        b = np.minimum(gti[:, 1], stop)
        return np.column_stack([a[b > a], b[b > a]])

    def check_gti(path, expected):
        actual = read_gti(path)
        expected = merge_intervals(expected)
        if actual.shape != expected.shape or not np.allclose(actual, expected, rtol=0, atol=1e-3):
            raise RuntimeError(
                f"GTI mismatch in {path}. Stop: do not sum overlapping or incorrectly timed cubes."
            )

    GTI = read_gti(EVENT_FILE)
    if not len(GTI):
        raise ValueError("Selected event file has no good-time intervals.")
    TSTART, TSTOP = float(GTI[0, 0]), float(GTI[-1, 1])
    with fits.open(EVENT_FILE) as f:
        h = f["EVENTS"].header
        dss = {k: h[k] for k in h if k.startswith(("DSTYP", "DSVAL", "DSREF", "DSUNI"))}
    print("Event selection keywords:", json.dumps(dss, indent=2))
    print("MET range:", TSTART, TSTOP, "; accepted time [days]:", np.diff(GTI).sum() / 86400)
    sc_intervals = []
    for p in sc_paths:
        with fits.open(p, memmap=True) as f:
            tab = f["SC_DATA"].data
            sc_intervals.extend(zip(tab["START"], tab["STOP"]))
    sc_coverage = merge_intervals(sc_intervals)
    j = 0
    for start, stop in GTI:
        while j < len(sc_coverage) and sc_coverage[j, 1] < start - 1e-3:
            j += 1
        if (
            j == len(sc_coverage)
            or sc_coverage[j, 0] > start + 1e-3
            or sc_coverage[j, 1] < stop - 1e-3
        ):
            raise ValueError(f"Spacecraft history does not fully cover GTI {start}:{stop}.")

    # ------------------------------------------------------------------
    # Energy grid of the templates: the analysis range widened by `padding`
    # analysis bins on each side (default |edisp_bins| + 2) so the templates
    # cover the true-energy margin needed for energy dispersion.
    # ------------------------------------------------------------------
    analysis_bpd = float(cfg.get("binning", {}).get("binsperdec", 8))
    # Fermipy rounds the number of observed bins. Use its actual logarithmic width.
    n_analysis_bins = int(np.round(analysis_bpd * np.log10(EMAX / EMIN)))
    if n_analysis_bins < 1:
        raise ValueError("No analysis energy bins.")
    analysis_bpd = n_analysis_bins / np.log10(EMAX / EMIN)
    for key, actual in [("tmin", TSTART), ("tmax", TSTOP)]:
        requested = selection.get(key)
        if requested is not None and (
            (key == "tmin" and actual < float(requested) - 1.0)
            or (key == "tmax" and actual > float(requested) + 1.0)
        ):
            raise ValueError("EVENT_FILE GTIs extend outside configured time range.")

    padding = ENERGY_PADDING_BINS
    if padding is None:
        padding = abs(int(gtlike.get("edisp_bins", -1))) + 2 if gtlike.get("edisp", False) else 0
    if padding < 0:
        raise ValueError("ENERGY_PADDING_BINS must be nonnegative.")
    ELO = EMIN / 10 ** (padding / analysis_bpd)
    EHI = EMAX * 10 ** (padding / analysis_bpd)
    NBINS = int(math.ceil(np.log10(EHI / ELO) * EXPOSURE_BINS_PER_DECADE))
    energy_pars = dict(
        ebinalg="LOG", emin=ELO, emax=EHI, enumbins=NBINS, bincalc="EDGE", ebinfile="NONE"
    )
    expected_energy = np.geomspace(ELO, EHI, NBINS + 1)
    if not np.isfinite(MAP_BINSZ) or MAP_BINSZ <= 0:
        raise ValueError("map-binsz must be positive and finite.")
    nx, ny = int(round(360 / MAP_BINSZ)), int(round(180 / MAP_BINSZ))
    if not np.isclose(nx * MAP_BINSZ, 360) or not np.isclose(ny * MAP_BINSZ, 180):
        raise ValueError("MAP_BINSZ must tile 360 and 180 degrees exactly.")
    map_pars = dict(
        cmap="none",
        nxpix=nx,
        nypix=ny,
        binsz=MAP_BINSZ,
        coordsys="GAL",
        xref=0.0,
        yref=0.0,
        axisrot=0.0,
        proj="CAR",
    )
    print("Template energy coverage [MeV]:", ELO, EHI, "; energy planes:", NBINS + 1)
    print("Full-sky map size:", nx, "x", ny)

    # ------------------------------------------------------------------
    # Emission profiles (FSSC baseline solar/lunar models, or user files)
    # ------------------------------------------------------------------
    PROFILE_ROOT = "https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/data/latSolarTemplate/"

    def sha256_file(path):
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
                h.update(block)
        return h.hexdigest()

    def get_profile(override, name):
        target = resolve_path(override) if override is not None else OUTDIR / "profiles" / name
        if not target.exists() and override is None:
            temp = target.with_suffix(".download")
            try:
                request = urllib.request.Request(
                    PROFILE_ROOT + name, headers={"User-Agent": "LAT-template-notebook"}
                )
                with (
                    urllib.request.urlopen(request, timeout=120) as response,
                    open(temp, "wb") as out,
                ):
                    shutil.copyfileobj(response, out)
                with fits.open(temp) as f:
                    f.verify("exception")
                temp.replace(target)
            except Exception:
                temp.unlink(missing_ok=True)
                raise RuntimeError(
                    f"Download failed. Download {PROFILE_ROOT + name} manually and pass it with --sun-profile/--moon-profile."
                )
        with fits.open(target) as f:
            for ext in ["ANGLES", "ENERGIES", "SST_PROFILE"]:
                if ext not in f:
                    raise ValueError(f"{target} is missing {ext}. Supply an SST intensity profile.")
            e = np.asarray(f["ENERGIES"].data.field(0), dtype=float).ravel()
            if ELO < e.min() * (1 - 1e-6) or EHI > e.max() * (1 + 1e-6):
                raise ValueError(
                    f"Requested energy range exceeds {target.name}: {e.min()} to {e.max()} MeV. "
                    "Adjust padding/range or supply a suitable profile."
                )
        print(target.name, "SHA256", sha256_file(target))
        return target

    profiles = {
        "SUN": get_profile(SUN_PROFILE, "solar_profile_v2r0.fits"),
        "MOON": get_profile(MOON_PROFILE, "lunar_profile_v2r0.fits"),
    }

    # ------------------------------------------------------------------
    # Run manifest and restartable tool runner. The manifest records all
    # inputs/settings; a changed manifest in an existing outdir is an error.
    # run_tool() skips a step whose output and receipt are already valid.
    # ------------------------------------------------------------------
    def file_stamp(p):
        p = Path(p).resolve()
        st = p.stat()
        return {"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}

    existing_lt = resolve_path(LIVETIME_CUBE) if LIVETIME_CUBE is not None else None
    if existing_lt is not None:
        check_gti(existing_lt, GTI)
    manifest = dict(
        schema=2,
        config=cfg,
        events=file_stamp(EVENT_FILE),
        spacecraft=[file_stamp(p) for p in sc_paths],
        existing_livetime=file_stamp(existing_lt) if existing_lt else None,
        irfs=IRFS,
        evtype=EVTYPE,
        zmax=ZMAX,
        energy=energy_pars,
        spatial=map_pars,
        versions=versions,
        lt_binsz=LIVETIME_BINSZ,
        dcostheta=DCOSTHETA,
        chunk_days={"SUN": SUN_CHUNK_DAYS, "MOON": MOON_CHUNK_DAYS},
        profiles={b: sha256_file(p) for b, p in profiles.items()},
    )
    manifest_text = json.dumps(manifest, sort_keys=True, indent=2, default=str)
    manifest_path = OUTDIR / "run_manifest.json"
    if manifest_path.exists() and manifest_path.read_text() != manifest_text:
        raise RuntimeError(
            "Inputs/settings changed. Choose a new OUTDIR to avoid reusing incompatible products."
        )
    manifest_path.write_text(manifest_text)
    run_key = hashlib.sha256(manifest_text.encode()).hexdigest()

    def run_tool(tool, outfile, **parameters):
        if stop_workers.is_set():
            raise RuntimeError("Generation cancelled after interruption or worker failure.")
        outfile = Path(outfile)
        tool_env = TOOL_ENV.copy()
        private_pfiles = OUTDIR / "pfiles" / outfile.name
        private_pfiles.mkdir(exist_ok=True)
        tool_env["PFILES"] = str(private_pfiles) + ";" + system_pfiles
        receipt = Path(str(outfile) + ".done.json")
        log = OUTDIR / "logs" / (outfile.name + ".log")
        parameters = dict(parameters, outfile=str(outfile), clobber="yes", chatter=2, mode="h")
        tool_args = [shutil.which(tool)] + [f"{key}={value}" for key, value in parameters.items()]
        task_key = hashlib.sha256(
            json.dumps([run_key, tool_args], sort_keys=True).encode()
        ).hexdigest()
        if outfile.exists() and receipt.exists():
            saved = json.loads(receipt.read_text())
            if saved.get("key") == task_key and saved.get("output") == file_stamp(outfile):
                print("Reuse:", outfile.name, flush=True)
                return outfile
        if args.finalize_array and tool == "gtltcubesun":
            raise RuntimeError(
                f"Missing or invalid completed array chunk: {outfile}; refusing to sum incomplete data."
            )
        receipt.unlink(missing_ok=True)
        print("Run:", tool, "->", outfile.name, flush=True)
        start = time.monotonic()
        with open(log, "w") as stream:
            stream.write(json.dumps(tool_args) + "\n")
            stream.flush()
            process = subprocess.Popen(
                tool_args, stdout=stream, stderr=subprocess.STDOUT, env=tool_env
            )
            try:
                next_heartbeat = time.monotonic() + 30
                while True:
                    if stop_workers.is_set():
                        raise RuntimeError("Stopping child tool after worker failure/interruption.")
                    try:
                        status = process.wait(timeout=1)
                        break
                    except subprocess.TimeoutExpired:
                        if time.monotonic() >= next_heartbeat:
                            print(
                                f"  {tool} [{outfile.name}]: {(time.monotonic() - start) / 60:.1f} min; log: {log.name}",
                                flush=True,
                            )
                            next_heartbeat = time.monotonic() + 30
            except BaseException:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                raise
        if status != 0 or not outfile.is_file():
            print(log.read_text(errors="replace")[-6000:])
            raise RuntimeError(f"{tool} failed (status {status}). See {log}.")
        with fits.open(outfile, memmap=True) as f:
            f.verify("exception")
        receipt.write_text(
            json.dumps(
                {"key": task_key, "output": file_stamp(outfile), "args": tool_args}, indent=2
            )
        )
        return outfile

    def time_chunks(gti, days):
        start, stop = gti[0, 0], gti[-1, 1]
        edges = np.arange(start, stop, days * 86400.0)
        edges = np.r_[edges, stop]
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            expected = clipped_gti(gti, lo, hi)
            if len(expected):
                yield i, float(lo), float(hi), expected

    # --prepare-array: only write the chunk plan for a Slurm array and stop.
    if args.prepare_array:
        from solar_array import chunk_specs

        planned = []
        for item in chunk_specs(GTI, SUN_CHUNK_DAYS, MOON_CHUNK_DAYS):
            path = OUTDIR / "chunks" / f"{item['body'].lower()}_lt_{item['number']:05d}.fits"
            parameters = dict(
                evfile=str(EVENT_FILE),
                scfile=scarg,
                body=item["body"],
                tmin=item["lo"],
                tmax=item["hi"],
                dcostheta=DCOSTHETA,
                thetasunmax=item["angle"],
                powerbinsun=item["power"],
                binsz=LIVETIME_BINSZ,
                phibins=0,
                zmax=ZMAX,
                outfile=str(path),
                clobber="yes",
                chatter=2,
                mode="h",
            )
            command = [shutil.which("gtltcubesun")] + [f"{k}={v}" for k, v in parameters.items()]
            planned.append(dict(outfile=str(path), command=command, gti=item["gti"]))
        plan = dict(outdir=str(OUTDIR), run_key=run_key, tasks=planned)
        (OUTDIR / "array_plan.json").write_text(json.dumps(plan, indent=2))
        (OUTDIR / "array_range.txt").write_text(f"0-{len(planned) - 1}%1000\n")
        print(
            f"Array preparation complete: {len(planned)} tasks; index range 0-{len(planned) - 1}",
            flush=True,
        )
        return

    # Step 1: ordinary livetime cube (reuse --ltcube if given) and all-sky exposure.
    ordinary_lt = existing_lt
    if ordinary_lt is None:
        ordinary_lt = run_tool(
            "gtltcube",
            OUTDIR / "ordinary_ltcube.fits",
            evfile=str(EVENT_FILE),
            scfile=scarg,
            dcostheta=DCOSTHETA,
            binsz=LIVETIME_BINSZ,
            zmax=ZMAX,
            phibins=0,
            tmin=0,
            tmax=0,
        )
    check_gti(ordinary_lt, GTI)
    average_exposure = run_tool(
        "gtexpcube2",
        OUTDIR / "average_exposure.fits",
        infile=str(ordinary_lt),
        irfs=IRFS,
        evtype=EVTYPE,
        edisp_bins=0,
        ignorephi="no",
        thmin=0,
        thmax=180,
        **map_pars,
        **energy_pars,
    )

    # Step 2: livetime in the Sun/Moon frame, chunked in time and then summed.
    # thetasunmax / powerbinsun follow the Solar System Tools recommendations.
    moving_lt = {}
    for body, days, angle, power in [
        ("SUN", SUN_CHUNK_DAYS, 180.0, 2.7),
        ("MOON", MOON_CHUNK_DAYS, 1.0, 2.0),
    ]:
        tasks = list(time_chunks(GTI, days))
        print(body, ":", len(tasks), "nonempty time chunks")
        files = []
        for number, lo, hi, expected in tasks:
            path = OUTDIR / "chunks" / f"{body.lower()}_lt_{number:05d}.fits"
            run_tool(
                "gtltcubesun",
                path,
                evfile=str(EVENT_FILE),
                scfile=scarg,
                body=body,
                tmin=lo,
                tmax=hi,
                dcostheta=DCOSTHETA,
                thetasunmax=angle,
                powerbinsun=power,
                binsz=LIVETIME_BINSZ,
                phibins=0,
                zmax=ZMAX,
            )
            check_gti(path, expected)
            files.append(path)
        if len(files) == 1:
            moving_lt[body] = files[0]
        else:
            listfile = OUTDIR / f"{body.lower()}_livetime_files.txt"
            listfile.write_text("".join(str(f) + "\n" for f in files))
            moving_lt[body] = run_tool(
                "gtltsumsun",
                OUTDIR / f"{body.lower()}_ltcube.fits",
                infile1="@" + str(listfile),
                infile2="none",
            )
        check_gti(moving_lt[body], GTI)

    # Steps 3-4: body-frame exposure and the final templates.
    templates = {}
    for body in ["SUN", "MOON"]:
        exposure = run_tool(
            "gtexphpsun",
            OUTDIR / f"{body.lower()}_moving_exposure.fits",
            infile=str(moving_lt[body]),
            irfs=IRFS,
            evtype=EVTYPE,
            binsz=0,
            ignorephi="no",
            thmin=0,
            thmax=180,
            **energy_pars,
        )
        templates[body] = run_tool(
            "gtsuntemp",
            OUTDIR / f"{body.lower()}_template.fits",
            expsun=str(exposure),
            avgexp=str(average_exposure),
            sunprof=str(profiles[body]),
            irfs=IRFS,
            evtype=EVTYPE,
            **map_pars,
            **energy_pars,
        )
    print("Final templates:", {k: str(v) for k, v in templates.items()})

    # Step 5: validate the templates and make preview/diagnostic plots.
    report = {}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
    from matplotlib.patches import Rectangle

    for ax, body in zip(axes, ["SUN", "MOON"]):
        with fits.open(templates[body], memmap=True) as f:
            f.verify("exception")
            cube = f[0].data
            if cube is None or cube.ndim != 3:
                raise ValueError(f"{body}: expected a 3D primary-HDU mapcube.")
            energy = np.asarray(f["ENERGIES"].data.field(0), dtype=float).ravel()
            if len(energy) != len(cube) or energy.shape != expected_energy.shape:
                raise ValueError(f"{body}: unexpected energy-plane count.")
            if not np.allclose(energy, expected_energy, rtol=1e-5, atol=0):
                raise ValueError(f"{body}: output energy grid differs from requested grid.")
            if not np.isfinite(cube).all() or np.any(cube < 0) or not np.any(cube > 0):
                raise ValueError(f"{body}: invalid intensity values; inspect tool logs.")
            k = int(np.argmin(abs(np.log(energy / PREVIEW_ENERGY_MEV))))
            plane = cube[k].copy()
            w = WCS(f[0].header).celestial
            unit = f[0].header.get("BUNIT", "ph cm^-2 s^-1 MeV^-1 sr^-1 (SST convention)")
            report[body] = dict(
                path=str(templates[body]),
                shape=list(cube.shape),
                energy_min_MeV=float(energy[0]),
                energy_max_MeV=float(energy[-1]),
                preview_energy_MeV=float(energy[k]),
                unit=unit,
            )
        # Nearest-pixel display sampled on a Galactic grid, without changing the FITS.
        ll, bb = np.meshgrid(np.linspace(30, -30, 481), np.linspace(-30, 30, 481))
        px, py = w.world_to_pixel(SkyCoord(l=ll * u.deg, b=bb * u.deg, frame="galactic"))
        ix, iy = np.rint(px).astype(int), np.rint(py).astype(int)
        valid = (ix >= 0) & (ix < plane.shape[1]) & (iy >= 0) & (iy < plane.shape[0])
        display = np.full(ll.shape, np.nan)
        display[valid] = plane[iy[valid], ix[valid]]
        im = ax.imshow(display, origin="lower", extent=[30, -30, -30, 30], cmap="magma", vmin=0)
        ax.add_patch(Rectangle((-20, -20), 40, 40, fill=False, edgecolor="cyan"))
        ax.set(
            xlabel="Galactic longitude [deg]",
            ylabel="Galactic latitude [deg]",
            title=f"{body}: {energy[k] / 1000:.3g} GeV",
        )
        fig.colorbar(im, ax=ax, label="ph cm⁻² s⁻¹ MeV⁻¹ sr⁻¹")
    fig.savefig(OUTDIR / "sun_moon_preview.png", dpi=160)
    plt.close(fig)
    (OUTDIR / "validation.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    from solar_diagnostics import plot_templates

    plot_templates(templates["SUN"], templates["MOON"], OUTDIR)

    # Example Fermitools XML model with both templates (normalizations fixed at 1).
    import xml.etree.ElementTree as ET

    root = ET.Element("source_library", title="Solar and lunar foregrounds")
    for body in ["SUN", "MOON"]:
        source = ET.SubElement(root, "source", name="LAT_" + body.lower(), type="DiffuseSource")
        spectrum = ET.SubElement(source, "spectrum", type="PowerLaw")
        for name, value, minimum, maximum, free in [
            ("Prefactor", "1", "0", "10", "1"),
            ("Index", "0", "-1", "1", "0"),
            ("Scale", "1000", "1", "1000000", "0"),
        ]:
            ET.SubElement(
                spectrum,
                "parameter",
                name=name,
                value=value,
                scale="1",
                min=minimum,
                max=maximum,
                free=free,
            )
        spatial = ET.SubElement(
            source, "spatialModel", type="MapCubeFunction", file=str(templates[body])
        )
        ET.SubElement(
            spatial,
            "parameter",
            name="Normalization",
            value="1",
            scale="1",
            min="0",
            max="10",
            free="0",
        )
    ET.indent(root)
    xmlpath = OUTDIR / "sun_moon_sources.xml"
    ET.ElementTree(root).write(xmlpath, encoding="utf-8", xml_declaration=True)
    print(xmlpath)


if __name__ == "__main__":
    main()
