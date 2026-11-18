"""Reuse expensive LAT products from an earlier preparation (``reuse_from``).

How reuse works:

1. :func:`stage_products` checks that the old run used identical data,
   selection, binning and response settings, then COPIES the selected events,
   counts cube, livetime cube and exposure map into the new ``fermipy/``
   folder. Any mismatch is an error (make a fresh preparation instead).
2. :func:`install_source_reuse` hooks into Fermipy's source-map step. For each
   source in the new model it compares a fingerprint (XML definition plus
   SHA-256 of every referenced file, see :func:`source_key`) with the old one.
   Unchanged sources are copied from the old source maps; new or changed ones
   are computed by Fermitools as usual.
3. Decisions are written to ``<outdir>/reuse_report.json``
   (``reused_sources``, ``regenerated_sources`` and the reason for each).

Old products are only read, never modified. The new output must be a
different directory. Only single-component (one event selection) analyses
are supported.
"""

import hashlib
import json
import os
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
import yaml


def checked(condition, message):
    """Raise ``ValueError("Reuse: " + message)`` unless ``condition`` holds."""
    if not condition:
        raise ValueError("Reuse: " + message)


def stage_products(previous, destination, cfg):
    """Validate an old output root and copy its reusable LAT products.

    Parameters
    ----------
    previous : path
        Old output root (``reuse_from``) containing ``fermipy_config.yaml``
        and ``fermipy/``.
    destination : path
        New ``<outdir>/fermipy`` folder (must not exist yet).
    cfg : dict
        New Fermipy configuration.

    Returns
    -------
    (old_fermipy_dir, report_dict)
        Passed on to :func:`install_source_reuse`.
    """
    old = Path(previous).resolve()
    dest = Path(destination).resolve()
    checked(old != dest.parent, "outdir must differ from reuse_from.")
    oldcfg_path = old / "fermipy_config.yaml"
    checked(oldcfg_path.is_file(), "reuse_from must contain fermipy_config.yaml.")
    oldcfg = yaml.safe_load(oldcfg_path.read_text())
    for section in ("data", "selection", "binning", "gtlike"):
        checked(
            oldcfg.get(section) == cfg.get(section),
            f"{section} differs; use fresh preparation for changed data/selection/grid/response.",
        )
    src = old / "fermipy"
    required = ["ft1_00.fits", "ccube_00.fits", "ltcube_00.fits", "bexpmap_00.fits"]
    for name in required:
        checked((src / name).is_file(), f"missing {src / name}")
    # Legacy outputs have no immutable input manifest. Refuse files (including
    # list members) modified after selection. Keep data inputs immutable.
    from astropy.io import fits

    def input_files(p):
        p = Path(p)
        yield p
        try:
            with fits.open(p):
                return
        except OSError:
            pass
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            q = Path(os.path.expandvars(line)).expanduser()
            checked(
                q.is_absolute(), "input-list members must use absolute paths for verified reuse."
            )
            checked(q.is_file(), f"missing input-list member {q}")
            yield q

    for key in ("evfile", "scfile"):
        for p in input_files(cfg["data"][key]):
            checked(
                p.stat().st_mtime_ns <= (src / "ft1_00.fits").stat().st_mtime_ns,
                f"{p} is newer than selected events; cannot verify old selection.",
            )
    checked(not dest.exists(), "destination fermipy directory already exists; choose a new outdir.")
    dest.mkdir(parents=True)
    for name in required:
        print("Reusing (copying):", src / name, flush=True)
        shutil.copy2(src / name, dest / name)  # No hardlinks: Fermitools may write files.
    report = {
        "reuse_from": str(old),
        "copied_products": required,
        "reused_sources": [],
        "regenerated_sources": [],
    }
    return src, report


def source_key(source, xmlpath, oldest_product_ns=None, hashes=None):
    """Fingerprint of one ``<source>`` XML element.

    The XML is canonicalized and every ``file=`` attribute is replaced by the
    SHA-256 of that file's contents, so moving an identical file does not
    change the key but editing it does. Returns None if a referenced file is
    missing (or, for legacy outputs, newer than the source map), meaning the
    source cannot safely be reused.
    """
    hashes = {} if hashes is None else hashes

    def walk(node):
        attrs = dict(node.attrib)
        for key, value in list(attrs.items()):
            if key == "file":
                p = Path(os.path.expandvars(value)).expanduser()
                if not p.is_absolute():
                    p = xmlpath.parent / p
                if not p.is_file():
                    raise OSError(str(p))
                if oldest_product_ns is not None and p.stat().st_mtime_ns > oldest_product_ns:
                    raise OSError("Reference file newer than source map: " + str(p))
                stat = p.stat()
                cache_key = (str(p.resolve()), stat.st_size, stat.st_mtime_ns)
                if cache_key not in hashes:
                    digest = hashlib.sha256()
                    with p.open("rb") as f:
                        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
                            digest.update(block)
                    hashes[cache_key] = digest.hexdigest()
                attrs[key] = hashes[cache_key]
        return node.tag, sorted(attrs.items()), [walk(child) for child in node]

    try:
        return json.dumps(walk(source), sort_keys=True)
    except OSError:
        return None


def map_files(base):
    """Main FITS plus exact numbered Fermitools chunks, in numeric order."""
    import re

    base = Path(base)
    pattern = re.compile(re.escape(base.name) + r"_(\d+)\.(?:fits|fit)$")
    chunks = []
    for path in base.parent.iterdir():
        match = pattern.fullmatch(path.name)
        if match and path.is_file():
            chunks.append((int(match.group(1)), path))
    return ([base] if base.is_file() else []) + [p for _, p in sorted(chunks)]


def hdu_records(path):
    """Read headers only; retain byte ranges including FITS block padding."""
    from astropy.io import fits

    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdus:
        for i, hdu in enumerate(hdus):
            info = hdu.fileinfo()
            start = info["hdrLoc"]
            size = info["datLoc"] + info["datSpan"] - start
            checked(start + size <= path.stat().st_size, f"Truncated FITS HDU in {path}")
            yield hdu.name.casefold(), (path, start, size), i


def source_index(base, names):
    """Index sources across all chunks; refuse ambiguous duplicate names."""
    wanted = {name.casefold() for name in names}
    result = {}
    for path in map_files(base):
        for name, record, i in hdu_records(path):
            if i and name in wanted:
                if name in result and result[name][0] == Path(base) and path != Path(base):
                    # Native loading uses the main file. saveSourceMaps may also
                    # emit the same source into a numbered chunk afterwards.
                    continue
                checked(name not in result, f"Duplicate source HDU {name} in {base} files")
                result[name] = record
    return result


def pixel_zero_rows(hdu):
    """Recognize the Fermitools unsigned-zero -> FITS-null writer bug.

    Likelihood FileUtils writes ipix as size_t; tip FitsColumn::setScalar
    uses fits_write_colnull with FitsPrimProps<unsigned long>::undefined()==0.
    Thus real PIX=0 is encoded as the signed J column's TNULL. The inverse is
    allowed only for the exact native schema and an unambiguous ordered map.
    No values or channels are changed, and no rows are dropped.
    """
    import numpy as np

    header = hdu.header
    if str(header.get("INDXSCHM", "")).strip().upper() != "SPARSE":
        return np.array([], dtype=np.int64)
    columns = {c.name.upper(): c for c in hdu.columns}
    if "PIX" not in columns:
        raise ValueError("Missing PIX column")
    pix = np.asarray(hdu.data[columns["PIX"].name])
    null = columns["PIX"].null
    rows = np.flatnonzero(pix == null) if null is not None else np.array([], dtype=np.int64)
    if not len(rows):
        return rows
    if (
        str(header.get("HPX_CONV", "")).strip() != "FGST_SRCMAP_SPARSE"
        or set(columns) != {"PIX", "CHANNEL", "VALUE"}
        or str(columns["PIX"].format) != "J"
        or str(columns["CHANNEL"].format) != "I"
        or str(columns["VALUE"].format) != "E"
        or null != -2147483648
        or any(c.bscale not in (None, 1) or c.bzero not in (None, 0) for c in columns.values())
    ):
        raise ValueError("PIX nulls do not match the native unsigned-zero writer bug schema")
    channel = np.asarray(hdu.data[columns["CHANNEL"].name])
    value = np.asarray(hdu.data[columns["VALUE"].name])
    npix = 12 * int(header["NSIDE"]) ** 2
    nenergy = int(header["N_ENERGY"])
    if (
        npix <= 0
        or nenergy <= 0
        or np.any(channel < 0)
        or np.any(channel >= nenergy)
        or not np.all(np.isfinite(value))
    ):
        raise ValueError("Invalid CHANNEL/VALUE or geometry; cannot repair PIX nulls")
    restored = pix.astype(np.int64)
    restored[rows] = 0
    if np.any(restored < 0) or np.any(restored >= npix):
        raise ValueError("Other invalid pixel indices; cannot repair PIX nulls")
    # Native SparseVector is written in increasing channel*npix+pixel order.
    # Strict ordering also rules out duplicate zero pixels within a channel.
    keys = channel.astype(np.int64) * npix + restored
    if len(keys) > 1 and np.any(keys[1:] <= keys[:-1]):
        raise ValueError("Restoring PIX=0 would violate native key ordering/uniqueness")
    return rows


def repair_pixel_zero_file(path, names):
    """Repair only a newly created output file; old cache files stay read-only."""
    from astropy.io import fits

    wanted = {name.casefold() for name in names}
    audit = []
    with fits.open(path, mode="update", memmap=True) as hdus:
        for hdu in hdus:
            if hdu.name.casefold() not in wanted:
                continue
            rows = pixel_zero_rows(hdu)
            if len(rows):
                channels = hdu.data["CHANNEL"][rows].astype(int).tolist()
                hdu.data["PIX"][rows] = 0
                if "CHECKSUM" in hdu.header:
                    hdu.add_checksum()
                elif "DATASUM" in hdu.header:
                    hdu.add_datasum()
                audit.append(
                    dict(source=hdu.name, restored_pixel=0, rows=int(len(rows)), channels=channels)
                )
                print(
                    f"Restored PIX=0 for {hdu.name}: {len(rows)} rows "
                    "(Fermitools unsigned-zero serialization bug).",
                    flush=True,
                )
        hdus.flush(output_verify="exception")
    return audit


def sparse_problems(index, allow_pixel_zero_repair=True):
    """Validate cached sparse HEALPix tables before handing them to Fermitools.

    Integer FITS nulls can remain finite in numpy; check TNULL explicitly.
    Only the verified unsigned-zero serialization pattern can be restored.
    """
    import numpy as np
    from astropy.io import fits

    grouped = {}
    for name, record in index.items():
        grouped.setdefault(record[0], set()).add(name)
    problems = {}
    for path, wanted in grouped.items():
        with fits.open(path, memmap=True) as hdus:
            for hdu in hdus:
                name = hdu.name.casefold()
                if (
                    name not in wanted
                    or str(hdu.header.get("INDXSCHM", "")).strip().upper() != "SPARSE"
                ):
                    continue
                try:
                    columns = {c.name.upper(): c for c in hdu.columns}
                    if not {"PIX", "CHANNEL", "VALUE"}.issubset(columns):
                        raise ValueError("missing PIX, CHANNEL, or VALUE column")
                    nside = int(hdu.header["NSIDE"])
                    nenergy = int(hdu.header["N_ENERGY"])
                    if nside <= 0 or nenergy <= 0:
                        raise ValueError("invalid NSIDE or N_ENERGY")
                    repair_rows = pixel_zero_rows(hdu) if allow_pixel_zero_repair else []
                    invalid = 0
                    for start in range(0, len(hdu.data), 1000000):
                        rows = hdu.data[start : start + 1000000]
                        arrays = {
                            key: np.asarray(rows[col.name])
                            for key, col in columns.items()
                            if key in {"PIX", "CHANNEL", "VALUE"}
                        }
                        if allow_pixel_zero_repair and len(repair_rows):
                            arrays["PIX"] = arrays["PIX"].copy()
                            arrays["PIX"][arrays["PIX"] == columns["PIX"].null] = 0
                        bad = np.zeros(len(rows), dtype=bool)
                        for key, a in arrays.items():
                            bad |= ~np.isfinite(a)
                            null = columns[key].null
                            if null is not None:
                                bad |= a == null
                        pix, channel = arrays["PIX"], arrays["CHANNEL"]
                        bad |= (pix < 0) | (pix >= 12 * nside * nside)
                        bad |= (channel < 0) | (channel >= nenergy)
                        invalid += int(np.count_nonzero(bad))
                    if invalid:
                        problems[name] = (
                            f"{invalid} invalid sparse rows (null, nonfinite, or out-of-range index)"
                        )
                except (ValueError, KeyError, TypeError, OSError) as exc:
                    problems[name] = str(exc)
    return problems


def assemble_maps(target, structural_file, records, names=()):
    """Stream exact HDU bytes into a standard FITS file, without loading cubes.

    Only structural HDUs and explicitly selected sources are copied. Old fixed
    aggregate caches are excluded. Original files are never changed.
    """
    import tempfile

    structural = []
    for name, record, i in hdu_records(Path(structural_file)):
        if i == 0 or name in {"skymap", "ebounds", "gti"}:
            structural.append(record)
    checked(bool(structural), "No structural FITS HDUs found")
    fd, temporary = tempfile.mkstemp(prefix=target.name + ".", suffix=".tmp", dir=target.parent)
    try:
        with os.fdopen(fd, "wb") as output:
            for path, offset, length in structural + list(records):
                with path.open("rb") as source:
                    source.seek(offset)
                    while length:
                        block = source.read(min(length, 8 * 1024 * 1024))
                        checked(bool(block), f"Unexpected end of {path}")
                        output.write(block)
                        length -= len(block)
        audit = repair_pixel_zero_file(Path(temporary), names)
        invalid = sparse_problems(
            source_index(Path(temporary), names), allow_pixel_zero_repair=False
        )
        checked(not invalid, f"Assembled maps remain invalid: {invalid}")
        checked(not target.exists(), f"Destination appeared during assembly: {target}")
        os.replace(temporary, target)
        return audit
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def install_source_reuse(component, oldwork, report, report_path):
    """Wrap Fermipy's source-map creation to reuse unchanged source maps.

    Replaces ``component._create_srcmaps`` with a version that builds a new
    ``srcmap`` file containing only the verified, unchanged sources (copied
    byte-for-byte from the old maps), then lets Fermitools compute all
    others. Cached sparse maps with invalid rows are regenerated instead.
    """
    oldwork, report_path = Path(oldwork), Path(report_path)
    original = component._create_srcmaps

    def create(overwrite=False, **kwargs):
        oldmap = oldwork / "srcmap_00.fits"
        oldxml = oldwork / "srcmdl_00.xml"
        newxml = Path(component.files["srcmdl"])
        target = Path(component.files["srcmap"])
        checked(
            not map_files(target),
            "Source-map products already exist in new output; choose a new outdir.",
        )
        if not map_files(oldmap) or not oldxml.is_file():
            report["source_map_note"] = "Old XML/maps absent; regenerating all source maps."
            original(overwrite=overwrite, **kwargs)
            report_path.write_text(json.dumps(report, indent=2))
            return
        tree = ET.parse(newxml)
        sources = tree.getroot().findall("source")
        names = [s.get("name") for s in sources]
        checked(len({n.casefold() for n in names}) == len(names), "Duplicate model source names")
        old_sources = {s.get("name"): s for s in ET.parse(oldxml).getroot().findall("source")}
        indices = source_index(oldmap, old_sources)
        print("Checking cached sparse source maps for invalid rows...", flush=True)
        invalid_old = sparse_problems(indices)
        hashes = {}
        receipt = oldwork.parent / "source_fingerprints.json"
        saved = json.loads(receipt.read_text()) if receipt.is_file() else {}
        fingerprints = {s.get("name"): source_key(s, newxml, hashes=hashes) for s in sources}
        reuse, reasons = [], {}
        for source in sources:
            name = source.get("name")
            if name not in old_sources:
                reasons[name] = "new_source"
            elif name.casefold() not in indices:
                reasons[name] = "missing_source_hdu"
            elif name.casefold() in invalid_old:
                reasons[name] = "invalid_cached_sparse_map: " + invalid_old[name.casefold()]
                print(f"Regenerating {name}: {reasons[name]}", flush=True)
            else:
                # Legacy referenced-file age is checked against the actual chunk.
                path = indices[name.casefold()][0]
                previous = (
                    saved.get(name)
                    if saved
                    else source_key(old_sources[name], oldxml, path.stat().st_mtime_ns, hashes)
                )
                if previous is None:
                    reasons[name] = "old_fingerprint_unavailable"
                elif fingerprints[name] is None:
                    reasons[name] = "new_fingerprint_unavailable"
                elif previous != fingerprints[name]:
                    reasons[name] = "source_definition_changed"
                else:
                    reuse.append(name)
        reuse_set = set(reuse)
        missing = [n for n in names if n not in reuse_set]
        report.update(
            reused_sources=reuse,
            regenerated_sources=missing,
            regeneration_reasons=reasons,
            source_map_files=[str(p) for p in map_files(oldmap)],
            reused_source_files={n: str(indices[n.casefold()][0]) for n in reuse},
        )
        report_path.write_text(json.dumps(report, indent=2))
        print(
            f"Reusing {len(reuse)} source maps; generating {len(missing)} new/changed maps.",
            flush=True,
        )
        if not reuse:
            # All-sky standard setup writes chunked products and creates maps
            # missing from the main file in native BinnedAnalysis.
            original(overwrite=overwrite, **kwargs)
            report["source_map_note"] = "No reusable maps; native Fermipy setup handles generation."
        else:
            checked(not overwrite, "Selective reuse requires setup(overwrite=False).")
            print("Preparing verified cached maps for native Fermipy loading...", flush=True)
            records = [indices[n.casefold()] for n in reuse]
            report["pixel_zero_repairs"] = assemble_maps(
                target, Path(component.files["ccube"]), records, reuse
            )
            actual = source_index(target, reuse)
            checked(len(actual) == len(reuse), "Assembled reused-map inventory incomplete")
            # Use the standard native loader for missing/new/changed sources:
            # SourceMapCache::getSourceMap loads a named main-file HDU if present,
            # otherwise createSourceMap computes it. No partial XML or separate
            # gtsrcmaps call/merge is needed. The original Fermipy hook sees the
            # prepared file and skips gtsrcmaps just as in a normal cached run.
            original(overwrite=False, **kwargs)
            report["source_map_note"] = (
                "Native BinnedAnalysis loads reused HDUs and computes omitted sources. "
                "No custom generation or merge of newly generated source maps."
            )
            print(
                f"Native likelihood initialization will compute {len(missing)} omitted source maps.",
                flush=True,
            )
        report["pixel_zero_repair_basis"] = (
            "Likelihood writes PIX as size_t; tip treats unsigned zero as null. "
            "Restored only exact native schema with valid ordered unique keys."
        )
        report_path.write_text(json.dumps(report, indent=2))
        # fit_isotropic writes full fingerprints only after native setup succeeds.

    component._create_srcmaps = create
