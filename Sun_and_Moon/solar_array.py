#!/usr/bin/env python3
"""Slurm array worker: compute one Sun or Moon livetime chunk (gtltcubesun).

Used between ``generate_sun_moon_templates.py --prepare-array`` (which writes
``array_plan.json``) and ``--finalize-array`` (which sums the chunks):

    python solar_array.py chunk <outdir>/array_plan.json <task_index>
    python solar_array.py count <events.fits> [--sun-chunk-days D] [--moon-chunk-days D]

``chunk`` runs task ``task_index`` of the plan, after checking that the inputs
and the Fermitools installation are unchanged, and writes a receipt when the
chunk's GTIs match the plan. Finished chunks are skipped, so a failed array
can simply be resubmitted. ``count`` prints the number of tasks.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time
import fcntl
import numpy as np
from astropy.io import fits


def merge_gti(values):
    """Sort and merge (start, stop) intervals; returns an (N, 2) array."""
    result = []
    for lo, hi in sorted(values):
        if hi <= lo:
            continue
        if result and lo <= result[-1][1] + 1e-6:
            result[-1][1] = max(result[-1][1], float(hi))
        else:
            result.append([float(lo), float(hi)])
    return np.asarray(result, dtype=float).reshape(-1, 2)


def read_gti(path):
    """Merged GTIs of a FITS file's GTI extension."""
    with fits.open(path, memmap=True) as h:
        return merge_gti(zip(h["GTI"].data["START"], h["GTI"].data["STOP"]))


def chunk_specs(gti, sun_days=30.0, moon_days=180.0):
    """Split the observation into Sun and Moon time chunks.

    Returns one dict per non-empty chunk with the body, chunk number, time
    range, gtltcubesun angle/power settings and the GTIs it must contain.
    Together the chunks of each body cover the GTIs exactly once.
    """
    if not len(gti):
        raise ValueError("Empty GTI")
    result = []
    for body, days, angle, power in [("SUN", sun_days, 180.0, 2.7), ("MOON", moon_days, 1.0, 2.0)]:
        if not np.isfinite(days) or days <= 0:
            raise ValueError("Invalid chunk days")
        edges = np.r_[np.arange(gti[0, 0], gti[-1, 1], days * 86400.0), gti[-1, 1]]
        for number, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            expected = merge_gti(zip(np.maximum(gti[:, 0], lo), np.minimum(gti[:, 1], hi)))
            if len(expected):
                result.append(
                    dict(
                        body=body,
                        number=number,
                        lo=float(lo),
                        hi=float(hi),
                        angle=angle,
                        power=power,
                        gti=expected.tolist(),
                    )
                )
    return result


def stamp(path):
    """Identity of a file: resolved path, size and modification time."""
    p = Path(path).resolve()
    s = p.stat()
    return dict(path=str(p), size=s.st_size, mtime_ns=s.st_mtime_ns)


def work(planfile, index):
    """Run array task ``index`` of ``planfile`` (see module docstring)."""
    planfile = Path(planfile).resolve()
    plan = json.loads(planfile.read_text())
    out = Path(plan["outdir"])
    run_lock = open(out / ".run.lock", "a")
    fcntl.flock(run_lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
    manifest_text = (out / "run_manifest.json").read_text()
    if hashlib.sha256(manifest_text.encode()).hexdigest() != plan["run_key"]:
        raise ValueError("Plan and run manifest differ")
    manifest = json.loads(manifest_text)
    for item in [manifest["events"], *manifest["spacecraft"]]:
        if stamp(item["path"]) != item:
            raise ValueError("Input changed: " + item["path"])
    if not 0 <= index < len(plan["tasks"]):
        raise ValueError("Array index outside task list")
    task = plan["tasks"][index]
    output = Path(task["outfile"])
    command = task["command"]
    receipt = Path(str(output) + ".done.json")
    if shutil.which("gtltcubesun") != command[0]:
        raise ValueError("Activate the same Fermitools environment as the preparation job")
    lock = open(str(output) + ".lock", "a")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    expected_key = hashlib.sha256(
        json.dumps([plan["run_key"], command], sort_keys=True).encode()
    ).hexdigest()
    if output.exists() and receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved.get("key") == expected_key and saved.get("output") == stamp(output):
            print("Already complete:", output, flush=True)
            return
    receipt.unlink(missing_ok=True)
    pfiles = os.environ.get("PFILES", "")
    if ";" not in pfiles:
        raise ValueError("Fermitools PFILES is not initialized")
    private = out / "pfiles" / output.name
    private.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["PFILES"] = str(private) + ";" + pfiles.split(";", 1)[1]
    for name in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        env[name] = "1"

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"Signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    log = out / "logs" / (output.name + ".log")
    print(f"Array task {index + 1}/{len(plan['tasks'])}: {output.name}", flush=True)
    with log.open("w") as f:
        f.write(json.dumps(command) + "\n")
        f.flush()
        proc = subprocess.Popen(command, stdout=f, stderr=subprocess.STDOUT, env=env)
        start = time.monotonic()
        try:
            while True:
                try:
                    status = proc.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    print(f"{output.name}: {(time.monotonic() - start) / 60:.1f} min", flush=True)
        except BaseException:
            proc.terminate()
            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
            raise
    if status != 0:
        raise RuntimeError(f"gtltcubesun failed ({status}): {log}")
    with fits.open(output) as h:
        h.verify("exception")
    actual = read_gti(output)
    expected = np.asarray(task["gti"])
    if actual.shape != expected.shape or not np.allclose(actual, expected, atol=1e-3, rtol=0):
        raise ValueError("Chunk GTIs do not match the assigned interval")
    temp = Path(str(receipt) + ".tmp")
    temp.write_text(
        json.dumps(dict(key=expected_key, output=stamp(output), args=command), indent=2)
    )
    temp.replace(receipt)
    print("Completed:", output.name, flush=True)


def main():
    """Command-line entry point."""
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="mode", required=True)
    c = sub.add_parser("count")
    c.add_argument("events")
    c.add_argument("--sun-chunk-days", type=float, default=30.0)
    c.add_argument("--moon-chunk-days", type=float, default=180.0)
    w = sub.add_parser("chunk")
    w.add_argument("plan")
    w.add_argument("index", type=int)
    a = p.parse_args()
    if a.mode == "count":
        print(len(chunk_specs(read_gti(a.events), a.sun_chunk_days, a.moon_chunk_days)))
    else:
        work(a.plan, a.index)


if __name__ == "__main__":
    main()
