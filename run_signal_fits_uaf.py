"""
Run DCB signal-shape fits for every UAF_training signal model group.

Reads signal_uaf_manifest.json (group -> [(mass, h5), ...]), and for each group
calls fit_signalshapes.py over its 19 mass points, writing
signal_fits_uaf/<group>/sig_fit_<mass>.json (+ per-mass fit plots).

No interpolation here (per request). Groups run GROUP_WORKERS at a time; each
group fits its masses sequentially inside fit_signalshapes.py.
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

MANIFEST = "signal_uaf_manifest.json"
OUTROOT = "signal_fits_uaf"
M_MIN, M_MAX, BIN = "5", "130", "0.2"
GROUP_WORKERS = int(os.environ.get("SIG_WORKERS", "2"))

man = json.load(open(MANIFEST))


def run_group(item):
    group, entries = item
    entries = sorted(entries)
    outdir = os.path.join(OUTROOT, group)
    os.makedirs(outdir, exist_ok=True)
    # Skip fits already produced (unless SIG_FORCE=1); refit only what's missing.
    if os.environ.get("SIG_FORCE", "0") != "1":
        entries = [(m, f) for m, f in entries
                   if not os.path.exists(os.path.join(outdir, "sig_fit_%d.json" % m))]
    if not entries:
        return group, 0, 19, 19
    masses = [str(m) for m, _ in entries]
    files = [f for _, f in entries]
    cmd = [sys.executable, "fit_signalshapes.py", "--dcb-model",
           "--m-min", M_MIN, "--m-max", M_MAX, "--bin-size", BIN,
           "--no-error-band", "-M"] + masses + ["-i"] + files + ["-o", outdir]
    log = os.path.join(outdir, "fit_log.txt")
    with open(log, "w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=f).returncode
    n = len([m for m, _ in sorted(man[group])
             if os.path.exists(os.path.join(outdir, "sig_fit_%d.json" % m))])
    return group, rc, n, len(man[group])


def main():
    items = sorted(man.items())
    print("Fitting %d groups (%d concurrent)..." % (len(items), GROUP_WORKERS), flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=GROUP_WORKERS) as ex:
        for group, rc, n, tot in ex.map(run_group, items):
            print("  %-34s rc=%d  %d/%d fits" % (group, rc, n, tot), flush=True)
            results.append((group, n, tot))
    print("\n%-34s %s" % ("group", "fits ok"))
    print("-" * 50)
    total_ok = total = 0
    for g, n, tot in sorted(results):
        total_ok += n; total += tot
        print("%-34s %d/%d" % (g, n, tot))
    print("-" * 50)
    print("TOTAL: %d/%d signal fits" % (total_ok, total))


if __name__ == "__main__":
    main()
