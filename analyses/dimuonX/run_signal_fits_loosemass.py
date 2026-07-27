"""
Run DCB signal-shape fits for every CATHODE_training_loosemass_v1 signal group.

Same as run_signal_fits_uaf.py but reads signal_loosemass_manifest.json and writes
to signal_fits_loosemass/<group>/sig_fit_<mass>.json (+ per-mass fit plots).
Interpolation is done separately by run_interpolation_loosemass.py.

Groups run GROUP_WORKERS at a time; each group fits its masses sequentially
inside fit_signalshapes.py. Set SIG_FORCE=1 to refit everything, SIG_WORKERS=N
to change concurrency.
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

# fit_signalshapes.py lives in the general framework, two levels up from here.
FRAMEWORK_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

MANIFEST = "signal_loosemass_manifest.json"
OUTROOT = "signal_fits_loosemass"
M_MIN, M_MAX, BIN = "5", "130", "0.2"
GROUP_WORKERS = int(os.environ.get("SIG_WORKERS", "4"))

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
        return group, 0, len(man[group]), len(man[group])
    masses = [str(m) for m, _ in entries]
    files = [f for _, f in entries]
    cmd = [sys.executable, os.path.join(FRAMEWORK_DIR, "fit_signalshapes.py"), "--dcb-model",
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
