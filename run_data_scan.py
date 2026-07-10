"""
Driver: run the adaptive fit for every (window, mass) in scan_plan.json.

Uses the validated configuration (bkg_mc_config func-forms, +/-7 sigma start,
0.5 sigma bins, adaptive shrink to 4 sigma on GoF failure), 2B signal shapes,
each window's local sideband-included spectrum (data_spectra/<window>.h5).

Runs NWORK adaptive fits concurrently. Results under data_fits/<window>/M<mass>/.
Writes data_scan_summary.json + prints a per-mass table.
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

CONFIG = "bkg_mc_config.json"
SIGDIR = "signal_fits/2B"
SPECTRA = "data_spectra"
OUTROOT = "data_fits"
NWORK = int(os.environ.get("DS_WORKERS", "3"))
NSIG_START = 7.0
NSIG_MIN = 4.0
BIN_FRAC = 0.5

plan = json.load(open("scan_plan.json"))


def run_one(job):
    tag, mass = job
    outdir = os.path.join(OUTROOT, tag, "M%s" % mass)
    cmd = [sys.executable, "run_fit_adaptive.py",
           "-M", str(mass), "-c", CONFIG, "-s", SIGDIR,
           "-o", outdir,
           "--n-sigma-start", str(NSIG_START), "--n-sigma-min", str(NSIG_MIN),
           "--bin-frac", str(BIN_FRAC),
           "-i", os.path.join(SPECTRA, "%s.h5" % tag)]
    log = os.path.join(outdir, "adaptive_log.txt")
    os.makedirs(outdir, exist_ok=True)
    with open(log, "w") as f:
        ret = subprocess.run(cmd, stdout=f, stderr=f)
    return tag, mass, ret.returncode


def collect(tag, mass):
    """Return (best_nsigma, sbprob, n_sig, n_sig_unc) from the best/ dir, or None."""
    best = os.path.join(OUTROOT, tag, "M%s" % mass, "best")
    rj = os.path.join(best, "fit_results_%s.json" % mass)
    if not os.path.exists(rj):
        return None
    r = json.load(open(rj))
    return r


def main():
    jobs = []
    for tag, info in plan.items():
        for m in info["masses"]:
            jobs.append((tag, m))
    jobs.sort()
    print("Running %d adaptive fits (%d concurrent)..." % (len(jobs), NWORK), flush=True)
    done = 0
    with ThreadPoolExecutor(max_workers=NWORK) as ex:
        for tag, mass, rc in ex.map(run_one, jobs):
            done += 1
            print("  [%d/%d] %s M%s  rc=%d" % (done, len(jobs), tag, mass, rc), flush=True)

    # summarize
    summary = {}
    print("\n" + "=" * 92)
    print("%-9s %-6s %-6s %-17s %-7s %-8s %-9s %s" %
          ("window", "mass", "stat", "N_sig +- unc", "signif", "sbProb", "obs_lim", "window[GeV]"))
    print("-" * 92)
    npass = 0
    for tag, info in plan.items():
        for m in info["masses"]:
            r = collect(tag, m)
            row = summary.setdefault(tag, {})
            if r is None:
                print("%-9s %-6s %-6s" % (tag, m, "FAIL"))
                row[str(m)] = {"status": "fail"}
                continue
            npass += 1
            ns = r.get("obs_excess_events", 0.0)
            nsu = r.get("obs_excess_events_unc", 0.0)
            sig = r.get("signif", float("nan"))
            p = r.get("sbfit_prob", -1)
            lim = r.get("obs_lim_events", float("nan"))
            mlo, mhi = r.get("m_min"), r.get("m_max")
            print("%-9s %-6s %-6s %7.1f +- %-6.1f %-7.2f %-8.3f %-9.1f [%.2f,%.2f]" %
                  (tag, m, "OK", ns, nsu, sig, p, lim, mlo, mhi))
            row[str(m)] = {"status": "ok", "center": info["center"],
                           "n_sig": ns, "n_sig_unc": nsu, "signif": sig,
                           "sbfit_prob": p, "obs_lim_events": lim,
                           "m_min": mlo, "m_max": mhi}
    print("-" * 92)
    print("passed %d / %d" % (npass, len(jobs)))
    json.dump(summary, open("data_scan_summary.json", "w"), indent=2)
    print("wrote data_scan_summary.json")


if __name__ == "__main__":
    main()
