"""
Scan window size (+-N sigma) x bin size (f * sigma) on the signal-injection
mocks, to find a configuration that minimizes the yield pull while keeping a
good goodness of fit and non-railed F-test orders.

Only the direct-fit signal shape is used (interp is ~identical) to halve runtime.
For each (N, f, mass) it records:
  - pull = (N_extr - N_inj_in_window) / sigma_extr
  - recovery = N_extr / N_inj
  - Z (discovery signif), sbfit_prob (S+B GOF)
  - selected F-test orders per family (bern/polyExp/exp/expPoly) and whether any
    rails against its configured ceiling.
"""
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

MASSES = [15, 35, 55, 75]
M_DATA_MIN, M_DATA_MAX = 10.0, 80.0
OUTROOT = "config_scan"

# (n_sigma, bin_fraction)
CONFIGS = [
    (6, 0.25), (6, 0.35), (6, 0.5),
    (7, 0.25), (7, 0.35), (7, 0.5),
]

# ceilings from bkg_mc_config.json func_forms (for rail detection)
CEIL = {"bern": 7, "polyExp": 4, "exp": 5, "expPoly": 4}
FAMILY_ORDER = ["bern", "polyExp", "exp", "expPoly"]  # order they appear in the log


def sigma_of(mass):
    with open(f"signal_fits/2B/case_interpolation_M{mass}.0.json") as f:
        return json.load(f)["sigma"]


def window(mass, nsig, binf):
    s = sigma_of(mass)
    lo = max(M_DATA_MIN, mass - nsig * s)
    hi = min(M_DATA_MAX, mass + nsig * s)
    return round(lo, 4), round(hi, 4), round(binf * s, 4)


def inj_in_window(mass, lo, hi):
    with h5py.File(f"mock_data/inj_M{mass}.h5") as f:
        inj = f["masses"][:]
    return int(np.count_nonzero((inj >= lo) & (inj <= hi)))


def tag_of(nsig, binf):
    return f"N{nsig}_b{binf:g}"


def run_one(job):
    nsig, binf, mass = job
    lo, hi, bs = window(mass, nsig, binf)
    n_inj = inj_in_window(mass, lo, hi)
    outdir = f"{OUTROOT}/{tag_of(nsig, binf)}/M{mass}"
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, "doFit.py", "-c", "bkg_mc_config.json", "-M", str(mass),
        "--m-min", str(lo), "--m-max", str(hi), "--bin-size", str(bs),
        "--sig_norm", str(n_inj), "-i", f"mock_data/mock_M{mass}.h5",
        "-s", f"signal_fits/2B/sig_fit_{mass}.json", "-o", outdir,
    ]
    with open(os.path.join(outdir, "fit_log.txt"), "w") as f:
        subprocess.run(cmd, stdout=f, stderr=f)
    return job


def parse_orders(logpath):
    if not os.path.exists(logpath):
        return []
    with open(logpath) as f:
        txt = f.read()
    return [int(m) for m in re.findall(r"Chose order (\d+) based on F-test", txt)]


def main():
    os.makedirs(OUTROOT, exist_ok=True)
    jobs = [(n, b, m) for (n, b) in CONFIGS for m in MASSES]
    print(f"Running {len(jobs)} fits (3 concurrent)...", flush=True)
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(run_one, jobs))

    rows = []
    for (nsig, binf) in CONFIGS:
        for mass in MASSES:
            lo, hi, bs = window(mass, nsig, binf)
            n_inj = inj_in_window(mass, lo, hi)
            outdir = f"{OUTROOT}/{tag_of(nsig, binf)}/M{mass}"
            rj = os.path.join(outdir, f"fit_results_{mass}.0.json")
            orders = parse_orders(os.path.join(outdir, "fit_log.txt"))
            om = dict(zip(FAMILY_ORDER, orders)) if len(orders) >= 4 else {}
            railed = [f for f in FAMILY_ORDER if om.get(f, -1) >= CEIL[f]]
            if not os.path.exists(rj):
                rows.append((nsig, binf, mass, n_inj, None, None, None, None, None, om, railed))
                continue
            r = json.load(open(rj))
            n_extr = r["obs_excess_events"]
            n_err = r["obs_excess_events_unc"]
            pull = (n_extr - n_inj) / n_err if n_err and n_err > 0 else float("nan")
            rec = n_extr / n_inj if n_inj else float("nan")
            rows.append((nsig, binf, mass, n_inj, n_extr, n_err, pull, rec,
                         r.get("sbfit_prob", -1), r.get("signif", -1), om, railed))

    # ---- report ----
    print("\n" + "=" * 100)
    hdr = f"{'win':<5}{'bin':<6}{'M':<5}{'N_inj':<8}{'N_extr':<9}{'pull':<8}{'rec':<7}{'sbfitP':<8}{'Z':<7}{'orders(be/pE/ex/eP)':<22}{'railed'}"
    print(hdr)
    print("-" * 100)
    summary = []
    for (nsig, binf, mass, n_inj, n_extr, n_err, pull, rec, sbp, Z, om, railed) in rows:
        ostr = "/".join(str(om.get(f, "?")) for f in FAMILY_ORDER)
        if n_extr is None:
            print(f"{nsig:<5}{binf:<6}{mass:<5}{n_inj:<8}{'FAILED':<9}")
            continue
        print(f"{nsig:<5}{binf:<6}{mass:<5}{n_inj:<8}{n_extr:<9.0f}"
              f"{pull:<8.2f}{rec:<7.2f}{sbp:<8.3f}{Z:<7.2f}{ostr:<22}{','.join(railed)}")
        summary.append(dict(nsig=nsig, binf=binf, mass=mass, n_inj=n_inj,
                            n_extr=n_extr, pull=pull, rec=rec, sbfit_prob=sbp,
                            Z=Z, orders=om, railed=railed))
    print("=" * 100)

    # per-config aggregate (max |pull|, min sbfitP, any railed)
    print("\nPer-config summary (want: small max|pull|, sbfitP>~0.05, no railed):")
    print(f"{'win':<5}{'bin':<6}{'max|pull|':<11}{'mean|pull|':<12}{'min sbfitP':<12}{'any railed?'}")
    print("-" * 60)
    for (nsig, binf) in CONFIGS:
        sub = [s for s in summary if s["nsig"] == nsig and s["binf"] == binf]
        if len(sub) < len(MASSES):
            print(f"{nsig:<5}{binf:<6}  incomplete ({len(sub)}/{len(MASSES)})")
            continue
        pulls = [abs(s["pull"]) for s in sub if s["pull"] == s["pull"]]
        maxp = max(pulls) if pulls else float("nan")
        meanp = np.mean(pulls) if pulls else float("nan")
        minsbp = min(s["sbfit_prob"] for s in sub)
        anyrail = any(s["railed"] for s in sub)
        print(f"{nsig:<5}{binf:<6}{maxp:<11.2f}{meanp:<12.2f}{minsbp:<12.3f}{anyrail}")

    with open(os.path.join(OUTROOT, "scan_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved {OUTROOT}/scan_summary.json")


if __name__ == "__main__":
    main()
