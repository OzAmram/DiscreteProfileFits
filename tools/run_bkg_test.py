"""
Background-only GOF test on a given background file (default: 3% subsample that
emulates the anomaly-score-cut statistics of the real search).

For each mass at a chosen +-N sigma window and f*sigma bin size, fit and report:
  - bkgfit chi2/ndof and p-value (background-only fit to the data)
  - sbfit  chi2/ndof and p-value (signal floated; on pure bkg r should be ~0)
  - N_spur (fitted signal on background-only = spurious signal) and z_spur
  - selected F-test orders per family and whether any rail against the ceiling
"""
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

MASSES = [int(x) for x in os.environ.get("BKG_MASSES", "15 35 55 75").split()]
M_DATA_MIN, M_DATA_MAX = 10.0, 80.0
N_SIGMA = float(os.environ.get("BKG_NSIGMA", "7"))
BIN_FACTOR = float(os.environ.get("BKG_BINFRAC", "0.25"))
BKG = os.environ.get("BKG_FILE", "bkg_mc_3pct.h5")
OUTROOT = os.environ.get("BKG_OUT", f"bkg3pct_N{N_SIGMA:g}_b{BIN_FACTOR:g}")

CEIL = {"bern": 7, "polyExp": 4, "exp": 5, "expPoly": 4}
FAM = ["bern", "polyExp", "exp", "expPoly"]


def sigma_of(mass):
    with open(f"signal_fits/2B/case_interpolation_M{mass}.0.json") as f:
        return json.load(f)["sigma"]


def window(mass):
    s = sigma_of(mass)
    lo = max(M_DATA_MIN, mass - N_SIGMA * s)
    hi = min(M_DATA_MAX, mass + N_SIGMA * s)
    return round(lo, 4), round(hi, 4), round(BIN_FACTOR * s, 4)


def run_one(mass):
    lo, hi, bs = window(mass)
    outdir = f"{OUTROOT}/M{mass}"
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, "doFit.py", "-c", "bkg_mc_config.json", "-M", str(mass),
        "--m-min", str(lo), "--m-max", str(hi), "--bin-size", str(bs),
        "--sig_norm", "100", "-i", BKG,
        "-s", f"signal_fits/2B/sig_fit_{mass}.json", "-o", outdir,
    ]
    with open(os.path.join(outdir, "fit_log.txt"), "w") as f:
        subprocess.run(cmd, stdout=f, stderr=f)
    return mass


def orders(logpath):
    if not os.path.exists(logpath):
        return {}
    txt = open(logpath).read()
    o = [int(m) for m in re.findall(r"Chose order (\d+) based on F-test", txt)]
    return dict(zip(FAM, o)) if len(o) >= 4 else {}


def main():
    os.makedirs(OUTROOT, exist_ok=True)
    print(f"BKG-only GOF: {BKG}  window=+-{N_SIGMA:g}sig  bin={BIN_FACTOR:g}sig", flush=True)
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(run_one, MASSES))

    print("\n" + "=" * 104)
    print(f"{'M':<5}{'nbins':<7}{'bkg chi2/ndof':<16}{'bkgP':<8}"
          f"{'sb chi2/ndof':<16}{'sbP':<8}{'N_spur':<15}{'z_spur':<8}"
          f"{'orders(be/pE/ex/eP)':<20}{'railed'}")
    print("-" * 104)
    for M in MASSES:
        rj = f"{OUTROOT}/M{M}/fit_results_{M}.0.json"
        om = orders(f"{OUTROOT}/M{M}/fit_log.txt")
        railed = ",".join(f for f in FAM if om.get(f, -1) >= CEIL[f])
        ostr = "/".join(str(om.get(f, "?")) for f in FAM)
        if not os.path.exists(rj):
            print(f"{M:<5}FAILED")
            continue
        r = json.load(open(rj))
        bc, bn = r.get("bkgfit_chi2", -1), r.get("bkgfit_ndof", -1)
        sc, sn = r.get("sbfit_chi2", -1), r.get("sbfit_ndof", -1)
        ns, nse = r.get("obs_excess_events", 0), r.get("obs_excess_events_unc", 0)
        z = ns / nse if nse and nse > 0 else float("nan")
        print(f"{M:<5}{bn+1:<7}{('%.1f/%d=%.2f' % (bc, bn, bc/bn if bn else 0)):<16}"
              f"{r.get('bkgfit_prob', -1):<8.3f}"
              f"{('%.1f/%d=%.2f' % (sc, sn, sc/sn if sn else 0)):<16}"
              f"{r.get('sbfit_prob', -1):<8.3f}"
              f"{('%.0f+/-%.0f' % (ns, nse)):<15}{z:<8.2f}{ostr:<20}{railed}")
    print("=" * 104)


if __name__ == "__main__":
    main()
