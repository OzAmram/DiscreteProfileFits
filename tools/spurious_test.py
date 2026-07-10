"""
Spurious-signal test: fit the *pure background* MC (no injected signal) with the
signal PDF floated, at a chosen window/bin config. The fitted signal yield is the
spurious signal -- the bias from the background envelope absorbing/creating a peak
where there is none. This is the reviewer-standard way to quantify the absorption
bias, with no injected signal to confuse the measurement.

Reports, per mass:
  N_spur      = fitted signal yield on background-only MC
  N_spur_unc  = its uncertainty
  z_spur      = N_spur / N_spur_unc  (spurious signal in units of its own error)
  frac_stat   = N_spur / sqrt(B_win) (spurious signal vs sqrt(background) -- the
                natural signal-size scale; want << 1)
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

MASSES = [15, 35, 55, 75]
M_DATA_MIN, M_DATA_MAX = 10.0, 80.0
N_SIGMA = float(os.environ.get("SPUR_NSIGMA", "6"))
BIN_FACTOR = float(os.environ.get("SPUR_BINFRAC", "0.5"))
BKG = "bkg_mc_masses.h5"
OUTROOT = f"spurious_test_N{N_SIGMA:g}_b{BIN_FACTOR:g}"


def sigma_of(mass):
    with open(f"signal_fits/2B/case_interpolation_M{mass}.0.json") as f:
        return json.load(f)["sigma"]


def window(mass):
    s = sigma_of(mass)
    lo = max(M_DATA_MIN, mass - N_SIGMA * s)
    hi = min(M_DATA_MAX, mass + N_SIGMA * s)
    return round(lo, 4), round(hi, 4), round(BIN_FACTOR * s, 4)


def bwin(mass, lo, hi):
    with h5py.File(BKG) as f:
        b = f["masses"][:]
    return int(np.count_nonzero((b >= lo) & (b <= hi)))


def run_one(mass):
    lo, hi, bs = window(mass)
    outdir = f"{OUTROOT}/M{mass}"
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, "doFit.py", "-c", "bkg_mc_config.json", "-M", str(mass),
        "--m-min", str(lo), "--m-max", str(hi), "--bin-size", str(bs),
        "--sig_norm", "1000", "-i", BKG,
        "-s", f"signal_fits/2B/sig_fit_{mass}.json", "-o", outdir,
    ]
    with open(os.path.join(outdir, "fit_log.txt"), "w") as f:
        subprocess.run(cmd, stdout=f, stderr=f)
    return mass


def main():
    os.makedirs(OUTROOT, exist_ok=True)
    print(f"Spurious-signal test at +-{N_SIGMA:g} sigma, bin={BIN_FACTOR:g} sigma "
          f"on pure background MC ({BKG})", flush=True)
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(run_one, MASSES))

    print("\n" + "=" * 84)
    print(f"{'M':<5}{'window':<18}{'B_win':<10}{'N_spur':<16}{'z_spur':<9}{'N_spur/sqrtB':<14}{'sbfitP'}")
    print("-" * 84)
    for M in MASSES:
        lo, hi, bs = window(M)
        B = bwin(M, lo, hi)
        rj = f"{OUTROOT}/M{M}/fit_results_{M}.0.json"
        if not os.path.exists(rj):
            print(f"{M:<5}[{lo},{hi}]  FAILED")
            continue
        r = json.load(open(rj))
        ns = r["obs_excess_events"]
        nse = r["obs_excess_events_unc"]
        z = ns / nse if nse and nse > 0 else float("nan")
        frac = ns / np.sqrt(B) if B else float("nan")
        print(f"{M:<5}{('[%.1f,%.1f]' % (lo, hi)):<18}{B:<10}"
              f"{('%.0f +/- %.0f' % (ns, nse)):<16}{z:<9.2f}{frac:<14.2f}"
              f"{r.get('sbfit_prob', -1):<.3f}")
    print("=" * 84)
    print("Interpretation: |z_spur| < ~0.5 and |N_spur/sqrtB| small => negligible "
          "absorption bias.\nN_spur/sqrtB is spurious yield in units of a 1-sigma "
          "counting signal.")


if __name__ == "__main__":
    main()
