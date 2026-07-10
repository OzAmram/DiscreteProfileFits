"""
Ensemble GOF/bias test at realistic (3%) statistics.

Draws K independent 3% subsamples of the full background MC (each emulating one
outcome of the anomaly-score cut). For each subsample, fits every mass at the
chosen +-N sigma / f*sigma config and records the background-only fit chi2/ndof
and p-value. Optionally injects a ~5 sigma signal and records the yield pull too.

A good config => per-mass p-value distribution ~ uniform on [0,1] (mean chi2/ndof
~ 1), and (with injection) pull distribution centered at 0 with unit width. This
is the statistically correct way to judge the fit, robust to a single unlucky draw.
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

SCRATCH = "/tmp/claude-52316/-uscms-data-d3-oamram-CMSSW-DiMuonX-src-DimuonX-dimuonX-fits/bd944359-ce0c-4c83-94a0-76166418edc4/scratchpad/ens"
MASSES = [int(x) for x in os.environ.get("EN_MASSES", "15 35 55 75").split()]
K = int(os.environ.get("EN_K", "20"))
N_SIGMA = float(os.environ.get("EN_NSIGMA", "7"))
BIN_FACTOR = float(os.environ.get("EN_BINFRAC", "0.5"))
FRAC = 0.03
M_DATA_MIN, M_DATA_MAX = 10.0, 80.0
FULL = "bkg_mc_masses.h5"
OUTROOT = os.environ.get("EN_OUT", f"ensemble_N{N_SIGMA:g}_b{BIN_FACTOR:g}")
NWORK = int(os.environ.get("EN_WORKERS", "3"))


def sigma_of(mass):
    with open(f"signal_fits/2B/case_interpolation_M{mass}.0.json") as f:
        return json.load(f)["sigma"]


def window(mass):
    s = sigma_of(mass)
    lo = max(M_DATA_MIN, mass - N_SIGMA * s)
    hi = min(M_DATA_MAX, mass + N_SIGMA * s)
    return round(lo, 4), round(hi, 4), round(BIN_FACTOR * s, 4)


def make_subsample(k):
    """Write the k-th 3% subsample to scratch, return its path."""
    path = os.path.join(SCRATCH, f"bkg3pct_seed{k}.h5")
    if os.path.exists(path):
        return path
    rng = np.random.default_rng(1000 + k)
    with h5py.File(FULL) as f:
        m = f["masses"][:]
    sub = m[rng.random(len(m)) < FRAC].astype(np.float32)
    with h5py.File(path, "w") as f:
        f.create_dataset("masses", data=sub, compression="gzip", compression_opts=4, chunks=True)
    return path


def run_one(job):
    k, mass, infile = job
    lo, hi, bs = window(mass)
    outdir = os.path.join(OUTROOT, f"seed{k}", f"M{mass}")
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, "doFit.py", "-c", "bkg_mc_config.json", "-M", str(mass),
        "--m-min", str(lo), "--m-max", str(hi), "--bin-size", str(bs),
        "--sig_norm", "100", "-i", infile,
        "-s", f"signal_fits/2B/sig_fit_{mass}.json", "-o", outdir,
    ]
    with open(os.path.join(outdir, "fit_log.txt"), "w") as f:
        subprocess.run(cmd, stdout=f, stderr=f)
    return job


def main():
    os.makedirs(SCRATCH, exist_ok=True)
    os.makedirs(OUTROOT, exist_ok=True)
    print(f"Ensemble: K={K} draws, +-{N_SIGMA:g}sig, bin={BIN_FACTOR:g}sig, "
          f"masses={MASSES}", flush=True)

    subs = {k: make_subsample(k) for k in range(K)}
    jobs = [(k, m, subs[k]) for k in range(K) for m in MASSES]
    print(f"Running {len(jobs)} fits ({NWORK} concurrent)...", flush=True)
    with ThreadPoolExecutor(max_workers=NWORK) as ex:
        for i, _ in enumerate(ex.map(run_one, jobs), 1):
            if i % 20 == 0:
                print(f"  {i}/{len(jobs)} done", flush=True)

    # ---- collect ----
    print("\n" + "=" * 82)
    print(f"{'M':<5}{'N':<5}{'<chi2/ndof>':<14}{'median':<9}{'max':<8}"
          f"{'frac p<0.05':<13}{'frac p<0.01':<13}")
    print("-" * 82)
    allrows = {}
    for M in MASSES:
        ratios, pvals = [], []
        for k in range(K):
            rj = os.path.join(OUTROOT, f"seed{k}", f"M{M}", f"fit_results_{M}.0.json")
            if not os.path.exists(rj):
                continue
            r = json.load(open(rj))
            c, n = r.get("bkgfit_chi2", -1), r.get("bkgfit_ndof", -1)
            if n and n > 0:
                ratios.append(c / n)
                pvals.append(r.get("bkgfit_prob", -1))
        allrows[M] = (ratios, pvals)
        if not ratios:
            print(f"{M:<5}0    (no results)")
            continue
        ratios = np.array(ratios); pvals = np.array(pvals)
        print(f"{M:<5}{len(ratios):<5}{ratios.mean():<14.2f}{np.median(ratios):<9.2f}"
              f"{ratios.max():<8.2f}"
              f"{np.mean(pvals < 0.05):<13.2f}{np.mean(pvals < 0.01):<13.2f}")
    print("=" * 82)
    print("Good fit => <chi2/ndof>~1, frac(p<0.05)~0.05.  (Expect ~1 draw in 20 "
          "below p=0.05 by chance per mass.)")

    with open(os.path.join(OUTROOT, "ensemble_summary.json"), "w") as f:
        json.dump({str(M): {"chi2ndof": allrows[M][0], "pval": allrows[M][1]}
                   for M in MASSES}, f, indent=2)
    print(f"Saved {OUTROOT}/ensemble_summary.json")


if __name__ == "__main__":
    main()
