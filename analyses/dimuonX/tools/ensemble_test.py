"""
Ensemble GOF/bias test at realistic (3%) statistics.

Draws K independent 3% subsamples of the full background MC (each emulating one
outcome of the anomaly-score cut). For each subsample, fits every mass through
`run_fit_adaptive.py` -- the SAME wrapper the data scan uses -- and records the
background-only fit chi2/ndof and p-value.

The point is to validate the deployed procedure end to end, so the adaptive
window/bin-size search is included rather than a single fixed window. Note the
consequence: the wrapper picks the attempt whose *S+B* GoF passes, so the
background-only p-values collected here are conditioned on that selection and
are not expected to be exactly uniform. Judge them as "no pathologies / mean
chi2/ndof ~ 1", not as a strict uniformity test.
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
# Adaptive-search settings -- these mirror run_data_scan.py's defaults so the
# ensemble exercises exactly the procedure applied to data.
N_SIGMA = float(os.environ.get("EN_NSIGMA", "7"))       # --n-sigma-start
N_SIGMA_MIN = float(os.environ.get("EN_NSIGMA_MIN", "4"))  # --n-sigma-min
BIN_FACTOR = float(os.environ.get("EN_BINFRAC", "0.5"))    # --bin-frac (start)
SIGDIR = os.environ.get("EN_SIGDIR", "signal_fits/2B_loosemass")
FRAC = 0.03
M_DATA_MAX = 80.0
FULL = "bkg_mc_masses.h5"
CONFIG = os.environ.get("EN_CONFIG", "dimuonX_config.json")
OUTROOT = os.environ.get("EN_OUT", f"ensemble_N{N_SIGMA:g}_b{BIN_FACTOR:g}")
# run_fit_adaptive.py lives in the general framework, three levels up from tools/.
FRAMEWORK_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
NWORK = int(os.environ.get("EN_WORKERS", "3"))


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
    outdir = os.path.join(OUTROOT, f"seed{k}", f"M{mass}")
    os.makedirs(outdir, exist_ok=True)
    # Same driver + flags as run_data_scan.py, so the adaptive window/bin-size
    # search is identical to the one applied to data.
    cmd = [
        sys.executable, os.path.join(FRAMEWORK_DIR, "run_fit_adaptive.py"), "-M", str(mass), "-c", CONFIG,
        "-s", SIGDIR, "-o", outdir,
        "--n-sigma-start", str(N_SIGMA), "--n-sigma-min", str(N_SIGMA_MIN),
        "--bin-frac", str(BIN_FACTOR), "--m-data-max", str(M_DATA_MAX),
        "-i", infile,
    ]
    with open(os.path.join(outdir, "adaptive_log.txt"), "w") as f:
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
    print("\n" + "=" * 96)
    print(f"{'M':<5}{'N':<5}{'<chi2/ndof>':<14}{'median':<9}{'max':<8}"
          f"{'frac p<0.05':<13}{'frac p<0.01':<13}{'<halfwidth/sig>':<17}{'coarsened':<10}")
    print("-" * 96)
    allrows = {}
    for M in MASSES:
        ratios, pvals, halfw, coarse = [], [], [], 0
        sigma = json.load(open(os.path.join(SIGDIR,
                          f"case_interpolation_M{float(M)}.json")))["sigma"]
        for k in range(K):
            rj = os.path.join(OUTROOT, f"seed{k}", f"M{M}", "best",
                              f"fit_results_{float(M)}.json")
            if not os.path.exists(rj):
                continue
            r = json.load(open(rj))
            c, n = r.get("bkgfit_chi2", -1), r.get("bkgfit_ndof", -1)
            if n and n > 0:
                ratios.append(c / n)
                pvals.append(r.get("bkgfit_prob", -1))
            # adaptive diagnostics: how far the window shrank / bins coarsened
            mn, mx = r.get("m_min"), r.get("m_max")
            if mn is not None and mx is not None:
                halfw.append(0.5 * (mx - mn) / sigma)
            bs = (r.get("script_options") or {}).get("bin_size")
            if bs and bs > (BIN_FACTOR * sigma) * 1.05:
                coarse += 1
        allrows[M] = {"chi2ndof": ratios, "pval": pvals,
                      "halfwidth_sig": halfw, "n_coarsened": coarse}
        if not ratios:
            print(f"{M:<5}0    (no results)")
            continue
        ratios = np.array(ratios); pvals = np.array(pvals)
        hw = np.mean(halfw) if halfw else float("nan")
        print(f"{M:<5}{len(ratios):<5}{ratios.mean():<14.2f}{np.median(ratios):<9.2f}"
              f"{ratios.max():<8.2f}"
              f"{np.mean(pvals < 0.05):<13.2f}{np.mean(pvals < 0.01):<13.2f}"
              f"{hw:<17.2f}{coarse:<10d}")
    print("=" * 96)
    print("Good fit => <chi2/ndof>~1, no pathological outliers. NOTE: p-values are")
    print("conditioned on the adaptive wrapper's S+B GoF selection, so frac(p<0.05)")
    print("is not expected to be exactly 0.05 -- see module docstring.")

    with open(os.path.join(OUTROOT, "ensemble_summary.json"), "w") as f:
        json.dump({str(M): allrows[M] for M in MASSES}, f, indent=2)
    print(f"Saved {OUTROOT}/ensemble_summary.json")


if __name__ == "__main__":
    main()
