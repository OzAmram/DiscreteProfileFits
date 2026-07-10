"""
Signal+background ensemble at realistic (3%) statistics.

For K independent 3% subsamples, inject a ~5 sigma signal (sized to each draw's
own in-window background), fit at the chosen +-N sigma / f*sigma config, and
record:
  - yield pull = (N_extr - N_inj_in_window) / sigma_extr
  - S+B fit chi2/ndof and p-value

Good config => pull distribution mean~0, width~1 (unbiased, correct errors) and
S+B chi2/ndof ~ 1 with ~uniform p-values.
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
TARGET_Z = float(os.environ.get("EN_TARGETZ", "5"))
SIGSRC = os.environ.get("EN_SIGSRC", "mc")  # "mc" = resampled MC events, "dcb" = sampled from DCB template
FRAC = 0.03
M_DATA_MIN, M_DATA_MAX = 10.0, 80.0
FULL = "bkg_mc_masses.h5"
OUTROOT = os.environ.get("EN_OUT", f"ensemble_sb_N{N_SIGMA:g}_b{BIN_FACTOR:g}")
NWORK = int(os.environ.get("EN_WORKERS", "3"))


def load_sig(mass):
    with h5py.File(f"signal_data/2B/signal_mass_{mass}.h5") as f:
        return f["masses"][:]


def sig_params(mass):
    with open(f"signal_fits/2B/sig_fit_{mass}.json") as f:
        p = json.load(f)
    return p["mean"], p["sigma"]


def _dcb(x, mean, sigma, alpha, n, alpha2, n2):
    t = (x - mean) / sigma
    Al = (n / alpha) ** n * np.exp(-0.5 * alpha ** 2); Bl = n / alpha - alpha
    Ar = (n2 / alpha2) ** n2 * np.exp(-0.5 * alpha2 ** 2); Br = n2 / alpha2 - alpha2
    return np.where(t < -alpha, Al * (Bl - t) ** (-n),
                    np.where(t > alpha2, Ar * (Br + t) ** (-n2), np.exp(-0.5 * t ** 2)))


def sample_dcb(mass, n, rng, lo, hi):
    """Accept-reject sample n masses from the DCB template over [lo,hi]."""
    with open(f"signal_fits/2B/sig_fit_{mass}.json") as f:
        p = json.load(f)
    args = (p["mean"], p["sigma"], p["alpha"], p["sign"], p["alpha2"], p["sign2"])
    xs = np.linspace(lo, hi, 4000)
    ymax = _dcb(xs, *args).max() * 1.05
    out = []
    while len(out) < n:
        xr = rng.uniform(lo, hi, size=2 * n)
        yr = rng.uniform(0, ymax, size=2 * n)
        out.extend(xr[yr < _dcb(xr, *args)].tolist())
    return np.array(out[:n], dtype=np.float32)


def window(mass, sigma):
    lo = max(M_DATA_MIN, mass - N_SIGMA * sigma)
    hi = min(M_DATA_MAX, mass + N_SIGMA * sigma)
    return round(lo, 4), round(hi, 4), round(BIN_FACTOR * sigma, 4)


def subsample(k):
    rng = np.random.default_rng(1000 + k)
    with h5py.File(FULL) as f:
        m = f["masses"][:]
    return m[rng.random(len(m)) < FRAC].astype(np.float32)


def build_mock(k, mass):
    """Return (mock_path, n_inj_in_window, sig_json_path, m_min, m_max, bs)."""
    rng = np.random.default_rng(5000 + k * 100 + mass)
    bkg = subsample(k)
    sig = load_sig(mass)
    mean, sigma = sig_params(mass)
    m_min, m_max, bs = window(mass, sigma)

    lo2, hi2 = mean - 2 * sigma, mean + 2 * sigma
    B_win = np.count_nonzero((bkg >= lo2) & (bkg <= hi2))
    if SIGSRC == "dcb":
        big = sample_dcb(mass, 20000, rng, m_min, m_max)
        f_win = np.count_nonzero((big >= lo2) & (big <= hi2)) / len(big)
        N_inject = int(round(TARGET_Z * np.sqrt(B_win) / f_win))
        inj = sample_dcb(mass, N_inject, rng, m_min, m_max)
    else:
        f_win = np.count_nonzero((sig >= lo2) & (sig <= hi2)) / len(sig)
        N_inject = int(round(TARGET_Z * np.sqrt(B_win) / f_win))
        idx = rng.choice(len(sig), size=N_inject, replace=N_inject > len(sig))
        inj = sig[idx]
    mock = np.concatenate([bkg, inj]).astype(np.float32)
    rng.shuffle(mock)

    n_inj = int(np.count_nonzero((inj >= m_min) & (inj <= m_max)))
    d = os.path.join(SCRATCH, f"sb_seed{k}_M{mass}")
    os.makedirs(d, exist_ok=True)
    mp = os.path.join(d, "mock.h5")
    with h5py.File(mp, "w") as f:
        f.create_dataset("masses", data=mock, compression="gzip", compression_opts=4, chunks=True)
    return mp, n_inj, f"signal_fits/2B/sig_fit_{mass}.json", m_min, m_max, bs


def run_one(job):
    k, mass = job
    mp, n_inj, sig_json, m_min, m_max, bs = build_mock(k, mass)
    outdir = os.path.join(OUTROOT, f"seed{k}", f"M{mass}")
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, "doFit.py", "-c", "bkg_mc_config.json", "-M", str(mass),
        "--m-min", str(m_min), "--m-max", str(m_max), "--bin-size", str(bs),
        "--sig_norm", str(max(n_inj, 1)), "-i", mp, "-s", sig_json, "-o", outdir,
    ]
    with open(os.path.join(outdir, "fit_log.txt"), "w") as f:
        subprocess.run(cmd, stdout=f, stderr=f)
    with open(os.path.join(outdir, "truth.json"), "w") as f:
        json.dump({"n_inj": n_inj}, f)
    return job


def main():
    os.makedirs(SCRATCH, exist_ok=True)
    os.makedirs(OUTROOT, exist_ok=True)
    print(f"S+B ensemble: K={K}, +-{N_SIGMA:g}sig, bin={BIN_FACTOR:g}sig, "
          f"targetZ={TARGET_Z:g}, masses={MASSES}", flush=True)
    jobs = [(k, m) for k in range(K) for m in MASSES]
    print(f"Running {len(jobs)} fits ({NWORK} concurrent)...", flush=True)
    with ThreadPoolExecutor(max_workers=NWORK) as ex:
        for i, _ in enumerate(ex.map(run_one, jobs), 1):
            if i % 20 == 0:
                print(f"  {i}/{len(jobs)} done", flush=True)

    print("\n" + "=" * 90)
    print(f"{'M':<5}{'N':<5}{'pull mean':<11}{'pull std':<10}{'<r>':<8}"
          f"{'<sb chi2/ndof>':<16}{'frac sbP<0.05':<14}")
    print("-" * 90)
    summary = {}
    for M in MASSES:
        pulls, rs, ratios, pvals = [], [], [], []
        for k in range(K):
            base = os.path.join(OUTROOT, f"seed{k}", f"M{M}")
            rj = os.path.join(base, f"fit_results_{M}.0.json")
            tj = os.path.join(base, "truth.json")
            if not (os.path.exists(rj) and os.path.exists(tj)):
                continue
            r = json.load(open(rj)); n_inj = json.load(open(tj))["n_inj"]
            ne, nee = r.get("obs_excess_events", 0), r.get("obs_excess_events_unc", 0)
            if nee and nee > 0:
                pulls.append((ne - n_inj) / nee)
                rs.append(ne / n_inj if n_inj else np.nan)
            c, n = r.get("sbfit_chi2", -1), r.get("sbfit_ndof", -1)
            if n and n > 0:
                ratios.append(c / n); pvals.append(r.get("sbfit_prob", -1))
        summary[M] = dict(pull=pulls, r=rs, chi2ndof=ratios, sbpval=pvals)
        if not pulls:
            print(f"{M:<5}0    (no results)"); continue
        pulls = np.array(pulls); ratios = np.array(ratios); pvals = np.array(pvals)
        print(f"{M:<5}{len(pulls):<5}{pulls.mean():<11.2f}{pulls.std():<10.2f}"
              f"{np.nanmean(rs):<8.2f}{ratios.mean():<16.2f}{np.mean(pvals<0.05):<14.2f}")
    print("=" * 90)
    print("Good => pull mean~0, pull std~1, <sb chi2/ndof>~1, frac(sbP<0.05)~0.05")
    with open(os.path.join(OUTROOT, "ensemble_sb_summary.json"), "w") as f:
        json.dump({str(M): summary[M] for M in MASSES}, f, indent=2)
    print(f"Saved {OUTROOT}/ensemble_sb_summary.json")


if __name__ == "__main__":
    main()
