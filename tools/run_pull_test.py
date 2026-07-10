"""
Signal-injection pull test.

For each mass and each signal-shape variant (direct DCB fit vs interpolated),
fit the signal-injected mock at a fixed +-8 sigma window with sig_norm set to the
injected in-window yield (so r ~ 1), then report the yield pull:
    pull = (N_extracted - N_injected_in_window) / sigma_extracted
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

MASSES = [15, 35, 55, 75]
SHAPES = {
    "fit":    "signal_fits/2B/sig_fit_{M}.json",
    "interp": "signal_fits/2B/case_interpolation_M{M}.0.json",
}
M_DATA_MIN, M_DATA_MAX = 10.0, 80.0
N_SIGMA = float(os.environ.get("PULL_NSIGMA", "8"))
BIN_FACTOR = float(os.environ.get("PULL_BINFRAC", "0.25"))
MOCK_DIR = os.environ.get("PULL_MOCKDIR", "mock_data")
_tag = os.environ.get("PULL_TAG", "")
OUTBASE = f"pull_test{_tag}_{int(N_SIGMA)}sig_bin{BIN_FACTOR:g}"


def window(mass):
    with open(f"signal_fits/2B/case_interpolation_M{mass}.0.json") as f:
        sigma = json.load(f)["sigma"]
    m_min = max(M_DATA_MIN, mass - N_SIGMA * sigma)
    m_max = min(M_DATA_MAX, mass + N_SIGMA * sigma)
    return round(m_min, 4), round(m_max, 4), round(BIN_FACTOR * sigma, 4)


def inj_in_window(mass, m_min, m_max):
    with h5py.File(f"{MOCK_DIR}/inj_M{mass}.h5") as f:
        inj = f["masses"][:]
    return int(np.count_nonzero((inj >= m_min) & (inj <= m_max)))


def run_one(job):
    mass, shape, sig_json, m_min, m_max, bs, sig_norm = job
    outdir = os.path.join(OUTBASE, f"M{mass}_{shape}")
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, "doFit.py", "-c", "bkg_mc_config.json",
        "-M", str(mass), "--m-min", str(m_min), "--m-max", str(m_max),
        "--bin-size", str(bs), "--sig_norm", str(sig_norm),
        "-i", f"{MOCK_DIR}/mock_M{mass}.h5", "-s", sig_json, "-o", outdir,
    ]
    with open(os.path.join(outdir, "fit_log.txt"), "w") as logf:
        subprocess.run(cmd, stdout=logf, stderr=logf)
    return outdir


def main():
    os.makedirs(OUTBASE, exist_ok=True)

    jobs = []
    truth = {}
    for M in MASSES:
        m_min, m_max, bs = window(M)
        n_inj = inj_in_window(M, m_min, m_max)
        truth[M] = (m_min, m_max, bs, n_inj)
        for shape, tmpl in SHAPES.items():
            jobs.append((M, shape, tmpl.format(M=M), m_min, m_max, bs, n_inj))

    print("Window / injected-in-window truth:")
    for M in MASSES:
        m_min, m_max, bs, n_inj = truth[M]
        print(f"  M{M}: [{m_min}, {m_max}]  bin={bs}  N_inj_in_window={n_inj}")

    print(f"\nLaunching {len(jobs)} fits (3 concurrent)...")
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(run_one, jobs))

    # ---- collect results ----
    print("\n" + "=" * 78)
    print(f"{'Mass':<6}{'shape':<8}{'N_inj':<8}{'N_extr':<14}{'r':<8}{'pull':<8}{'Z':<7}{'sbfit p':<8}")
    print("-" * 78)
    for M in MASSES:
        m_min, m_max, bs, n_inj = truth[M]
        for shape in SHAPES:
            outdir = os.path.join(OUTBASE, f"M{M}_{shape}")
            rj = os.path.join(outdir, f"fit_results_{M}.0.json")
            if not os.path.exists(rj):
                print(f"{('M%d'%M):<6}{shape:<8}{n_inj:<8}{'FIT FAILED'}")
                continue
            r = json.load(open(rj))
            n_extr = r["obs_excess_events"]
            n_err  = r["obs_excess_events_unc"]
            sig_norm = n_inj  # by construction
            rval = n_extr / sig_norm if sig_norm else float("nan")
            pull = (n_extr - n_inj) / n_err if n_err > 0 else float("nan")
            print(f"{('M%d'%M):<6}{shape:<8}{n_inj:<8}"
                  f"{('%.0f +/- %.0f' % (n_extr, n_err)):<14}"
                  f"{rval:<8.2f}{pull:<8.2f}{r.get('signif',-1):<7.2f}"
                  f"{r.get('sbfit_prob',-1):<8.3f}")
    print("=" * 78)


if __name__ == "__main__":
    main()
