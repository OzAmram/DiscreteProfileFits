"""
Adaptive-window fit wrapper.

Fits a single mass hypothesis starting with a ±8σ window and 0.25σ bin size.
If the s+b postfit chi2 p-value is below threshold, shrinks the window by
0.5σ on each side and retries. Gives up and reports failure at ±4σ.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys


def get_sigma(sig_dir, mass):
    sig_json = os.path.join(sig_dir, f"case_interpolation_M{float(mass)}.json")
    if not os.path.exists(sig_json):
        sys.exit(f"ERROR: Signal JSON not found: {sig_json}")
    with open(sig_json) as f:
        return json.load(f)["sigma"], sig_json


def compute_window(mass, sigma, n_sigma, m_data_min, m_data_max, bin_frac=0.25):
    bin_size = round(bin_frac * sigma, 4)
    m_min    = mass - n_sigma * sigma
    m_max    = mass + n_sigma * sigma
    if m_data_min is not None:
        m_min = max(m_min, m_data_min)
    if m_data_max is not None:
        m_max = min(m_max, m_data_max)
    return round(m_min, 4), round(m_max, 4), bin_size


def run_fit(mass, m_min, m_max, bin_size, config, sig_json, outdir, input_file=None):
    os.makedirs(outdir, exist_ok=True)

    cmd = [
        sys.executable, "doFit.py",
        "-c", config,
        "-M", str(mass),
        "--m-min",    str(m_min),
        "--m-max",    str(m_max),
        "--bin-size", str(bin_size),
        "-s", sig_json,
        "-o", outdir,
    ]
    if input_file is not None:
        cmd += ["-i", input_file]

    log_path = os.path.join(outdir, "fit_log.txt")
    print(f"    cmd: {' '.join(cmd)}")
    with open(log_path, "w") as logf:
        ret = subprocess.run(cmd, stdout=logf, stderr=logf)

    if ret.returncode != 0:
        print(f"    doFit.py exited with code {ret.returncode}")
        return None

    results_json = os.path.join(outdir, f"fit_results_{mass}.json")
    if not os.path.exists(results_json):
        print(f"    fit_results JSON not found — fit likely crashed")
        return None

    with open(results_json) as f:
        results = json.load(f)

    pval = results.get("sbfit_prob")
    if pval is None:
        print(f"    sbfit_prob not found in results JSON")
    return pval


def main():
    parser = argparse.ArgumentParser(
        description="Adaptive-window background fit for a single mass hypothesis")
    parser.add_argument("-M", "--mass", type=float, required=True,
                        help="Signal mass hypothesis (GeV)")
    parser.add_argument("-c", "--config", default="bkg_mc_config.json",
                        help="JSON config file")
    parser.add_argument("-s", "--sig-dir", default="signal_fits/2G",
                        help="Directory containing case_interpolation_M{mass}.json files")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Output directory (default: bkg_mc_fits/M{mass})")
    parser.add_argument("--pval-thresh", type=float, default=0.05,
                        help="Minimum acceptable s+b chi2 p-value (default 0.05)")
    parser.add_argument("--n-sigma-start", type=float, default=8.0,
                        help="Starting half-width in sigma units (default 8)")
    parser.add_argument("--n-sigma-min", type=float, default=4.0,
                        help="Minimum half-width before declaring failure (default 4)")
    parser.add_argument("--n-sigma-step", type=float, default=0.5,
                        help="Step size for shrinking the window (default 0.5)")
    parser.add_argument("--bin-frac", type=float, default=0.25,
                        help="Bin size as a fraction of sigma (default 0.25)")
    parser.add_argument("--m-data-min", type=float, default=None,
                        help="Lower data boundary; fit window is floored here (default: no floor)")
    parser.add_argument("--m-data-max", type=float, default=None,
                        help="Upper data boundary; fit window is capped here (default: no cap)")
    parser.add_argument("-i", "--input", default=None,
                        help="Override the config's input h5 file (e.g. a mock signal-injected dataset)")
    args = parser.parse_args()

    mass = args.mass
    sigma, sig_json = get_sigma(args.sig_dir, mass)
    outbase = args.outdir or f"bkg_mc_fits/M{int(mass)}"

    print(f"\nM = {mass} GeV   sigma = {sigma:.3f} GeV")
    print(f"Window scan: ±{args.n_sigma_start}σ → ±{args.n_sigma_min}σ "
          f"(step {args.n_sigma_step}σ),  p-value threshold = {args.pval_thresh}")
    print("="*60)

    n_sigma   = args.n_sigma_start
    success   = False
    final_dir = None

    while n_sigma >= args.n_sigma_min - 1e-9:
        m_min, m_max, bin_size = compute_window(
            mass, sigma, n_sigma, args.m_data_min, args.m_data_max, args.bin_frac)

        print(f"\n  ±{n_sigma:.1f}σ  →  [{m_min}, {m_max}] GeV,  bin_size={bin_size} GeV")

        attempt_dir = os.path.join(outbase, f"window_{n_sigma:.1f}sig")
        pval = run_fit(mass, m_min, m_max, bin_size, args.config, sig_json,
                       attempt_dir, input_file=args.input)

        if pval is None:
            print(f"  --> fit failed, shrinking window")
        else:
            print(f"  --> sbfit chi2 p-value = {pval:.4f}", end="")
            if pval >= args.pval_thresh:
                print(f"  PASS ✓")
                success   = True
                final_dir = attempt_dir
                break
            else:
                print(f"  FAIL (< {args.pval_thresh}), shrinking window")

        n_sigma -= args.n_sigma_step

    print("\n" + "="*60)
    if success:
        # Copy successful attempt to top-level outdir for easy access
        best_dir = os.path.join(outbase, "best")
        if os.path.exists(best_dir):
            shutil.rmtree(best_dir)
        shutil.copytree(final_dir, best_dir)

        win_min, win_max, _ = compute_window(
            mass, sigma, n_sigma, args.m_data_min, args.m_data_max, args.bin_frac)
        print(f"SUCCESS  M={mass} GeV  ±{n_sigma:.1f}σ  p={pval:.4f}")
        print(f"  Winning window : [{win_min}, {win_max}] GeV")
        print(f"  Results        : {best_dir}")
    else:
        print(f"FAILURE  M={mass} GeV  could not achieve p >= {args.pval_thresh} "
              f"at ±{args.n_sigma_min:.1f}σ")
    print("="*60)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
