"""
Adaptive-window fit wrapper.

Fits a single mass hypothesis starting with a +/-8sigma window and 0.25sigma
bin size. If the s+b postfit chi2 p-value is below threshold, shrinks the
window down to +/-4sigma (0.5sigma steps). If no window size at the starting
bin fraction passes, coarsens the bin size (0.1sigma steps, up to 1.0sigma by
default) and repeats the window scan at each. Reports the highest-p attempt
across the whole (bin size, window) grid if nothing passes outright.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

# doFit.py is always a sibling of this file, regardless of the caller's cwd.
DOFIT_PY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "doFit.py")


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
        sys.executable, DOFIT_PY,
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
    parser.add_argument("-c", "--config", default="dimuonX_config.json",
                        help="JSON config file")
    parser.add_argument("-s", "--sig-dir", required=True,
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
                        help="Starting bin size as a fraction of sigma (default 0.25)")
    parser.add_argument("--bin-frac-max", type=float, default=1.0,
                        help="Max bin size (in sigma) to escalate to when GoF fails at "
                             "every window; coarser bins are less sensitive to narrow "
                             "fluctuations landing in one bin (default 1.0)")
    parser.add_argument("--bin-frac-step", type=float, default=0.1,
                        help="Bin-size escalation step in sigma (default 0.1)")
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

    # bin sizes (in sigma) to try, escalating from the default only on failure
    bin_fracs = []
    bf = args.bin_frac
    while bf <= max(args.bin_frac_max, args.bin_frac) + 1e-9:
        bin_fracs.append(round(bf, 4))
        bf += args.bin_frac_step

    print(f"\nM = {mass} GeV   sigma = {sigma:.3f} GeV")
    print(f"Window scan: ±{args.n_sigma_start}σ → ±{args.n_sigma_min}σ "
          f"(step {args.n_sigma_step}σ);  bin size {args.bin_frac}σ → "
          f"{bin_fracs[-1]}σ (step {args.bin_frac_step}σ) on failure;  "
          f"p-value threshold = {args.pval_thresh}")
    print("="*60)

    success    = False
    final_dir  = final_nsig = final_bf = final_pval = None
    # Fallback: the highest-p attempt across the whole (bin size, window) grid. If
    # GoF is never achieved we keep this as 'best' and report its poor p-value,
    # rather than dropping the mass point from the scan entirely.
    best_dir = best_pval = best_nsig = best_bf = None

    for bf in bin_fracs:
        n_sigma = args.n_sigma_start
        while n_sigma >= args.n_sigma_min - 1e-9:
            m_min, m_max, bin_size = compute_window(
                mass, sigma, n_sigma, args.m_data_min, args.m_data_max, bf)

            print(f"\n  bin={bf:.1f}σ  ±{n_sigma:.1f}σ  →  [{m_min}, {m_max}] GeV,  "
                  f"bin_size={bin_size} GeV")

            attempt_dir = os.path.join(outbase, f"bin{bf:.1f}_window{n_sigma:.1f}sig")
            pval = run_fit(mass, m_min, m_max, bin_size, args.config, sig_json,
                           attempt_dir, input_file=args.input)

            if pval is None:
                print(f"  --> fit failed, shrinking window")
            else:
                if best_pval is None or pval > best_pval:
                    best_dir, best_pval, best_nsig, best_bf = attempt_dir, pval, n_sigma, bf
                print(f"  --> sbfit chi2 p-value = {pval:.4f}", end="")
                if pval >= args.pval_thresh:
                    print(f"  PASS ✓")
                    success = True
                    final_dir, final_nsig, final_bf, final_pval = attempt_dir, n_sigma, bf, pval
                    break
                else:
                    print(f"  FAIL (< {args.pval_thresh}), shrinking window")

            n_sigma -= args.n_sigma_step

        if success:
            break
        if bf != bin_fracs[-1]:
            print(f"\n  --> window scan failed at bin={bf:.1f}σ; increasing bin size to "
                  f"{round(bf + args.bin_frac_step, 4)}σ")

    print("\n" + "="*60)
    best_out = os.path.join(outbase, "best")
    if not success and best_dir is not None:
        # GoF never reached threshold: keep the best (highest-p) attempt overall.
        final_dir, final_nsig, final_bf, final_pval = best_dir, best_nsig, best_bf, best_pval

    if final_dir is not None:
        if os.path.exists(best_out):
            shutil.rmtree(best_out)
        shutil.copytree(final_dir, best_out)
        win_min, win_max, _ = compute_window(
            mass, sigma, final_nsig, args.m_data_min, args.m_data_max, final_bf)
        tag = "SUCCESS" if success else "POOR-GOF (kept)"
        print(f"{tag}  M={mass} GeV  bin={final_bf:.1f}σ  ±{final_nsig:.1f}σ  p={final_pval:.4f}")
        print(f"  Winning window : [{win_min}, {win_max}] GeV")
        print(f"  Results        : {best_out}")
        if not success:
            print(f"  NOTE: p < {args.pval_thresh} even at bin={bin_fracs[-1]}σ / "
                  f"±{args.n_sigma_min:.1f}σ; kept the highest-p attempt and reporting its p-value.")
    else:
        print(f"FAILURE  M={mass} GeV  no usable fit produced (all attempts crashed)")
    print("="*60)

    # exit 0 whenever a usable best/ exists (passing OR poor-GoF-but-kept);
    # exit 1 only if no fit could be produced at all.
    return 0 if final_dir is not None else 1


if __name__ == "__main__":
    sys.exit(main())
