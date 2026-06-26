"""
Closure test: exclude MS20 from the 2G interpolation training, regenerate
the signal shape at 20 GeV, and compare to the true fitted DCB.
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# ---------------------------------------------------------------------------
# Double Crystal Ball PDF (analytic, matches RooDoubleCB convention)
# ---------------------------------------------------------------------------
def dcb(x, mean, sigma, alpha, n, alpha2, n2):
    t = (x - mean) / sigma
    A_l = (n / alpha) ** n * np.exp(-0.5 * alpha ** 2)
    B_l = n / alpha - alpha
    A_r = (n2 / alpha2) ** n2 * np.exp(-0.5 * alpha2 ** 2)
    B_r = n2 / alpha2 - alpha2
    return np.where(
        t < -alpha,
        A_l * (B_l - t) ** (-n),
        np.where(t > alpha2, A_r * (B_r + t) ** (-n2), np.exp(-0.5 * t ** 2)),
    )


def dcb_norm(mean, sigma, alpha, n, alpha2, n2, lo, hi, npts=5000):
    """Numerical normalization of DCB over [lo, hi]."""
    xs = np.linspace(lo, hi, npts)
    ys = dcb(xs, mean, sigma, alpha, n, alpha2, n2)
    return np.trapz(ys, xs)


def eval_dcb_pdf(xs, p):
    vals = dcb(xs, p["mean"], p["sigma"], p["alpha"], p["sign"], p["alpha2"], p["sign2"])
    lo, hi = xs[0], xs[-1]
    norm = dcb_norm(p["mean"], p["sigma"], p["alpha"], p["sign"],
                    p["alpha2"], p["sign2"], lo, hi)
    return vals / norm


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def run_interpolation(json_dir, test_mass, out_dir, fit_range=(15, 75)):
    masses_arg = str(float(test_mass))
    cmd = [
        sys.executable, "interpolation.py",
        "-s", "case",
        "--json-dir", json_dir,
        "--masses", masses_arg,
        "-m", str(fit_range[0]),
        "-M", str(fit_range[1]),
        "-o", out_dir,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr[-2000:])
        raise RuntimeError("interpolation.py failed")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", default="2G")
    parser.add_argument("--exclude-mass", type=int, default=20)
    parser.add_argument("--out", default="signal_comparison/closure_test_M{mass}.png")
    args = parser.parse_args()

    fit_dir = f"signal_fits/{args.variant}/"
    data_h5 = f"signal_data/{args.variant}/signal_mass_{args.exclude_mass}.h5"
    true_json = os.path.join(fit_dir, f"sig_fit_{args.exclude_mass}.json")
    mass = args.exclude_mass
    out_png = args.out.format(mass=mass)

    # -----------------------------------------------------------------------
    # Step 1: build a temp dir with all sig_fits except the excluded mass
    # -----------------------------------------------------------------------
    tmpdir = tempfile.mkdtemp(prefix="closure_")
    try:
        for jf in glob.glob(os.path.join(fit_dir, "sig_fit_*.json")):
            import re
            m = re.search(r"sig_fit_(\d+)\.json", os.path.basename(jf))
            if m and int(m.group(1)) == mass:
                continue
            shutil.copy(jf, tmpdir)

        all_masses = sorted(
            int(re.search(r"sig_fit_(\d+)\.json", os.path.basename(f)).group(1))
            for f in glob.glob(os.path.join(tmpdir, "sig_fit_*.json"))
        )
        print(f"Training masses (excl. {mass} GeV): {all_masses}")

        # -----------------------------------------------------------------------
        # Step 2: run interpolation without the excluded mass
        # -----------------------------------------------------------------------
        run_interpolation(
            json_dir=tmpdir,
            test_mass=mass,
            out_dir=tmpdir,
            fit_range=(min(all_masses), max(all_masses)),
        )
        interp_json = os.path.join(tmpdir, f"case_interpolation_M{float(mass)}.json")
        with open(interp_json) as f:
            p_interp = json.load(f)

    finally:
        # keep tmpdir alive until after we read the json
        pass

    with open(true_json) as f:
        p_true = json.load(f)

    print("\nTrue parameters:")
    for k in ["mean", "sigma", "alpha", "sign", "alpha2", "sign2"]:
        print(f"  {k:8s}: true={p_true[k]:.4f}  interp={p_interp[k]:.4f}  "
              f"diff={100*(p_interp[k]-p_true[k])/p_true[k]:+.1f}%")

    # -----------------------------------------------------------------------
    # Step 3: make comparison plot
    # -----------------------------------------------------------------------
    fit_lo = 0.8 * mass
    fit_hi = 1.2 * mass

    # load data
    with h5py.File(data_h5) as hf:
        all_masses_data = hf["masses"][:]
    mask = (all_masses_data >= fit_lo) & (all_masses_data <= fit_hi)
    data_in_window = all_masses_data[mask]

    nbins = 60
    counts, edges = np.histogram(data_in_window, bins=nbins, range=(fit_lo, fit_hi))
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]
    n_total = len(data_in_window)

    xs = np.linspace(fit_lo, fit_hi, 2000)
    pdf_true = eval_dcb_pdf(xs, p_true)
    pdf_interp = eval_dcb_pdf(xs, p_interp)

    # scale PDFs to match data (events per bin)
    scale = n_total * bin_width

    fig = plt.figure(figsize=(8, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    ax_r = fig.add_subplot(gs[1], sharex=ax)

    ax.errorbar(centers, counts, yerr=np.sqrt(np.maximum(counts, 1)),
                fmt="ko", ms=3, capsize=2, label="MC data")
    ax.plot(xs, scale * pdf_true, "b-", lw=2, label="True DCB fit")
    ax.plot(xs, scale * pdf_interp, "r--", lw=2, label="Interpolated DCB")
    ax.set_ylabel("Events / bin", fontsize=12)
    ax.legend(fontsize=11)
    ax.set_title(f"Interpolation closure test: {args.variant} MS{mass} excluded from interpolation",
                 fontsize=11)
    ax.set_xlim(fit_lo, fit_hi)

    # ratio: interp / true
    ratio = pdf_interp / np.maximum(pdf_true, 1e-30)
    ax_r.plot(xs, ratio, "r-", lw=1.5)
    ax_r.axhline(1, color="k", ls="--", lw=1)
    ax_r.set_ylabel("Interp / True", fontsize=11)
    ax_r.set_xlabel(r"$m_{\mu\mu}$ [GeV]", fontsize=12)
    ax_r.set_ylim(0.85, 1.15)
    ax_r.set_xlim(fit_lo, fit_hi)
    plt.setp(ax.get_xticklabels(), visible=False)

    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out_png}")

    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
