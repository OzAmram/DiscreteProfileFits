import argparse
import glob
import json
import os
import re

import matplotlib.pyplot as plt
import numpy as np

DCB_PARAMS = ["mean", "sigma", "alpha", "sign", "alpha2", "sign2"]

PARAM_LABELS = {
    "mean":   r"Mean $\mu$ [GeV]",
    "sigma":  r"Width $\sigma$ [GeV]",
    "alpha":  r"$\alpha$ (left tail)",
    "sign":   r"$n$ (left power)",
    "alpha2": r"$\alpha_2$ (right tail)",
    "sign2":  r"$n_2$ (right power)",
}

COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
          "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]


def load_fit_dir(fit_dir, mass_min=None, mass_max=None):
    json_files = sorted(glob.glob(os.path.join(fit_dir, "sig_fit_*.json")))
    data = {}
    for jf in json_files:
        m = re.search(r"sig_fit_(\d+(?:\.\d+)?)\.json", os.path.basename(jf))
        if m is None:
            continue
        mass = float(m.group(1))
        if mass_min is not None and mass < mass_min:
            continue
        if mass_max is not None and mass > mass_max:
            continue
        with open(jf) as f:
            data[mass] = json.load(f)
    return data


def load_interp_dir(fit_dir, mass_min=None, mass_max=None):
    """Load the fine-grid interpolation curve (case_interpolation_M*.json)."""
    grid = {}
    for jf in glob.glob(os.path.join(fit_dir, "case_interpolation_M*.json")):
        m = re.search(r"_M([0-9.]+)\.json", os.path.basename(jf))
        if m is None:
            continue
        mass = float(m.group(1))
        if mass_min is not None and mass < mass_min:
            continue
        if mass_max is not None and mass > mass_max:
            continue
        with open(jf) as f:
            grid[mass] = json.load(f)
    return grid


def main():
    parser = argparse.ArgumentParser(
        description="Compare DCB signal shape parameters across variants"
    )
    parser.add_argument("--fit-dirs", nargs="+", required=True,
                        help="Directories containing sig_fit_{mass}.json files (one per variant)")
    parser.add_argument("--labels", nargs="+", required=True,
                        help="Legend labels for each fit directory (same order)")
    parser.add_argument("--out", default="signal_comparison",
                        help="Output directory for comparison plots")
    parser.add_argument("--params", nargs="+", default=DCB_PARAMS,
                        help="Which DCB parameters to plot")
    parser.add_argument("--mass-min", type=float, default=None,
                        help="Exclude fit points below this mass")
    parser.add_argument("--mass-max", type=float, default=None,
                        help="Exclude fit points above this mass")
    parser.add_argument("--interp", action="store_true",
                        help="Overlay the fine-grid interpolation curve "
                             "(case_interpolation_M*.json) as a smooth line and draw "
                             "the discrete fits as markers only")
    args = parser.parse_args()

    if len(args.fit_dirs) != len(args.labels):
        parser.error("--fit-dirs and --labels must have the same length")

    os.makedirs(args.out, exist_ok=True)

    datasets = {}
    interp = {}
    for label, fit_dir in zip(args.labels, args.fit_dirs):
        datasets[label] = load_fit_dir(fit_dir, mass_min=args.mass_min, mass_max=args.mass_max)
        masses = sorted(datasets[label].keys())
        print(f"  {label}: {len(masses)} mass points  {masses[0]:.0f}–{masses[-1]:.0f} GeV")
        if args.interp:
            interp[label] = load_interp_dir(fit_dir, mass_min=args.mass_min,
                                            mass_max=args.mass_max)

    for param in args.params:
        ylabel = PARAM_LABELS.get(param, param)
        fig, axes = plt.subplots(2, 1, figsize=(8, 6),
                                 gridspec_kw={"height_ratios": [3, 1]},
                                 sharex=True)
        ax, ax_ratio = axes

        ref_label = args.labels[0]
        ref_masses = np.array(sorted(datasets[ref_label].keys()))
        ref_vals = np.array([datasets[ref_label][m].get(param, np.nan) for m in ref_masses])

        for i, label in enumerate(args.labels):
            color = COLORS[i % len(COLORS)]
            marker = MARKERS[i % len(MARKERS)]
            d = datasets[label]
            masses = np.array(sorted(d.keys()))
            vals = np.array([d[m].get(param, np.nan) for m in masses])
            errs = np.array([d[m].get(param + "-err", 0.0) for m in masses])

            if args.interp and interp.get(label):
                # discrete fits as markers only; interpolation as the smooth line
                ax.errorbar(masses, vals, yerr=errs, fmt=marker, ls="none",
                            color=color, label=label, capsize=3, markersize=5)
                gm = np.array(sorted(interp[label].keys()))
                gv = np.array([interp[label][m].get(param, np.nan) for m in gm])
                ax.plot(gm, gv, "-", color=color, lw=1.3)
            else:
                ax.errorbar(masses, vals, yerr=errs, fmt=marker + "-",
                            color=color, label=label, capsize=3, markersize=5)

            if i > 0:
                ref_interp = np.interp(masses, ref_masses, ref_vals,
                                       left=np.nan, right=np.nan)
                with np.errstate(invalid="ignore"):
                    ratio = vals / ref_interp
                ax_ratio.plot(masses, ratio, marker + "-", color=color,
                              markersize=4, label=f"{label}/{ref_label}")

        ax.set_ylabel(ylabel, fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_title(f"Signal shape parameter: {param}", fontsize=12)

        ax_ratio.axhline(1.0, color="black", lw=1, ls="--")
        ax_ratio.set_ylabel(f"Ratio to {ref_label}", fontsize=10)
        ax_ratio.set_xlabel(r"$m_S$ [GeV]", fontsize=12)
        ax_ratio.set_ylim(0.7, 1.3)
        ax_ratio.grid(True, alpha=0.3)
        if len(args.labels) > 1:
            ax_ratio.legend(fontsize=9)

        fig.tight_layout()
        out_path = os.path.join(args.out, f"{param}_vs_mass.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"  Saved {out_path}")

    print(f"\nDone — {len(args.params)} comparison plots written to {args.out}/")


if __name__ == "__main__":
    main()
