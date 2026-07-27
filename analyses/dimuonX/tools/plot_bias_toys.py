"""
Plot the combine-toy bias test (tools/bias_test_toys.py).

Writes note/Figures/8Systematics_figures/Fitting/bias_toys.png -- the median pull
for every (generating family, injected signal) combination at each mass. The
conventional acceptance criterion for discrete profiling is |median pull| < 0.5,
which is drawn as a band.

Run from the repo root:  python3 tools/plot_bias_toys.py
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FAM_COLOR = {"bern": "#9467bd", "polyExp": "#d62728",
             "exp": "#2ca02c", "expPoly": "#1f77b4"}
FAM_LABEL = {"bern": "Bernstein", "polyExp": r"Poly $\times$ exp",
             "exp": "Sum of exp", "expPoly": "Exp of poly"}
THRESH = 0.5


def median_err(p):
    """Uncertainty on the median (~1.253 sigma/sqrt(N) for a normal)."""
    return 1.253 * np.std(p) / np.sqrt(len(p)) if len(p) > 1 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="bias_toys/bias_toys_summary.json")
    ap.add_argument("--out", default="note/Figures/8Systematics_figures/Fitting/bias_toys.png")
    args = ap.parse_args()

    d = json.load(open(args.summary))
    res = d["results"]
    masses = sorted({r["mass"] for r in res})
    zs = sorted({r["z"] for r in res})
    fams = [f for f in ["bern", "polyExp", "exp", "expPoly"]
            if any(r["gen_fam"] == f for r in res)]

    fig, axes = plt.subplots(1, len(zs), figsize=(4.6 * len(zs), 4.4), sharey=True)
    if len(zs) == 1:
        axes = [axes]

    worst = 0.0
    for ax, z in zip(axes, zs):
        for i, fam in enumerate(fams):
            xs, ys, es = [], [], []
            for m in masses:
                sel = [r for r in res if r["mass"] == m and r["z"] == z
                       and r["gen_fam"] == fam]
                if not sel or not sel[0]["pull"]:
                    continue
                p = np.array(sel[0]["pull"])
                # small horizontal offset per family so points don't overlap
                xs.append(m + (i - 1.5) * 1.2)
                ys.append(np.median(p))
                es.append(median_err(p))
                worst = max(worst, abs(np.median(p)))
            ax.errorbar(xs, ys, yerr=es, fmt="o", ms=6, capsize=3,
                        color=FAM_COLOR[fam], label=FAM_LABEL[fam])
        ax.axhspan(-THRESH, THRESH, color="gray", alpha=0.18, zorder=0)
        ax.axhline(0, color="k", lw=1)
        ax.axhline(THRESH, color="gray", ls="--", lw=1)
        ax.axhline(-THRESH, color="gray", ls="--", lw=1)
        ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
        ax.set_title("Injected signal: %s" % ("none" if z == 0 else r"$\sim%g\sigma$" % z),
                     fontsize=11)
        ax.set_ylim(-1.2, 1.2)
        ax.grid(alpha=0.3, axis="y")
    axes[0].set_ylabel("Median pull  " + r"$(r_{\mathrm{fit}} - r_{\mathrm{gen}})/\sigma_r$")
    axes[0].legend(fontsize=9, ncol=2, title="Generating function", title_fontsize=9)

    fig.suptitle("Bias test: toys generated from each background function, "
                 "fit with the full envelope", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print("wrote", args.out)
    print("largest |median pull| = %.3f  (criterion: < %.1f)" % (worst, THRESH))

    # printed table for the note text
    print("\n%-6s %-10s %-6s %6s %8s %8s" % ("mass", "gen", "z", "ntoys", "median", "width"))
    for m in masses:
        for fam in fams:
            for z in zs:
                sel = [r for r in res if r["mass"] == m and r["z"] == z
                       and r["gen_fam"] == fam]
                if not sel or not sel[0]["pull"]:
                    continue
                p = np.array(sel[0]["pull"])
                print("%-6d %-10s %-6g %6d %+8.3f %8.3f"
                      % (m, fam, z, len(p), np.median(p), np.std(p)))


if __name__ == "__main__":
    main()
