"""
Histogram the signal-yield pull (fitted - true)/unc, SEPARATED BY MASS, for both
the background-only and S+B scans, for a chosen bin size.

  - background-only : true signal = 0        -> pull = N_fit / sigma
  - S+B injection   : true = injected count  -> pull = (N_fit - N_inj) / sigma

Layout: one figure per bin size, rows = mass, cols = [bkg-only, S+B].
Usage:  python plot_signal_pulls.py <binfrac>   (e.g. 0.5 or 0.25)
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BINF = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
BKG_DIR = "ensemble_N7_b%g" % BINF
SB_DIR = "ensemble_sb_N7_b%g" % BINF
MASSES = [15, 35, 55, 75]
SEEDS = range(20)
OUT = "signal_pulls_bymass_N7_b%g.png" % BINF


def pulls(root, mass, sb):
    ps = []
    for k in SEEDS:
        base = os.path.join(root, "seed%d" % k, "M%d" % mass)
        rj = os.path.join(base, "fit_results_%.1f.json" % mass)
        if not os.path.exists(rj):
            continue
        r = json.load(open(rj))
        nf, u = r.get("obs_excess_events"), r.get("obs_excess_events_unc")
        if nf is None or not u or u <= 0:
            continue
        if sb:
            tj = os.path.join(base, "truth.json")
            if not os.path.exists(tj):
                continue
            nt = json.load(open(tj))["n_inj"]
        else:
            nt = 0.0
        ps.append((nf - nt) / u)
    return np.array(ps)


def draw(ax, arr, title, color):
    bins = np.linspace(-4, 4, 17)
    if len(arr):
        ax.hist(arr, bins=bins, color=color, edgecolor="white", linewidth=0.4)
    x = np.linspace(-4, 4, 200)
    binw = bins[1] - bins[0]
    ax.plot(x, len(arr) * binw * np.exp(-0.5 * x ** 2) / np.sqrt(2 * np.pi),
            "r--", lw=1.5)
    ax.axvline(0, color="k", lw=0.8, alpha=0.5)
    if len(arr):
        ax.axvline(arr.mean(), color=color, lw=1.5, ls=":")
        ax.set_title("%s   N=%d, mean=%+.2f, std=%.2f"
                     % (title, len(arr), arr.mean(), arr.std()), fontsize=9)
    else:
        ax.set_title("%s   (no data)" % title, fontsize=9)
    ax.set_xlim(-4, 4)


def main():
    fig, axes = plt.subplots(len(MASSES), 2, figsize=(11, 3.0 * len(MASSES)),
                             sharex=True)
    have_sb = os.path.isdir(SB_DIR)
    for i, m in enumerate(MASSES):
        pb = pulls(BKG_DIR, m, sb=False)
        draw(axes[i][0], pb, "M%d  bkg-only" % m, "#3b6ea5")
        if have_sb:
            ps = pulls(SB_DIR, m, sb=True)
            draw(axes[i][1], ps, "M%d  S+B" % m, "#2c8c5a")
        else:
            axes[i][1].set_title("M%d  S+B (ensemble not ready)" % m, fontsize=9)
        axes[i][0].set_ylabel("fits")
    for j in range(2):
        axes[-1][j].set_xlabel("signal pull  (fitted $-$ true)/$\\sigma$")
    fig.suptitle("Signal-yield pulls by mass  ($\\pm$7$\\sigma$ window, %g$\\sigma$ bins)"
                 "   red = unit Gaussian" % BINF, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT, dpi=120)
    print("Saved %s" % OUT)
    # numeric summary
    for m in MASSES:
        pb = pulls(BKG_DIR, m, sb=False)
        row = "  M%-3d bkg: mean=%+.2f std=%.2f (N=%d)" % (m, pb.mean(), pb.std(), len(pb)) if len(pb) else "  M%-3d bkg: --" % m
        if have_sb:
            ps = pulls(SB_DIR, m, sb=True)
            if len(ps):
                row += "   | S+B: mean=%+.2f std=%.2f (N=%d)" % (ps.mean(), ps.std(), len(ps))
        print(row)


if __name__ == "__main__":
    main()
