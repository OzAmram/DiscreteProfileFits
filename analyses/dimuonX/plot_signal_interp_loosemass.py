"""
Summary plots of the interpolated DCB parameters vs mass for the loosemass signals.

For each group under signal_fits_loosemass/<group>/ this overlays the interpolated
parameter curve (from case_interpolation_M<mass>.json, the 0.1 GeV grid) on top of
the discrete fit points (from sig_fit_<M>.json) for all six DCB parameters:
  <group>/interp_summary.png   (2x3 panel: curve + fit markers, one per group)
It also writes one combined figure overlaying every group's interpolated curve:
  signal_fits_loosemass/interp_summary_allgroups.png
"""
import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTROOT = "signal_fits_loosemass"
PARAMS = ["mean", "sigma", "alpha", "sign", "alpha2", "sign2"]
PLABEL = {"mean": "mean [GeV]", "sigma": r"$\sigma$ [GeV]", "alpha": r"$\alpha$",
          "sign": "n", "alpha2": r"$\alpha_2$", "sign2": r"$n_2$"}

groups = sorted(d for d in glob.glob(os.path.join(OUTROOT, "*"))
                if os.path.isdir(d) and glob.glob(os.path.join(d, "sig_fit_*.json")))


def load_fit_points(gdir):
    pts = {}
    for jf in glob.glob(os.path.join(gdir, "sig_fit_*.json")):
        d = json.load(open(jf))
        pts[float(d["mass"])] = d
    ms = np.array(sorted(pts))
    return ms, {p: np.array([pts[m][p] for m in ms]) for p in PARAMS}, \
        {p: np.array([pts[m].get(p + "-err", 0.0) for m in ms]) for p in PARAMS}


def load_interp_curve(gdir):
    grid = {}
    for jf in glob.glob(os.path.join(gdir, "case_interpolation_M*.json")):
        m = re.search(r"_M([0-9.]+)\.json", os.path.basename(jf))
        if not m:
            continue
        grid[float(m.group(1))] = json.load(open(jf))
    ms = np.array(sorted(grid))
    return ms, {p: np.array([grid[m][p] for m in ms]) for p in PARAMS}


cmap = plt.get_cmap("tab10")
allcurves = {}
for gi, gdir in enumerate(groups):
    g = os.path.basename(gdir)
    fm, fv, fe = load_fit_points(gdir)
    im, iv = load_interp_curve(gdir)
    allcurves[g] = (im, iv)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, p in zip(axes.ravel(), PARAMS):
        ax.plot(im, iv[p], "-", color="C0", lw=1.5, label="interpolation", zorder=1)
        ax.errorbar(fm, fv[p], yerr=fe[p], fmt="o", ms=5, color="k",
                    capsize=2, label="fit points", zorder=2)
        ax.set_xlabel("resonance mass [GeV]")
        ax.set_ylabel(PLABEL[p])
        ax.grid(alpha=0.3)
    axes.ravel()[0].legend(loc="best", fontsize=9)
    fig.suptitle(g, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(gdir, "interp_summary.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("wrote %s" % out)

# combined overlay: one panel per param, one curve per group
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
for ax, p in zip(axes.ravel(), PARAMS):
    for gi, g in enumerate(sorted(allcurves)):
        im, iv = allcurves[g]
        ax.plot(im, iv[p], "-", color=cmap(gi % 10), lw=1.3, label=g)
    ax.set_xlabel("resonance mass [GeV]")
    ax.set_ylabel(PLABEL[p])
    ax.grid(alpha=0.3)
axes.ravel()[0].legend(loc="upper left", fontsize=6, ncol=1)
fig.suptitle("Interpolated DCB parameters vs mass — all loosemass groups", fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(OUTROOT, "interp_summary_allgroups.png")
fig.savefig(out, dpi=120)
plt.close(fig)
print("wrote %s" % out)
