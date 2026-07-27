"""
Leave-one-out interpolation closure test for the loosemass DCB signal shapes.

For every group and every *interior* fit mass point, remove that point, rebuild the
interpolation from the remaining points using the exact 'case' recipe
(mean,sigma -> TSpline3; alpha,sign,alpha2,sign2 -> pol1 fit), evaluate the
parameters at the held-out mass, and compare to the directly fitted values there.
Endpoints (lowest/highest mass) are skipped because holding them out would be
extrapolation, not interpolation, for the spline.

Outputs:
  signal_loo_interp/loo_residuals.png   (per-parameter residual vs mass, all groups)
  signal_loo_interp/loo_summary.json    (per-parameter median/max residuals)
and prints a summary table. The physically important closure metric is the peak
position (mean) and width (sigma): a good interpolation reproduces the held-out
sigma to well within its fit uncertainty.
"""
import glob
import json
import os
from array import array

import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

OUTROOT = "signal_fits_loosemass"
OUTDIR = "signal_loo_interp"
os.makedirs(OUTDIR, exist_ok=True)

SPLINE = ["mean", "sigma"]
POL1 = ["alpha", "sign", "alpha2", "sign2"]
PARAMS = SPLINE + POL1
PLABEL = {"mean": "mean [GeV]", "sigma": r"$\sigma$ [GeV]", "alpha": r"$\alpha$",
          "sign": "n", "alpha2": r"$\alpha_2$", "sign2": r"$n_2$"}
FIT_LO, FIT_HI = 10.0, 100.0

groups = sorted(d for d in glob.glob(os.path.join(OUTROOT, "*"))
                if os.path.isdir(d) and glob.glob(os.path.join(d, "sig_fit_*.json")))


def load_group(gdir):
    pts = {}
    for jf in glob.glob(os.path.join(gdir, "sig_fit_*.json")):
        d = json.load(open(jf))
        pts[float(d["mass"])] = d
    ms = sorted(pts)
    return ms, pts


def interp_at(masses, vals, errs, param, target):
    """Rebuild the 'case' interpolation from (masses, vals) and eval at target."""
    x = array("d", masses)
    if param in SPLINE:
        g = ROOT.TGraph(len(masses), x, array("d", vals))
        sp = ROOT.TSpline3("sp_%s" % param, g)
        return sp.Eval(target)
    else:
        ge = ROOT.TGraphErrors(len(masses), x, array("d", vals),
                               array("d", [0.0] * len(masses)), array("d", errs))
        f = ROOT.TF1("f_%s" % param, "pol1", FIT_LO, FIT_HI)
        for _ in range(3):
            ge.Fit(f, "QN0", "", FIT_LO, FIT_HI)
        return f.Eval(target)


# residuals[param] = list of (mass, true, pred, err, group)
residuals = {p: [] for p in PARAMS}
counter = 0
for gdir in groups:
    g = os.path.basename(gdir)
    ms, pts = load_group(gdir)
    interior = ms[1:-1]  # skip endpoints (extrapolation for spline)
    for held in interior:
        keep = [m for m in ms if m != held]
        for p in PARAMS:
            vals = [pts[m][p] for m in keep]
            errs = [pts[m].get(p + "-err", 0.0) for m in keep]
            counter += 1
            # unique ROOT object names per call to avoid clashes
            pred = interp_at(keep, vals, errs, p, held)
            true = pts[held][p]
            err = pts[held].get(p + "-err", 0.0)
            residuals[p].append((held, true, pred, err, g))

# ---- summary ----
summary = {}
print("\n%-8s %10s %10s %10s %10s" %
      ("param", "med|d|", "max|d|", "med rel%", "max rel%"))
print("-" * 52)
for p in PARAMS:
    arr = residuals[p]
    d = np.array([pr - tr for (_, tr, pr, _, _) in arr])
    rel = np.array([100.0 * (pr - tr) / tr if tr != 0 else 0.0
                    for (_, tr, pr, _, _) in arr])
    summary[p] = {"median_abs": float(np.median(np.abs(d))),
                  "max_abs": float(np.max(np.abs(d))),
                  "median_rel_pct": float(np.median(np.abs(rel))),
                  "max_rel_pct": float(np.max(np.abs(rel))),
                  "n": len(arr)}
    print("%-8s %10.4f %10.4f %10.3f %10.3f" %
          (p, summary[p]["median_abs"], summary[p]["max_abs"],
           summary[p]["median_rel_pct"], summary[p]["max_rel_pct"]))
print("-" * 52)
print("total held-out points per parameter: %d" % len(residuals["mean"]))
json.dump(summary, open(os.path.join(OUTDIR, "loo_summary.json"), "w"), indent=2)

# ---- plot: relative residual vs mass, per parameter, all groups overlaid ----
cmap = plt.get_cmap("tab10")
gnames = sorted({r[4] for r in residuals["mean"]})
gcolor = {g: cmap(i % 10) for i, g in enumerate(gnames)}

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
for ax, p in zip(axes.ravel(), PARAMS):
    for (m, tr, pr, er, g) in residuals[p]:
        rel = 100.0 * (pr - tr) / tr if tr != 0 else 0.0
        ax.plot(m, rel, "o", ms=4, color=gcolor[g])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("held-out mass [GeV]")
    ax.set_ylabel("(interp - fit) / fit  [%]  " + PLABEL[p])
    ax.grid(alpha=0.3)
handles = [plt.Line2D([], [], marker="o", ls="", color=gcolor[g], label=g)
           for g in gnames]
axes.ravel()[0].legend(handles=handles, fontsize=6, loc="best")
fig.suptitle("Leave-one-out interpolation closure — relative residuals", fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(os.path.join(OUTDIR, "loo_residuals.png"), dpi=120)
plt.close(fig)
print("wrote %s" % os.path.join(OUTDIR, "loo_residuals.png"))

# ---- extra: mean/sigma residual in absolute (physical) units ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, p, unit in zip(axes, ["mean", "sigma"], ["GeV", "GeV"]):
    for (m, tr, pr, er, g) in residuals[p]:
        ax.errorbar(m, pr - tr, yerr=er, fmt="o", ms=4, color=gcolor[g],
                    capsize=1.5, elinewidth=0.7)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("held-out mass [GeV]")
    ax.set_ylabel("interp - fit  [%s]  (%s)" % (unit, p))
    ax.grid(alpha=0.3)
fig.suptitle("Leave-one-out closure: mean & sigma (error bars = fit uncertainty)",
             fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(OUTDIR, "loo_mean_sigma_abs.png"), dpi=120)
plt.close(fig)
print("wrote %s" % os.path.join(OUTDIR, "loo_mean_sigma_abs.png"))
