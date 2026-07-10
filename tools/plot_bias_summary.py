"""
Summary: the low-mass 'spurious signal' is a parent-MC fluctuation, not method bias.

Left  : per-mass spurious pull mean for the 3% subsample ensemble (every draw
        inherits the SAME parent MC, so a local MC fluctuation biases all of them)
        vs smooth-truth toys (bias of the fit method alone).
Right : full-MC (3.76M) background-only residual at M35, showing the ~1.5-sigma
        dip near 34-35 GeV that the smooth PDF cannot absorb.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

MASSES = [15, 35, 55]


def pull_mean(root, m, sb=False):
    ps = []
    for k in range(20):
        base = os.path.join(root, "seed%d" % k, "M%d" % m)
        rj = os.path.join(base, "fit_results_%d.0.json" % m)
        if not os.path.exists(rj):
            continue
        r = json.load(open(rj))
        u = r.get("obs_excess_events_unc")
        if not u or u <= 0:
            continue
        nt = 0.0
        if sb:
            tj = os.path.join(base, "truth.json")
            nt = json.load(open(tj))["n_inj"] if os.path.exists(tj) else 0.0
        ps.append((r["obs_excess_events"] - nt) / u)
    a = np.array(ps)
    return (a.mean(), a.std() / np.sqrt(len(a))) if len(a) else (np.nan, np.nan)


fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

# ---- Left: grouped bars ----
sub = [pull_mean("ensemble_N7_b0.5", m) for m in MASSES]
smo = [pull_mean("smooth_toy_N7_b0.5", m) for m in MASSES]
x = np.arange(len(MASSES))
w = 0.36
axL.axhspan(-0.2, 0.2, color="0.85", zorder=0, label="|mean| < 0.2 (negligible)")
axL.bar(x - w / 2, [s[0] for s in sub], w, yerr=[s[1] for s in sub],
        capsize=4, color="#3b6ea5", label="3% subsamples (inherit parent MC)")
axL.bar(x + w / 2, [s[0] for s in smo], w, yerr=[s[1] for s in smo],
        capsize=4, color="#c26b3e", label="smooth-truth toys (method only)")
axL.axhline(0, color="k", lw=0.8)
axL.set_xticks(x); axL.set_xticklabels(["M%d" % m for m in MASSES])
axL.set_ylabel("spurious signal pull  (N$_{fit}$/$\\sigma$), mean $\\pm$ err")
axL.set_title("Low-mass bias vanishes with a smooth truth\n($\\pm$7$\\sigma$, 0.5$\\sigma$ bins)")
axL.legend(fontsize=9, frameon=False, loc="lower left")
for i, (s, m2) in enumerate(zip(sub, smo)):
    axL.text(x[i] - w / 2, s[0] - 0.06, "%+.2f" % s[0], ha="center", va="top", fontsize=8)
    axL.text(x[i] + w / 2, m2[0] - 0.06, "%+.2f" % m2[0], ha="center", va="top", fontsize=8)

# ---- Right: full-MC M35 residual ----
allm = h5py.File("bkg_mc_masses.h5")["masses"][:]
lo, hi = 30.9693, 39.0307
sel = allm[(allm >= lo) & (allm <= hi)]
nb = 28
h, edges = np.histogram(sel, bins=nb, range=(lo, hi))
ctr = 0.5 * (edges[:-1] + edges[1:])
# smooth reference from a wide fit
wlo, whi = 35 - 18, 35 + 18
wsel = allm[(allm >= wlo) & (allm <= whi)]
hw, ew = np.histogram(wsel, bins=140, range=(wlo, whi))
cw = 0.5 * (ew[:-1] + ew[1:]); mask = hw > 0
xs = (cw - 35) / 18
coef = np.polyfit(xs[mask], np.log(hw[mask]), 4)
binw_wide = ew[1] - ew[0]
ref = np.exp(np.polyval(coef, (ctr - 35) / 18)) / binw_wide * (edges[1] - edges[0])
resid = h - ref
axR.errorbar(ctr, resid, yerr=np.sqrt(h), fmt="o", color="k", ms=4, capsize=2)
axR.axhline(0, color="0.5", lw=1)
sig = 0.576
axR.axvspan(35 - 2 * sig, 35 + 2 * sig, color="#c26b3e", alpha=0.18,
            label="signal $\\pm2\\sigma$ window")
axR.set_xlabel("m$_{\\mu\\mu}$ [GeV]")
axR.set_ylabel("full-MC counts $-$ smooth reference")
axR.set_title("Parent MC (3.76M) has a $\\sim$1.5$\\sigma$ dip at 34-35 GeV\n"
              "$\\to$ N$_{spur}$ = $-$500$\\pm$171 (z=$-$2.9), inherited by all subsamples")
axR.legend(fontsize=9, frameon=False)

fig.tight_layout()
fig.savefig("bias_summary.png", dpi=130)
print("Saved bias_summary.png")
print("subsample:", [(m, "%.2f" % s[0]) for m, s in zip(MASSES, sub)])
print("smooth   :", [(m, "%.2f" % s[0]) for m, s in zip(MASSES, smo)])
