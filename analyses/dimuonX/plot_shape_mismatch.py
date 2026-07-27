#!/usr/bin/env python3
"""
Aggregate + plot the signal-shape mismatch sensitivity test.

Reads shape_mismatch_test/results_ensemble.json (K toys/cell) if present,
else results.json (single realization).  Per (injected-truth, template, mass)
cell it aggregates the toys -> mean and standard error on the mean (error bar).

Outputs in shape_mismatch_test/:
  shape_mismatch_Z.png       - 2x2 per-mass 3x3 heatmaps of mean recovered Z (+-sem)
  shape_mismatch_rbias.png   - 2x2 per-mass 3x3 heatmaps of mean recovered r (+-sem)
  shape_mismatch_Z_errbars.png   - 2x2 per-mass errorbar plots: Z vs template, per injectant
  shape_mismatch_r_errbars.png   - 2x2 per-mass errorbar plots: r vs template, per injectant
  shape_mismatch_summary.txt - text tables with mean +- sem + headline loss numbers
"""
import json
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(REPO, "shape_mismatch_test")
CLASSES = ["narrow", "medium", "wide"]
LABEL = {"narrow": "narrow\n(VLL)", "medium": "medium\n(TP)", "wide": "wide\n(H2XH3)"}
SHORT = {"narrow": "narrow", "medium": "medium", "wide": "wide"}
COLOR = {"narrow": "tab:blue", "medium": "tab:orange", "wide": "tab:red"}
MASSES = [15, 35, 55, 75]


def load():
    ens = os.path.join(OUT, "results_ensemble.json")
    single = os.path.join(OUT, "results.json")
    path = ens if os.path.exists(ens) else single
    with open(path) as f:
        d = json.load(f)
    # skip calibration sub-fits if any slipped in
    recs = [r for r in d["results"] if r.get("sub", "") in ("", None) or str(r.get("sub", "")).startswith("k")]
    return d, recs, path


def aggregate(recs, key):
    """(inj,tmpl,mass) -> (mean, sem, std, n) over toys, ignoring None/failed."""
    buckets = defaultdict(list)
    for r in recs:
        v = r.get(key)
        if r.get("ok") and v is not None and np.isfinite(v):
            buckets[(r["inj"], r["template"], r["mass"])].append(v)
    agg = {}
    for k, vals in buckets.items():
        a = np.array(vals, float)
        n = len(a)
        agg[k] = (a.mean(), a.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0,
                  a.std(ddof=1) if n > 1 else 0.0, n)
    return agg


def grid(agg, mass, idx=0):
    M = np.full((3, 3), np.nan)
    for i, inj in enumerate(CLASSES):
        for j, tmpl in enumerate(CLASSES):
            if (inj, tmpl, mass) in agg:
                M[i, j] = agg[(inj, tmpl, mass)][idx]
    return M


def _dark(im, val):
    r, g, b, _ = im.cmap(im.norm(val))
    return (0.299 * r + 0.587 * g + 0.114 * b) < 0.5


def heatmap_grid(agg, title, cmap, fname, fmt="{:.2f}", center=None, vmin=None, vmax=None):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, mass in zip(axes.flat, MASSES):
        Mv, Me = grid(agg, mass, 0), grid(agg, mass, 1)
        if center is not None:
            vabs = np.nanmax(np.abs(Mv - center)) or 0.1
            im = ax.imshow(Mv, cmap=cmap, vmin=center - vabs, vmax=center + vabs)
        else:
            im = ax.imshow(Mv, cmap=cmap, vmin=vmin, vmax=vmax)
        for i in range(3):
            for j in range(3):
                if np.isnan(Mv[i, j]):
                    continue
                txt = fmt.format(Mv[i, j]) + f"\n$\\pm${Me[i,j]:.2f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=9,
                        color="white" if _dark(im, Mv[i, j]) else "black",
                        fontweight="bold" if i == j else "normal")
                if i == j:
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                               edgecolor="lime", lw=2.4))
        ax.set_xticks(range(3)); ax.set_xticklabels([LABEL[c] for c in CLASSES], fontsize=8)
        ax.set_yticks(range(3)); ax.set_yticklabels([LABEL[c] for c in CLASSES], fontsize=8)
        ax.set_xlabel("template shape", fontsize=9)
        ax.set_ylabel("injected true shape", fontsize=9)
        ax.set_title(f"$m_S$ = {mass} GeV", fontsize=11)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = os.path.join(OUT, fname)
    fig.savefig(p, dpi=150); plt.close(fig)
    print("saved", p)


def errbar_grid(agg, ylabel, title, fname, hline=None):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True)
    x = np.arange(3)
    for ax, mass in zip(axes.flat, MASSES):
        for inj in CLASSES:
            y = [agg.get((inj, t, mass), (np.nan,))[0] for t in CLASSES]
            e = [agg.get((inj, t, mass), (np.nan, 0))[1] for t in CLASSES]
            ax.errorbar(x, y, yerr=e, marker="o", ms=6, capsize=4, lw=1.6,
                        color=COLOR[inj], label=f"inject {SHORT[inj]}")
            # mark the matched (correct-template) point with a ring
            j = CLASSES.index(inj)
            ax.scatter([x[j]], [y[j]], s=170, facecolors="none",
                       edgecolors=COLOR[inj], lw=2.2, zorder=5)
        if hline is not None:
            ax.axhline(hline, color="gray", ls="--", lw=1)
        ax.set_xticks(x); ax.set_xticklabels([SHORT[c] for c in CLASSES])
        ax.set_xlabel("template shape used in fit")
        ax.set_ylabel(ylabel)
        ax.set_title(f"$m_S$ = {mass} GeV", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(title + "  (open ring = matched template)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = os.path.join(OUT, fname)
    fig.savefig(p, dpi=150); plt.close(fig)
    print("saved", p)


def main():
    d, recs, path = load()
    ntoys = d.get("ntoys", 1)
    zagg = aggregate(recs, "signif")
    ragg = aggregate(recs, "r_fit")

    heatmap_grid(zagg, f"Mean recovered discovery Z  (K={ntoys} toys/cell, green box = matched)",
                 "viridis", "shape_mismatch_Z.png", fmt="{:.2f}", vmin=0)
    heatmap_grid(ragg, f"Mean recovered signal strength r  (1.0 = correct yield; K={ntoys})",
                 "RdBu_r", "shape_mismatch_rbias.png", fmt="{:.2f}", center=1.0)
    errbar_grid(zagg, "recovered Z",
                f"Discovery significance vs template shape  (K={ntoys} toys, error = SEM)",
                "shape_mismatch_Z_errbars.png", hline=d.get("target_Z"))
    errbar_grid(ragg, "recovered r",
                f"Signal-strength recovery vs template shape  (K={ntoys} toys, error = SEM)",
                "shape_mismatch_r_errbars.png", hline=1.0)

    # ---- text summary ----
    L = ["SIGNAL-SHAPE MISMATCH SENSITIVITY TEST",
         f"source: {os.path.basename(path)}   K={ntoys} toys/cell   "
         f"target Z={d.get('target_Z')}   calibrated={d.get('calibrated')}",
         "classes: narrow=VLL  medium=TP  wide=H2XH3   (error = std/sqrt(K))", ""]
    for mass in MASSES:
        L.append(f"=== m_S = {mass} GeV ===")
        for name, agg in (("Z", zagg), ("r", ragg)):
            L.append(f"  recovered {name}  (rows=injected truth, cols=template):")
            L.append("             " + "".join(f"{c:>14}" for c in CLASSES))
            for inj in CLASSES:
                cells = []
                for t in CLASSES:
                    if (inj, t, mass) in agg:
                        m_, e_, *_ = agg[(inj, t, mass)]
                        cells.append(f"{m_:6.2f}+-{e_:4.2f}")
                    else:
                        cells.append("     n/a    ")
                L.append(f"   {inj:>7}  " + "  ".join(cells))
        for inj in CLASSES:
            if (inj, inj, mass) not in zagg:
                continue
            zdiag = zagg[(inj, inj, mass)][0]
            zoff = {t: zagg[(inj, t, mass)][0] for t in CLASSES if (inj, t, mass) in zagg}
            wt = min(zoff, key=zoff.get)
            roff = {t: ragg[(inj, t, mass)][0] for t in CLASSES if (inj, t, mass) in ragg}
            wr = min(roff, key=roff.get)
            L.append(f"   inject {inj:>7}: matched Z={zdiag:.2f}; worst-template={wt} "
                     f"Z={zoff[wt]:.2f} ({100*(1-zoff[wt]/zdiag):+.0f}% Z); "
                     f"worst yield r={roff[wr]:.2f} (tmpl {wr})")
        L.append("")
    txt = "\n".join(L)
    with open(os.path.join(OUT, "shape_mismatch_summary.txt"), "w") as f:
        f.write(txt)
    print("\n" + txt)


if __name__ == "__main__":
    main()
