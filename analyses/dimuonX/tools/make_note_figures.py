"""
Build the figures for the Fitting section of the analysis note.

Writes into note/Figures/8Systematics_figures/Fitting/ so the note directory
mirrors the AN layout and the files can be copied across verbatim.

Figures produced here:
  interp_params_vs_mass.png -- DCB parameters vs mass with the interpolation
                               overlaid (the scheme is smooth)
  gof_ensemble.png          -- background-only chi2/ndof over the 3% MC ensemble

The note's bias figure is NOT made here: it comes from tools/plot_bias_toys.py,
which reads the combine-toy study (tools/bias_test_toys.py). The older
subsample/injection pull plot was a proof of concept and is no longer in the
note, so pull_figure() below writes to the scratch area rather than into note/.

Run from the repo root:  python3 tools/make_note_figures.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "note/Figures/8Systematics_figures/Fitting"
SIGDIR = "signal_fits/2B_loosemass"
# The 2B interpolation is built (and the scan runs) only over this range. Signal
# MC exists up to 95 GeV, but those fits sit on the Z peak and are degenerate
# (mu plateaus near 79 GeV, tail params unconstrained), so they are neither used
# nor shown.
M_FIT_MIN, M_FIT_MAX = 10.0, 75.0
ENS_BKG = "ensemble_capped_adaptive/ensemble_summary.json"
ENS_SB = "ensemble_sb_capped_adaptive/ensemble_sb_summary.json"

# DCB parameters, with the labels used in the note's Eq. (dcb)
PARAMS = [
    ("mean", r"$\mu$ [GeV]"),
    ("sigma", r"$\sigma$ [GeV]"),
    ("alpha", r"$\alpha_1$"),
    ("alpha2", r"$\alpha_2$"),
    ("sign", r"$n_1$"),
    ("sign2", r"$n_2$"),
]


def interp_params_figure():
    """DCB params from the discrete mass points, with the interpolation overlaid."""
    # fitted points (the 5 GeV grid the signal MC was generated at)
    fits = {}
    for f in sorted(os.listdir(SIGDIR)):
        if f.startswith("sig_fit_") and f.endswith(".json"):
            m = float(f[len("sig_fit_"):-len(".json")])
            if M_FIT_MIN - 1e-9 <= m <= M_FIT_MAX + 1e-9:
                fits[m] = json.load(open(os.path.join(SIGDIR, f)))
    fit_m = np.array(sorted(fits))

    # interpolated templates (the fine grid actually used in the scan)
    interp = {}
    for f in os.listdir(SIGDIR):
        if f.startswith("case_interpolation_M") and f.endswith(".json"):
            m = float(f[len("case_interpolation_M"):-len(".json")])
            interp[m] = json.load(open(os.path.join(SIGDIR, f)))
    int_m = np.array(sorted(interp))

    fig, axes = plt.subplots(3, 2, figsize=(11, 10))
    for ax, (key, label) in zip(axes.flat, PARAMS):
        y = np.array([fits[m][key] for m in fit_m])
        err_key = key + "-err"
        yerr = np.array([fits[m].get(err_key, 0.0) for m in fit_m])
        ax.errorbar(fit_m, y, yerr=yerr, fmt="ko", ms=4, capsize=2,
                    label="Fits to signal MC", zorder=3)
        ax.plot(int_m, [interp[m][key] for m in int_m], "r-", lw=1.5,
                label="Interpolation", zorder=2)
        ax.set_xlabel(r"$m_{\mu\mu}$ [GeV]")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        # The lowest mass point has very loose tail parameters; let the y-range
        # follow the interpolation rather than those error bars.
        yi = np.array([interp[m][key] for m in int_m])
        lo, hi = yi.min(), yi.max()
        pad = max(0.35 * (hi - lo), 0.15 * abs(hi) if hi else 1.0)
        ax.set_ylim(lo - pad, hi + pad)
    axes.flat[0].legend(fontsize=10)
    fig.suptitle("Double Crystal Ball parameters vs. resonance mass (2B signal)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    p = os.path.join(OUT, "interp_params_vs_mass.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


def scan_grid_figure(plan_json="scan_plan.json"):
    """The scan grid: spacing vs mass, next to the local resolution it tracks."""
    plan = json.load(open(plan_json))
    masses = np.array(sorted(m for w in plan.values() for m in w["masses"]))
    steps = np.diff(masses)
    mid = masses[:-1]

    nodes = {}
    for f in os.listdir(SIGDIR):
        if f.startswith("case_interpolation_M") and f.endswith(".json"):
            m = float(f[len("case_interpolation_M"):-len(".json")])
            nodes[m] = json.load(open(os.path.join(SIGDIR, f)))["sigma"]
    nm = np.array(sorted(nodes))
    sig_at = np.array([nodes[nm[np.argmin(abs(nm - m))]] for m in mid])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(nm, [nodes[m] for m in nm], "-", color="#4c72b0", lw=2,
                 label=r"Signal resolution $\sigma(m)$")
    axes[0].plot(mid, steps, "k.", ms=5, label="Grid spacing")
    axes[0].set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    axes[0].set_ylabel("[GeV]")
    axes[0].set_title("Spacing tracks the resolution", fontsize=11)
    axes[0].legend(fontsize=10)
    axes[0].grid(alpha=0.3)

    axes[1].plot(mid, steps / sig_at, "k.", ms=5)
    axes[1].axhline(1.0, color="r", ls="--", lw=1.5, label=r"$1\,\sigma$")
    axes[1].set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    axes[1].set_ylabel(r"Grid spacing / $\sigma$")
    axes[1].set_ylim(0, 2)
    axes[1].set_title("Spacing in units of the local resolution", fontsize=11)
    axes[1].legend(fontsize=10)
    axes[1].grid(alpha=0.3)

    fig.suptitle("Mass scan grid: %d hypotheses from %.1f to %.1f GeV"
                 % (len(masses), masses[0], masses[-1]), fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(OUT, "scan_grid.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


def gof_figure():
    d = json.load(open(ENS_BKG))
    masses = sorted(d, key=float)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    # pooled chi2/ndof
    allr = np.concatenate([np.array(d[m]["chi2ndof"]) for m in masses])
    axes[0].hist(allr, bins=np.linspace(0, 2, 21), color="#4c72b0",
                 edgecolor="k", alpha=0.85)
    axes[0].axvline(1.0, color="r", ls="--", lw=1.5, label=r"$\chi^2/$ndof $= 1$")
    axes[0].set_xlabel(r"$\chi^2/$ndof (background-only)")
    axes[0].set_ylabel("Fits")
    axes[0].set_title(f"All masses pooled (N = {len(allr)})", fontsize=11)
    axes[0].legend(fontsize=10)

    # per mass
    means = [np.mean(d[m]["chi2ndof"]) for m in masses]
    stds = [np.std(d[m]["chi2ndof"]) for m in masses]
    x = [float(m) for m in masses]
    axes[1].errorbar(x, means, yerr=stds, fmt="ko", ms=6, capsize=4)
    axes[1].axhline(1.0, color="r", ls="--", lw=1.5)
    axes[1].set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    axes[1].set_ylabel(r"$\langle \chi^2/$ndof$\rangle$")
    axes[1].set_ylim(0, 2)
    axes[1].set_title("Mean per mass (bars = RMS over draws)", fontsize=11)
    axes[1].grid(alpha=0.3)

    fig.suptitle("Background-only goodness of fit, 3% background MC ensemble",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUT, "gof_ensemble.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


def pull_figure():
    d = json.load(open(ENS_SB))
    masses = sorted(d, key=float)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    allp = np.concatenate([np.array(d[m]["pull"]) for m in masses])
    axes[0].hist(allp, bins=np.linspace(-3, 3, 25), color="#55a868",
                 edgecolor="k", alpha=0.85)
    axes[0].axvline(0.0, color="r", ls="--", lw=1.5, label="Unbiased")
    axes[0].set_xlabel(r"Yield pull $(N_{\mathrm{fit}} - N_{\mathrm{inj}})/\sigma_N$")
    axes[0].set_ylabel("Fits")
    axes[0].set_title(f"All masses pooled (N = {len(allp)})", fontsize=11)
    axes[0].legend(fontsize=10)

    means = [np.mean(d[m]["pull"]) for m in masses]
    errs = [np.std(d[m]["pull"]) / np.sqrt(len(d[m]["pull"])) for m in masses]
    x = [float(m) for m in masses]
    axes[1].errorbar(x, means, yerr=errs, fmt="ko", ms=6, capsize=4)
    axes[1].axhline(0.0, color="r", ls="--", lw=1.5)
    axes[1].set_xlabel(r"$m_{\mu\mu}$ [GeV]")
    axes[1].set_ylabel("Mean pull")
    axes[1].set_ylim(-1.5, 1.5)
    axes[1].set_title("Mean per mass (bars = error on the mean)", fontsize=11)
    axes[1].grid(alpha=0.3)

    fig.suptitle(r"Signal+background injection: extracted-yield pulls ($\sim 5\sigma$ signal)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(OUT, "sb_pulls.png")
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    interp_params_figure()
    gof_figure()
    # Kept as an internal cross-check only -- superseded in the note by the
    # combine-toy bias test (tools/bias_test_toys.py + tools/plot_bias_toys.py).
    if os.environ.get("MAKE_PULL_FIGURE"):
        pull_figure()
