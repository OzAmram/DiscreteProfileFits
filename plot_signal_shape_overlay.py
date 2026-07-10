"""
Overlay the fitted DCB signal shapes of all signal models at fixed dimuon mass,
to check shape consistency across production modes. Two panels: M=20 and M=40 GeV.
Each curve is the DCB from that model's sig_fit_<M>.json, normalized to unit area.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FITROOT = "signal_fits_uaf"
MASSES = [20, 40]

# short, readable labels; grouped by the resolution class the fits reveal
LABEL = {
    "VLLVLLToZHTo2MuInv_MVLL-250": "VLL→ZH (M=250)",
    "VLLVLLToZHTo2MuInv_MVLL-500": "VLL→ZH (M=500)",
    "TpTpTo2T2STo2Mu2B_MTp1000":   "T'→2T2S 2B",
    "TpTpTo2T2STo2Mu2G_MTp1000":   "T'→2T2S 2G",
    "TpTpTo2T2STo2MuInv_MTp1000":  "T'→2T2S Inv",
    "TTH2to2Mu":                   "TTH2→2μ",
    "H2toH1H3to2Mu_MH2-250":       "H2→H1H3 (M=250)",
    "H2toH1H3to2Mu_MH2-500":       "H2→H1H3 (M=500)",
    "H2toH1toInvH3to2Mu_MH2-250":  "H2→H1(inv)H3 (M=250)",
    "H2toH1toInvH3to2Mu_MH2-500":  "H2→H1(inv)H3 (M=500)",
}
# color families: blues=VLL(narrow), greens=TpTp(medium), reds/oranges=H2/TTH(wide)
COLOR = {
    "VLLVLLToZHTo2MuInv_MVLL-250": "#1f4e9c",
    "VLLVLLToZHTo2MuInv_MVLL-500": "#5b8bd6",
    "TpTpTo2T2STo2Mu2B_MTp1000":   "#0f7d3b",
    "TpTpTo2T2STo2Mu2G_MTp1000":   "#3faa62",
    "TpTpTo2T2STo2MuInv_MTp1000":  "#79c98f",
    "TTH2to2Mu":                   "#8c3b00",
    "H2toH1H3to2Mu_MH2-250":       "#c0392b",
    "H2toH1H3to2Mu_MH2-500":       "#e0685c",
    "H2toH1toInvH3to2Mu_MH2-250":  "#d97a00",
    "H2toH1toInvH3to2Mu_MH2-500":  "#f0a640",
}
ORDER = list(LABEL.keys())


def dcb(x, mean, sigma, a1, n1, a2, n2):
    """RooFit DoubleCB (Gaussian core, power-law tails at -a1 and +a2)."""
    t = (x - mean) / sigma
    out = np.empty_like(t)
    core = (t > -a1) & (t < a2)
    out[core] = np.exp(-0.5 * t[core] ** 2)
    L = t <= -a1
    A1 = (n1 / abs(a1)) ** n1 * np.exp(-0.5 * a1 ** 2)
    B1 = n1 / abs(a1) - abs(a1)
    out[L] = A1 * (B1 - t[L]) ** (-n1)
    R = t >= a2
    A2 = (n2 / abs(a2)) ** n2 * np.exp(-0.5 * a2 ** 2)
    B2 = n2 / abs(a2) - abs(a2)
    out[R] = A2 * (B2 + t[R]) ** (-n2)
    return out


fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
for ax, M in zip(axes, MASSES):
    half = 2.5 if M == 20 else 3.5
    x = np.linspace(M - half, M + half, 2000)
    for g in ORDER:
        jf = os.path.join(FITROOT, g, "sig_fit_%d.json" % M)
        if not os.path.exists(jf):
            continue
        r = json.load(open(jf))
        y = dcb(x, r["mean"], r["sigma"], r["alpha"], r["sign"], r["alpha2"], r["sign2"])
        y /= np.trapz(y, x)  # unit area
        ax.plot(x, y, color=COLOR[g], lw=1.8,
                label="%-22s σ=%.2f" % (LABEL[g], r["sigma"]))
    ax.set_title("Signal DCB shapes at m$_{\\mu\\mu}$ = %d GeV" % M, fontsize=13)
    ax.set_xlabel("m$_{\\mu\\mu}$ [GeV]")
    ax.set_ylabel("normalized density [1/GeV]")
    ax.set_xlim(M - half, M + half)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8.5, frameon=False, loc="upper right",
              prop={"family": "monospace", "size": 8})

fig.suptitle("Signal-shape consistency across models  "
             "(unit-area DCB; widths cluster: VLL < T' < H2/TTH)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig("signal_shape_overlay.png", dpi=140)
print("Saved signal_shape_overlay.png")
