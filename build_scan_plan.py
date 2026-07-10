"""
Build the data-fit scan plan:
  - global mass grid at ~1 sigma fidelity (step = local signal sigma),
    snapped to the 0.2 GeV interpolation-template nodes
  - assign each mass to the window whose 5 GeV center is closest
    (each window owns its central +/-1.25 GeV -> contiguous, non-overlapping)
Prints a per-window summary and writes scan_plan.json.
"""
import glob
import json
import os

import numpy as np

BASE = "/eos/uscms/store/user/tvami/CATHODE/Flow_Learning/Results_vr"
SERIES = "nsplit_20260708_104130"
SIGDIR = "signal_fits/2B"

# window edges (GeV*10) -> center
WIN_TAGS = ["125_175", "150_200", "175_225", "200_250", "225_275", "250_300",
            "275_325", "300_350", "325_375", "350_400", "375_425", "400_450",
            "425_475", "450_500"]


def win_center(tag):
    lo, hi = [int(x) / 10.0 for x in tag.split("_")]
    return 0.5 * (lo + hi)


def win_file(tag):
    return os.path.join(BASE, "%s_%s" % (SERIES, tag),
                        "stage3_classifier_masked", "global_anomalous_masses.h5")


# sigma(m) from the 2B interpolation nodes (0.2 GeV grid)
NODES = {}
for f in glob.glob(os.path.join(SIGDIR, "case_interpolation_M*.json")):
    m = float(f.split("_M")[-1][:-5])
    NODES[m] = json.load(open(f))["sigma"]
node_masses = np.array(sorted(NODES))


def sigma(m):
    return NODES[node_masses[np.argmin(np.abs(node_masses - m))]]


def snap(m):
    return float(node_masses[np.argmin(np.abs(node_masses - m))])


centers = [win_center(t) for t in WIN_TAGS]
OWN = 1.25  # each window owns +/-1.25 GeV around its center
m_lo = min(centers) - OWN
m_hi = max(centers) + OWN

# global adaptive grid, snapped + de-duped
grid = []
m = m_lo
while m <= m_hi + 1e-6:
    grid.append(snap(m))
    m += sigma(m)
grid = sorted(set(grid))

# assign each mass to nearest center (ties -> lower center)
plan = {t: [] for t in WIN_TAGS}
for mass in grid:
    ci = int(np.argmin([abs(mass - c) + 1e-6 * c for c in centers]))
    plan[WIN_TAGS[ci]].append(mass)

# ---- report ----
print("=" * 70)
print("DATA-FIT SCAN PLAN  (2B signal, closest-center assignment)")
print("=" * 70)
print("%-9s %-8s %-9s %-7s  %s" % ("window", "center", "owns", "Nfit", "masses (GeV)"))
print("-" * 70)
total = 0
out = {}
for t in WIN_TAGS:
    c = win_center(t)
    ms = plan[t]
    total += len(ms)
    out[t] = {"center": c, "file": win_file(t), "masses": ms}
    owns = "[%.2f,%.2f]" % (c - OWN, c + OWN)
    mstr = " ".join("%.1f" % x for x in ms)
    print("%-9s %-8.1f %-9s %-7d  %s" % (t, c, owns, len(ms), mstr))
print("-" * 70)
print("total mass points: %d   range [%.1f, %.1f] GeV" % (total, grid[0], grid[-1]))
json.dump(out, open("scan_plan.json", "w"), indent=2)
print("wrote scan_plan.json")
