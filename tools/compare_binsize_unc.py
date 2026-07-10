"""
Compare signal-yield pull AND yield uncertainty across bin sizes for the
smooth-toy test. The point: finer bins should shrink sigma(N) while keeping the
pull (bias) small.
"""
import json
import os

import numpy as np

MASSES = [15, 35, 55]
CONFIGS = [0.5, 0.25, 0.1]  # bin fraction


def gather(binf, m):
    root = "smooth_toy_N7_b%g" % binf
    pulls, uncs = [], []
    for k in range(20):
        rj = os.path.join(root, "seed%d" % k, "M%d" % m, "fit_results_%d.0.json" % m)
        if not os.path.exists(rj):
            continue
        r = json.load(open(rj))
        u = r.get("obs_excess_events_unc")
        if not u or u <= 0:
            continue
        pulls.append(r["obs_excess_events"] / u)
        uncs.append(u)
    return np.array(pulls), np.array(uncs)


print("%-6s %-6s %-22s %-16s" % ("mass", "bin", "pull mean(std)  N", "median sigma(N)"))
print("-" * 60)
for m in MASSES:
    base_unc = None
    for binf in CONFIGS:
        p, u = gather(binf, m)
        if not len(p):
            print("M%-4d %-6.2f  (no data)" % (m, binf))
            continue
        med = np.median(u)
        if binf == 0.5:
            base_unc = med
        rel = "" if base_unc is None else "  (%.1f%% of 0.5sig)" % (100 * med / base_unc)
        print("M%-4d %-6.2f  %+.2f(%.2f)  N=%-3d   %8.1f  mean=%.1f%s"
              % (m, binf, p.mean(), p.std(), len(p), med, u.mean(), rel))
    print()
