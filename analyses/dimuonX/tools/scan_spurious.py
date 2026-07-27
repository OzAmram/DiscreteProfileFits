"""
Scan window size and bin size to minimise the spurious-signal bias.

Runs the background-only ensemble (ensemble_test.py) for several (Nsigma, binfrac)
configs, then reports, per mass, the mean/median signal PULL = N_fit/sigma (true
signal = 0). A config with all per-mass means ~0 has negligible spurious signal.

Configs are run sequentially (each already uses 3 workers). Subsamples are cached
in SCRATCH and shared across configs (seeded, reproducible).
"""
import json
import os
import subprocess
import sys

import numpy as np

MASSES = [15, 35, 55, 75]
SEEDS = range(20)
# (Nsigma window, bin fraction)
CONFIGS = [
    (4, 0.5),
    (5, 0.5),
    (10, 0.5),
    (5, 0.25),
    (10, 0.25),
]
# already on disk: (7,0.5)->ensemble_N7_b0.5 , (7,0.25)->ensemble_N7_b0.25


def outdir(nsig, binf):
    return "ensemble_N%g_b%g" % (nsig, binf)


def run_config(nsig, binf):
    d = outdir(nsig, binf)
    env = dict(os.environ)
    env["EN_NSIGMA"] = str(nsig)
    env["EN_BINFRAC"] = str(binf)
    env["EN_MASSES"] = " ".join(str(m) for m in MASSES)
    env["EN_K"] = "20"
    env["EN_WORKERS"] = "3"
    env["EN_OUT"] = d
    print("\n########## running %s ##########" % d, flush=True)
    # ensemble_test.py lives beside this script (tools/); reference it by that path.
    _here = os.path.dirname(os.path.abspath(__file__))
    subprocess.run([sys.executable, os.path.join(_here, "ensemble_test.py")], env=env)
    return d


def pulls(d):
    out = {}
    for m in MASSES:
        ps = []
        for k in SEEDS:
            rj = os.path.join(d, "seed%d" % k, "M%d" % m, "fit_results_%.1f.json" % m)
            if not os.path.exists(rj):
                continue
            r = json.load(open(rj))
            u = r.get("obs_excess_events_unc")
            if not u or u <= 0:
                continue
            ps.append(r["obs_excess_events"] / u)
        out[m] = np.array(ps)
    return out


def report(configs):
    print("\n" + "=" * 78)
    print("SPURIOUS-SIGNAL PULL (N_fit/sigma, true=0) by config")
    print("=" * 78)
    print("%-14s" % "config" + "".join("%-15s" % ("M%d" % m) for m in MASSES) + "|mean|")
    print("-" * 78)
    for d in configs:
        pl = pulls(d)
        cells = []
        absmeans = []
        for m in MASSES:
            a = pl[m]
            if len(a):
                cells.append("%+.2f(%.2f)" % (a.mean(), a.std()))
                absmeans.append(abs(a.mean()))
            else:
                cells.append("--")
        worst = max(absmeans) if absmeans else float("nan")
        print("%-14s" % d + "".join("%-15s" % c for c in cells) + "  %.2f" % worst)
    print("=" * 78)
    print("cell = mean(std) of pull;  last col = worst |mean| across masses (smaller=better)")


def main():
    dirs = []
    for nsig, binf in CONFIGS:
        dirs.append(run_config(nsig, binf))
    # include the two existing configs in the report
    alldirs = ["ensemble_N7_b0.5", "ensemble_N7_b0.25"] + dirs
    alldirs = [d for d in alldirs if os.path.isdir(d)]
    report(alldirs)


if __name__ == "__main__":
    main()
