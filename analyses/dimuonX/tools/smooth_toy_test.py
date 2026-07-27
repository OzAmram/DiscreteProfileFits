"""
Decisive test of METHOD bias, free of parent-MC fluctuations.

The subsample ensemble draws every toy from the SAME parent MC, so any local
statistical fluctuation in the parent (e.g. the ~1.5-sigma dip at 33.5-35 GeV)
is inherited by all draws and shows up as a consistent "spurious signal". That is
an artifact of the finite MC, not of the fit method.

Here we instead generate background toys from a SMOOTH truth: fit the full MC over
a WIDE range (so a local wiggle is averaged out), then accept-reject sample toys
at the realistic 3% in-window statistics, inject nothing, fit with the same S+B
model, and measure the pull = N_fit/sigma. If the method is unbiased the pull
distribution is centred at 0 -- proving the subsample bias was the MC, not us.
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

# doFit.py lives in the general framework, three levels up from tools/.
FRAMEWORK_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))

SCRATCH = "/tmp/claude-52316/-uscms-data-d3-oamram-CMSSW-DiMuonX-src-DimuonX-dimuonX-fits/bd944359-ce0c-4c83-94a0-76166418edc4/scratchpad/smooth"
FULL = "bkg_mc_masses.h5"
MASSES = [int(x) for x in os.environ.get("SM_MASSES", "15 35 55").split()]
K = int(os.environ.get("SM_K", "20"))
NSIG = float(os.environ.get("SM_NSIGMA", "7"))
BINF = float(os.environ.get("SM_BINFRAC", "0.5"))
FRAC = 0.03
WIDE = float(os.environ.get("SM_WIDE", "18"))   # +-GeV wide range for smooth fit
DEG = int(os.environ.get("SM_DEG", "4"))
OUT = os.environ.get("SM_OUT", "smooth_toy_N%g_b%g" % (NSIG, BINF))
NWORK = int(os.environ.get("SM_WORKERS", "3"))


def win(mass):
    s = json.load(open("signal_fits/2B_loosemass/sig_fit_%d.json" % mass))["sigma"]
    lo = max(10.0, mass - NSIG * s)
    hi = min(80.0, mass + NSIG * s)
    return round(lo, 4), round(hi, 4), round(BINF * s, 4)


def smooth_shape(allm, mass, lo, hi):
    """Fit log(counts) vs m over a wide range; return a callable density on [lo,hi]
    and the expected 3% in-window yield."""
    wlo, whi = max(10.0, mass - WIDE), min(80.0, mass + WIDE)
    sel = allm[(allm >= wlo) & (allm <= whi)]
    nb = 120
    h, edges = np.histogram(sel, bins=nb, range=(wlo, whi))
    ctr = 0.5 * (edges[:-1] + edges[1:])
    m = h > 0
    # scale mass to [-1,1] for conditioning
    xs = (ctr - 0.5 * (wlo + whi)) / (0.5 * (whi - wlo))
    coef = np.polyfit(xs[m], np.log(h[m]), DEG)
    binw = edges[1] - edges[0]

    def dens(marr):
        x = (marr - 0.5 * (wlo + whi)) / (0.5 * (whi - wlo))
        return np.exp(np.polyval(coef, x)) / binw   # per-GeV density (full stats)

    # expected 3% yield in [lo,hi]
    grid = np.linspace(lo, hi, 2000)
    n_full = np.trapz(dens(grid), grid)
    return dens, FRAC * n_full


def sample(dens, n, lo, hi, rng):
    grid = np.linspace(lo, hi, 4000)
    ymax = dens(grid).max() * 1.05
    out = []
    while len(out) < n:
        xr = rng.uniform(lo, hi, 3 * n)
        yr = rng.uniform(0, ymax, 3 * n)
        out.extend(xr[yr < dens(xr)].tolist())
    return np.array(out[:n], dtype=np.float32)


def run_one(job):
    k, mass, lo, hi, bs, nexp = job
    rng = np.random.default_rng(70000 + k * 100 + mass)
    n = rng.poisson(nexp)
    toy = sample(SHAPES[mass][0], n, lo, hi, rng)
    d = os.path.join(SCRATCH, "M%d_seed%d" % (mass, k))
    os.makedirs(d, exist_ok=True)
    tp = os.path.join(d, "toy.h5")
    with h5py.File(tp, "w") as f:
        f.create_dataset("masses", data=toy, compression="gzip", chunks=True)
    outdir = os.path.join(OUT, "seed%d" % k, "M%d" % mass)
    os.makedirs(outdir, exist_ok=True)
    cmd = [sys.executable, os.path.join(FRAMEWORK_DIR, "doFit.py"), "-c", "dimuonX_config.json", "-M", str(mass),
           "--m-min", str(lo), "--m-max", str(hi), "--bin-size", str(bs),
           "--sig_norm", "100", "-i", tp,
           "-s", "signal_fits/2B_loosemass/sig_fit_%d.json" % mass, "-o", outdir]
    with open(os.path.join(outdir, "log.txt"), "w") as f:
        subprocess.run(cmd, stdout=f, stderr=f)
    return job


SHAPES = {}


def main():
    os.makedirs(SCRATCH, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    allm = h5py.File(FULL)["masses"][:]
    jobs = []
    for mass in MASSES:
        lo, hi, bs = win(mass)
        dens, nexp = smooth_shape(allm, mass, lo, hi)
        SHAPES[mass] = (dens, nexp)
        print("M%d: window [%.2f,%.2f] bin=%.3f  expected 3%% yield=%.0f"
              % (mass, lo, hi, bs, nexp), flush=True)
        for k in range(K):
            jobs.append((k, mass, lo, hi, bs, nexp))
    print("Running %d smooth-toy fits (%d concurrent)..." % (len(jobs), NWORK), flush=True)
    with ThreadPoolExecutor(max_workers=NWORK) as ex:
        for i, _ in enumerate(ex.map(run_one, jobs), 1):
            if i % 20 == 0:
                print("  %d/%d" % (i, len(jobs)), flush=True)

    print("\n=== SMOOTH-TOY spurious pull (method bias, no MC fluctuation) ===")
    for mass in MASSES:
        ps = []
        for k in range(K):
            rj = os.path.join(OUT, "seed%d" % k, "M%d" % mass, "fit_results_%d.0.json" % mass)
            if not os.path.exists(rj):
                continue
            r = json.load(open(rj))
            u = r.get("obs_excess_events_unc")
            if u and u > 0:
                ps.append(r["obs_excess_events"] / u)
        a = np.array(ps)
        if len(a):
            print("  M%d: N=%d  pull mean=%+.2f  std=%.2f" % (mass, len(a), a.mean(), a.std()))


if __name__ == "__main__":
    main()
