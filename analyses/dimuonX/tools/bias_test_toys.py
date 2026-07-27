"""
Bias test with combine-generated toys.

The earlier bias/pull tests drew subsamples of the background MC, so every toy
inherited the parent sample's own statistical fluctuations (notably a ~1.5 sigma
deficit at 34-35 GeV), which showed up as a fake negative "bias". This test
removes that entirely: the truth is a SMOOTH FITTED SHAPE, not the MC itself.

Procedure, per mass hypothesis:
  1. Fit the FULL background MC (done separately, see --truth-dir) to pin down a
     smooth background shape at high statistics. The resulting combine workspace
     holds all four functional-form families in a RooMultiPdf, each already at
     its F-test-selected order with its parameters at the full-MC best fit.
  2. For each generating family in turn, generate toys with `combine -M
     GenerateOnly`, freezing pdf_index to that family so the toy truth is exactly
     that shape, and setting the background yield to 3% of the full-MC yield in
     the window so the toys carry realistic search statistics.
  3. Fit every toy back with `combine -M FitDiagnostics` letting pdf_index FLOAT,
     i.e. through the same discrete-profiling envelope the real fit uses.
  4. Record the pull (r_fit - r_gen)/sigma_r.

Generating from each family in turn and fitting with the full envelope is the
standard discrete-profiling bias test: it measures the bias induced by the
function choice itself, and subsumes the parent-MC question.

A signal is injected at r_gen = 0 (does the fit invent a signal from nothing?)
and at ~2 and ~5 sigma equivalents (does it recover a real signal unbiased?).

Run from the repo root, AFTER the truth workspaces exist:
    python3 tools/bias_test_toys.py --masses 15 35 55 75 --ntoys 100
"""
import argparse
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import h5py

# Import ROOT up front: its lazy attribute loading does not work from inside the
# worker threads, so touching ROOT.TFile there raises AttributeError.
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

FULL_MC = "bkg_mc_masses.h5"
SIGDIR = "signal_fits/2B_loosemass"
CONFIG = "dimuonX_config.json"
FRAC = 0.03            # anomaly-selection efficiency the search operates at
SIG_NORM = 100.0       # r is in units of this many signal events (config sig_norm)
# combine renames the background yield when building the workspace; this is the
# name that actually exists post-text2workspace (multi_pdf_norm does NOT).
NORM_PAR = "shapeBkg_background_mass_mumu__norm"


def sig_params(mass):
    with open(os.path.join(SIGDIR, "case_interpolation_M%s.json" % float(mass))) as f:
        return json.load(f)


def dcb(x, p):
    t = (x - p["mean"]) / p["sigma"]
    a1, n1, a2, n2 = p["alpha"], p["sign"], p["alpha2"], p["sign2"]
    A1 = (n1 / a1) ** n1 * np.exp(-0.5 * a1 ** 2); B1 = n1 / a1 - a1
    A2 = (n2 / a2) ** n2 * np.exp(-0.5 * a2 ** 2); B2 = n2 / a2 - a2
    return np.where(t < -a1, A1 * (B1 - t) ** (-n1),
                    np.where(t > a2, A2 * (B2 + t) ** (-n2), np.exp(-0.5 * t ** 2)))


def gen_strengths(mass, m_min, m_max, zs):
    """r_gen values giving roughly the requested significances at 3% stats."""
    p = sig_params(mass)
    mean, sigma = p["mean"], p["sigma"]
    lo2, hi2 = mean - 2 * sigma, mean + 2 * sigma

    with h5py.File(FULL_MC) as f:
        m = f["masses"][:]
    # background in +-2sigma, scaled to the 3% the search actually sees
    B_win = FRAC * np.count_nonzero((m >= lo2) & (m <= hi2))

    # fraction of the signal template inside +-2sigma (over the fit window)
    xs = np.linspace(m_min, m_max, 4000)
    ys = dcb(xs, p)
    tot = np.trapz(ys, xs)
    core = np.trapz(ys[(xs >= lo2) & (xs <= hi2)], xs[(xs >= lo2) & (xs <= hi2)])
    f_win = core / tot

    out = {}
    for z in zs:
        n_inj = 0.0 if z == 0 else z * np.sqrt(B_win) / f_win
        out[z] = round(n_inj / SIG_NORM, 5)
    return out, B_win


def truth_info(truth_dir, mass):
    """Window + F-test orders + pdf_index order, from the full-MC truth fit."""
    rj = os.path.join(truth_dir, "fit_results_%s.json" % float(mass))
    r = json.load(open(rj))
    forms = r["func_forms"]
    # pdf_index follows the order the families were added, i.e. config key order
    fams = list(forms.keys())
    orders = {fam: forms[fam][r["final_func_forms"][fam]] for fam in fams}
    return r["m_min"], r["m_max"], fams, orders


def pull_of(r, r_gen, r_err, r_lo, r_hi):
    """(r - r_gen) divided by the uncertainty pointing back towards r_gen.

    --robustFit + cminDefaultMinimizerStrategy 0 routinely leaves the symmetric
    rErr at 0 while still filling the MINOS errors (doFit.py works around the
    same thing), so the asymmetric errors are the primary handle here. Using the
    directional error is also the right choice for a pull: if the fit landed
    above the truth, the relevant uncertainty is the downward one.
    """
    lo, hi, sym = abs(r_lo), abs(r_hi), abs(r_err)
    err = (lo if r >= r_gen else hi) or sym or (0.5 * (lo + hi))
    if not err or err <= 0:
        return None
    return (r - r_gen) / err


def run(cmd, log, cwd):
    with open(log, "w") as f:
        return subprocess.run(cmd, shell=True, stdout=f, stderr=f, cwd=cwd).returncode


def one_config(job):
    """Generate toys from family `k`, fit them back with the floating envelope.

    Runs subprocesses ONLY -- no ROOT here. ROOT is not thread safe, and reading
    the outputs from inside the worker threads deadlocks the pool; the results
    are collected serially by collect_config() once the pool has drained.
    """
    mass, k, fam, z, r_gen, n3pct, truth_dir, ntoys, seed = job
    tag = "M%s_gen%s_z%g" % (mass, fam, z)
    wdir = os.path.abspath(truth_dir)

    gen_n = "_gen_%s_z%g" % (fam, z)
    fit_n = "_fit_%s_z%g" % (fam, z)

    # 1) generate: pdf_index frozen to this family, yield fixed to 3% stats
    gen = (
        "combine -M GenerateOnly workspace_test_mumu.root -m {m} -t {t} --saveToys "
        "--expectSignal {r} --setParameters pdf_index={k},{norm}={n} "
        "--freezeParameters pdf_index,{norm} -n {n_} -s {s} --rMin -5 --rMax 10"
    ).format(m=float(mass), t=ntoys, r=r_gen, k=k, norm=NORM_PAR, n=n3pct,
             n_=gen_n, s=seed)
    rc = run(gen, os.path.join(wdir, "log%s.txt" % gen_n), wdir)
    if rc != 0:
        return tag, None, "GenerateOnly rc=%d" % rc

    toyfile = "higgsCombine%s.GenerateOnly.mH%g.%s.root" % (gen_n, float(mass), seed)
    if not os.path.exists(os.path.join(wdir, toyfile)):
        return tag, None, "toy file missing: %s" % toyfile

    # 2) fit back: pdf_index FLOATS -> full discrete profiling, as in the real fit.
    # Flags match the ones doFit.py uses for the real fit.
    fit = (
        "combine -M FitDiagnostics workspace_test_mumu.root -m {m} -t {t} "
        "--toysFile {tf} --rMin -5 --rMax 10 --robustFit 1 "
        "--cminDefaultMinimizerStrategy 0 -n {n_} -s {s}"
    ).format(m=float(mass), t=ntoys, tf=toyfile, n_=fit_n, s=seed)
    rc = run(fit, os.path.join(wdir, "log%s.txt" % fit_n), wdir)
    if rc != 0:
        return tag, None, "FitDiagnostics rc=%d" % rc

    res = os.path.join(wdir, "fitDiagnostics%s.root" % fit_n)
    if not os.path.exists(res):
        return tag, None, "fitDiagnostics output missing"
    return tag, res, None


def collect_config(job, res_path):
    """Read one fitDiagnostics output. Serial: ROOT is not thread safe."""
    mass, k, fam, z, r_gen, n3pct, truth_dir, ntoys, seed = job
    f = ROOT.TFile(res_path)
    t = f.Get("tree_fit_sb")
    if not t:
        f.Close()
        return None, "tree_fit_sb missing"
    pulls, rs, errs = [], [], []
    for e in t:
        if getattr(e, "fit_status", 0) < 0:      # fit failed outright
            continue
        p = pull_of(e.r, r_gen, e.rErr, e.rLoErr, e.rHiErr)
        if p is None:
            continue
        rs.append(e.r)
        errs.append(abs(e.r - r_gen) / abs(p) if p else 0.0)
        pulls.append(p)
    f.Close()
    return {"mass": mass, "gen_fam": fam, "z": z, "r_gen": r_gen,
            "n_toys": len(pulls), "pull": pulls, "r": rs, "rErr": errs}, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masses", nargs="+", type=int, default=[15, 35, 55, 75])
    ap.add_argument("--zs", nargs="+", type=float, default=[0.0, 2.0, 5.0])
    ap.add_argument("--ntoys", type=int, default=100)
    ap.add_argument("--truth-root", default="bias_toys")
    ap.add_argument("--out", default="bias_toys/bias_toys_summary.json")
    ap.add_argument("-j", "--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    jobs = []
    for mass in args.masses:
        tdir = os.path.join(args.truth_root, "truth_M%d" % mass)
        m_min, m_max, fams, orders = truth_info(tdir, mass)
        rgens, B_win = gen_strengths(mass, m_min, m_max, args.zs)
        with h5py.File(FULL_MC) as f:
            m = f["masses"][:]
        n3pct = round(FRAC * np.count_nonzero((m >= m_min) & (m <= m_max)), 2)
        print("M%-3d window=[%.2f,%.2f] orders=%s  N_bkg(3%%)=%.0f  r_gen=%s"
              % (mass, m_min, m_max, orders, n3pct,
                 {z: rgens[z] for z in args.zs}), flush=True)
        for k, fam in enumerate(fams):
            for z in args.zs:
                # Distinct seed per config: a shared seed makes the toy sets
                # correlated across generating families, so their pulls would
                # not be independent measurements.
                seed = args.seed + 1000 * mass + 10 * k + int(z)
                jobs.append((mass, k, fam, z, rgens[z], n3pct, tdir,
                             args.ntoys, seed))

    print("\nRunning %d configs x %d toys (%d workers)...\n"
          % (len(jobs), args.ntoys, args.workers), flush=True)
    results, fails = [], []

    # Phase 1: all the combine work, in parallel (subprocesses only).
    done = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        for i, (job, (tag, res, err)) in enumerate(
                zip(jobs, ex.map(one_config, jobs)), 1):
            if err:
                fails.append((tag, err))
                print("  [%d/%d] %-22s FAILED: %s" % (i, len(jobs), tag, err), flush=True)
            else:
                done.append((job, tag, res))
                print("  [%d/%d] %-22s toys fit" % (i, len(jobs), tag), flush=True)

    # Phase 2: read the outputs serially.
    print("\nCollecting...\n", flush=True)
    for job, tag, res_path in done:
        r, err = collect_config(job, res_path)
        if err:
            fails.append((tag, err))
            print("  %-22s FAILED: %s" % (tag, err), flush=True)
            continue
        results.append(r)
        p = np.array(r["pull"])
        print("  %-22s n=%3d  median pull = %+.3f" % (
            tag, len(p), np.median(p) if len(p) else float("nan")), flush=True)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump({"results": results, "failures": fails}, f, indent=2)
    print("\nSaved %s  (%d ok, %d failed)" % (args.out, len(results), len(fails)))


if __name__ == "__main__":
    main()
