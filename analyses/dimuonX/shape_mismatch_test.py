#!/usr/bin/env python3
"""
Signal-shape mismatch sensitivity test.

Inject ~5 sigma of a "true" signal (REAL signal MC events) into the 3% anomaly-cut
background MC, then fit with DCB templates from the SAME and from DIFFERENT width
classes.  Quantify how much discovery sensitivity is lost when the template does
not match the true signal shape.

Three width classes (representative loosemass families):
    narrow -> VLL   (VLLVLLToZHTo2MuInv_MVLL-250)
    medium -> TP    (TpTpTo2T2STo2Mu2G_MTp1000)
    wide   -> H2XH3 (H2toH1H3to2Mu_MH2-250)

At each of M = 15, 35, 55, 75 GeV we build one mock per injected class and fit it
with all three class templates -> full 3x3 (injectant x template) matrix per mass.

Injection sizing (per injected class): counting estimate Z ~ S_win / sqrt(B_win)
in a +-2 sigma window of the TRUE shape, solved for N_inject at target Z = 5.
So the DIAGONAL (correct template) is ~5 sigma by construction; off-diagonal
(mismatched template) recovers less -- the shortfall is the sensitivity loss.

Fit window is COMMON across the three templates at a given mass (+-8 * sigma_max,
0.25*sigma_max bins, capped to [10, 80]).  Since the mock and window are identical
across the three template fits, the ONLY thing that varies is the signal shape
JSON -> clean isolation of the template-mismatch effect.

Primary metric: recovered discovery significance Z (combine dLL).  Also recorded:
fitted signal strength r (=1 at the injected truth, since sig_norm = N_inject) and
the expected upper limit.

Run under cmsenv (combine on PATH).  See run_shape_mismatch.sh.
"""
import argparse
import concurrent.futures as cf
import json
import os
import subprocess
import sys

import h5py
import numpy as np

REPO    = os.path.dirname(os.path.abspath(__file__))
FRAMEWORK_DIR = os.path.abspath(os.path.join(REPO, "..", ".."))  # doFit.py etc. live here
BKG     = os.path.join(REPO, "bkg_mc_3pct.h5")
SIGDATA = os.path.join(REPO, "signal_data_loosemass")
FITS    = os.path.join(REPO, "signal_fits_loosemass")
CONFIG  = os.path.join(REPO, "dimuonX_config.json")
OUT     = os.path.join(REPO, "shape_mismatch_test")

FAMILIES = {
    "narrow": "VLLVLLToZHTo2MuInv_MVLL-250",
    "medium": "TpTpTo2T2STo2Mu2G_MTp1000",
    "wide":   "H2toH1H3to2Mu_MH2-250",
}
CLASS_LABEL = {"narrow": "VLL (narrow)", "medium": "TP (medium)", "wide": "H2XH3 (wide)"}
CLASSES = ["narrow", "medium", "wide"]
MASSES  = [15, 35, 55, 75]
TARGET_Z = 5.0
SEED = 42


def interp_params(fam, mass):
    p = os.path.join(FITS, FAMILIES[fam], f"case_interpolation_M{float(mass)}.json")
    with open(p) as f:
        return json.load(f)


def template_json(fam, mass):
    return os.path.join(FITS, FAMILIES[fam], f"case_interpolation_M{float(mass)}.json")


def fit_window(mass):
    """Common fit window at this mass: +-8 sigma_max, 0.25 sigma_max bins, capped [10,80]."""
    smax = max(interp_params(c, mass)["sigma"] for c in CLASSES)
    lo = max(10.0, round(mass - 8.0 * smax, 3))
    hi = min(80.0, round(mass + 8.0 * smax, 3))
    bs = round(0.25 * smax, 4)
    return lo, hi, bs


def counting_N(inj, mass, bkg, sig):
    """Counting-estimate N_inject for a target Z=5 discovery in a +-2 sigma window."""
    sp = interp_params(inj, mass)
    mean, sigma = sp["mean"], sp["sigma"]
    lo, hi = mean - 2 * sigma, mean + 2 * sigma
    B_win = int(np.count_nonzero((bkg >= lo) & (bkg <= hi)))
    f_win = np.count_nonzero((sig >= lo) & (sig <= hi)) / len(sig)
    S_win = TARGET_Z * np.sqrt(B_win)
    N = int(round(S_win / f_win))
    info = dict(mean=float(mean), sigma=float(sigma), B_win=B_win,
                f_win=float(f_win), S_win=float(S_win),
                counting_N=N, counting_Z=float(S_win / np.sqrt(B_win)))
    return N, info


def write_dataset(path, masses):
    with h5py.File(path, "w") as f:
        f.create_dataset("masses", data=masses.astype(np.float32),
                         compression="gzip", compression_opts=4, chunks=True)


def make_toy(bkg, sig, N_inject, rng, bootstrap_bkg=True, poisson_sig=True):
    """A toy pseudo-dataset: (bootstrapped) background + Poisson(N) injected signal."""
    if bootstrap_bkg:
        b = bkg[rng.choice(len(bkg), size=len(bkg), replace=True)]
    else:
        b = bkg
    n = int(rng.poisson(N_inject)) if poisson_sig else int(N_inject)
    n = max(1, n)
    s = sig[rng.choice(len(sig), size=n, replace=(n > len(sig)))]
    d = np.concatenate([b, s]).astype(np.float32)
    rng.shuffle(d)
    return d


def _clean_light(outdir):
    """Drop heavy per-fit artefacts, keep fit_results_*.json + fit_log.txt."""
    for fn in os.listdir(outdir):
        if fn.startswith("fit_results_") or fn == "fit_log.txt":
            continue
        try:
            os.remove(os.path.join(outdir, fn))
        except (IsADirectoryError, OSError):
            pass


def run_one_fit(inj, template, mass, mock, sig_norm, sub="", light=False):
    lo, hi, bs = fit_window(mass)
    name = f"inj-{inj}_tmpl-{template}" + (f"_{sub}" if sub else "")
    outdir = os.path.join(OUT, f"M{mass}", name)
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        sys.executable, os.path.join(FRAMEWORK_DIR, "doFit.py"),
        "-c", CONFIG,
        "-M", str(mass),
        "--m-min", str(lo), "--m-max", str(hi), "--bin-size", str(bs),
        "-i", mock,
        "-s", template_json(template, mass),
        "--sig_norm", str(sig_norm),
        "-l", f"inj{inj}_t{template}",
        "-o", outdir,
    ]
    logp = os.path.join(outdir, "fit_log.txt")
    with open(logp, "w") as log:
        subprocess.run(cmd, cwd=REPO, stdout=log, stderr=subprocess.STDOUT)

    rj = os.path.join(outdir, f"fit_results_{float(mass)}.json")
    rec = dict(inj=inj, template=template, mass=mass, match=(inj == template),
               sub=sub, window=[lo, hi, bs], sig_norm=sig_norm, ok=False)
    if os.path.exists(rj):
        with open(rj) as f:
            r = json.load(f)
        rec.update(ok=True,
                   signif=r.get("signif"),
                   pval=r.get("pval"),
                   r_fit=r.get("obs_excess_events", 0.0) / sig_norm if sig_norm else None,
                   excess_events=r.get("obs_excess_events"),
                   excess_events_unc=r.get("obs_excess_events_unc"),
                   exp_lim_events=r.get("exp_lim_events"),
                   obs_lim_events=r.get("obs_lim_events"))
    if light:
        _clean_light(outdir)
    return rec


def calibrate_cell(inj, mass, calibrate):
    """Size the injection N for one cell so the matched-template Z ~= TARGET_Z."""
    rng = np.random.default_rng(SEED + 1000 * mass + (hash(inj) % 1000))
    with h5py.File(BKG) as f:
        bkg = f["masses"][:]
    with h5py.File(os.path.join(SIGDATA, FAMILIES[inj], f"mass_{mass}.h5")) as f:
        sig = f["masses"][:]
    N, truth = counting_N(inj, mass, bkg, sig)
    truth.update(inj=inj, mass=mass)
    if calibrate:
        toy = make_toy(bkg, sig, N, rng, bootstrap_bkg=False, poisson_sig=False)
        mpath = os.path.join(OUT, "mocks", f"calib_{inj}_M{mass}.h5")
        os.makedirs(os.path.dirname(mpath), exist_ok=True)
        write_dataset(mpath, toy)
        cal = run_one_fit(inj, inj, mass, mpath, float(N), sub="calib", light=True)
        Z0 = cal.get("signif") or 0.0
        truth["calib_N"], truth["calib_Z"] = N, Z0
        if Z0 > 0.5:
            N = int(round(N * TARGET_Z / Z0))
        N = max(10, N)
    truth["N_inject"] = N
    print(f"[calib] inj={inj:6s} M{mass}: N={N} "
          f"(counting {truth['counting_N']}, calibZ={truth.get('calib_Z')})", flush=True)
    return truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="Build+fit a single combo (medium/medium M35) to validate.")
    ap.add_argument("--no-calib", action="store_true",
                    help="Skip Z-calibration; size injection by the counting estimate only.")
    ap.add_argument("--ntoys", type=int, default=1,
                    help="Toys per (injectant,template,mass) cell (ensemble error bars).")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--masses", type=int, nargs="+", default=MASSES)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(os.path.join(OUT, "mocks"), exist_ok=True)

    masses = args.masses
    injectants = CLASSES
    templates = CLASSES
    if args.smoke:
        masses, injectants, templates, args.ntoys = [35], ["medium"], ["medium"], 1

    calibrate = (not args.no_calib) and (not args.smoke)
    cells = [(inj, m) for m in masses for inj in injectants]

    # ---- phase 1: calibrate injection size per cell (parallel) ----
    print(f"Phase 1: calibrating {len(cells)} cells "
          f"(calibrate={calibrate}) ...\n", flush=True)
    truths = {}
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(calibrate_cell, inj, m, calibrate): (inj, m) for inj, m in cells}
        for fut in cf.as_completed(futs):
            t = fut.result()
            truths[f"{t['inj']}_M{t['mass']}"] = t

    # ---- phase 2: build toy datasets per cell ----
    toy_paths = {}   # (inj, mass, k) -> path
    for inj, m in cells:
        rng = np.random.default_rng(SEED + 7 + 1000 * m + (hash(inj) % 1000))
        N = truths[f"{inj}_M{m}"]["N_inject"]
        with h5py.File(BKG) as f:
            bkg = f["masses"][:]
        with h5py.File(os.path.join(SIGDATA, FAMILIES[inj], f"mass_{m}.h5")) as f:
            sig = f["masses"][:]
        for k in range(args.ntoys):
            # a single toy (k=0, ntoys=1) uses the fixed dataset, no resampling,
            # so it reproduces the single-realization run.
            bs_bkg = args.ntoys > 1
            po_sig = args.ntoys > 1
            toy = make_toy(bkg, sig, N, rng, bootstrap_bkg=bs_bkg, poisson_sig=po_sig)
            p = os.path.join(OUT, "mocks", f"toy_{inj}_M{m}_k{k}.h5")
            write_dataset(p, toy)
            toy_paths[(inj, m, k)] = p

    # ---- phase 3: fit every (injectant, template, mass, toy) ----
    jobs = []
    for inj, m in cells:
        N = truths[f"{inj}_M{m}"]["N_inject"]
        for tmpl in templates:
            for k in range(args.ntoys):
                jobs.append((inj, tmpl, m, k, toy_paths[(inj, m, k)], float(N)))
    print(f"\nPhase 3: {len(jobs)} fits ({args.ntoys} toys/cell) "
          f"with {args.workers} workers ...\n", flush=True)

    light = args.ntoys > 1
    results = []
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {}
        for inj, tmpl, m, k, mock, sig_norm in jobs:
            fut = ex.submit(run_one_fit, inj, tmpl, m, mock, sig_norm,
                            sub=f"k{k}" if light else "", light=light)
            futs[fut] = (inj, tmpl, m, k)
        done = 0
        for fut in cf.as_completed(futs):
            inj, tmpl, m, k = futs[fut]
            rec = fut.result()
            rec["toy"] = k
            results.append(rec)
            done += 1
            z = rec.get("signif")
            print(f"[fit {done:3d}/{len(jobs)}] M{m} inj={inj:6s} tmpl={tmpl:6s} k{k} "
                  f"{'OK' if rec['ok'] else 'FAIL'}  "
                  f"Z={z if z is None else round(z, 2)}  r={rec.get('r_fit')}", flush=True)

    out = dict(target_Z=TARGET_Z, seed=SEED, families=FAMILIES, ntoys=args.ntoys,
               calibrated=calibrate, truths=truths, results=results)
    outjson = os.path.join(OUT, "results_smoke.json" if args.smoke
                           else ("results.json" if args.ntoys == 1 else "results_ensemble.json"))
    with open(outjson, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved {outjson}")


if __name__ == "__main__":
    main()
