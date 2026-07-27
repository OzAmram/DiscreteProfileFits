"""
Re-plot a background-only ensemble from its saved combine outputs (does NOT re-fit).

Same idea as replot_ensemble.py but for ensembles with NO signal injection: it
refreshes each fit's postfit plots + fit_results json (so the corrected
makePostfitPlot chi2 is recorded) and re-aggregates the background-only GOF the
same way ensemble_test.py does.

Usage:  python replot_bkg_ensemble.py <ensemble_dir>
"""
import io
import json
import os
import sys
import types
# Archived under tools/; add repo root to the path so `import makePostfitPlot`
# resolves when run as tools/replot_bkg_ensemble.py from the repo root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from contextlib import redirect_stdout, redirect_stderr

import numpy as np

import makePostfitPlot as M

ROOT_DIR = sys.argv[1] if len(sys.argv) > 1 else "ensemble_N7_b0.5"
MASSES = [15, 35, 55, 75]
SEEDS = range(20)


def replot_one(d, mass):
    rj = os.path.join(d, "fit_results_%.1f.json" % mass)
    fd = os.path.join(d, "fitDiagnostics_test_mumu.root")
    ws = os.path.join(d, "datacardInputs_mass_mumu.root")
    if not (os.path.exists(rj) and os.path.exists(fd) and os.path.exists(ws)):
        return False
    r = json.load(open(rj))
    so = r["script_options"]
    nbins_fine = int((r["m_max"] - r["m_min"]) / so["bin_size"])
    common = dict(
        inputWSFile=ws, fitDiagFile=fd, cat="mass_mumu", mass=mass,
        mMin=r["m_min"], mMax=r["m_max"], nBins=nbins_fine, pdfNBins=600,
        sigNorm=so["sig_norm"], poiName="m", outDir=d + "/", ext="",
        lumi=so.get("lumi", ""), sqrts=so.get("sqrts", 13.6),
        drawSignal=True, jsonFile=rj,
    )
    buf = io.StringIO()
    with redirect_stdout(buf), redirect_stderr(buf):
        M.make_postfit_plot(types.SimpleNamespace(**common))
        ob = types.SimpleNamespace(**common)
        ob.bkgOnly = True
        ob.drawSignal = False
        M.make_postfit_plot(ob)
    return True


def main():
    print("Re-plotting %s (using saved combine outputs, no re-fit)" % ROOT_DIR, flush=True)
    ok = fail = 0
    for k in SEEDS:
        for m in MASSES:
            d = os.path.join(ROOT_DIR, "seed%d" % k, "M%d" % m)
            if replot_one(d, float(m)):
                ok += 1
            else:
                fail += 1
    print("  re-plotted %d dirs (%d skipped: missing outputs)" % (ok, fail), flush=True)

    # ---- aggregate background-only GOF (same metrics as ensemble_test.py) ----
    print("\n" + "=" * 82)
    print("%-5s%-5s%-14s%-9s%-8s%-13s%-13s"
          % ("M", "N", "<chi2/ndof>", "median", "max", "frac p<0.05", "frac p<0.01"))
    print("-" * 82)
    summary = {}
    for M_ in MASSES:
        ratios, pvals = [], []
        for k in SEEDS:
            rj = os.path.join(ROOT_DIR, "seed%d" % k, "M%d" % M_,
                              "fit_results_%.1f.json" % M_)
            if not os.path.exists(rj):
                continue
            r = json.load(open(rj))
            c, n = r.get("bkgfit_chi2", -1), r.get("bkgfit_ndof", -1)
            if n and n > 0:
                ratios.append(c / n)
                pvals.append(r.get("bkgfit_prob", -1))
        summary[M_] = dict(chi2ndof=ratios, pval=pvals)
        if not ratios:
            print("%-5d0    (no results)" % M_); continue
        ra, pa = np.array(ratios), np.array(pvals)
        print("%-5d%-5d%-14.2f%-9.2f%-8.2f%-13.2f%-13.2f"
              % (M_, len(ra), ra.mean(), np.median(ra), ra.max(),
                 np.mean(pa < 0.05), np.mean(pa < 0.01)))
    print("=" * 82)
    print("Good fit => <chi2/ndof>~1, frac(p<0.05)~0.05  "
          "(~1 draw in 20 below p=0.05 by chance per mass).")
    out = os.path.join(ROOT_DIR, "ensemble_summary_fixed.json")
    with open(out, "w") as f:
        json.dump({str(m): summary[m] for m in MASSES}, f, indent=2)
    print("Saved " + out)


if __name__ == "__main__":
    main()
