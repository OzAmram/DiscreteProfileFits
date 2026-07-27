"""
Re-plot an existing S+B ensemble from its saved combine outputs (does NOT re-fit).

The combine fits (yields, pulls, fitDiagnostics) are unchanged; only the postfit
chi2 computation in makePostfitPlot changed (it now uses combine's actual
pdf_index instead of a chi2 scan). So we just re-run make_postfit_plot on each
existing fit dir to refresh sbfit_chi2/bkgfit_chi2 in its fit_results json, then
re-aggregate the same way ensemble_sb.py does.

Usage:  python replot_ensemble.py <ensemble_dir>
"""
import io
import json
import os
import sys
# Archived under tools/; add repo root to the path so `import makePostfitPlot`
# resolves when run as tools/replot_ensemble.py from the repo root.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import types
from contextlib import redirect_stdout, redirect_stderr

import numpy as np

import makePostfitPlot as M

ROOT_DIR = sys.argv[1] if len(sys.argv) > 1 else "ensemble_sb_N7_b0.5"
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

    # ---- aggregate (same metrics as ensemble_sb.py) ----
    print("\n" + "=" * 96)
    print("%-5s%-5s%-11s%-10s%-8s%-16s%-14s%-14s"
          % ("M", "N", "pull mean", "pull std", "<r>",
             "<sb chi2/ndof>", "frac sbP<0.05", "frac bkgP<0.05"))
    print("-" * 96)
    summary = {}
    for M_ in MASSES:
        pulls, rs, sbr, sbp, bkp = [], [], [], [], []
        for k in SEEDS:
            base = os.path.join(ROOT_DIR, "seed%d" % k, "M%d" % M_)
            rj = os.path.join(base, "fit_results_%.1f.json" % M_)
            tj = os.path.join(base, "truth.json")
            if not (os.path.exists(rj) and os.path.exists(tj)):
                continue
            r = json.load(open(rj))
            n_inj = json.load(open(tj))["n_inj"]
            ne, nee = r.get("obs_excess_events", 0), r.get("obs_excess_events_unc", 0)
            if nee and nee > 0:
                pulls.append((ne - n_inj) / nee)
                rs.append(ne / n_inj if n_inj else np.nan)
            c, n = r.get("sbfit_chi2", -1), r.get("sbfit_ndof", -1)
            if n and n > 0:
                sbr.append(c / n); sbp.append(r.get("sbfit_prob", -1))
            bc, bn = r.get("bkgfit_chi2", -1), r.get("bkgfit_ndof", -1)
            if bn and bn > 0:
                bkp.append(r.get("bkgfit_prob", -1))
        summary[M_] = dict(pull=pulls, r=rs, sb_chi2ndof=sbr, sbpval=sbp, bkgpval=bkp)
        if not sbr:
            print("%-5d0    (no results)" % M_); continue
        pa, sa, spa, bpa = map(np.array, (pulls, sbr, sbp, bkp))
        print("%-5d%-5d%-11.2f%-10.2f%-8.2f%-16.2f%-14.2f%-14.2f"
              % (M_, len(sbr), pa.mean(), pa.std(), np.nanmean(rs),
                 sa.mean(), np.mean(spa < 0.05), np.mean(bpa < 0.05)))
    print("=" * 96)
    print("Good => pull mean~0, pull std~1, <sb chi2/ndof>~1, "
          "frac(sbP<0.05)~0.05, frac(bkgP<0.05)~0.05")
    out = os.path.join(ROOT_DIR, "ensemble_sb_summary_fixed.json")
    with open(out, "w") as f:
        json.dump({str(m): summary[m] for m in MASSES}, f, indent=2)
    print("Saved " + out)


if __name__ == "__main__":
    main()
