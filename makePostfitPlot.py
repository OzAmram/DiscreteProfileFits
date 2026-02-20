#!/usr/bin/env python3
"""
makePostfitPlot.py
==================
Plot the postfit signal + background from CMS combine FitDiagnostics,
for a workspace built with a RooMultiPdf background model.

Follows the style of makeMultipdfPlot.py (elfontan/flashggFinalFit).

Usage (run inside CMSSW environment):
    python makePostfitPlot.py [options]

Default inputs (relative to output dir):
    datacardInputs_mass_mumu.root  -- workspace with multi_pdf, data_obs, signal
    fitDiagnostics_test_mumu.root  -- postfit parameter values from combine
"""

import os
import sys
import json
import math
from optparse import OptionParser
from collections import OrderedDict as od

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


# ---------------------------------------------------------------------------
# Command-line options
# ---------------------------------------------------------------------------

def get_options():
    parser = OptionParser()
    parser.add_option(
        "--inputWSFile", dest="inputWSFile",
        default="datacardInputs_mass_mumu.root",
        help="Input workspace ROOT file (produced by DataCardMaker)")
    parser.add_option(
        "--fitDiagFile", dest="fitDiagFile",
        default="fitDiagnostics_test_mumu.root",
        help="FitDiagnostics ROOT file from combine -M FitDiagnostics")
    parser.add_option(
        "--cat", dest="cat", default="mass_mumu",
        help="Channel / category name (used in PDF and data names)")
    parser.add_option(
        "--mass", dest="mass", default=15.0, type="float",
        help="Signal mass hypothesis in GeV")
    parser.add_option(
        "--mMin", dest="mMin", default=12.0, type="float",
        help="Lower edge of physical mass range [GeV]")
    parser.add_option(
        "--mMax", dest="mMax", default=18.0, type="float",
        help="Upper edge of physical mass range [GeV]")
    parser.add_option(
        "--nBins", dest="nBins", default=60, type="int",
        help="Number of display bins across the mass range")
    parser.add_option(
        "--pdfNBins", dest="pdfNBins", default=600, type="int",
        help="Number of fine bins for smooth PDF curves")
    parser.add_option(
        "--sigNorm", dest="sigNorm", default=100.0, type="float",
        help="Signal yield at mu=1 (events); must match DataCardMaker value)")
    parser.add_option(
        "--poiName", dest="poiName", default="m",
        help="Name of the mass observable in the workspace")
    parser.add_option(
        "--outDir", dest="outDir", default="plots_postfit",
        help="Directory for output plots")
    parser.add_option(
        "--ext", dest="ext", default="",
        help="Optional extra string appended to output file names")
    parser.add_option(
        "--lumi", dest="lumi", default="",
        help="Luminosity string for the upper-right label, e.g. '27.0 fb^{-1}'")
    parser.add_option(
        "--sqrts", dest="sqrts", default="13.6",
        help="Centre-of-mass energy label (TeV)")
    parser.add_option(
        "--noSignal", dest="drawSignal", action="store_false", default=True,
        help="Suppress the signal curve on the plot")
    parser.add_option(
        "--jsonFile", dest="jsonFile", default="",
        help="Json file to store results in (optional)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Function name for legend from PDF object name
# ---------------------------------------------------------------------------
def pdf_legend_label(pname, is_best_fit):
    if "bern" in pname.lower():
        fname = "Bern. poly"
    elif "exp" in pname.lower():
        fname = "Sum of exp."
    elif "pow" in pname.lower():
        fname = "Power law"
    elif "lau" in pname.lower():
        fname = "Laurent series"
    else:
        fname = pname
    if is_best_fit:
        return "%s (best fit)" % fname
    return fname


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def make_postfit_plot(opt):
    """
    Produce the postfit signal+background plot.

    Parameters
    ----------
    opt : namespace / object with attributes matching the command-line options
          (inputWSFile, fitDiagFile, cat, mass, mMin, mMax, nBins, pdfNBins,
           sigNorm, poiName, outDir, ext, lumi, sqrts, drawSignal, jsonFile)
    """
    print(opt.__dict__)

    # ---------------------------------------------------------------------------
    # Physical mass parameters
    # ---------------------------------------------------------------------------
    xmin      = opt.mMin
    xmax      = opt.mMax
    m_range   = xmax - xmin                          # GeV
    bin_size  = m_range / opt.nBins                  # GeV per display bin

    # ---------------------------------------------------------------------------
    # Helper: copy a normalized [0,1] TH1 into a physical [xmin,xmax] TH1
    # Since the bins are uniformly spaced this is just a relabeling of the x-axis.
    # ---------------------------------------------------------------------------
    def norm_to_phys(h_norm, name):
        N = h_norm.GetNbinsX()
        h_phys = ROOT.TH1D(name, name, N, xmin, xmax)
        for i in range(1, N + 1):
            h_phys.SetBinContent(i, h_norm.GetBinContent(i))
            h_phys.SetBinError(i, h_norm.GetBinError(i))
        h_phys.SetDirectory(0)
        return h_phys

    # ---------------------------------------------------------------------------
    # 1.  Open workspace
    # ---------------------------------------------------------------------------
    print("\n --> Opening workspace : %s" % opt.inputWSFile)
    if not os.path.exists(opt.inputWSFile):
        print("ERROR: workspace file not found: %s" % opt.inputWSFile)
        sys.exit(1)

    f_ws = ROOT.TFile(opt.inputWSFile)
    w    = f_ws.Get("w")
    if not w:
        print("ERROR: workspace 'w' not found in %s" % opt.inputWSFile)
        sys.exit(1)

    # Mass observable (lives in [0, 1] normalized space)
    m = w.var(opt.poiName)
    if not m:
        print("ERROR: variable '%s' not found in workspace" % opt.poiName)
        sys.exit(1)

    m_arglist = ROOT.RooArgList(m)

    # Data
    data = w.data("data_obs")
    if not data:
        print("ERROR: 'data_obs' not found in workspace")
        sys.exit(1)
    print("   data_obs: %.0f events" % data.sumEntries())

    # MultiPdf and index category
    multipdf = w.pdf("multi_pdf")
    if not multipdf:
        print("ERROR: 'multi_pdf' not found in workspace")
        sys.exit(1)

    pdf_index_cat = w.cat("pdf_index")
    if not pdf_index_cat:
        # fallback name used by older DataCardMaker versions
        pdf_index_cat = w.cat("pdf_index_mass_mumu")

    pdfindex_initial = pdf_index_cat.getIndex() if pdf_index_cat else 0

    # Background normalization
    bkg_norm_var = w.var("multi_pdf_norm")
    bkg_norm_fit = bkg_norm_var.getVal() if bkg_norm_var else data.sumEntries()

    # Signal PDF  (may not be present for bkg-only workspaces)
    sig_pdf_name = "model_signal_m_%s" % opt.cat
    sig_pdf      = w.pdf(sig_pdf_name)
    if not sig_pdf:
        print("   WARNING: signal PDF '%s' not found; signal will not be drawn" % sig_pdf_name)
        opt.drawSignal = False


    # ---------------------------------------------------------------------------
    # 2.  Read postfit parameter values from FitDiagnostics tree_fit_sb
    # ---------------------------------------------------------------------------
    r_fit       = 0.0
    pdfindex_fit = pdfindex_initial

    print(" --> Reading postfit values from: %s" % opt.fitDiagFile)
    if not os.path.exists(opt.fitDiagFile):
        print("   WARNING: fitDiagnostics file not found; using pre-fit workspace values")
    else:
        f_diag = ROOT.TFile(opt.fitDiagFile)
        tree_fit_sb = f_diag.Get("tree_fit_sb")

        if tree_fit_sb and tree_fit_sb.GetEntries() > 0:
            tree_fit_sb.GetEntry(0)
            r_fit = float(tree_fit_sb.r)
            print("   Postfit signal strength  r = %.4f" % r_fit)

            # Use fit_s RooFitResult to set ALL floating parameters (including
            # background shape parameters like bern_p0, exp_p0, etc.).
            # tree_fit_sb only stores explicit datacard nuisances, so background
            # shape parameters -- which are freely floating workspace variables but
            # not declared as flatParam/param in the datacard -- are missing from
            # the tree and would silently keep their pre-fit workspace values.
            fit_s = f_diag.Get("fit_s")
            if fit_s:
                final_pars = fit_s.floatParsFinal()
                par_iter   = final_pars.createIterator()
                par = par_iter.Next()
                n_set = 0
                while par:
                    wvar = w.var(par.GetName())
                    if wvar:
                        wvar.setVal(par.getVal())
                        n_set += 1
                    par = par_iter.Next()
                print("   Set %d postfit parameters from fit_s RooFitResult" % n_set)
            else:
                # Fallback: iterate workspace vars and match to tree leaves.
                # This will miss freely-floating background shape parameters.
                print("   WARNING: fit_s not found in fitDiagnostics; "
                      "background shape parameters may use pre-fit values")
                all_vars  = w.allVars()
                var_iter  = all_vars.createIterator()
                v = var_iter.Next()
                while v:
                    vname = v.GetName()
                    leaf  = tree_fit_sb.GetLeaf(vname)
                    if leaf and not v.isConstant():
                        v.setVal(leaf.GetValue())
                    v = var_iter.Next()

            # Refresh bkg norm after parameter update
            if bkg_norm_var:
                bkg_norm_fit = bkg_norm_var.getVal()

            # Update pdf index from tree (the discrete cat is stored as an int branch)
            leaf_idx = tree_fit_sb.GetLeaf("pdf_index")
            if leaf_idx:
                pdfindex_fit = int(leaf_idx.GetValue())
                if pdf_index_cat:
                    pdf_index_cat.setIndex(pdfindex_fit)

            f_diag.Close()

        else:
            print("   WARNING: tree_fit_sb empty; using pre-fit workspace values")
            r_fit = 0.0

    print("   Background norm (postfit) = %.1f events" % bkg_norm_fit)
    print("   Best-fit pdf index        = %d"           % pdfindex_fit)
    sig_events = r_fit * opt.sigNorm
    print("   Signal events (r*norm)    = %.1f"         % sig_events)


    # ---------------------------------------------------------------------------
    # 3.  Build histograms from the individual component PDFs in multiPdf
    # ---------------------------------------------------------------------------
    n_pdfs  = multipdf.getNumPdfs()
    print("\n --> MultiPdf contains %d component PDFs" % n_pdfs)

    bpdfs         = od()   # ordered dict: name -> pdf object
    bpdf_bf_name  = None

    for ipdf in range(n_pdfs):
        pdf   = multipdf.getPdf(ipdf)
        pname = pdf.GetName()
        bpdfs[pname] = pdf
        print("   [%d] %s" % (ipdf, pname))
        if ipdf == pdfindex_fit:
            bpdf_bf_name = pname

    if bpdf_bf_name is None:
        bpdf_bf_name = list(bpdfs.keys())[0]
    print("   Best-fit PDF: %s" % bpdf_bf_name)

    hists = od()

    # Fine-binned PDF histograms (smooth curves) and coarse-binned (for ratio)
    for pname, bpdf in bpdfs.items():

        # Fine-binned: for smooth curve in the upper pad
        h_fine_norm = bpdf.createHistogram(
            "htmp_%s_fine" % pname, m,
            ROOT.RooFit.Binning(opt.pdfNBins, 0., 1.))
        h_fine_norm.Scale(bkg_norm_fit * float(opt.pdfNBins) / float(opt.nBins))
        hists[pname] = norm_to_phys(h_fine_norm, "h_%s" % pname)

        # Coarse-binned: only needed for the best-fit PDF to compute residuals
        if pname == bpdf_bf_name:
            h_coarse_norm = bpdf.createHistogram(
                "htmp_%s_coarse" % pname, m,
                ROOT.RooFit.Binning(opt.nBins, 0., 1.))
            h_coarse_norm.Scale(bkg_norm_fit)
            hists["%s_nBins" % pname] = norm_to_phys(
                h_coarse_norm, "h_%s_nBins" % pname)

    # Data histogram
    h_data_norm = m.createHistogram("htmp_data", ROOT.RooFit.Binning(opt.nBins, 0., 1.))
    h_data_norm.SetBinErrorOption(ROOT.TH1.kPoisson)
    data.fillHistogram(h_data_norm, m_arglist)
    hists["data"] = norm_to_phys(h_data_norm, "h_data")
    hists["data"].SetBinErrorOption(ROOT.TH1.kPoisson)

    # Signal histogram (fine-binned for drawing, coarse-binned for chi2)
    if opt.drawSignal and sig_pdf:
        h_sig_norm = sig_pdf.createHistogram(
            "htmp_sig", m,
            ROOT.RooFit.Binning(opt.pdfNBins, 0., 1.))
        h_sig_norm.Scale(sig_events * float(opt.pdfNBins) / float(opt.nBins))
        hists["signal"] = norm_to_phys(h_sig_norm, "h_signal")

        # Zero out fine bins where expected signal < 0.1 events (per display bin)
        # so the curve doesn't appear as a flat line across the tails
        h_sig = hists["signal"]
        for ibin in range(1, h_sig.GetNbinsX() + 1):
            if h_sig.GetBinContent(ibin) < 0.1:
                h_sig.SetBinContent(ibin, 0.)

        h_sig_coarse_norm = sig_pdf.createHistogram(
            "htmp_sig_coarse", m,
            ROOT.RooFit.Binning(opt.nBins, 0., 1.))
        h_sig_coarse_norm.Scale(sig_events)
        hists["signal_nBins"] = norm_to_phys(h_sig_coarse_norm, "h_signal_nBins")


    # ---------------------------------------------------------------------------
    # 4.  Build ratio histograms  (quantity - best-fit bkg)
    # ---------------------------------------------------------------------------
    hists_ratio = od()

    h_bkg_bf_coarse = hists["%s_nBins" % bpdf_bf_name]

    # Alternative PDFs relative to best-fit
    for pname in bpdfs:
        h_r = hists[pname].Clone("h_ratio_%s" % pname)
        h_r.Add(hists[bpdf_bf_name], -1.)
        h_r.SetDirectory(0)
        hists_ratio[pname] = h_r

    # Data minus best-fit background
    h_ratio_data = hists["data"].Clone("h_ratio_data")
    h_ratio_data.Reset()
    for ibin in range(1, opt.nBins + 1):
        bval   = hists["data"].GetBinContent(ibin)
        berr   = hists["data"].GetBinError(ibin)
        bkgval = h_bkg_bf_coarse.GetBinContent(ibin)
        h_ratio_data.SetBinContent(ibin, bval - bkgval)
        h_ratio_data.SetBinError(ibin, berr)
    h_ratio_data.SetDirectory(0)
    hists_ratio["data"] = h_ratio_data

    # Signal in the ratio pad (same fine histogram, already above zero)
    if "signal" in hists:
        hists_ratio["signal"] = hists["signal"].Clone("h_ratio_signal")
        hists_ratio["signal"].SetDirectory(0)


    # ---------------------------------------------------------------------------
    # 5.  Chi2 / ndof and p-value  (signal + best-fit background vs data)
    # ---------------------------------------------------------------------------
    chi2_val    = 0.0
    n_bins_chi2 = 0
    for ibin in range(1, opt.nBins + 1):
        d  = hists["data"].GetBinContent(ibin)
        b  = h_bkg_bf_coarse.GetBinContent(ibin)
        s  = hists["signal_nBins"].GetBinContent(ibin) if "signal_nBins" in hists else 0.0
        sb = b + s
        if sb > 0:
            chi2_val    += (d - sb) ** 2 / sb
            n_bins_chi2 += 1

    # Free parameters: bkg shape + bkg norm + signal strength r (if signal present)
    bkg_par_set  = bpdfs[bpdf_bf_name].getParameters(ROOT.RooArgSet(m))
    n_free_shape = bkg_par_set.selectByAttrib("Constant", False).getSize()
    n_free_total = n_free_shape + 1                           # +1 for multi_pdf_norm
    if "signal_nBins" in hists:
        n_free_total += 1                                     # +1 for signal strength r
    ndof         = n_bins_chi2 - n_free_total
    chi2_prob    = ROOT.TMath.Prob(chi2_val, ndof) if ndof > 0 else -1.0

    print("\n --> Chi2 / ndof = %.2f / %d  (p-value = %.4f)" % (chi2_val, ndof, chi2_prob))


    # ---------------------------------------------------------------------------
    # 6.  Canvas and pad setup  (matching reference script layout)
    # ---------------------------------------------------------------------------
    ROOT.TGaxis.SetMaxDigits(4)
    ROOT.TGaxis.SetExponentOffset(-0.05, 0.00, "y")

    canv = ROOT.TCanvas("canv_%s" % opt.cat, "canv_%s" % opt.cat, 700, 700)

    pad1 = ROOT.TPad("pad1_%s" % opt.cat, "pad1_%s" % opt.cat, 0, 0.35, 1, 1)
    pad1.SetTickx()
    pad1.SetTicky()
    pad1.SetBottomMargin(0.04)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.05)
    pad1.Draw()

    pad2 = ROOT.TPad("pad2_%s" % opt.cat, "pad2_%s" % opt.cat, 0, 0, 1, 0.35)
    pad2.SetTickx()
    pad2.SetTicky()
    pad2.SetTopMargin(0.03)
    pad2.SetBottomMargin(0.30)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.05)
    pad2.Draw()

    padSizeRatio = 0.65 / 0.35   # ratio of pad1 height to pad2 height

    # Colour palette for background PDFs (index 0 = best-fit, drawn last / on top)
    colour_opts = [
        ROOT.kAzure   + 1,   # best-fit (first entry in bpdfs)
        ROOT.kGreen   + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange  + 1,
        ROOT.kCyan    + 2,
    ]

    xaxis_title = "m_{#mu#mu} [GeV]"
    yaxis_title = "Events / %.1f GeV" % bin_size


    # ===========================================================================
    # Upper pad:  data + background PDFs + signal
    # ===========================================================================
    pad1.cd()

    # Invisible axis frame driven by data range
    h_axes = hists["data"].Clone("h_axes_%s" % opt.cat)
    h_axes.Reset()
    data_max    = hists["data"].GetMaximum()
    data_max_err = hists["data"].GetBinError(hists["data"].GetMaximumBin())
    h_axes.SetMaximum((data_max + data_max_err) * 2.2)
    h_axes.SetMinimum(0.)
    h_axes.SetTitle("")
    h_axes.GetXaxis().SetTitle("")
    h_axes.GetXaxis().SetLabelSize(0)
    h_axes.GetYaxis().SetTitle(yaxis_title)
    h_axes.GetYaxis().SetTitleSize(0.062)
    h_axes.GetYaxis().SetTitleOffset(0.95)
    h_axes.GetYaxis().SetLabelSize(0.050)
    h_axes.GetYaxis().SetLabelOffset(0.007)
    h_axes.Draw()

    # ---- Legends side-by-side in the upper headroom ----
    #leg of data + signal
    leg0 = ROOT.TLegend(0.13, 0.60, 0.54, 0.8)
    leg0.SetFillStyle(0)
    leg0.SetLineColor(0)
    leg0.SetTextSize(0.045)
    leg0.AddEntry(hists["data"], "Data", "ep")
    if "signal" in hists:
        leg0.AddEntry(
            hists["signal"],
            "Signal (#mu=%.2f), m=%.0f GeV" % (r_fit, opt.mass),
            "fl")

    # leg for backgrounds
    leg = ROOT.TLegend(0.54, 0.6, 0.93, 0.80)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetTextSize(0.05)

    # ---- Signal ----
    if "signal" in hists:
        hists["signal"].SetLineWidth(2)
        hists["signal"].SetLineColor(ROOT.kRed)
        #hists["signal"].SetFillColor(ROOT.kRed - 9)
        #hists["signal"].SetFillStyle(3001)
        hists["signal"].Draw("HIST same")

    # ---- Background PDFs ----
    bnames = list(bpdfs.keys())
    for bidx, pname in enumerate(bnames):
        col = colour_opts[bidx % len(colour_opts)]
        is_bf = (pname == bpdf_bf_name)
        hists[pname].SetLineWidth(3 if is_bf else 2)
        hists[pname].SetLineColor(col)
        hists[pname].SetLineStyle(ROOT.kSolid if is_bf else ROOT.kDashed)
        hists[pname].Draw("HIST same C")
        leg.AddEntry(hists[pname], pdf_legend_label(pname, is_bf), "l")

    # ---- Data points ----
    hists["data"].SetMarkerStyle(20)
    hists["data"].SetMarkerColor(ROOT.kBlack)
    hists["data"].SetLineColor(ROOT.kBlack)
    hists["data"].SetMarkerSize(0.85)
    hists["data"].Draw("Same PE")

    leg0.Draw("Same")
    leg.Draw("Same")

    # ---- Chi2 / ndof label ----
    lat_chi2 = ROOT.TLatex()
    lat_chi2.SetTextFont(42)
    lat_chi2.SetTextAlign(11)
    lat_chi2.SetNDC()
    lat_chi2.SetTextSize(0.05)
    lat_chi2.DrawLatex(
        0.17, 0.54,
        "#chi^{2}/ndof = %.1f/%d  (p = %.2f)" % (chi2_val, ndof, chi2_prob))

    # ---- CMS / lumi labels ----
    lat_cms = ROOT.TLatex()
    lat_cms.SetTextFont(61)
    lat_cms.SetTextAlign(11)
    lat_cms.SetNDC()
    lat_cms.SetTextSize(0.075)
    lat_cms.DrawLatex(0.14, 0.92, "CMS")

    lat_prelim = ROOT.TLatex()
    lat_prelim.SetTextFont(52)
    lat_prelim.SetTextAlign(11)
    lat_prelim.SetNDC()
    lat_prelim.SetTextSize(0.058)
    lat_prelim.DrawLatex(0.14, 0.83, "Preliminary")

    if opt.lumi:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetNDC()
        lat_lumi.SetTextSize(0.050)
        lat_lumi.DrawLatex(
            0.95, 0.92,
            "%s (%s TeV)" % (opt.lumi, opt.sqrts))
    elif opt.sqrts:
        lat_lumi = ROOT.TLatex()
        lat_lumi.SetTextFont(42)
        lat_lumi.SetTextAlign(31)
        lat_lumi.SetNDC()
        lat_lumi.SetTextSize(0.050)
        lat_lumi.DrawLatex(0.95, 0.92, "(%s TeV)" % opt.sqrts)


    # ===========================================================================
    # Lower pad:  data − best-fit background  (+ alternative PDFs)
    # ===========================================================================
    pad2.cd()

    # Y-range driven by data residuals
    res_vals = [abs(hists_ratio["data"].GetBinContent(i))
                + hists_ratio["data"].GetBinError(i)
                for i in range(1, opt.nBins + 1)]
    res_max  = max(res_vals) if res_vals else 1.0
    res_max  = max(res_max * 1.5, 5.0)

    h_axes_ratio = hists_ratio["data"].Clone("h_axes_ratio_%s" % opt.cat)
    h_axes_ratio.Reset()
    h_axes_ratio.SetMaximum( res_max)
    h_axes_ratio.SetMinimum(-res_max)
    h_axes_ratio.SetTitle("")
    h_axes_ratio.GetXaxis().SetTitle(xaxis_title)
    h_axes_ratio.GetXaxis().SetTitleSize(0.08 * padSizeRatio)
    h_axes_ratio.GetXaxis().SetTitleOffset(0.60)
    h_axes_ratio.GetXaxis().SetLabelSize(0.046 * padSizeRatio)
    h_axes_ratio.GetXaxis().SetLabelOffset(0.007)
    h_axes_ratio.GetXaxis().SetTickLength(0.03 * padSizeRatio)
    h_axes_ratio.GetYaxis().SetTitle("Data #minus Bkg")
    h_axes_ratio.GetYaxis().SetTitleSize(0.065 * padSizeRatio)
    h_axes_ratio.GetYaxis().SetTitleOffset(0.40)
    h_axes_ratio.GetYaxis().SetLabelSize(0.042 * padSizeRatio)
    h_axes_ratio.GetYaxis().SetLabelOffset(0.007)
    h_axes_ratio.GetYaxis().SetNdivisions(505)
    h_axes_ratio.Draw()

    # Zero line
    zero_line = ROOT.TLine(xmin, 0., xmax, 0.)
    zero_line.SetLineColor(ROOT.kBlack)
    zero_line.SetLineStyle(ROOT.kDashed)
    zero_line.SetLineWidth(1)
    zero_line.Draw()

    # Signal overlay in ratio pad
    if "signal" in hists_ratio:
        hists_ratio["signal"].SetLineWidth(2)
        hists_ratio["signal"].SetLineColor(ROOT.kRed)
        hists_ratio["signal"].SetFillColor(ROOT.kRed - 9)
        hists_ratio["signal"].SetFillStyle(3001)
        hists_ratio["signal"].Draw("HIST same")

    # Alternative background PDFs (show deviation from best-fit)
    for bidx, pname in enumerate(bnames):
        if pname == bpdf_bf_name:
            continue
        col = colour_opts[bidx % len(colour_opts)]
        hists_ratio[pname].SetLineWidth(2)
        hists_ratio[pname].SetLineColor(col)
        hists_ratio[pname].SetLineStyle(ROOT.kDashed)
        hists_ratio[pname].Draw("HIST same C")

    # Data residuals
    hists_ratio["data"].SetMarkerStyle(20)
    hists_ratio["data"].SetMarkerColor(ROOT.kBlack)
    hists_ratio["data"].SetLineColor(ROOT.kBlack)
    hists_ratio["data"].SetMarkerSize(0.85)
    hists_ratio["data"].Draw("Same PE")

    # Ratio pad sub-label
    lat_ratio = ROOT.TLatex()
    lat_ratio.SetTextFont(42)
    lat_ratio.SetTextAlign(33)
    lat_ratio.SetNDC()
    lat_ratio.SetTextSize(0.042 * padSizeRatio)
    lat_ratio.DrawLatex(0.89, 0.93, "Best fit bkg subtracted")

    canv.Update()


    # ---------------------------------------------------------------------------
    # 7.  Save output
    # ---------------------------------------------------------------------------
    if not os.path.isdir(opt.outDir):
        os.makedirs(opt.outDir)

    ext_str  = ("_" + opt.ext) if opt.ext else ""
    out_base = os.path.join(opt.outDir, "postfit_multipdf_%s%s" % (opt.cat, ext_str))

    canv.SaveAs(out_base + ".pdf")
    canv.SaveAs(out_base + ".png")
    print("\n --> Saved: %s.pdf" % out_base)
    print(" --> Saved: %s.png" % out_base)

    f_ws.Close()


    # ---------------------------------------------------------------------------
    # 8.  Edit output json with chi2 (if provided)
    # ---------------------------------------------------------------------------
    if len(opt.jsonFile) > 0 and os.path.exists(opt.jsonFile):
        with open(opt.jsonFile, 'r') as f:
            results = json.load(f)

        results['sbfit_chi2'] = chi2_val
        results['sbfit_ndof'] = ndof
        results['sbfit_prob'] = chi2_prob

        with open(opt.jsonFile, 'w') as f:
            json.dump(results, f, indent=4)


if __name__ == "__main__":
    (opt, args) = get_options()
    make_postfit_plot(opt)
