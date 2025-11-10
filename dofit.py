import os
import pickle
import json
import random
import optparse
from Fitter import Fitter
from DataCardMaker import DataCardMaker
from Utils import *
from array import array
from fit_signalshapes import  fit_signalmodel

def dofit(options):

    label = options.label
    plot_dir = options.plotDir
    if(plot_dir[-1] != '/'):
        plot_dir += '/'

    if(not os.path.exists(plot_dir)):
        os.system("mkdir %s" % plot_dir)

    if(os.path.exists(plot_dir + "fit_results_{}.json".format(options.mass))):
        #remove old results
        os.system("rm %s" % plot_dir + "fit_results_{}.json".format(options.mass))

    fine_bin_size = 0.05
    mass = options.mass 
    #binsx = list(np.arange(11, 30, 0.5))
    binsx = list(np.arange(11, 20, 0.5))

    if(options.m_max < 0. and options.rebin): 
        options.m_max = get_m_max(options.inputFile) + 5.0
        options.m_max = max(1.2*options.mass, options.m_max)
    #print("m MAX %.2f" % options.m_max)
    
    if(options.m_min > 0 and options.m_min < binsx[-1]):
        start_idx = 0
        while(binsx[start_idx] < options.m_min):
            start_idx +=1
        binsx = binsx[start_idx:]

        if(abs(options.m_min - binsx[0]) < 50.):
            binsx[0] = options.m_min
        else:
            binsx.insert(0, options.m_min)
        print("Will start fit from %.0f GeV" % binsx[0])

    if(options.m_max > 0 and options.m_max < binsx[-1]):
        print("rebinning with max m %.2f" % options.m_max)
        end_idx = len(binsx)-1
        while(binsx[end_idx]   > options.m_max and end_idx > 0): 
            end_idx -=1
        binsx = binsx[:end_idx]

        if(abs(options.m_max - binsx[-1]) < 50.):
            binsx[-1] = options.m_max
        else:
            binsx.append(options.m_max)
        print("Will end fit at %.0f GeV" % binsx[-1])
        print(binsx)


    # round to smallest precision we are storing mass values with, otherwise
    # get weird effects related to bin size
    #roundTo(binsx, fine_bin_size)



    nbins_fine = int((binsx[-1] - binsx[0])/fine_bin_size)
    print('nbins_fine', nbins_fine)

    histos_sb = ROOT.TH1F("m_sb", "m_sb" ,nbins_fine, binsx[0], binsx[-1])
    
    
    load_h5_sb(options.inputFile, histos_sb)
    print("************ Found %i total events \n" % histos_sb.GetEntries())
    print(histos_sb.Integral())

    

    if(options.rebin):
        bins_nonzero = get_rebinning(binsx, histos_sb)
        print("Rebinning to avoid zero bins!")
        print("old", binsx)
        print("new", bins_nonzero)
        bins = bins_nonzero
    else:
        bins = binsx
    roobins = ROOT.RooBinning(len(bins)-1, array('d', bins), "mbins")

    if(options.refit_sig):
        print ("########## FIT SIGNAL AND SAVE PARAMETERS ############")
        sig_file_name = "sig_fit.root"

        fit_signalmodel(options.inputFile, sig_file_name, mass, binsx, nbins_fine, plot_dir,return_fit=False,
                        dcb_model=options.dcb_model)

    else:  # use precomputed signal shape
        if(not os.path.exists(options.sig_shape)):
            print("Sig file %s doesn't exist" % options.sig_shape)
            exit(1)
        sig_file_name = options.sig_shape



    print("\n\n ############# FIT BACKGROUND AND SAVE PARAMETERS ###########")
    #orderToTry = [2, 3, 4, 5]

    #orderToTry = [2, 3, 4]
    orderToTry = [3,4]
    chi2s = [0]*len(orderToTry)
    fit_params = [0] * len(orderToTry)
    ndofs = [0]*len(orderToTry)
    probs = [0]*len(orderToTry)
    fit_errs = [0]*len(orderToTry)
    bkg_fnames = [""]*len(orderToTry)

    #polyExp, bern, exp, 
    #func_form = "polyExp"
    func_form = "exp"
    #func_form = "bern"


    if options.blinded:
        print("BLIND FIT TO DO ")
        exit(1)

    fitting_histogram = histos_sb
    data_name = "data_bkg"

    for i, order in enumerate(orderToTry):
        print("Trying %i parameter background fit" % order)
        bkg_fnames[i] = str(order) + 'par_bkg_fit%i.root' % i
        bkg_outfile = ROOT.TFile(bkg_fnames[i], 'RECREATE')

        nPars = get_nPars(order, func_form)

        model_name = "model_b" + str(i)
        fitter_bkg = Fitter(['m_fine'], debug = False)
        fitter_bkg.bkgShape(name=model_name, poi='m_fine', order=order, func_form=func_form )
        fitter_bkg.importBinnedData(fitting_histogram, ['m_fine'], data_name)
        
        fres = fitter_bkg.fit(model_name, data_name, options=[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0),  ROOT.RooFit.Minos(1), ROOT.RooFit.Minimizer("Minuit2")])
        #Running fit two times seems to improve things sometimes (better initial guesses for params?)
        #fres = fitter_bkg.fit(model_name, data_name, options=[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0),  ROOT.RooFit.Minos(1), ROOT.RooFit.Minimizer("Minuit2")])

        chi2_fine = fitter_bkg.projection(
            model=model_name, data=data_name, poi="m_fine",
            filename=plot_dir + str(order) + "par_bkg_fit.png", binning=0, logy=False)

        #chi2_binned = fitter_bkg.projection(
        #     model_name, data_name, "m_fine",
        #     plot_dir + str(order) + "par_bkg_fit_binnedx.png",roobins,True)


        bkg_outfile.cd()

        m = fitter_bkg.getVar('m_fine')
        m.setBins(nbins_fine)
        model = fitter_bkg.getFunc(model_name)
        dataset = fitter_bkg.getData(data_name)

        #rescale so pdfs are in evts per 0.5 GeV
        fit_range_low = roobins.lowBound()
        fit_range_high = roobins.highBound()
        n = roobins.numBoundaries() - 1
        #RootFit default normalization is full range divided by number of bins
        default_norm = (fit_range_high - fit_range_low)/ n
        rescale = 0.5/ default_norm
        fit_norm = ROOT.RooFit.Normalization(rescale,ROOT.RooAbsReal.Relative)

        #use toys to sample errors rather than linear method, 
        #needed b/c fn's usually have strong correlation of params
        linear_errors = False


        frame = m.frame()
        dataset.plotOn(frame, ROOT.RooFit.Name(data_name), ROOT.RooFit.Invisible(), ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), 
                ROOT.RooFit.Rescale(rescale))

        model.plotOn(frame, ROOT.RooFit.VisualizeError(fres, 1, linear_errors), ROOT.RooFit.FillColor(ROOT.kRed - 7), ROOT.RooFit.LineColor(ROOT.kRed - 7), ROOT.RooFit.Name(fres.GetName()), 
                       fit_norm)

        model.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed + 1), ROOT.RooFit.Name(model_name),  fit_norm)
        

        useBinAverage = True
        hpull = frame.pullHist(data_name, model_name, useBinAverage)
        hresid = frame.residHist(data_name, model_name, False, useBinAverage)
        dhist = ROOT.RooHist(frame.findObject(data_name, ROOT.RooHist.Class()))



        #redraw data (so on top of model curves)
        if(options.rebin):
            dataset.plotOn(frame, ROOT.RooFit.Name(data_name),   ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
                       ROOT.RooFit.Rescale(rescale))
        else:
            dataset.plotOn(frame, ROOT.RooFit.Name(data_name),  ROOT.RooFit.XErrorSize(0), ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
                       ROOT.RooFit.Rescale(rescale))

        useBinAverage = True
        hpull = frame.pullHist(data_name, model_name, useBinAverage)
        hresid = frame.residHist(data_name, model_name, False, useBinAverage)
        dhist = ROOT.RooHist(frame.findObject(data_name, ROOT.RooHist.Class()))



        #get fractional error on fit
        central = frame.getCurve(model_name);
        curve =  frame.getCurve("fitresults");
        upBound = ROOT.TGraph(central.GetN());
        loBound = ROOT.TGraph(central.GetN());
        norm = get_roohist_sum(dhist)

        for j in range(curve.GetN()):
            if( j < central.GetN() ): upBound.SetPoint(j, curve.GetX()[j], curve.GetY()[j]);
            else: loBound.SetPoint( 2*central.GetN() - j, curve.GetX()[j], curve.GetY()[j]);


        fit_hist = model.createHistogram("h_model_fit", m, ROOT.RooFit.Binning(roobins))
        print("Hist fit norm", fit_hist.Integral())
        fit_hist.Scale(norm / fit_hist.Integral())

        #Get hist of pulls:  (data - fit) / tot_unc
        hresid_norm = get_pull_hist(model, frame, central, curve, hresid, fit_hist,  bins)


        #abs because somtimes order is reversed
        err_on_sig = abs(upBound.Eval(options.mass) - loBound.Eval(options.mass))/2.
        frac_err_on_sig = err_on_sig / central.Eval(options.mass)
        bkg_fit_frac_err = frac_err_on_sig


        my_chi2, my_ndof = calculateChi2(hpull, nPars, excludeZeros = True, dataHist = dhist)
        my_prob = ROOT.TMath.Prob(my_chi2, my_ndof)

        PlotFitResults(frame, fres.GetName(), nPars, hresid_norm, data_name,
                       [model_name], my_chi2, my_ndof,
                       f"{func_form}_order{order}_bkg_fit_binned",
                       plot_dir, plot_label = label)


        #TODO
        graphs = {}
        for parName in fitter_bkg.par_names:
            graphs[pa] = ROOT.TGraphErrors()

        #largest_frac_err = 0.
        bkg_fit_params = dict()
        for var, graph in graphs.items():
            print(var)
            value, error = fitter_bkg.fetch(var)
            bkg_fit_params[var] = (value, error)
            graph.SetPoint(0, mass, value)
            graph.SetPointError(0, 0.0, error)
            #frac_err = abs(error/value)
            #largest_frac_err = max(frac_err, largest_frac_err)
        bkg_fit_params['cov'] = convert_matrix(fres.covarianceMatrix())
        print(bkg_fit_params['cov'])

        bkg_outfile.cd()
        for name, graph in graphs.items():
            graph.Write(name)
        bkg_outfile.Close()

        print("#############################")
        print("% Order %i results: " % order)
        print("bkg fit chi2/nbins (fine binning) ", chi2_fine)
        print("My chi2, ndof, prob", my_chi2, my_ndof, my_prob)
        print("My chi/ndof, chi2/nbins", my_chi2/my_ndof, my_chi2/(my_ndof + nPars))
        print("Fit func fractional unc at sig mass ", bkg_fit_frac_err)
        print("#############################")

        chi2s[i] = my_chi2
        ndofs[i] = my_ndof
        probs[i] = my_prob
        fit_params[i] = bkg_fit_params
        fit_errs[i] = bkg_fit_frac_err
        fitter_bkg.delete()

    exit(1)

    best_i = f_test(nParsToTry, ndofs, chi2s, fit_errs, thresh = options.ftest_thresh, err_thresh = options.err_thresh)
    nPars_bkg = nParsToTry[best_i]
    bkg_fname = bkg_fnames[best_i]
    print("\n Chose %i parameters based on F-test ! \n" % nPars_bkg)

    # Fit to total data
    #histos_sb.Print("range")
    #histos_bkg.Print("range")
    #histos_sig.Print("range")

    sb_fname = "sb_fit.root"
    sb_outfile = ROOT.TFile(sb_fname, 'RECREATE')
    sb_outfile.cd()
    histos_sb.Write("m_sb")
    sb_outfile.Close()
    sig_data_name = 'm_sb'
    sb_label = "raw"

    card = DataCardMaker(sb_label)
    if options.dcb_model:
        card.addDCBSignalShape('model_signal_m', 'm', sig_file_name,
                               {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})
    else:
        card.addSignalShape('model_signal_m', 'm', sig_file_name,
                            {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})



    sig_norm = card.addFixedYieldFromFile('model_signal_m', 0, sig_file_name,
                                          "m_sig", norm = options.sig_norm)
    #sig_norm = card.addFloatingYield('model_signal_m', 0, sig_file_name,
    #                                 "m_sig", constant=False)
    card.addSystematic("CMS_scale_j", "param", [0.0, options.scale_j_unc])
    card.addSystematic("CMS_res_j", "param", [0.0, options.res_j_unc])

    exit(1)
    card.addBkgShapeNoTag('model_bkg_m', 'm', bkg_fname, order=order_bkg)
    card.addFloatingYield('model_bkg_m', 1, sb_fname, "m_sb")

    for i in range(0, nPars_bkg):
        card.addSystematic("CMS_mumu_p%i" % i, "flatParam", [])

    card.addSystematic("model_bkg_m_mumu_norm", "flatParam", [])
    card.importBinnedData("sb_fit.root", sig_data_name,
                          ["m"], 'data_obs', 1.0)

    if(options.sig_norm_unc > 0):
        card.addSystematic("SigEff", "lnN", values = {"model_signal_m" : 1. + options.sig_norm_unc})
    card.makeCard()
    card.delete()



    cmd = (
        "text2workspace.py datacard_mumu_{l2}.txt "
        + "-o workspace_mumu_{l1}_{l2}.root "
        + "&& combine -M FitDiagnostics workspace_mumu_{l1}_{l2}.root --cminPreFit 1 "
        + "-m {mass} -n _{l1}_{l2} --robustFit 1"
        + "&& combine -M Significance workspace_mumu_{l1}_{l2}.root --usePLC "
        + "-m {mass} -n significance_{l1}_{l2} "
        + "&& combine -M Significance workspace_mumu_{l1}_{l2}.root --usePLC "
        + "-m {mass} --pvalue -n pvalue_{l1}_{l2} "
        + "&& combine -M AsymptoticLimits workspace_mumu_{l1}_{l2}.root "
        + "-m {mass} -n lim_{l1}_{l2} "
        ).format(mass=mass, l1=label, l2=sb_label)
    print(cmd)
    os.system(cmd)
    workspace_name = 'workspace_mumu_{l1}_{l2}.root'.format(l1=label, l2=sb_label)
    sbfit_chi2, sbfit_ndof = checkSBFit(workspace_name, sb_label, bins, label + "_" + sb_label, nPars_bkg, 
            plot_dir = plot_dir, draw_sig = options.draw_sig, plot_label = label)

    sbfit_prob = ROOT.TMath.Prob(sbfit_chi2, sbfit_ndof)

    f_signif_name = ('higgsCombinesignificance_{l1}_{l2}.'
                     + 'Significance.mH{mass:.0f}.root'
                     ).format(mass=mass, l1=label, l2=sb_label)
    f_exp_signif_name = ('higgsCombine_exp_significance_{l1}_{l2}.'
                     + 'Significance.mH{mass:.0f}.root'
                     ).format(mass=mass, l1=label, l2=sb_label)
    f_limit_name = ('higgsCombinelim_{l1}_{l2}.'
                    + 'AsymptoticLimits.mH{mass:.0f}.root'
                    ).format(mass=mass, l1=label, l2=sb_label)
    f_pval_name = ('higgsCombinepvalue_{l1}_{l2}.'
                   + 'Significance.mH{mass:.0f}.root'
                   ).format(mass=mass, l1=label, l2=sb_label)
    f_diagnostics_name = ('fitDiagnostics_{l1}_{l2}.root'
                   ).format(l1=label, l2=sb_label)

    f_signif = ROOT.TFile(f_signif_name, "READ")
    res1 = f_signif.Get("limit")
    res1.GetEntry(0)
    signif = res1.limit
    print("Significance is %.3f \n" % signif)

    f_diagnostics = ROOT.TFile(f_diagnostics_name, "READ")
    params = f_diagnostics.Get("tree_fit_sb")
    params.GetEntry(0)
    sig_strength = params.r
    sig_strength_unc = params.rErr


    #expected significance
    print('sig_norm %.3f' % sig_norm)



    true_sig_strength = get_sig_in_window(options.inputFile, binsx[0], binsx[-1]) /  sig_norm
    print("True sig strength %.3f" % true_sig_strength)

    cmd = ("combine -M Significance workspace_mumu_{l1}_{l2}.root -t -1 --expectSignal %.3f --toysFreq " % (true_sig_strength)
        + "-m {mass} -n _exp_significance_{l1}_{l2} ").format(mass = mass, l1 = label, l2 = sb_label)
    print(cmd)
    os.system(cmd)

    f_exp_signif = ROOT.TFile(f_exp_signif_name, "READ")
    res_e = f_exp_signif.Get("limit")
    res_e.GetEntry(0)
    exp_signif = res_e.limit

    exp_pval = 0.5-(0.5*(1+ROOT.Math.erf(exp_signif/np.sqrt(2)))-0.5*(1+ROOT.Math.erf(0/np.sqrt(2))))
    print("Asimov significance is %.3f \n" % exp_signif)


    f_limit = ROOT.TFile(f_limit_name, "READ")
    res2 = f_limit.Get("limit")
    eps = 0.01
    obs_limit = -1
    exp_limit = exp_low = exp_high = exp_two_low = exp_two_high = -1

    
    for i in range(6):
        res2.GetEntry(i)
        if(res2.quantileExpected == -1):  # obs limit
            obs_limit = res2.limit
        elif(abs(res2.quantileExpected - 0.5) < eps):  # exp limit
            exp_limit = res2.limit
        elif(abs(res2.quantileExpected - 0.025) < eps):  # 2sigma, low
            exp_two_low = res2.limit
        elif(abs(res2.quantileExpected - 0.16) < eps):  # 1sigma, low
            exp_low = res2.limit
        elif(abs(res2.quantileExpected - 0.84) < eps):  # 1sigma, high
            exp_high = res2.limit
        elif(abs(res2.quantileExpected - 0.975) < eps):  # 2sigma, high
            exp_two_high = res2.limit

    print("Obs limit is %.3f (%.1f events)" % (obs_limit, obs_limit*sig_norm))
    print("Expected was %.3f (%.1f events)" % (exp_limit, exp_limit*sig_norm))
    print("Expected range %.1f-%.1f (one sigma), %.1f-%.1f (two sigma)" % (exp_low * sig_norm, exp_high*sig_norm, exp_two_low * sig_norm, exp_two_high * sig_norm))

    f_pval = ROOT.TFile(f_pval_name, "READ")
    res3 = f_pval.Get("limit")
    res3.GetEntry(0)
    pval = res3.limit
    print("p-value is %.3f \n" % pval)

    #signal yeild in +/- 2 sigma
    check_rough_sig(options.inputFile, options.mass*0.9, options.mass*1.1)

    f_diagnostics = ROOT.TFile(f_diagnostics_name, "READ")
    f_diagnostics.ls()
    params = f_diagnostics.Get("tree_fit_sb")
    params.GetEntry(0)
    sig_strength = params.r
    sig_strength_unc = params.rErr
    print('r ', sig_strength, 'unc', sig_strength_unc)


    f_signif.Close()
    f_limit.Close()
    f_pval.Close()
    f_diagnostics.Close()



    results = dict()

    # bkg fit results
    results['bkgfit_chi2'] = chi2s[best_i]
    results['bkgfit_ndof'] = ndofs[best_i]
    results['bkgfit_prob'] = probs[best_i]
    results['bkgfit_frac_err'] = fit_errs[best_i]
    results['bkg_fit_params'] = fit_params[best_i]
    results['sbfit_chi2'] = sbfit_chi2
    results['sbfit_ndof'] = sbfit_ndof
    results['sbfit_prob'] = sbfit_prob
    results['nPars_bkg'] = nPars_bkg
    results['signif'] = signif
    results['asimov_signif'] = exp_signif
    results['asimov_pval'] = exp_pval
    results['pval'] = pval
    results['obs_excess_events'] = sig_strength*sig_norm
    results['obs_excess_events_unc'] = sig_strength_unc*sig_norm
    results['obs_lim_events'] = obs_limit*sig_norm
    results['exp_lim_events'] = exp_limit*sig_norm
    results['exp_lim_1sig_low'] = exp_low * sig_norm
    results['exp_lim_2sig_low'] = exp_two_low * sig_norm
    results['exp_lim_1sig_high'] = exp_high * sig_norm
    results['exp_lim_2sig_high'] = exp_two_high* sig_norm
    results['sig_norm_unc'] = options.sig_norm_unc
    results['mass'] = options.mass
    results['m_min'] = options.m_min
    results['m_max'] = options.m_max
    results['script_options'] = vars(options)

    print("Saving fit results to %s" % plot_dir + "fit_results_{}.json".format(options.mass))
        
    with open(plot_dir + "fit_results_{}.json".format(options.mass), "w") as jsonfile:
        json.dump(results, jsonfile, indent=4)

    return results


def fitting_options():
    parser = optparse.OptionParser()
    parser.add_option("--scale_j_unc", type=float, default=0.01,
                      help="Uncertainty on signal mean from JES")
    parser.add_option("--res_j_unc", type=float, default=0.035,
                      help="Uncertainty on signal width from JER")

    parser.add_option("--m_min", type=float, default=-1.0,
                      help="Minimum m for the fit")
    parser.add_option("--m_max", type=float, default=-1.0,
                      help="Maximum m for the fit")
    parser.add_option("--sig_norm", type=float, default=100.0,
                      help="Signal normalization (definition of mu=1)")
    parser.add_option("--ftest_thresh", type=float, default=0.05,
                      help="Threshold to prefer a function in the f-test")
    parser.add_option("--err_thresh", type=float, default=0.5,
                      help="Threshold on fit unc to be included in f-test")
    parser.add_option("-s", "--sig_shape", default="sig_shape_M15.root",
                      help="Pre-saved signal shape")
    parser.add_option("--refit_sig", default=False, action="store_true",
                      help="""Fit the signal events (using truth labels)
                      to get signal shape""")
    parser.add_option("--rebin", default=False, action="store_true",
                      help="""Rebin  to make sure no bins less than 5 evts""")
    parser.add_option("-M", "-M", dest="mass", type=float, default=15.,
                      help="Signal mass hypothesis")
    parser.add_option("-i", "--inputFile", dest="inputFile",
                      default='fit_inputs/no_selection_03p.h5',
                      help="input h5 file")
    parser.add_option("-p", "--plotDir", dest="plotDir", default='plots/',
                      help="Where to put the plots")
    parser.add_option("-l", "--label", dest="label", default='test',
                      help="Label for file names")
    parser.add_option("--no_draw_sig", dest="draw_sig", action = 'store_false', help="Don't draw separate signal and bkg contribution on S+B fit plots")
    parser.set_defaults(draw_sig = True)
    parser.add_option("-b", "--blinded", dest="blinded", action="store_true",
                      default=False,
                      help="Blinding the signal region for the fit.")
    parser.add_option("--dcb-model", dest="dcb_model", action="store_true",
                      default=False,
                      help="""Whether to use double crystal ball model for signal shape instead
                      of default model (gaussian core with single crystal ball)""")
    parser.add_option("--sig_norm_unc", dest="sig_norm_unc", type=float, default= -1.0, help="Fractional uncertainty on signal normalization (for limits)")
    return parser


if __name__ == "__main__":
    parser = fitting_options()
    (options, args) = parser.parse_args()
    dofit(options)

