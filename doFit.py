import os
import pickle
import json
import random
import optparse
import math
from Fitter import Fitter, shape_map
from DataCardMaker import DataCardMaker
from Utils import *
from array import array
from fit_signalshapes import  fit_signalmodel

def dofit(options):

    label = options.label
    plot_dir = options.outDir
    if(plot_dir[-1] != '/'):
        plot_dir += '/'

    if(not os.path.exists(plot_dir)):
        os.system("mkdir %s" % plot_dir)

    if(os.path.exists(plot_dir + "fit_results_{}.json".format(options.mass))):
        #remove old results
        os.system("rm %s" % plot_dir + "fit_results_{}.json".format(options.mass))

    fine_bin_size = 0.1
    mass = options.mass 
    binsx = list(np.arange(options.m_min, options.m_max + options.bin_size, options.bin_size))


    nbins_fine = int((binsx[-1] - binsx[0])/fine_bin_size)
    print('nbins_fine', nbins_fine)

    #transform to 0 to 1 range
    xmin, xmax = binsx[0], binsx[-1]

    histos_sb = ROOT.TH1F("m_sb", "m_sb" ,nbins_fine, 0., 1.)
    
    
    load_h5_sb(options.inputFile, histos_sb, xmin=xmin, xmax=xmax)
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

    sb_fname = plot_dir + "sb_fit.root"
    sb_outfile = ROOT.TFile(sb_fname, 'RECREATE')
    sb_outfile.cd()
    histos_sb.Write("m_sb")
    sb_outfile.Close()
    sig_data_name = 'm_sb'
    fit_label = "mumu"
    poi_name = "m"

    print(" CREATE WORKSPACE ") 
    card = DataCardMaker(fit_label, outdir=plot_dir)

    card.importBinnedData(sb_fname, sig_data_name,
                          poi_name, 'data_obs', 1.0)

    if options.dcb_model:
        card.addDCBSignalShape('model_signal_m', sig_file_name,
                               {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0}, 
                               xmin=xmin, xmax=xmax)
    else:
        card.addSignalShape('model_signal_m',  sig_file_name,
                            {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0},
                            xmin=xmin, xmax=xmax)


    sig_norm = card.addFixedYieldFromFile('model_signal_m', 0, sig_file_name,
                                          "m_sig", norm = options.sig_norm)
    card.addSystematic("CMS_scale_j", "param", [0.0, options.scale_j_unc])
    card.addSystematic("CMS_res_j", "param", [0.0, options.res_j_unc])

    if options.blinded:
        print("BLIND FIT TO DO ")
        exit(1)

    # Normalized bin edges in 0-1 space (workspace variable range) — used for post-fit plot
    bins_norm = [(b - xmin) / (xmax - xmin) for b in binsx]

    fitting_histogram = histos_sb
    data_name = "data_bkg"

    func_forms = {
            "exp": [1, 2, 3, 4],
            "bern": [2, 3, 4, 5], 
            }
    for func_form, orderToTry in func_forms.items():
        print("\n \n Fitting with functional form %s " % func_form)

        chi2s = [0]*len(orderToTry)
        fit_params = [0] * len(orderToTry)
        ndofs = [0]*len(orderToTry)
        probs = [0]*len(orderToTry)
        fit_errs = [0]*len(orderToTry)
        bkg_fnames = [""]*len(orderToTry)

        for i, order in enumerate(orderToTry):
            print("Trying %i parameter background fit" % order)
            bkg_fnames[i] = plot_dir + func_form + "_" + str(order) + 'par_bkg_fit%i.json' % i

            nPars = get_nPars(order, func_form)

            model_name = "model_b" + str(i)
            fitter_bkg = Fitter(['m_fine'], debug = False, outdir=plot_dir)
            fitter_bkg.importBinnedData(fitting_histogram, 'm_fine', data_name)
            fitter_bkg.bkgShape(name=model_name, poi='m_fine', order=order, func_form=func_form )
            
            fres = fitter_bkg.fit(model_name, data_name, options=[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0),  ROOT.RooFit.Minos(1), ROOT.RooFit.Minimizer("Minuit2")])
            #Running fit two times sometimes seems to improve things sometimes (better initial guesses for params)
            #fres = fitter_bkg.fit(model_name, data_name, options=[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0),  ROOT.RooFit.Minos(1), ROOT.RooFit.Minimizer("Minuit2")])

            chi2_fine, ndof_fine = fitter_bkg.projection(
                model=model_name, data=data_name, poi="m_fine",
                filename=plot_dir +func_form + "_" + str(order) + "par_bkg_fit.png", binning=0, logy=False)

            chi2_prob = ROOT.TMath.Prob(chi2_fine, ndof_fine)


            bkg_fit_params = dict()
            for parName in fitter_bkg.par_names:
                value, error = fitter_bkg.fetch(parName)
                bkg_fit_params[parName] = (value, error)

            bkg_fit_params['par_names'] = fitter_bkg.par_names
            bkg_fit_params['cov'] = convert_matrix(fres.covarianceMatrix())


            with open(bkg_fnames[i], "w") as jsonfile:
                json.dump(bkg_fit_params, jsonfile, indent=4)

            print("#############################")
            print("Order %i results: " % order)
            print("bkg fit chi2/nbins (fine binning) ", chi2_fine, ndof_fine, chi2_fine/ndof_fine, chi2_prob)
            print("#############################")

            chi2s[i] = chi2_fine
            ndofs[i] = ndof_fine
            probs[i] = chi2_prob
            fit_params[i] = bkg_fit_params
            fit_errs[i] = 0. # Deprecated
            fitter_bkg.delete()

        #F-test on this functional for to determine best num of parameters
        best_i = f_test(orderToTry, ndofs, chi2s, fit_errs, thresh = options.ftest_thresh, err_thresh = options.err_thresh)
        best_order = orderToTry[best_i]
        print("\n Chose order %i based on F-test ! \n" % best_order)

        #Add this functional form to workspace of the final fit
        shape_builder = shape_map[func_form]
        bkg_model,_,bkg_pars = shape_builder(func_form, card.poi, order=best_order)

        card.bkg_shapes.append(bkg_model)
        card.bkg_pars.extend(bkg_pars)
        func_name = func_form + "_order" + str(best_order)
        card.bkg_shape_names.append(func_name)


    print("#### Building multi pdf ### ")
    card.buildBkgShape()

    if(options.sig_norm_unc > 0):
        card.addSystematic("SigEff", "lnN", values = {"model_signal_m" : 1. + options.sig_norm_unc})

    print("making card")
    card.makeCard()
    #card.delete()

    cmd = (
        " cd {plot_dir} ; "
        + "text2workspace.py datacard_mass_{l2}.txt -o workspace_{l1}_{l2}.root; "
        + "combine -M FitDiagnostics workspace_{l1}_{l2}.root -m {mass} -n _{l1}_{l2} --robustFit 1 --cminDefaultMinimizerStrategy 0 --saveWorkspace; "
        + "combine -M AsymptoticLimits workspace_{l1}_{l2}.root -m {mass} -n lim_{l1}_{l2}; "
        ).format(plot_dir=plot_dir, mass=mass, l1=label, l2=fit_label)
    print(cmd)
    os.system(cmd)
    workspace_name = 'workspace_{l1}_{l2}.root'.format(l1=label, l2=fit_label)

    f_limit_name = (plot_dir + 'higgsCombinelim_{l1}_{l2}.'
                    + 'AsymptoticLimits.mH{mass:.0f}.root'
                    ).format(mass=mass, l1=label, l2=fit_label)
    f_diagnostics_name = (plot_dir + 'fitDiagnostics_{l1}_{l2}.root'
                   ).format(l1=label, l2=fit_label)


    f_diagnostics = ROOT.TFile(f_diagnostics_name, "READ")
    params_sb = f_diagnostics.Get("tree_fit_sb")
    params_sb.GetEntry(0)
    sig_strength = params_sb.r
    sig_strength_unc = params_sb.rErr
    nll_sb = params_sb.nll_min


    print("Signal strength is %.3f +/- %.3f" % (sig_strength, sig_strength_unc))

    params_b = f_diagnostics.Get("tree_fit_b")
    params_b.GetEntry(0)
    nll_b = params_b.nll_min

    delta_ll = nll_b - nll_sb
    print("dLL is %.2e" % delta_ll)

    delta_ll = max(delta_ll, 0.)
    signif = np.sqrt(2 * delta_ll)
    pval = 0.5 * math.erfc(signif / math.sqrt(2))

    print("Asymptotic significance (from dLL) is %.2f" % (signif))


    print('sig_norm %.3f' % sig_norm)



    true_sig_strength = get_sig_in_window(options.inputFile, binsx[0], binsx[-1]) /  sig_norm
    print("True sig strength %.3f" % true_sig_strength)

    #expected significance TODO
    exp_signif = 0.5
    exp_pval = 1.0


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


    #signal yield in +/- 2 sigma
    check_rough_sig(options.inputFile, options.mass*0.9, options.mass*1.1)


    f_limit.Close()
    f_diagnostics.Close()

    results = dict()

    # bkg fit results
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
    parser.add_option("-M", "-M", dest="mass", type=float, default=15.,
                      help="Signal mass hypothesis")
    parser.add_option("-i", "--inputFile", dest="inputFile",
                      default='fit_inputs/no_selection_03p.h5',
                      help="input h5 file")
    parser.add_option("-o", "--outDir", dest="outDir", default='plots/',
                      help="Where to put the output")
    parser.add_option("-s", "--sig_shape", default="sig_shape_M15.json",
                      help="Pre-saved signal shape")

    parser.add_option("--m-min", type=float, default=11.0,
                      help="Minimum m for the fit")
    parser.add_option("--m-max", type=float, default=20.0,
                      help="Maximum m for the fit")
    parser.add_option("--bin-size", type=float, default=0.2,
                      help="Size of bins")

    parser.add_option("--scale_j_unc", type=float, default=0.01,
                      help="Uncertainty on signal mean from JES")
    parser.add_option("--res_j_unc", type=float, default=0.035,
                      help="Uncertainty on signal width from JER")
    parser.add_option("--sig_norm", type=float, default=100.0,
                      help="Signal normalization (definition of mu=1)")
    parser.add_option("--ftest_thresh", type=float, default=0.05,
                      help="Threshold to prefer a function in the f-test")
    parser.add_option("--err_thresh", type=float, default=0.5,
                      help="Threshold on fit unc to be included in f-test")
    parser.add_option("--refit_sig", default=False, action="store_true",
                      help="""Fit the signal events (using truth labels)
                      to get signal shape""")
    parser.add_option("--rebin", default=False, action="store_true",
                      help="""Rebin  to make sure no bins less than 5 evts""")
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

