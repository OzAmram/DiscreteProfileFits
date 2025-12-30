import argparse
import json
import os
from Utils import *
from Fitter import Fitter

def fit_signalmodel(input_file, sig_file_name, mass, x_bins, fine_bins,
                    plot_dir, return_fit=False, dcb_model=False, fit_range = 0.2, plot_label =""):

    fine_bin_size = 0.2
    bins_fine = int((x_bins[-1]-x_bins[0])/fine_bin_size)

    mlow = (1.0 - fit_range) * mass
    mhigh = (1.0 + fit_range) * mass

    bins_sig_fit = array(
        'f', truncate(
             [x_bins[0] + ib*fine_bin_size for ib in range(bins_fine + 1)],
             mlow, mhigh)
        )

    large_bins_sig_fit = array('f', truncate(x_bins, mlow, mhigh))
    roobins_sig_fit = ROOT.RooBinning(len(large_bins_sig_fit) - 1,
                                      array('d', large_bins_sig_fit),
                                      "mbins_sig")

    histos_sig = ROOT.TH1F("m_sig", "m_sig",
                           len(bins_sig_fit) - 1, bins_sig_fit)

    load_h5_sb(input_file, histos_sig)
    sig_outfile = ROOT.TFile(sig_file_name, "RECREATE")
    fitter = Fitter(['m_fine'])

    if dcb_model:
        fitter.signalDCB('model_s', "m_fine", mass)
    else:
        fitter.signalResonance('model_s', "m_fine", mass)

    fitter.w.var("MH").setVal(mass)
    fitter.importBinnedData(histos_sig, ['m_fine'], 'data')
    fres = fitter.fit('model_s', 'data', [ROOT.RooFit.Save(1)])
    if fres:
        fname = sig_file_name.replace('.root', '.txt')
        with open(fname, 'w') as f:
            f.write('%d\n' % fres.status())
    m_fine = fitter.getVar('m_fine')
    m_fine.setBins(len(bins_sig_fit))

    chi2_fine, ndof = fitter.projection("model_s", "data", "m_fine",
                                  plot_dir + plot_label + "signal_fit.png")

    fitter.projection("model_s", "data", "m_fine",
                      plot_dir + plot_label +  "signal_fit_log.png", 0, True)

    chi2, ndof = fitter.projection("model_s", "data", "m_fine",
                             plot_dir + plot_label + "signal_fit_binned.png",
                             roobins_sig_fit)

    fitter.projection("model_s", "data", "m_fine",
                      plot_dir + plot_label + "signal_fit_log_binned.png",
                      roobins_sig_fit, logy=True)

    print("Fit done")
    sig_outfile.cd()
    histos_sig.Write()
    print("cd, write")

    if dcb_model:
        graphs = {'mean': ROOT.TGraphErrors(),
                  'sigma': ROOT.TGraphErrors(),
                  'alpha': ROOT.TGraphErrors(),
                  'alpha2': ROOT.TGraphErrors(),
                  'sign': ROOT.TGraphErrors(),
                  'sign2': ROOT.TGraphErrors()}
    else:
        graphs = {'mean': ROOT.TGraphErrors(),
                  'sigma': ROOT.TGraphErrors(),
                  'alpha': ROOT.TGraphErrors(),
                  'sign': ROOT.TGraphErrors(),
                  'scalesigma': ROOT.TGraphErrors(),
                  'sigfrac': ROOT.TGraphErrors()}

    for var, graph in graphs.items():
        value, error = fitter.fetch(var)
        graph.SetPoint(0, mass, value)
        graph.SetPointError(0, 0.0, error)

    sig_outfile.cd()
    for name, graph in graphs.items():
        graph.Write(name)
        graph.Print()

    sig_outfile.Close()

    print (
        """
        #############################
        signal fit chi2 (fine binning), %.3f
        signal fit chi2 (large binning), %.3f
        #############################
        """ % (chi2_fine, chi2)
        )

    if return_fit:
        return fitter
    else:
        return None




def fit_signals(options):

    out_dir = os.path.abspath(options.outDir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fine_bin_size = 0.2
    masses = options.masses

    binsx = list(np.arange(11, 40, 0.2))

    # round to smallest precision we are storing mass values with, otherwise
    # get weird effects related to bin size
    #roundTo(binsx, fine_bin_size)
    bins_fine = int((binsx[-1]-binsx[0])/fine_bin_size)

    if options.dcbModel:
        full_graphs = {'mean': ROOT.TGraphErrors(),
                       'sigma': ROOT.TGraphErrors(),
                       'alpha': ROOT.TGraphErrors(),
                       'alpha2': ROOT.TGraphErrors(),
                       'sign': ROOT.TGraphErrors(),
                       'sign2': ROOT.TGraphErrors()}
    else:
        full_graphs = {'mean': ROOT.TGraphErrors(),
                       'sigma': ROOT.TGraphErrors(),
                       'alpha': ROOT.TGraphErrors(),
                       'sign': ROOT.TGraphErrors(),
                       'scalesigma': ROOT.TGraphErrors(),
                       'sigfrac': ROOT.TGraphErrors()}

    for i, mass in enumerate(masses):
        print("########## FIT SIGNAL %.2f AND SAVE PARAMETERS ############" % mass)
        sig_file_name = os.path.join(out_dir, "sig_fit_{}.root".format(mass))
        plot_label = "M%i_" % mass
        current_fit = fit_signalmodel(options.inputFiles[i], sig_file_name,
                                      mass, binsx, bins_fine, out_dir + "/",
                                      return_fit=True, dcb_model=options.dcbModel, 
                                      fit_range = options.fitRange, plot_label = plot_label)

        print("return")
        parameters = dict()
        print(i, mass)
        for var, graph in full_graphs.items():
            val, err = current_fit.fetch(var)
            parameters[var] = val
            parameters['%s-err' % var] = err
            graph.SetPoint(i, mass, val)
            graph.SetPointError(i, 0.0, err)
            graph.Print()

        #parameters["script_options"] = vars(options)
        #sig_file_name = os.path.join(out_dir, "sig_fit_{}.json".format(mass))
        #with open(sig_file_name, 'w') as f:
        #    json.dump(parameters, f, indent=4)

    print("Done with loop!")

    full_file_name = os.path.join(out_dir, "full_fit.root")
    full_outfile = ROOT.TFile(full_file_name, "RECREATE")
    full_outfile.cd()

    for name, graph in full_graphs.items():
        graph.Write(name)

    full_outfile.Close()
    print("Done!")
    return


def fitting_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-M", "-M", dest="masses", type=int, nargs="+",
                        help="List of injected signal masses")
    parser.add_argument("-i", "--inputFiles", dest="inputFiles", nargs="+",
                        help="""List of input h5 files.
                        Must have same length as mass list""")
    parser.add_argument("-o", "--outDir", dest="outDir", default='plots/',
                        help="Where to store output files")
    parser.add_argument("--dcb-model", "--dcbModel", dest="dcbModel", action="store_true",
                        default=False,
                        help="Whether or not to use double crystal ball model")
    parser.add_argument("--fitRange", dest="fitRange", type = float, default=0.1,
                        help="What mass range to perform fit to signal shape over (in terms of frac. of signal mass)")
    return parser


if __name__ == "__main__":
    parser = fitting_options()
    options = parser.parse_args()
    fit_signals(options)
