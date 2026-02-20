import argparse
import json
import os
from Utils import *
from Fitter import Fitter

def fit_signalmodel(input_file, sig_file_name, mass, x_bins,
                    plot_dir, dcb_model=False, fit_range = 0.2, plot_label =""):

    mlow = (1.0 - fit_range) * mass
    mhigh = (1.0 + fit_range) * mass

    bins_sig_fit = array(
        'f', truncate( x_bins, mlow, mhigh))


    histos_sig = ROOT.TH1F("m_sig", "m_sig",
                           len(bins_sig_fit) - 1, bins_sig_fit)

    load_h5_sb(input_file, histos_sig)
    histos_sig.Print("range")

    fitter = Fitter(['m_fine'])

    fitter.signalShape('model_s', "m_fine", mass, options.sig_shape)

    fitter.w.var("MH").setVal(mass)
    fitter.importBinnedData(histos_sig, 'm_fine', 'data')
    fres = fitter.fit('model_s', 'data', [ROOT.RooFit.Save(1)])

    if fres:
        fname = sig_file_name.replace('.json', '.txt')
        with open(fname, 'w') as f:
            f.write('%d\n' % fres.status())
    m_fine = fitter.getVar('m_fine')
    m_fine.setBins(len(bins_sig_fit))

    chi2, ndof = fitter.projection("model_s", "data", "m_fine",
                                  plot_dir + plot_label + "signal_fit.png")

    fitter.projection("model_s", "data", "m_fine",
                      plot_dir + plot_label +  "signal_fit_log.png", logy=True)

    print("Fit done")

    print (
        """
        #############################
        signal fit chi2 / ndof: %.3f / %.3f =   %.3f
        """ % (chi2, ndof, chi2/ndof)
        )

    params = {}
    for par in fitter.sig_shape_params:
        value, error = fitter.fetch(par)
        params[par] = value
        params[par + "-err"] = error

    params['mass'] = mass
    print("Results:")
    print(params)

    with open(sig_file_name, 'w') as f:
        json.dump(params, f, indent=4)


    return fitter




def fit_signals(options):

    out_dir = os.path.abspath(options.outDir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    masses = options.masses

    binsx = list(np.arange(options.m_min, options.m_max + options.bin_size, options.bin_size))

    # round to smallest precision we are storing mass values with, otherwise
    # get weird effects related to bin size
    #roundTo(binsx, fine_bin_size)

    for i, mass in enumerate(masses):
        print("########## FIT SIGNAL %.2f AND SAVE PARAMETERS ############" % mass)
        sig_file_name = os.path.join(out_dir, "sig_fit_{}.json".format(mass))
        plot_label = "M%i_" % mass
        current_fit = fit_signalmodel(options.inputFiles[i], sig_file_name,
                                      mass, binsx, out_dir + "/",
                                      dcb_model=options.dcbModel, 
                                      fit_range = options.fitRange, plot_label = plot_label)

    print("Done with loop!")
    return


def fitting_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("--m-min", type=float, default=11.0,
                      help="Minimum m for the fit")
    parser.add_argument("--m-max", type=float, default=40.0,
                      help="Maximum m for the fit")
    parser.add_argument("--bin-size", type=float, default=0.2,
                      help="Size of bins")
    parser.add_argument("--sig-shape", default='DCB', 
            help="Functional form for signal shape. Options: gaus, CB, CBgaus, DCB")

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
    parser.add_argument("--fitRange", dest="fitRange", type = float, default=0.2,
                        help="What mass range to perform fit to signal shape over (in terms of frac. of signal mass)")
    return parser


if __name__ == "__main__":
    parser = fitting_options()
    options = parser.parse_args()
    fit_signals(options)
