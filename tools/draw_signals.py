import ROOT
import json
from numpy import random
from array import array
import os,sys,commands
from Fitter import *
from DataCardMaker import *
ROOT.gROOT.SetBatch(True)

f_shape = ""


f_sig1 = "XToYY_shapes/sig_fit_3000.root"
#f_sig1 = "XToYY_shapes/sig_fit_5000.root"
#f_sig2 = "XToYY_TNT_shapes/sig_fit_3000.root"
#f_sig3 = "XToYY_TNT_lowsig_shapes/sig_fit_3000.root"
f_sig2 = 'old/interp_test/graviton_interpolation_M3000.0.root'
#f_sig2 = 'interp_test2/case_interpolation_M5000.0.root'

#fs = [f_sig1, f_sig2]
d = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/TagNTrain/plots/sig_shape_uncs_XYY_M3000/"
fs = [d + 'sig_shape_nom/sig_fit_3000.root', d + 'JES_up/sig_fit_3000.root', d + 'JES_down/sig_fit_3000.root']
labels = ["Nominal", "JES Up", "JES Down"]
colors = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue, ROOT.kGreen]

fout = "XYY_sig_shape_variation_JES.png"
xlow = 2500
xhigh = 3500

#labels = ["No Cut", "TNT Cut (5 #sigma inj.)", "TNT Cut (0 #sigma inj.)"]


var = ROOT.RooRealVar("mjj", "mjj", xlow, xhigh)
frame = var.frame()
#dcb1 = get_signal_shape(f_sig)
card = DataCardMaker('sig_test')
legend = ROOT.TLegend(0.2, 0.2)
frame.SetTitle("")
ls = []

for i in range(len(fs)):
    card.addDCBSignalShape('sig%i' %i, 'mjj', fs[i], {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})
    sig = card.w.pdf('sig%i_JJ_sig_test' %i )
    sig.plotOn(frame, ROOT.RooFit.LineColor(colors[i]), ROOT.RooFit.Name(labels[i]) )
    sig.Print()
    l = ROOT.TLine(0,0,0,0)
    l.SetLineWidth(2)
    l.SetLineColor(colors[i])
    ls.append(l)
    legend.AddEntry(l, labels[i],  "l") 

c1 =ROOT.TCanvas("c1","",800,800)
c1.cd()
frame.Draw()
legend.Draw()
c1.Print(fout)

