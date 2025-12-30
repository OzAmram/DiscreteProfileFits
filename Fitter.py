import ROOT
import json
from numpy import random
from array import array
import os,sys
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

def polyExpShape(name = 'model', poi = None, order=4, start_vals=None):
    #Sum of exponentials

    #initial values
    p_init = 0.0
    pmin = -10.0
    pmax = 10.0

    c_init = 0.2
    cmin = 0.
    cmax = 1.

    p_list = ROOT.RooArgList()
    par_names = []

    exp_par = ROOT.RooRealVar(f"polyExp_e", f"polyExp_e", p_init, pmin, pmax)
    exp = ROOT.RooExponential(f"polyExp_exp", f"polyExp_exp", poi, exp_par)

    for i in range(order):
        poly_par = ROOT.RooRealVar(f"polyExp_p{i}", f"polyExp_p{i}", p_init, pmin, pmax)
        p_list.add(poly_par)
        par_names.append(f"polyExp_p{i}")

    poly = ROOT.RooPolynomial("polyExp_poly", "polyExp_poly", poi, p_list)
    p_objs = [poly, exp]

    shape = ROOT.RooProdPdf(name+"_shape", name+"_shape", exp, poly)

    return shape, par_names, p_objs

def bernShape(name = 'model', poi = None, order=4, start_vals=None):
    #sum of bernstein polynomials

    #initial values
    if(start_vals is None):
        start_vals = [0.2, 1.5, 2.0, 2.0, 2.0]
    pmin = 0.0
    pmax = 20.0

    par_list = ROOT.RooArgList(poi)
    p_objs = []
    par_names = []
    for i in range(order):
        poly_par = ROOT.RooRealVar(f"bern_p{i}", f"bern_p{i}", start_vals[i], pmin, pmax)
        par_list.add(poly_par)
        par_names.append(f"bern_p{i}")
        p_objs.append(poly_par)

    shape = ROOT.RooBernstein(name+"_shape", name+"_shape", poi, par_list)

    return shape, par_names, p_objs


def expShape(name, poi = None, order=4, start_vals=None):
    #Sum of exponentials

    #initial values
    p_init = -0.5
    pmin = -10.0
    pmax = 10.0

    c_init = 0.2
    cmin = 0.
    cmax = 1.

    e_list = ROOT.RooArgList()
    c_list = ROOT.RooArgList()
    par_names = []
    objs = []

    for i in range(order):
        exp_par = ROOT.RooRealVar(f"exp_p{i}", f"exp_p{i}", p_init, pmin, pmax)
        exp = ROOT.RooExponential(f"exp_{i}", f"exp_{i}", poi, exp_par)
        coef = ROOT.RooRealVar(f"exp_c{i}", f"exp_c{i}", c_init, cmin, cmax)

        e_list.add(exp)
        par_names.append(f"exp_p{i}")

        if(i<(order-1)): # don't include last coeff because redundant
            c_list.add(coef)
            par_names.append(f"exp_c{i}")

        objs.extend([exp_par, exp, coef])

    recursiveFraction=True
    shape = ROOT.RooAddPdf(name+"_shape", name+"_shape", e_list, c_list, recursiveFraction)

    return shape, par_names, objs




class Fitter(object):
    def __init__(self,poi = ['x'], debug = False):
        self.cache_name = "cache%i.root"%(random.randint(0, 1e+6))
        print("Making cache %s "  % self.cache_name)
        self.cache=ROOT.TFile(self.cache_name,"RECREATE")
        self.cache.cd()
        self.debug = debug
        self.cleanedup = False
        self.objs = []
        self.par_names = []

        self.w=ROOT.RooWorkspace("w","w")
        self.dimensions = len(poi)
        self.poi=poi
        for v in poi:
            self.w.factory(v+"[1,161]")


    def __del__(self):
        if(not self.cleanedup): self.delete()

    def delete(self):
        os.system("rm %s" % self.cache_name)
        self.cleanedup = True

    def importBinnedData(self,histogram,poi = ["x"],name = "data", regions=[]):
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            var = self.w.var(p)
            cList.add(var)
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print ('Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead')
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            binningx =[]
            for i in range(1,bins+2):
                #v = mmin + i * (mmax-mmin)/float(N)
                binningx.append(axis.GetBinLowEdge(i))
            if(len(regions) == 0):
                var.setMax(maxi)
                var.setMin(mini)
            else:
                for reg_name,reg_low,reg_high in regions:
                    var.setRange(reg_name, reg_low, reg_high)
            if(self.debug): 
                print (" set binning "+str(binningx)) 
            var.setBinning(ROOT.RooBinning(len(binningx)-1,array("d",binningx)))
            #a = self.w.var(p).getBinning()
            #for b in range(0,a.numBins()+1):
                #print a.binLow(b)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)
        getattr(self.w,'import')(dataHist)

        self.objs.append(dataHist)

        self.w.Print()
        return dataHist

    def fetch(self,var):
        #import pdb; pdb.set_trace()
        self.w.var(var).Print()
        print("Fetching value " ,self.w.var(var).getVal()  )
        print("Fetching error " ,self.w.var(var).getError())
        return (self.w.var(var).getVal(), self.w.var(var).getError())

    def getFunc(self,model = "model"):
        return self.w.pdf(model)

    def getData(self,data = "data"):
        return self.w.data(data)

    def getVar(self,var = "mjj"):
        return self.w.var(var)

    def getW(self):
        return self.w

    def fit(self,model = "model",data="data",options=[]):
        fit_data = self.w.data(data)
        print(fit_data)
        print(options, options[0])
        self.w.Print()
        self.w.pdf(model).Print()

        if len(options)==0:
            fitresults = self.w.pdf(model).fitTo(fit_data)
        if len(options)==1:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0])	    
        if len(options)==2:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0],options[1])
        if len(options)==3:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0],options[1],options[2])
        if len(options)==4:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0],options[1],options[2],options[3])

        if fitresults:
            fitresults.Print() 
            f = ROOT.TFile.Open('fitresults.root','RECREATE')
            fitresults.SetName("fitresults")
            fitresults.Write()
            f.Close()	 
        return fitresults 

    def getLegend(self):
        self.legend = ROOT.TLegend(0.7510112,0.7183362,0.8502143,0.919833)
        self.legend.SetTextSize(0.032)
        self.legend.SetLineColor(0)
        self.legend.SetShadowColor(0)
        self.legend.SetLineStyle(1)
        self.legend.SetLineWidth(1)
        self.legend.SetFillColor(0)
        self.legend.SetFillStyle(0)
        self.legend.SetMargin(0.35)
        return self.legend

    def projection(self,model = "model",data="data",poi="x",filename="fit.root",binning=0,logy=False,xtitle='x',mass=1000):


        f = ROOT.TFile.Open("fitresults.root",'READ')

        linear_errors = False
        ndata = self.w.data(data).numEntries()
        self.frame=self.w.var(poi).frame(ROOT.RooFit.Bins(ndata))
        nbins = self.w.var(poi).getBinning().numBins()

        if f: 
            fr = f.Get('fitresults')
        else:
            fr = 0
            print("No fit result found (fitresults.root), plotting model only")

        if binning:
            self.w.data(data).plotOn(self.frame,ROOT.RooFit.Binning(binning),ROOT.RooFit.Invisible())
            if fr: 
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.VisualizeError(fr,1, linear_errors),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fr.GetName()))
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.LineColor(ROOT.kRed+1))	 
            self.w.data(data).plotOn(self.frame,ROOT.RooFit.Binning(binning))
        else: 
            self.w.data(data).plotOn(self.frame,
                                     ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),
                                     ROOT.RooFit.Binning(nbins),
                                     ROOT.RooFit.Invisible()
                                    )
            if fr: 
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.VisualizeError(fr,1, linear_errors),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fr.GetName()))
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.LineColor(ROOT.kRed+1))

            self.w.data(data).plotOn(self.frame,
                                     ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),
                                     ROOT.RooFit.Binning(nbins)
                                    )
        self.legend = self.getLegend()
        self.legend.AddEntry( self.w.pdf(model)," Full PDF","l")

        self.c=ROOT.TCanvas("c","c")
        if logy:
            self.frame.SetMinimum(1)
            self.frame.SetMaximum(1e+7)
            self.c.SetLogy()
        self.c.cd()
        self.frame.Draw()
        self.frame.GetYaxis().SetTitle('events')
        self.frame.GetXaxis().SetTitle(xtitle)
        self.frame.SetTitle('')
        self.c.Draw()

        self.legend.Draw("same")	    
        self.c.SaveAs(filename)

        chi2 = ROOT.RooChi2Var("chi2", "chi2", self.w.pdf(model), self.w.data(data))
        chi2_val = chi2.getVal()

        pullDist = self.frame.pullHist()
        nfloat = self.w.pdf(model).getParameters(self.w.data(data)).selectByAttrib("Constant", False).getSize()

        ndof = nbins - nfloat

        print(auto_chi2, auto_chi2_v2, chi2_val, ndof, chi2_val/ndof)

        return chi2_val, ndof


    def signalResonance(self,name = 'model',poi="MVV",mass=0):

        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        self.w.factory("MH[1000]")
        self.w.factory("mean[%.1f,%.1f,%.1f]"%(mass,0.8*mass,1.2*mass))
        self.w.factory("sigma[%.1f,%.1f,%.1f]"%(mass*0.02,mass*0.005,mass*0.10))
        self.w.factory("alpha[0.85,0.60,1.20]")
        self.w.factory("sign[6,0.1,150]")
        self.w.factory("scalesigma[2.0,1.2,3.6]")
        gsigma = ROOT.RooFormulaVar("gsigma","@0*@1", ROOT.RooArgList(self.w.var("sigma"),self.w.var("scalesigma")))
        getattr(self.w,'import')(gsigma)
        self.w.factory("Gaussian::gauss(%s,mean,gsigma)"%poi)
        self.w.factory("CBShape::cb(%s,mean,sigma,alpha,sign)"%poi)
        self.w.factory('SUM::'+name+'(sigfrac[0.0,0.0,0.850]*gauss,cb)')
    
    def signalDCB(self, name='model', poi="MVV", mass=0):
    
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        
        self.w.factory("MH[1000]")
        self.w.factory("mean[%.1f,%.1f,%.1f]"%(mass,0.8*mass,1.2*mass))
        self.w.factory("sigma[%.1f,%.1f,%.1f]"%(mass*0.02,mass*0.005,mass*0.10))
        self.w.factory("alpha[1.2,0.0,18]")
        self.w.factory("alpha2[1.2,0.0,10]")
        self.w.factory("sign[5,0,600]")
        self.w.factory("sign2[5,0,50]")
        self.w.factory("DoubleCB::"+name+"(%s,mean,sigma,alpha,sign,alpha2,sign2)"%poi)




    def bkgShape(self, func_form='bern', name = 'model',poi="m",start_vals = None,order=2):

        if(type(poi) == str): poi = self.w.var(poi)


        if(func_form == 'bern'):
            shape, par_names, objs = bernShape(name = name, poi=poi, order=order, start_vals=start_vals)
        elif(func_form == 'polyExp'):
            shape, par_names, objs = polyExpShape(name = name, poi=poi, order=order, start_vals=start_vals)
        elif(func_form == 'exp'):
            shape, par_names, objs = expShape(name = name, poi=poi, order=order, start_vals=start_vals)
        else:
            print("Shape %s not implemented!" % shape)
            exit(1)
        
        self.objs.extend(objs)
        self.par_names.extend(par_names)

        norm = 10000.
        norm_name = name + "_norm"
        norm_var = ROOT.RooRealVar(norm_name, norm_name, norm, 0., 1e8)

        self.w.Print()
        model = ROOT.RooExtendPdf(name, name, shape, norm_var)

        print("bkg shape")
        self.w.Print()
        getattr(self.w,'import')(model)
        self.objs.append((shape, norm_var, model))

        return model


    def addDCBSignalShape(self, name, poi, jsonFile, scale={},
                          resolution={}):

        
        pdfName="_".join([name,])
        
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.items():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.items():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
            
        if(type(poi) == str): 
            poi = self.w.var(poi)
    
        f = ROOT.TFile(jsonFile,'READ')
        meanG = f.Get('mean')
        sigmaG = f.Get('sigma')
        alpha_oneG = f.Get('alpha')
        sign_oneG = f.Get('sign')
        alpha_twoG = f.Get('alpha2')
        sign_twoG = f.Get('sign2')

        x = ROOT.Double(0.)
        mean = ROOT.Double(0.)
        meanG.GetPoint(0,x,mean)
        sigma = ROOT.Double(0.)
        sigmaG.GetPoint(0,x,sigma)
        alpha_one = ROOT.Double(0.)
        alpha_oneG.GetPoint(0, x, alpha_one)
        alpha_two = ROOT.Double(0.)
        alpha_twoG.GetPoint(0, x, alpha_two)
        sign_one = ROOT.Double(0.)
        sign_oneG.GetPoint(0,x,sign_one)
        sign_two = ROOT.Double(0.)
        sign_twoG.GetPoint(0,x,sign_two)

        meanVar = "_".join(["MEAN", name, ])
        self.w.factory(
            "expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(
                name=meanVar, param=mean, vv_syst=scaleStr,
                vv_systs=','.join(scaleSysts)))

        sigmaVar = "_".join(["SIGMA", name, ])
        self.w.factory(
            "expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(
                name=sigmaVar, param=sigma, vv_syst=resolutionStr,
                vv_systs=','.join(resolutionSysts)))

        alphaOneVar = "_".join(["ALPHAONE", name, ])        
        alpha_one = ROOT.RooRealVar(alphaOneVar,alphaOneVar,alpha_one)
        getattr(self.w, 'import')(alpha_one)

        alphaTwoVar = "_".join(["ALPHATWO", name, ])        
        alpha_two = ROOT.RooRealVar(alphaTwoVar, alphaTwoVar, alpha_two)
        getattr(self.w, 'import')(alpha_two)
        
        signOneVar = "_".join(["SIGNONE", name, ])
        sign_one = ROOT.RooRealVar(signOneVar, signOneVar, sign_one)    
        getattr(self.w, 'import')(sign_one)
        
        signTwoVar = "_".join(["SIGNTWO", name, ])
        sign_two = ROOT.RooRealVar(signTwoVar, signTwoVar, sign_two)    
        getattr(self.w, 'import')(sign_two)
        
        #dcbFunc = "_".join(["dcb", name, ])
        #import pdb; pdb.set_trace()
        dcb = ROOT.RooDoubleCB(pdfName, pdfName, poi,
                                      self.w.function(meanVar),
                                      self.w.function(sigmaVar), alpha_one,
                                      sign_one, alpha_two, sign_two)

        getattr(self.w, 'import')(dcb)

        self.objs.append(dcb)

        return dcb, [alpha_one, alpha_two, sign_one, sign_two]


