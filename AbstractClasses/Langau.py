from AbstractClasses.Elementary import Elementary
from ROOT import TMath, TF1, gROOT
from array import array


class Langau(Elementary):

    def __init__(self, verbose=False):
        Elementary.__init__(self, verbose=verbose)
        self.FitResults = {
            "FitFunction": None,
            "Peak": None,
            "FWHM": None,
            "Chi2": None,
            "NDF": None,
            "FitParameters": array('d', [0]*4),
            "FitParamErrors": array('d', [0]*4)
        }
    
    def __del__(self):
        ffitold = gROOT.GetListOfFunctions().FindObject(self.FunName)
        if ffitold: del ffitold
    
    def LoadConfig(self):
        self.ShowAndWait = False
    
    @staticmethod
    def LangauFun(x, par):
        # Fit parameters:
        # par[0]=Width (scale) parameter of Landau density
        # par[1]=Most Probable (MP, location) parameter of Landau density
        # par[2]=Total area (integral -inf to inf, normalization constant)
        # par[3]=Width (sigma) of convoluted Gaussian function
        #
        # In the Landau distribution (represented by the CERNLIB approximation),
        # the maximum is located at x=-0.22278298 with the location parameter=0.
        # This shift is corrected within this function, so that the actual
        # maximum is identical to the MP parameter.
    
        # Numeric constants
        invsq2pi = 0.3989422804014  # (2 pi)^(-1/2)
        mpshift  = -0.22278298   # Landau maximum location
    
        # Control constants
        np = 100.0  # number of convolution steps
        sc =   5.0  # convolution extends to +-sc Gaussian sigmas
    
        # MP shift correction
        mpc = par[1] - mpshift * par[0]
    
        # Range of convolution integral
        xlow = x[0] - sc * par[3]
        xupp = x[0] + sc * par[3]
    
        step = (xupp-xlow) / np
    
        # Variables
        sum = 0.0
        i=1.0
        
        # Convolution integral of Landau and Gaussian by sum
        while i<=np/2:
            xx = xlow + (i-.5) * step
            fland = TMath.Landau(xx,mpc,par[0]) / par[0]
            sum += fland * TMath.Gaus(x[0],xx,par[3])
    
            xx = xupp - (i-.5) * step
            fland = TMath.Landau(xx,mpc,par[0]) / par[0]
            sum += fland * TMath.Gaus(x[0],xx,par[3])
    
            i += 1
    
        return (par[2] * step * sum * invsq2pi / par[3])
    
    def LangauFit(self, his, fitrange=None, startvalues=None, parlimitslo=None, parlimitshi=None):
        # Once again, here are the Landau * Gaussian parameters:
        #   par[0]=Width (scale) parameter of Landau density
        #   par[1]=Most Probable (MP, location) parameter of Landau density
        #   par[2]=Total area (integral -inf to inf, normalization constant)
        #   par[3]=Width (sigma) of convoluted Gaussian function
        #
        # Variables for langaufit call:
        #   his             histogram to fit
        #   fitrange[2]     lo and hi boundaries of fit range
        #   startvalues[4]  reasonable start values for the fit
        #   parlimitslo[4]  lower parameter limits
        #   parlimitshi[4]  upper parameter limits

        # ------------------------------------------------
        # Setting fit range and start values

        # fitrange for fit:
        if fitrange is None:
            fitrange = array("d", [0,0])
            fitrange[0]=0.3*his.GetMean()
            fitrange[1]=3.0*his.GetMean()
        else:
            assert(len(fitrange)==2), "Wrong Dimensions of fitrange"

        # lower parameter limits
        if parlimitslo is None:
            parlimitslo = array("d", [0,0,0,0])
            parlimitslo[0]=2.0 # par[0]=Width (scale) parameter of Landau density
            parlimitslo[1]=20.0 # par[1]=Most Probable (MP, location) parameter of Landau density
            parlimitslo[2]=1.0 # par[2]=Total area (integral -inf to inf, normalization constant)
            parlimitslo[3]=1.0 # par[3]=Width (sigma) of convoluted Gaussian function
        else:
            assert(len(fitrange)==4), "Wrong Dimensions of parlimitslo"

        # upper parameter limits
        if parlimitshi is None:
            parlimitshi = array("d", [0,0,0,0])
            parlimitshi[0]=50.0         # par[0]=Width (scale) parameter of Landau density
            parlimitshi[1]=250.0        # par[1]=Most Probable (MP, location) parameter of Landau density
            parlimitshi[2]=(his.GetXaxis().GetLast()-his.GetXaxis().GetFirst())*his.GetMaximum()   # par[2]=Total area (integral -inf to inf, normalization constant)
            parlimitshi[3]=50.0         # par[3]=Width (sigma) of convoluted Gaussian function
        else:
            assert(len(fitrange)==4), "Wrong Dimensions of parlimitshi"

        # Startvalues for fit:
        if startvalues is None:
            startvalues = array("d", [0,0,0,0])
            startvalues[0]=11.0       # par[0]=Width (scale) parameter of Landau density
            startvalues[1]=his.GetBinCenter(his.GetMaximumBin())      # par[1]=Most Probable (MP, location) parameter of Landau density
            startvalues[2]=40.0       # par[2]=Total area (integral -inf to inf, normalization constant)
            startvalues[3]=11.0       # par[3]=Width (sigma) of convoluted Gaussian function
        else:
            assert(len(fitrange)==4), "Wrong Dimensions of startvalues"
        
        self.FunName = "Fitfcn_%s" % his.GetName()
    
        ffitold = gROOT.GetListOfFunctions().FindObject(self.FunName)
        if ffitold: del ffitold
    
        ffit = TF1(self.FunName,self.LangauFun,fitrange[0],fitrange[1],4)
        ffit.SetParameters(startvalues)
        ffit.SetParNames("Width","MP","Area","GSigma")
        
        i = 0
        while i<4: 
            ffit.SetParLimits(i, parlimitslo[i], parlimitshi[i])
            i += 1
        
        his.Fit(self.FunName,"RB0")   # fit within specified range, use ParLimits, do not plot
        
        ffit.GetParameters(self.FitResults["FitParameters"])    # obtain fit parameters
        i = 0
        while i<4:
            self.FitResults["FitParamErrors"][i] = ffit.GetParError(i)     # obtain fit parameter errors
            i += 1
    
        self.FitResults["Chi2"] = ffit.GetChisquare()  # obtain chi^2
        self.FitResults["NDF"] = ffit.GetNDF()           # obtain ndf
        self.FitResults["FitFunction"] = ffit           # Fit function
        
        # Get Peak and FWHM: 
        self.LangauProperties(self.FitResults["FitParameters"])

        print "MPV from fit: {0:0.3f} +- {1:0.03f}".format(self.FitResults["FitParameters"][1], self.FitResults["FitParamErrors"][1])
        print "MPV from func: ", self.FitResults["Peak"]
    
        return ffit              # return fit function
    
    def LangauProperties(self, params):
    
       #  Seaches for the location (x value) at the maximum of the
       #  Landau-Gaussian convolute and its full width at half-maximum.
       # 
       #  The search is probably not very efficient, but it's a first try.
    
        i = 0
        MAXCALLS = 10000
        
        #  Search for maximum
        
        p = params[1] - 0.1 * params[0]
        step = 0.05 * params[0]
        lold = -2.0
        l    = -1.0
    
        x = array("d", [0])
    
        while ( (l != lold) and (i < MAXCALLS) ):
            i += 1
            lold = l
            x[0] = p + step
            l = self.LangauFun(x,params)
            if (l < lold):
                step = -step/10
            p += step
    
        if (i == MAXCALLS):
            return (-1)
        
        self.FitResults["Peak"] = x[0]
        
        fy = l/2
        
        
        #  Search for right x location of fy
        
        p = self.FitResults["Peak"] + params[0]
        step = params[0]
        lold = -2.0
        l    = -1e300
        i    = 0
        
        
        while ( (l != lold) and (i < MAXCALLS) ):
            i += 1
            lold = l
            x[0] = p + step
            l = TMath.Abs(self.LangauFun(x,params) - fy)
            if (l > lold):
                step = -step/10
            p += step
    
        
        if (i == MAXCALLS):
            return (-2)
        fxr = x[0]
        
        
        #  Search for left x location of fy
        
        p = self.FitResults["Peak"] - 0.5 * params[0]
        step = -params[0]
        lold = -2.0
        l    = -1e300
        i    = 0
        
        while ( (l != lold) and (i < MAXCALLS) ):
            i += 1
            lold = l
            x[0] = p + step
            l = TMath.Abs(self.LangauFun(x,params) - fy)
            if (l > lold):
                step = -step/10
            p += step
        
        if (i == MAXCALLS):
            return (-3)
        
        fxl = x[0]
        
        self.FitResults["FWHM"] = fxr - fxl
        return (0)
    
    def GetMPV(self): # difference to "Peak"?
        return self.FitResults["FitParameters"][1]

    def GetMPVError(self):
        return self.FitResults["FitParamErrors"][1]

    def GetFWHM(self):
        assert(self.FitResults["FWHM"] is not None), "No fit made. Execute self.LangauFit() before calling self.GetFWHM()"
        return self.FitResults["FWHM"]
    
    