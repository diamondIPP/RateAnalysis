from ROOT import TMath, TF1, gROOT, TH1F, gStyle
from array import array

#-----------------------------------------------------------------------
#
# Convoluted Landau and Gaussian Fitting Function
#         (using ROOT's Landau and Gauss functions)
#
#  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
#  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
#   Markus Friedl (Markus.Friedl@cern.ch)
#
#  Adapted for pyROOT by M.Seeli
#
#-----------------------------------------------------------------------

def langaufun(x, par):

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

def langaufit(his, fitrange, startvalues, parlimitslo, parlimitshi, fitparams, fiterrors, ChiSqr, NDF):
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
    #   fitparams[4]    returns the final fit parameters
    #   fiterrors[4]    returns the final fit errors
    #   ChiSqr          returns the chi square
    #   NDF             returns ndf
    
    FunName = "Fitfcn_%s" % his.GetName()

    ffitold = gROOT.GetListOfFunctions().FindObject(FunName)
    if ffitold: del ffitold

    ffit = TF1(FunName,langaufun,fitrange[0],fitrange[1],4)
    ffit.SetParameters(startvalues)
    ffit.SetParNames("Width","MP","Area","GSigma")
    
    i = 0
    while i<4: 
        ffit.SetParLimits(i, parlimitslo[i], parlimitshi[i])
        i += 1
    
    his.Fit(FunName,"RB0")   # fit within specified range, use ParLimits, do not plot

    ffit.GetParameters(fitparams)    # obtain fit parameters
    i = 0
    while i<4:
        fiterrors[i] = ffit.GetParError(i)     # obtain fit parameter errors
        i += 1

    ChiSqr[0] = ffit.GetChisquare()  # obtain chi^2
    NDF[0] = ffit.GetNDF()           # obtain ndf

    return ffit              # return fit function

def langaupro(params, maxx, FWHM):

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
        l = langaufun(x,params)

        if (l < lold):
            step = -step/10

        p += step

    if (i == MAXCALLS):
        return (-1)
    
    maxx = x[0]
    
    fy = l/2
    
    
    #  Search for right x location of fy
    
    p = maxx + params[0]
    step = params[0]
    lold = -2.0
    l    = -1e300
    i    = 0
    
    
    while ( (l != lold) and (i < MAXCALLS) ):
        i += 1

        lold = l
        x[0] = p + step
        l = TMath.Abs(langaufun(x,params) - fy)

        if (l > lold):
            step = -step/10

        p += step

    
    if (i == MAXCALLS):
        return (-2)
    
    fxr = x[0]
    
    
    #  Search for left x location of fy
    
    p = maxx - 0.5 * params[0]
    step = -params[0]
    lold = -2.0
    l    = -1e300
    i    = 0
    
    while ( (l != lold) and (i < MAXCALLS) ):
        i += 1

        lold = l
        x[0] = p + step
        l = TMath.Abs(langaufun(x,params) - fy)

        if (l > lold):
            step = -step/10

        p += step

    
    if (i == MAXCALLS):
        return (-3)
    
    
    fxl = x[0]
    
    FWHM = fxr - fxl
    return (0)

if __name__ == "__main__":
    # Fill Histogram
    data = [0,0,0,0,0,0,2,6,11,18,18,55,90,141,255,323,454,563,681,
            737,821,796,832,720,637,558,519,460,357,291,279,241,212,
            153,164,139,106,95,91,76,80,80,59,58,51,30,49,23,35,28,
            23,22,27,27,24,20,16,17,14,20,12,12,13,10,17,7,6,12,6,
            12,4,9,9,10,3,4,5,2,4,1,5,5,1,7,1,6,3,3,3,4,5,4,4,2,2,7,2,4]
    hSNR = TH1F("snr","Signal-to-noise",400,0,400)

    for i in xrange(98):
        hSNR.Fill(i,data[i])

    # Fitting SNR histo
    print "Fitting...\n"

    # Setting fit range and start values
    fr = array("d", [0,0])
    sv = array("d", [0,0,0,0])
    pllo = array("d", [0,0,0,0])
    plhi = array("d", [0,0,0,0])
    fp = array("d", [0,0,0,0])
    fpe = array("d", [0,0,0,0])

    # fitrange for fit:
    fr[0]=0.3*hSNR.GetMean()
    fr[1]=3.0*hSNR.GetMean()

    # lower parameter limits
    pllo[0]=0.5 # par[0]=Width (scale) parameter of Landau density
    pllo[1]=5.0 # par[1]=Most Probable (MP, location) parameter of Landau density
    pllo[2]=1.0 # par[2]=Total area (integral -inf to inf, normalization constant)
    pllo[3]=0.4 # par[3]=Width (sigma) of convoluted Gaussian function

    # upper parameter limits
    plhi[0]=5.0         # par[0]=Width (scale) parameter of Landau density
    plhi[1]=50.0        # par[1]=Most Probable (MP, location) parameter of Landau density
    plhi[2]=1000000.0   # par[2]=Total area (integral -inf to inf, normalization constant)
    plhi[3]=5.0         # par[3]=Width (sigma) of convoluted Gaussian function

    # Startvalues for fit:
    sv[0]=1.8       # par[0]=Width (scale) parameter of Landau density
    sv[1]=20.0      # par[1]=Most Probable (MP, location) parameter of Landau density
    sv[2]=50000.0   # par[2]=Total area (integral -inf to inf, normalization constant)
    sv[3]=3.0       # par[3]=Width (sigma) of convoluted Gaussian function

    chisqr = array("d", [0]) # returns the chi square
    ndf = array("d", [0]) # returns ndf

    fitsnr = langaufit(hSNR,fr,sv,pllo,plhi,fp,fpe,chisqr,ndf)

    SNRPeak = array("d", [0])
    SNRFWHM = array("d", [0])
    langaupro(fp,SNRPeak,SNRFWHM)

    print "Fitting done\nPlotting results...\n"

    # Global style settings
    gStyle.SetOptStat(1111)
    gStyle.SetOptFit(111)
    gStyle.SetLabelSize(0.03,"x")
    gStyle.SetLabelSize(0.03,"y")

    hSNR.GetXaxis().SetRange(0,70)
    hSNR.Draw()
    fitsnr.Draw("lsame")
    raw_input("wait")
