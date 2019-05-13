# -----------------------------------------------------------------------
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
# -----------------------------------------------------------------------
from ROOT import TMath, TF1, gROOT, TH1F, gStyle
from array import array
from math import pi, sqrt
from numpy import array, zeros


class Langau:
    def __init__(self, histo, nconv=100, fit_range=None):
        self.Histo = histo
        self.Max = histo.GetMaximum()
        self.NConvolutions = nconv
        self.NSigma = 5.
        self.FitRange = [k * histo.GetMean() for k in [.1, 3]] if fit_range is None else fit_range
        self.ParLimits = self.init_par_limits()
        self.StartValues = self.init_start_values()

        self.Fit = None
        self.Parameters = zeros(4)
        self.ParErrors = zeros(4)
        self.Chi2 = None
        self.NDF = None

    def init_par_limits(self):
        # todo: set limits based on histogram
        h = self.Histo
        limits = [(0, 0)] * 4
        limits[0] = [self.estimate_sigma() / 5 * k for k in [0, 3]]    # par[0]=Width (scale) parameter of Landau density
        limits[1] = [self.get_x_max() * k for k in [0.5, 1.5]]          # par[1]=Most Probable (MP, location) parameter of Landau density
        limits[2] = [h.Integral() * 500 * k for k in [.001, 10]]          # par[2]=Total area (integral -inf to inf, normalization constant)
        limits[3] = [self.estimate_sigma() * k for k in [.5, 3]]        # par[3]=Width (sigma) of convoluted Gaussian function
        # print limits
        return limits

    def estimate_sigma(self):
        fit = self.Histo.Fit('gaus', 'qs0', '', *[self.get_x_max() * k for k in [.7, 1.3]])
        return fit.Parameter(2)

    def get_x_max(self):
        return self.Histo.GetBinCenter(self.Histo.GetMaximumBin())

    def init_start_values(self):
        h = self.Histo
        values = zeros(4)
        values[0] = self.estimate_sigma() / 5
        values[1] = self.get_x_max()
        values[2] = h.Integral() * 500
        values[3] = self.estimate_sigma()
        return values

    def Mean(self, xmin, xmax):
        return self.Fit.Mean(xmax, xmin)

    def langau(self, x, par):

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

        x = x[0]

        # MP shift correction
        mpshift = -0.22278298  # Landau maximum location
        mpc = par[1] - mpshift * par[0]

        # Range of convolution integral
        xmin, xmax = [x + i * self.NSigma * par[3] for i in [-1, 1]]
        step = (xmax - xmin) / self.NConvolutions

        # Variables
        sum_int = 0.
        i = 1.

        # Convolution integral of Landau and Gaussian by sum
        while i <= self.NConvolutions / 2:
            xx = xmin + (i - .5) * step
            fland = TMath.Landau(xx, mpc, par[0]) / par[0]
            sum_int += fland * TMath.Gaus(x, xx, par[3])

            xx = xmax - (i - .5) * step
            fland = TMath.Landau(xx, mpc, par[0]) / par[0]
            sum_int += fland * TMath.Gaus(x, xx, par[3])

            i += 1

        return par[2] * step * sum_int / sqrt(2 * pi) / par[3]

    def langaufit(self):

        his = self.Histo
        name = 'Fitfcn_{n}'.format(n=his.GetName())

        ffitold = gROOT.GetListOfFunctions().FindObject(name)
        if ffitold:
            ffitold.Delete()

        ffit = TF1(name, self.langau, 0, his.GetXaxis().GetXmax(), 4)
        ffit.SetNpx(1000)
        ffit.SetParameters(self.StartValues)
        ffit.SetParNames('Width', 'MP', 'Area', 'GSigma')

        for i in xrange(4):
            ffit.SetParLimits(i, *self.ParLimits[i])

        his.Fit(name, 'QRB0', '', *self.FitRange)  # fit within specified range, use ParLimits, do not plot

        ffit.GetParameters(self.Parameters)  # obtain fit parameters
        self.ParErrors = [ffit.GetParError(i) for i in xrange(4)]

        self.Chi2 = ffit.GetChisquare()  # obtain chi^2
        self.NDF = ffit.GetNDF()  # obtain ndf
        self.Fit = ffit
        return ffit  # return fit function

    def langaupro(self, params, maxx, FWHM):

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
        l = -1.0

        x = array("d", [0])

        while ((l != lold) and (i < MAXCALLS)):
            i += 1

            lold = l
            x[0] = p + step
            l = self.langau(x, params)

            if (l < lold):
                step = -step / 10

            p += step

        if (i == MAXCALLS):
            return (-1)

        maxx = x[0]

        fy = l / 2

        #  Search for right x location of fy

        p = maxx + params[0]
        step = params[0]
        lold = -2.0
        l = -1e300
        i = 0

        while ((l != lold) and (i < MAXCALLS)):
            i += 1

            lold = l
            x[0] = p + step
            l = TMath.Abs(self.langau(x, params) - fy)

            if (l > lold):
                step = -step / 10

            p += step

        if (i == MAXCALLS):
            return (-2)

        fxr = x[0]

        #  Search for left x location of fy

        p = maxx - 0.5 * params[0]
        step = -params[0]
        lold = -2.0
        l = -1e300
        i = 0

        while ((l != lold) and (i < MAXCALLS)):
            i += 1

            lold = l
            x[0] = p + step
            l = TMath.Abs(self.langau(x, params) - fy)

            if (l > lold):
                step = -step / 10

            p += step

        if (i == MAXCALLS):
            return (-3)

        fxl = x[0]

        FWHM = fxr - fxl
        return (0)


if __name__ == "__main__":
    # Fill Histogram
    data = [0, 0, 0, 0, 0, 0, 2, 6, 11, 18, 18, 55, 90, 141, 255, 323, 454, 563, 681,
            737, 821, 796, 832, 720, 637, 558, 519, 460, 357, 291, 279, 241, 212,
            153, 164, 139, 106, 95, 91, 76, 80, 80, 59, 58, 51, 30, 49, 23, 35, 28,
            23, 22, 27, 27, 24, 20, 16, 17, 14, 20, 12, 12, 13, 10, 17, 7, 6, 12, 6,
            12, 4, 9, 9, 10, 3, 4, 5, 2, 4, 1, 5, 5, 1, 7, 1, 6, 3, 3, 3, 4, 5, 4, 4, 2, 2, 7, 2, 4]
    hSNR = TH1F("snr", "Signal-to-noise", 400, 0, 400)
    for n in xrange(98):
        hSNR.Fill(n, data[n])

    z = Langau(hSNR)

    # Fitting SNR histo
    print "Fitting...\n"
    fitsnr = z.langaufit()

    # z.langaupro(fp, SNRPeak, SNRFWHM)

    print "Fitting done\nPlotting results...\n"

    # Global style settings
    gStyle.SetOptStat(1111)
    gStyle.SetOptFit(111)
    gStyle.SetLabelSize(0.03, "x")
    gStyle.SetLabelSize(0.03, "y")

    hSNR.GetXaxis().SetRange(0, 70)
    hSNR.Draw()
    fitsnr.Draw("lsame")
    raw_input("wait")
