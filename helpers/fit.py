#!/usr/bin/env python
from ROOT import TF1, Math, TMath
from scipy.special import erf
from helpers.draw import *


class Fit(object):
    """ general class to perform fits on a histgram"""
    def __init__(self, name='fit', h=None, fit_range=None, npx=1000, invert=False):
        self.Name = name
        self.Histo = h

        # Range and Values
        self.XMin, self.XMax = [0, 1000] if fit_range is None else fit_range
        if h is not None:
            self.Values = get_h_values(h)
            self.X = get_hist_args(h)

        # Fit
        self.ParNames = self.set_par_names()
        self.NPars = len(self.ParNames)
        self.Invert = invert
        self.clear_old()
        self.Fit = self.init_fit()
        self.Fit.SetNpx(npx)
        self.set_par_limits()
        self.set_start_values()

    def __call__(self, *args, **kwargs):
        return self.fit(*args, **kwargs)

    def clear_old(self):
        old = get_object(self.Name)
        if old:
            old.Delete()

    def init_fit(self):
        return TF1()

    def set_par_names(self):
        pass

    def set_par_limits(self):
        pass

    def set_start_values(self):
        pass

    def get_chi2(self):
        return self.Fit.GetChisquare() / self.Fit.GetNDF()

    def get_parameter(self, i):
        return ufloat(self.Fit.GetParameter(i), self.Fit.GetParError(i))

    def get_parameters(self):
        for i in range(self.NPars):
            print('{}: {:2.1f}'.format(self.Fit.GetParName(i), self.Fit.GetParameter(i)))

    def _get_rise_time(self, p=.1, show=False, off_par=6):
        maxval = self.Fit.GetMaximum() - self.Fit.GetParameter(off_par)
        t1, t0 = self.Fit.GetX((1 - p) * maxval, 0, self.Fit.GetX(maxval)), self.Fit.GetX(p * maxval)
        if show:
            Draw.vertical_line(t1, -100, 1e5)
            Draw.vertical_line(t0, -100, 1e5)
        return t1 - t0

    def fit(self, n=1, show=True, minuit=True):
        if minuit:
            Math.MinimizerOptions.SetDefaultMinimizer('Minuit2', 'Migrad')
        for _ in range(n):
            self.Histo.Fit(self.Fit, 'qs0', '', self.XMin, self.XMax)
        if show:
            self.Fit.Draw('same')

    def draw(self, *args, **kwargs):
        pass


class Crystalball(Fit):
    def __init__(self, h=None, fit_range=None, inv=False, npx=1000):
        """ :parameter:  0 - scale, 1 - alpha, 2 - n, 3 - mean, 4 - sigma, 5 - offset """
        Fit.__init__(self, 'cystalball', h, fit_range, npx, inv)

    def set_par_names(self):
        return ['c', 'alpha', 'n', 'mean', 'sigma', 'offset']

    def init_fit(self):
        return TF1(self.Name, partial(crystalball, inv=self.Invert), self.XMin, self.XMax, self.NPars)

    def draw(self, c=1, alpha=1, n=1, m=20, sigma=2, off=0):
        self.Fit.SetParameters(c, alpha, n, m, sigma, off)
        Draw.histo(self.Fit)

    def set_par_limits(self):
        if self.Histo is not None:
            maxval = max(self.Values).n
            max_x = self.X[where(self.Values == max(self.Values))][0].n
            self.Fit.SetParLimits(0, 1, 2 * maxval)
            self.Fit.SetParLimits(1, .1, 10)
            self.Fit.SetParLimits(2, 1, 50)
            self.Fit.SetParLimits(3, .9 * max_x, 1.1 * max_x)
            self.Fit.SetParLimits(4, 1e-2, self.XMax - self.XMin)
            self.Fit.SetParLimits(5, -.1 * maxval, .1 * maxval)
            self.Fit.SetParameters(maxval, .5, 1, max_x, (self.XMax - self.XMin) / 4., 0)


class ErfLand(Fit):
    def __init__(self, h=None, fit_range=None, npx=1000):
        """ :parameter:  0 - erf-scale, 1 - alpha, 2 - n, 3 - mean, 4 - sigma, 5 - offset """
        Fit.__init__(self, 'erfland', h, fit_range, npx)

    def set_par_names(self):
        return ['landau-scale', 'mpv', 'sigma', 'erf-scale', 'xoff', 'width', 'offset', 'x0']

    def init_fit(self):
        return TF1(self.Name, erfland, self.XMin, self.XMax, 8)

    def draw(self, c0=1, mpv=7, sigma=2, c1=1, xoff=3, w=2, yoff=0):
        self.Fit.SetParameters(c0, mpv, sigma, c1, xoff, w, yoff)
        Draw.histo(self.Fit)

    def get_rise_time(self, p=.1, show=False):
        return self._get_rise_time(p, show, off_par=6)

    def set_par_limits(self):
        if self.Histo is not None:
            maxval = max(self.Values).n
            max_x = self.X[where(self.Values == max(self.Values))][0].n
            w = self.XMax - self.XMin
            self.Fit.SetParLimits(0, 1, 10 * maxval)
            self.Fit.SetParLimits(1, .9 * max_x, 1.1 * max_x)
            self.Fit.SetParLimits(2, 1e-2, w)
            self.Fit.SetParLimits(3, 1, 10 * maxval)
            self.Fit.SetParLimits(4, 1, 1.5 * max_x)
            self.Fit.SetParLimits(5, .1, 1)
            self.Fit.SetParLimits(6, -10, 10)
            w1 = max_x - self.Histo.GetBinCenter(self.Histo.FindFirstBinAbove(.1 * maxval))
            self.Fit.FixParameter(7, max_x - .4 * w1)
            self.Fit.SetParameters(maxval * 5, max_x, 3, maxval / 2., .5, max_x - 10, 0, max_x - .5 * w1)


class Langau(Fit):
    def __init__(self, h=None, nconv=100, fit_range=None, npx=1000):

        self.NConvolutions = nconv
        self.NSigma = 5.
        Fit.__init__(self, 'langau', h, fit_range, npx)
        self.XMin, self.XMax = [k * self.Histo.GetMean() for k in [.1, 3]] if fit_range is None else fit_range
        print(self.XMin, self.XMax)

    def init_fit(self):
        return TF1(self.Name, partial(langau, nconv=self.NConvolutions, nsigma=self.NSigma), 0, self.get_x_max() * 3, self.NPars)

    def set_par_names(self):
        return ['Width', 'MPV', 'Area', 'GSigma']

    def set_par_limits(self):
        sigma = self.estimate_sigma()
        self.Fit.SetParLimits(0, 0, .6 * sigma)                                 # Width (scale) parameter of Landau density
        self.Fit.SetParLimits(1, *array([.5, 1.5]) * self.get_x_max())          # Most Probable (MPV, location) parameter of Landau density
        self.Fit.SetParLimits(2, *array([.5, 5000]) * self.Histo.Integral())    # Total area (integral -inf to inf, normalization constant)
        self.Fit.SetParLimits(3, *array([.5, 3]) * sigma)                       # Width (sigma) of convoluted Gaussian function
        self.Fit.SetParameters(sigma / 5, self.get_x_max(), self.Histo.Integral() * 500, sigma)

    def estimate_sigma(self):
        fit = self.Histo.Fit('gaus', 'qs0', '', *array([.7, 1.3]) * self.get_x_max())
        return fit.Parameter(2)

    def get_x_max(self):
        return 1000 if self.Histo is None else self.Histo.GetBinCenter(self.Histo.GetMaximumBin())

    def get_mpv(self):
        return self.get_parameter(1)


class NLandau(Fit):

    def __init__(self, h=None, fit_range=None, npx=100, n=3):
        self.MPV = find_mpv_fwhm(h)[0].n
        self.Max = h.GetMaximum()
        self.W = get_fwhm(h).n / 2
        self.N = n
        super().__init__('TripleLandau', h, fit_range, npx)

    def init_fit(self):
        return TF1('TripelLandau', ' + '.join('landau({})'.format(3 * i) for i in range(0, self.N)), self.XMin, self.XMax)

    def set_par_names(self):
        return [n for i in range(self.N) for n in 'C{0} MPV{0} Sigma{0}'.format(i).split()]

    def set_par_limits(self):
        means = linspace(-self.W, self.W, self.N + 1) + self.MPV
        for j, i in enumerate(range(0, self.N * 3, 3)):
            self.Fit.SetParLimits(i, .2 * self.Max, self.Max * 5)
            self.Fit.SetParLimits(i + 1, *means[j:j + 2])
            self.Fit.SetParLimits(i + 2, .3 * self.W, 1.5 * self.W)

    def set_start_values(self):
        for i in range(0, self.N * 3, 3):
            self.Fit.SetParameter(i, self.Max)
            self.Fit.SetParameter(i + 1, self.MPV + (i - 1) * self.W / 2)
            self.Fit.SetParameter(i + 2, self.W)


def erfland(x, pars):
    c0, mpv, sigma, c1, xoff, w, yoff, x0 = [float(p) for p in pars]
    return yoff + (c0 * TMath.Landau(x[0], mpv, sigma) if x[0] > x0 else c1 * (erf(w * (x[0] - xoff)) + 1))


def crystalball(x, pars, inv=False):
    scale, alpha, n, m, sigma, off = [float(p) for p in pars]
    x[0] *= -1 if inv else 1
    m *= -1 if inv else 1
    if (x[0] - m) / sigma > -alpha:
        return gauss(x[0], scale, m, sigma, off)
    else:
        a = (n / abs(alpha)) ** n * exp(-abs(alpha) ** 2 / 2)
        b = n / abs(alpha) - abs(alpha)
        return scale * a * (b - (x[0] - m) / sigma) ** -n + off


def langau(x, pars, nconv, nsigma=5):
    # Convoluted Landau and Gaussian Fitting Function
    #         (using ROOT's Landau and Gauss functions)
    # Fit parameters:
    # par[0]=Width (scale) parameter of Landau density
    # par[1]=Most Probable (MP, location) parameter of Landau density
    # par[2]=Total area (integral -inf to inf, normalization constant)
    # par[3]=Width (sigma) of convoluted Gaussian function

    x = x[0]

    # MP shift correction
    mpshift = -0.22278298  # Landau maximum location
    mpc = pars[1] - mpshift * pars[0]

    # Range of convolution integral
    xmin, xmax = [x + i * nsigma * pars[3] for i in [-1, 1]]
    step = (xmax - xmin) / nconv

    # Variables
    sum_int = 0.
    i = 1.

    # Convolution integral of Landau and Gaussian by sum
    while i <= nconv / 2:
        xx = xmin + (i - .5) * step
        fland = TMath.Landau(xx, mpc, pars[0]) / pars[0]
        sum_int += fland * TMath.Gaus(x, xx, pars[3])

        xx = xmax - (i - .5) * step
        fland = TMath.Landau(xx, mpc, pars[0]) / pars[0]
        sum_int += fland * TMath.Gaus(x, xx, pars[3])

        i += 1

    return pars[2] * step * sum_int / sqrt(2 * pi) / pars[3]


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
    l = -1.0

    x = array("d", [0])

    while ((l != lold) and (i < MAXCALLS)):
        i += 1

        lold = l
        x[0] = p + step
        l = langau(x, params)

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
        l = TMath.Abs(langau(x, params) - fy)

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
        l = TMath.Abs(langau(x, params) - fy)

        if (l > lold):
            step = -step / 10

        p += step

    if (i == MAXCALLS):
        return (-3)

    fxl = x[0]

    FWHM = fxr - fxl
    return (0)


if __name__ == '__main__':
    z = Crystalball(fit_range=[-10, 20])
