#!/usr/bin/env python
# --------------------------------------------------------
#       class for the pulse height calibration of the digital CMS pixel chips
# created on August 2nd 2021 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from helpers.utils import join, file_exists, save_pickle, choose, prep_kw, deepcopy, warning, is_iter
from plotting.draw import FitRes, Draw
from plotting.fit import Erf
from numpy import genfromtxt, arange, array, concatenate, full
from src.sub_analysis import SubAnalysis
from src.dut import Plane


class Calibration(SubAnalysis):

    def __init__(self, pix_analysis):
        super().__init__(pix_analysis, pickle_dir='Calibration')

        self.TelescopeID = self.Run.Converter.TelescopeID
        self.Dir = join(self.Run.Converter.TrackingDir, 'data', 'calibrations', f'telescope{self.TelescopeID}')
        self.N = self.Ana.N

        self.Fit = self.Draw.make_f('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -5000, 255 * 7, npx=100)
        if self.check_fit_files():
            self.Parameters = self.load_fitpars()
        self.HasPoints = self.check_files()
        if self.HasPoints:
            self.Vcals = self.load_vcals()
            self.Points = self.load_calibration_points()

    # ----------------------------------------
    # region INIT
    def load_fit_files(self):
        return [join(self.Dir, f'ROC{i}.txt') for i in range(self.Ana.NRocs)]

    def check_fit_files(self):
        return True if all(file_exists(f, warn=False) for f in self.load_fit_files()) else False

    def load_files(self):
        return [join(self.Dir, f'cal{i}.txt') for i in range(self.Ana.NRocs)]

    def check_files(self):
        return True if all(file_exists(f, warn=False) for f in self.load_files()) else False

    def verify_fits(self):
        """ check if fits represent the data. """

    @save_pickle('FitPars', run='TelescopeID', dut='', print_dur=True)
    def load_fitpars(self, _redo=False):
        return array([genfromtxt(f, skip_header=3, usecols=arange(4)).reshape((Plane.NCols, Plane.NRows, 4)) for f in self.load_fit_files()])

    @save_pickle('Vcals', run='TelescopeID', dut='')
    def load_vcals(self, _redo=False):
        return [concatenate([genfromtxt(f, 'i2', skip_header=1, max_rows=1)[2:], 7 * genfromtxt(f, 'i2', skip_header=2, max_rows=1)[2:]]) for f in self.load_files()]

    @save_pickle('Points', run='TelescopeID', dut='', print_dur=True)
    def load_calibration_points(self, _redo=False):
        s = [arr.size for arr in self.Vcals]
        return [genfromtxt(f, 'i2', skip_header=4, usecols=arange(s[i])).reshape((Plane.NCols, Plane.NRows, s[i])) for i, f in enumerate(self.load_files())]
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region SET
    def set_parameters(self, col, row):
        self.Fit.SetParameters(*self.Parameters[self.N][col][row])
    # endregion SET
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_x(self, plane=None):
        return self.Vcals[choose(plane, self.N)]

    def get_y(self, col, row, plane=None):
        return self.Points[choose(plane, self.N)][col][row]

    def get_points(self, col, row, zero_sup=False):
        x, y = self.get_x(), self.get_y(col, row)
        return (x[y != 0], y[y != 0]) if zero_sup else (x, y)

    def get_vcal(self, col, row, adc):
        self.set_parameters(col, row)
        return -999 if self.get_min(col, row) > adc - 1 else self.Fit.GetX(adc)

    def get_adc(self, col, row, vcal):
        self.set_parameters(col, row)
        return self.Fit(vcal)

    def get_chi2(self, col, row):
        x, y = self.get_points(col, row, zero_sup=True)
        self.set_parameters(col, row)
        return 999 if y.size < 4 else sum((y - array([self.Fit(ix) for ix in x])) ** 2) / (1 if y.size == 4 else y.size - 4)  # DOF = array size - 4 fit pars

    @save_pickle('Chi2s', run='TelescopeID')
    def get_chi2s(self, _redo=False):
        return array([[self.get_chi2(col, row) for row in range(Plane.NRows)] for col in range(Plane.NCols)])

    def get_threshold(self, col, row, vcal=True):
        thresh = self.get_vcal(col, row, 0)
        return thresh if thresh == -999 else thresh * (self.Bins.Vcal2Ke if not vcal else 1)

    def get_thresholds(self, cols=None, rows=None, pix=None, vcal=True):
        cols, rows = self.Cut.get_fid_lines(cols, rows, pix)
        return array([[col, row, self.get_threshold(col, row, vcal)] for col in cols for row in rows if self.get_threshold(col, row, vcal) != -999])

    def get_vcals(self, col, row, adc):
        return array([self.get_vcal(col[i], row[i], adc[i]) for i in range(col.size)])

    def get_adcs(self, col, row, vcal):
        vcal = vcal if is_iter(vcal) else full(col.size, vcal)
        return array([self.get_vcal(col[i], row[i], vcal[i]) for i in range(col.size)])

    def get_min(self, col, row):
        _, _, off, scale = self.Parameters[self.N, col, row]
        return (off - 1) * scale
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_fit(self, col, row, plane=None, **kwargs):
        if self.HasPoints:
            self.Draw.graph(self.get_x(plane), self.get_y(col, row, plane), x_tit='vcal', y_tit='adc')
        self.Fit.SetParameters(*self.Parameters[choose(plane, self.N)][col][row])
        self.Draw.function(self.Fit, f'Calibration Fit for Pix {col} {row}', **prep_kw(kwargs, x_tit='vcal', y_tit='adc', lw=2, color=632, draw_opt='same' if self.HasPoints else None))

    def draw_erf_fit(self, col, row, **kwargs):
        if not self.HasPoints:
            return warning(f'no calibration points found for telescope {self.TelescopeID} ...')
        g = self.Draw.graph(*self.get_points(col, row, zero_sup=True), f'Calibration Fit for Pix {col} {row}', **prep_kw(kwargs, x_tit='vcal', y_tit='adc'))
        self.Fit.SetParameters(309.2062, 112.8961, 1.022439, 35.89524)
        fit = deepcopy(self.Fit) if g.GetN() > 3 else self.Draw.make_f('pol1', 0, 3000)
        g.Fit(fit, 'q')
        return FitRes(fit)

    def draw_chi2(self, **dkw):
        x = self.get_chi2s()
        self.Draw.distribution(x[x != 999], **prep_kw(dkw, lf=.1, q=.001, x_tit='#chi^{2} / DOF', title='#chi^{2} distribution of the calibration fits'))

    @staticmethod
    def draw_erf():
        e0 = Erf(fit_range=[0, 100], pars=[20, 200, 20, 20], npx=1000)
        e0.draw(line_style=7, x_tit='Charge [au]', y_tit='Pulse Height [adc]', **Draw.mode(2), gridy=True, ndivy=3, l_off_x=5, tick_size=0, center_tit=True)
        Draw.make_tf1('b', lambda x: 0 if x < e0.Fit.GetX(0) else e0.Fit(x), 0, 100, npx=1000).Draw('same')
        # Erf(fit_range=[0, 100], pars=[20, 150, 20, 30]).draw(draw_opt='same', line_style=7, color=1)
        # Erf(fit_range=[0, 100], pars=[50, 200, 20, 20]).draw(draw_opt='same', line_style=7, color=2)
    # endregion DRAW
    # ----------------------------------------


if __name__ == '__main__':
    pass
