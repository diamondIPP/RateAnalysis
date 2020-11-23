#!/usr/bin/env python
# --------------------------------------------------------
#       Pedestal analysis of the pad waveforms
# created on May 23rd 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import gPad
from numpy import log, diff

from src.sub_analysis import PadSubAnalysis
from helpers.save_plots import *


class PedestalAnalysis(PadSubAnalysis):
    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='Pedestal')
        self.SignalName = self.get_signal_name()
        self.RawName = self.get_signal_name(peak_int=1)

        self.Histogram = None

    # ----------------------------------------
    # region GET
    def get_bins(self, bin_size=.1):
        return self.Bins.make(-150, 150, bin_size)

    def get_signal_name(self, region=None, peak_int=None):
        return self.Ana.get_signal_name(region=region, peak_int=peak_int, sig_type='pedestal')

    def get_signal_var(self, name=None, peak_int=None):
        return self.Ana.get_signal_var(choose(name, self.get_signal_name(peak_int=peak_int)), off_corr=False, evnt_corr=False)

    def get_raw_var(self):
        return self.get_signal_var(self.RawName)

    def get_all_signal_names(self):
        return self.Ana.get_all_signal_names('pedestal')

    def get_short_name(self, name=None):
        return self.get_all_signal_names()[choose(name, self.SignalName)]

    def get_par(self, par=1, name=None, cut=None, redo=False):
        suffix = '{}_{}'.format(self.Cut(cut).GetName(), self.get_short_name(name))
        return fit2u(do_pickle(self.make_simple_pickle_path(suf=suffix), partial(self.draw_disto_fit, name=name, cut=self.Cut(cut), show=False, redo=redo), redo=redo), par=par)

    def get_mean(self, name=None, cut=None, redo=False, raw=False):
        return self.get_par(1, name, cut, redo) * (self.Polarity if raw else 1)

    def get_noise(self, name=None, cut=None, redo=False):
        return self.get_par(2, name, cut, redo)

    def get_fwhm(self):
        return self.get_noise() * 2 * sqrt(2 * log(2))

    def get_raw_mean(self, cut=None, redo=False):
        return self.get_mean(self.RawName, cut, redo)

    def get_raw_noise(self, cut=None, redo=False):
        return self.get_noise(self.RawName, cut, redo)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_under_signal(self, name=None, show=True):
        x = self.get_tree_vec(var=choose(name, self.Ana.get_signal_var(off_corr=False, evnt_corr=False)), cut=self.Cut.get_pulser().Value)
        self.Draw.distribution(x, self.get_bins(), 'Pedestal under Signal', x_tit='Pedestal [mV]', y_off=1.8, show=show, lm=.13, x_range=ax_range(x, 0, .1, .1, thresh=5))

    def draw_distribution(self, name=None, cut=None, logy=False, show=True, save=True, redo=False, prnt=True, normalise=None):
        def f():
            info('Drawing pedestal distribution for {d} of run {r}'.format(d=self.DUT.Name, r=self.Run.Number), prnt=prnt)
            x = self.get_tree_vec(var=self.get_signal_var(name), cut=self.Cut(cut))
            return self.Draw.distribution(x, self.get_bins(max(.1, 30 / sqrt(x.size))), 'Pedestal', x_tit='Pedestal [mV]', show=False, x_range=ax_range(x, 0, .1, .1, thresh=5), y_off=1.8)
        h = do_pickle(self.make_simple_pickle_path('Disto', '{}_{}'.format(self.Cut(cut).GetName(), self.get_short_name(name))), f, redo=redo)
        format_histo(h, normalise=normalise, sumw2=False)
        self.Draw(h, 'PedestalDistribution{}'.format(self.Cut(cut).GetName()), show, save=save, logy=logy, prnt=prnt, lm=.13, stats=None)
        return h

    def draw_disto_fit(self, name=None, cut=None, logy=False, show=True, save=True, redo=False, prnt=True, draw_cut=False, normalise=None):
        cut = self.Cut.generate_custom(exclude='ped sigma') if draw_cut else self.Cut(cut)
        h = self.draw_distribution(name, cut, logy, show=show, save=save, redo=redo, prnt=prnt, normalise=normalise)
        fit_pars = do_pickle(self.make_simple_pickle_path(suf='{}_fwhm_{}'.format(cut.GetName(), self.get_short_name(name))), fit_fwhm, redo=True, h=h, show=True)
        f = deepcopy(h.GetFunction('gaus'))
        f.SetNpx(1000)
        f.SetRange(h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        f.SetLineStyle(2)
        h.GetListOfFunctions().Add(f)
        format_statbox(h, fit=True)
        if draw_cut:
            b = self.__draw_cut(h)
            h.Draw('same')
            b.Draw('l')
            gPad.RedrawAxis()
        self.Draw.save_plots('PedestalDistributionFit{}'.format(cut.GetName()), save=save, prnt=prnt, show=show)
        SaveDraw.server_pickle(self.make_simple_pickle_path(suf='{}_fwhm_{}'.format(cut.GetName(), self.get_short_name(name))), fit_pars)
        return fit_pars

    @staticmethod
    def __draw_cut(h):
        fit = h.GetListOfFunctions()[2]
        xmin, xmax = fit.GetParameter(1) - 3 * fit.GetParameter(2), fit.GetParameter(1) + 3 * fit.GetParameter(2)
        b = Draw.box(xmin, -10, xmax, 1e7, line_color=2, width=2, fillstyle=3001, name='ped', style=7)
        legend = Draw.make_legend(.65, y2=.64, nentries=1, margin=.45, scale=1)
        legend.AddEntry(b, 'cut (3 sigma)', 'lf')
        legend.Draw()
        return b

    def draw_pulse_height(self, name=None, rel_t=False, show=True):
        x, y = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var(name)], cut=self.Cut())
        self.Draw.profile(x, y, self.Bins.get_time(), 'Pedestal Trend', y_tit='Pedestal [mV]', show=show, **self.get_t_args(rel_t), stats=set_statbox(entries=True, w=.2))

    def draw_signal_time(self, name=None, rel_t=False, show=True):
        x, y = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var(name)], cut=self.Cut())
        y_range = ax_range(y, 0, .2, .2)
        return self.Draw.histo_2d(x, y, self.Bins.get_time() + self.get_bins(.5), 'Pedestal vs. Time', y_tit='Pedestal [mV]', pal=53, show=show, y_range=y_range, **self.get_t_args(rel_t))

    def test(self):
        x, y = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var()], cut=self.Cut())
        return binned_stats(x, y, mean_sigma, self.Bins.get_time()[-1])

    def draw_trend(self, signal_name=None, bin_size=None, sigma=False, rel_t=False, redo=False, show=True):
        def f():
            (x, y), bins = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var(signal_name)], cut=self.Cut()), self.Bins.get_time(bin_size)[-1]
            x, y = bins[:-1] + diff(bins) / 2, binned_stats(x, y, mean_sigma, bins)[:, 1 if sigma else 0]
            name = 'Sigma' if sigma else 'Pedestal'
            return self.Draw.graph(x, y, title='{} Trend'.format(name), y_tit='{} [mV]'.format(name), show=False)
        g = do_pickle(self.make_simple_pickle_path('Trend', '{}_{}{}'.format(int(sigma), self.get_short_name(signal_name), self.Bins(bin_size))), f, redo=redo)
        format_histo(g, **self.get_t_args(rel_t))
        self.Draw(g, show=show)
        return g

    def fit_trend(self, signal_name=None, show=True, sigma=False):
        g = self.draw_trend(signal_name, show=False, sigma=sigma)
        fit = g.Fit('pol0', 'qs')
        format_statbox(g, fit=True)
        self.Draw(g, '{n}Fit'.format(n=g.GetTitle()), show, lm=.12)
        return FitRes(fit)

    def draw_sigma_trend(self, signal_name=None, bin_size=None, show=True):
        return self.draw_trend(signal_name, bin_size, sigma=True, show=show)

    def compare(self, sigma=False):
        values = [self.get_par(par=2 if sigma else 1, name=name) for name in self.get_all_signal_names()]
        g = self.Draw.graph(arange(len(values)), values, title='Pedestal Comparison', x_tit='Region Integral', y_tit='Mean {}'.format('Sigma' if sigma else 'Pedestal'), x_off=1.5, bm=.2,
                            x_range=ax_range(0, len(values) - 1, .2, .2))
        set_bin_labels(g, list(self.get_all_signal_names().values()))
    # endregion DRAW
    # ----------------------------------------
