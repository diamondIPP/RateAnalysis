#!/usr/bin/env python
# --------------------------------------------------------
#       Pedestal analysis of the pad waveforms
# created on May 23rd 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import gPad
from numpy import log

from src.sub_analysis import PadSubAnalysis
from helpers.save_plots import *


class PedestalAnalysis(PadSubAnalysis):
    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='Pedestal')
        self.SignalName = self.get_signal_name()
        self.RawName = self.get_signal_name(peak_int=1)
        self.Region = self.Ana.get_region('pedestal')[0]

        self.Histogram = None
        
    def __call__(self, err=False):
        return (self.get_mean(), self.get_noise()) if err else (self.get_mean().n, self.get_noise().n)

    # ----------------------------------------
    # region GET
    def get_bins(self, bin_size=.1, ph_max=150):
        return self.Bins.make(-ph_max, ph_max, bin_size)

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

    def get_under_signal(self, err=True, redo=False):
        def f():
            return FitRes(self.draw_under_signal(show=False).Fit('gaus', 'qs'))
        return do_pickle('USig', f, redo=redo).get_pars(err)[1:]

    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_map(self, res=None, fid=False, cut=None, show=True):
        cut = self.Cut.generate_custom(exclude=['fiducial'], prnt=False) if not fid and cut is None else self.Cut(cut)
        h0 = self.Draw.prof2d(*self.get_tree_vec(self.Ana.get_track_vars() + [self.get_signal_var()], self.Cut()), self.Bins.get_global(res), show=False)
        fid_vals = [h0.GetBinContent(ibin) for ibin in range(h0.GetNbinsY() * h0.GetNbinsX()) if h0.GetBinEntries(ibin) > 20]
        v = self.get_tree_vec(self.Ana.get_track_vars() + [self.get_signal_var()], cut)
        self.Draw.prof2d(*v, self.Bins.get_global(res), 'Pedestal Map', show=show, z_range=ax_range(fid_vals, 0))
        self.Cut.draw_fid()
        update_canvas()
        return fid_vals

    def draw_under_signal(self, name=None, cut=None, show=True):
        x = self.get_tree_vec(var=choose(name, self.Ana.get_signal_var(off_corr=False, evnt_corr=False)), cut=choose(cut, self.Cut.get_pulser().Value))
        return self.Draw.distribution(x, self.get_bins(), 'Pedestal under Signal', x_tit='Pedestal [mV]', y_off=1.8, show=show, lm=.13, x_range=ax_range(x, 0, .1, .1, thresh=5))

    def compare_under_signal(self, cut=None, bin_size=None, x_range=None):
        cut = choose(cut, self.Cut.get_pulser().Value)
        histos = [self.draw_distribution(name=n, cut=cut, show=False, redo=True, bin_size=bin_size) for n in [None, self.get_signal_name('aa'), self.get_signal_name('ad')]]
        s = self.Draw.stack(histos, 'Pedestal Comparison', ['Bucket 0', 'Bucket 1', 'Bucket 6'], scale=True)
        format_histo(s, x_range=choose(x_range, [-20, 20]))
        print([h.GetRMS() for h in histos])

    def draw_diffs(self, cut=None, bin_size=None, x_range=None):
        cut = choose(cut, self.Cut.get_pulser().Value)
        h = [self.draw_distribution(name=n, cut=cut, show=False, redo=True, bin_size=bin_size) for n in [None, None, self.get_signal_name('aa'), self.get_signal_name('ad')]]
        for ih in h:
            ih.Sumw2(True)
            ih.Scale(1 / ih.GetMaximum())
        [h[i].Add(h[i + 2], -1) for i in range(2)]
        graphs = [self.Draw.make_graph_from_profile(h[i]) for i in range(2)]
        mg = self.Draw.multigraph(graphs, 'Pedestal differences', ['#Delta 0-1', '#Delta 0-6'])
        format_histo(mg, x_range=choose(x_range, [-20, 20]))

    def draw_correlation(self, r0=None, r1='aa', cut=None, bin_size=.1, rnge=None, show=True):
        cut = choose(cut, self.Cut.get_pulser().Value)
        x, y = [self.get_tree_vec(var=self.get_signal_var(name), cut=self.Cut(cut)) for name in [self.get_signal_name(r) for r in [r0, r1]]]
        rnge = choose(rnge, [-15, 15])
        h = self.Draw.histo_2d(x, y, self.get_bins(bin_size, 300) * 2, 'Pedstal Correlation', x_tit='Bucket 5', y_tit='Bucket 6', x_range=rnge, y_range=rnge, show=show, grid=True)
        Draw.info('Correlation Factor: {:.2f}'.format(h.GetCorrelationFactor()))
        return h

    def draw_distribution(self, name=None, bin_size=None, cut=None, logy=False, show=True, save=True, redo=False, prnt=True, **kwargs):
        def f():
            info('Drawing pedestal distribution for {d} of run {r}'.format(d=self.DUT.Name, r=self.Run.Number), prnt=prnt)
            x = self.get_tree_vec(var=self.get_signal_var(name), cut=self.Cut(cut))
            bs = choose(bin_size, max(.1, 30 / sqrt(x.size)))
            return self.Draw.distribution(x, self.get_bins(bs), 'Pedestal', x_tit='Pedestal [mV]', show=False, x_range=ax_range(x, 0, .1, .1, thresh=5), y_off=1.8)
        h = do_pickle(self.make_simple_pickle_path('Disto', '{}_{}_{}'.format(self.Cut(cut).GetName(), self.get_short_name(name), bin_size)), f, redo=redo)
        format_histo(h, **kwargs)
        self.Draw(h, 'PedestalDistribution{}'.format(self.Cut(cut).GetName()), show, save=save, logy=logy, prnt=prnt, lm=.13, stats=None)
        return h

    def draw_disto_fit(self, name=None, bin_size=None, cut=None, logy=False, show=True, save=True, redo=False, prnt=True, draw_cut=False, **kwargs):
        cut = self.Cut.generate_custom(exclude='ped sigma') if draw_cut else self.Cut(cut)
        h = self.draw_distribution(name, bin_size, cut, logy, show=show, save=save, redo=redo, prnt=prnt, **kwargs)
        fit_pars = do_pickle(self.make_simple_pickle_path(suf='{}_fwhm_{}'.format(cut.GetName(), self.get_short_name(name))), fit_fwhm, redo=True, h=h, show=True)
        f = Draw.make_f('f', 'gaus', -100, 100, pars=fit_pars[...], npx=1000, line_style=2)
        h.GetFunction('gaus').Draw('same')
        f.Draw('same')
        format_statbox(h, fit=True, all_stat=True)
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
        g = do_pickle(self.make_simple_pickle_path('Trend', '{}_{}{}'.format(int(sigma), self.get_short_name(signal_name), self.Bins.w(bin_size))), f, redo=redo)
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
