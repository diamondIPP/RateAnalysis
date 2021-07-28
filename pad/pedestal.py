#!/usr/bin/env python
# --------------------------------------------------------
#       Pedestal analysis of the pad waveforms
# created on May 23rd 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import log

from src.sub_analysis import PadSubAnalysis
from helpers.save_plots import *


class PedestalAnalysis(PadSubAnalysis):
    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='Pedestal')
        self.SignalName = self.get_signal_name()
        self.RawName = self.get_signal_name(peak_int=1)
        self.Region = self.Ana.get_region('pedestal')[0]

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

    @save_pickle(suf_args='all')
    def get_fit(self, name=None, bin_scale=None, cut=None, _redo=False):
        return self.draw_disto_fit(name, bin_scale, cut, _redo, show=False, prnt=False)

    def get(self, par=1, name=None, bin_scale=None, cut=None, redo=False):
        return self.get_fit(name, bin_scale, cut, _redo=redo)[par]

    def get_mean(self, name=None, bin_scale=None, cut=None, redo=False, raw=False):
        return self.get(1, name, bin_scale, cut, redo) * (self.Polarity if raw else 1)

    def get_noise(self, name=None, bin_scale=None, cut=None, redo=False):
        return self.get(2, name, bin_scale, cut, redo)

    def get_fwhm(self, raw=False, redo=False):
        return (self.get_raw_noise if raw else self.get_noise)(redo=redo) * 2 * sqrt(2 * log(2))

    def get_raw_mean(self, cut=None, redo=False):
        return self.get_mean(self.RawName, 5, cut, redo)

    def get_raw_noise(self, cut=None, redo=False):
        return self.get_noise(self.RawName, 5, cut, redo)

    def get_under_signal(self, err=True, redo=False):
        def f():
            return FitRes(self.draw_under_signal(show=False).Fit('gaus', 'qs'))
        return do_pickle(self.make_simple_pickle_path('USig'), f, redo=redo).get_pars(err)[1:]
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
        histos = [self.draw_distribution(name=n, cut=cut, show=False, redo=True, bin_scale=bin_size) for n in [None, self.get_signal_name('aa'), self.get_signal_name('ad')]]
        s = self.Draw.stack(histos, 'Pedestal Comparison', ['Bucket 0', 'Bucket 1', 'Bucket 6'], scale=True)
        format_histo(s, x_range=choose(x_range, [-20, 20]))
        print([h.GetRMS() for h in histos])

    def draw_diffs(self, cut=None, bin_size=None, x_range=None):
        cut = choose(cut, self.Cut.get_pulser().Value)
        h = [self.draw_distribution(name=n, cut=cut, show=False, redo=True, bin_scale=bin_size) for n in [None, None, self.get_signal_name('aa'), self.get_signal_name('ad')]]
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

    @save_pickle('Disto', print_dur=True, suf_args='all')
    def get_distribution(self, sig=None, bin_scale=None, cut=None, _redo=False):
        x = self.get_tree_vec(var=self.get_signal_var(sig), cut=self.Cut(cut))
        bins = self.get_bins(max(.1, self.Bins.find_width(x)) * choose(bin_scale, 1))
        return self.Draw.distribution(x, bins, 'Pedestal Distribution', x_tit='Pedestal [mV]', show=False, x_range=ax_range(x, 0, .2, .4, thresh=5))

    def draw_distribution(self, sig=None, bin_scale=None, cut=None, redo=False, prefix='', **kwargs):
        h = self.get_distribution(sig, bin_scale, cut, _redo=redo)
        return self.Draw.distribution(h, **prep_kw(kwargs, file_name=f'{prefix}PedestalDistribution'))

    def draw_disto_fit(self, name=None, bin_scale=None, cut=None, redo=False, draw_cut=False, prefix='', **kwargs):
        cut = self.Cut.generate_custom(exclude='ped sigma') if draw_cut and cut is None else self.Cut(cut)
        h = self.draw_distribution(name, bin_scale, cut, redo, prefix, **kwargs)
        fit = fit_fwhm(h, show=True)
        Draw.make_f('f', 'gaus', -100, 100, pars=fit.Pars, npx=1000, line_style=2).Draw('same')
        h.GetFunction('gaus').Draw('same')
        self.__draw_cut(fit, draw_cut)
        format_statbox(h, fit=True, all_stat=True, form='.2f')
        self.Draw.save_plots('PedestalDistributionFit{}'.format(cut.GetName()), **kwargs)
        return fit

    def __draw_cut(self, fit, show=True):
        if show:
            xmin, xmax = fit[1] + fit[2] * self.Cut.get_ped_sigma() * array([-1, 1])
            b = Draw.box(xmin.n, -10, xmax.n, 1e7, line_color=2, width=2, fillcolor=2, style=7, opacity=.2)
            self.Draw.legend([b], [f'cut ({self.Cut.get_ped_sigma()} sigma)'], 'lf', y2=.54, margin=.45, w=.3)

    def draw_pulse_height(self, name=None, rel_t=False, show=True):
        x, y = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var(name)], cut=self.Cut())
        self.Draw.profile(x, y, self.Bins.get_time(), 'Pedestal Trend', y_tit='Pedestal [mV]', show=show, **self.get_t_args(rel_t), stats=set_statbox(entries=True, w=.2))

    def draw_signal_time(self, name=None, rel_t=False, show=True):
        x, y = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var(name)], cut=self.Cut())
        y_range = ax_range(y, 0, .2, .2)
        return self.Draw.histo_2d(x, y, self.Bins.get_time() + self.get_bins(.5), 'Pedestal vs. Time', y_tit='Pedestal [mV]', pal=53, show=show, y_range=y_range, **self.get_t_args(rel_t))

    @save_pickle('Trend', suf_args='[0, 1, 2]')
    def get_trend(self, signal_name=None, bin_size=None, sigma=False, _redo=False):
        (x, y), bins = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var(signal_name)], cut=self.Cut()), self.Bins.get_time(bin_size)[-1]
        x, y = bins[:-1] + diff(bins) / 2, binned_stats(x, y, mean_sigma, bins)[:, 1 if sigma else 0]
        name = 'Sigma' if sigma else 'Pedestal'
        return self.Draw.graph(x, y, title=f'{name} Trend', y_tit=f'{name} [mV]', show=False)

    def draw_trend(self, signal_name=None, bin_size=None, sigma=False, fit=False, rel_t=False, redo=False, **kwargs):
        g = self.get_trend(signal_name, bin_size, sigma, _redo=redo)
        g.Fit('pol0', f'qs{"" if fit else "0"}')
        return self.Draw.graph(g, **self.get_t_args(rel_t), **kwargs, stats=set_statbox(fit=fit, form='.2f'))

    def fit_trend(self, signal_name=None, show=True, sigma=False):
        g = self.draw_trend(signal_name, show=False, sigma=sigma)
        fit = g.Fit('pol0', 'qs')
        format_statbox(g, fit=True)
        self.Draw(g, '{n}Fit'.format(n=g.GetTitle()), show, lm=.12)
        return FitRes(fit)

    def draw_sigma_trend(self, signal_name=None, bin_size=None, show=True):
        return self.draw_trend(signal_name, bin_size, sigma=True, show=show)

    def compare(self, sigma=False):
        values = [self.get(par=2 if sigma else 1, name=name) for name in self.get_all_signal_names()]
        g = self.Draw.graph(arange(len(values)), values, title='Pedestal Comparison', x_tit='Region Integral', y_tit='Mean {}'.format('Sigma' if sigma else 'Pedestal'), x_off=1.5, bm=.2,
                            x_range=ax_range(0, len(values) - 1, .2, .2))
        set_bin_labels(g, list(self.get_all_signal_names().values()))
    # endregion DRAW
    # ----------------------------------------
