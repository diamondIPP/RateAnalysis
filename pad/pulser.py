# --------------------------------------------------------
#       PULSER ANALYSIS
# created on June 16th 2016 by M. Reichmann
# --------------------------------------------------------
from pad.waveform import Waveform
from src.sub_analysis import PadSubAnalysis
from plotting.save import *
from helpers.utils import do_pickle, partial, fit2u, save_pickle


class PulserAnalysis(PadSubAnalysis):

    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='Pulser')
        self.Cut = self.Ana.Cut.get_pulser()
        self.Polarity = self.load_polarity()
        self.SignalName = self.get_signal_name()
        self.SignalRegion = array(self.Run.IntegralRegions[self.Ana.DUT.Number - 1]['pulser']) * self.Ana.DigitiserBinWidth
        self.PedestalName = self.load_pedestal_name()

        self.Type = str(self.Ana.Run.Info['pulser']) if 'pulser' in self.Ana.Run.Info else None

        self.DigitiserBinWidth = self.Ana.DigitiserBinWidth
        self.Waveform = Waveform(self)
        self.Pedestal = self.Ana.Pedestal
        self.PeakName = self.Ana.get_peak_name(type_='pulser')
        self.Peaks = self.Ana.Peaks

    def __repr__(self):
        return f'{super().__repr__()} ({self.Type})'

    # ----------------------------------------
    # region INIT
    def load_pedestal_name(self, peak_int=None):
        region = self.Config.get_value('BASIC', 'pulser pedestal', default='ac')
        return self.Ana.get_pedestal_name(region, self.Ana.PeakIntegral if None else peak_int)

    def load_polarity(self):
        return int(self.Run.TreeConfig.get('General', 'polarities').split()[self.Channel])
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_events(self, cut=None, redo=False):  # for Waveform
        return self.Ana.get_events(self.Cut(cut), redo)

    def get_event_cut(self, cut=None, redo=False):
        cut = self.make_event_cut(self.get_events(cut, redo))
        return cut & (self.Peaks.get_n() == (1 if self.Type == 'extern' else 2))

    def get_signal_range(self, lfac=0, rfac=0, t_corr=True):
        return ax_range(self.SignalRegion / (1 if t_corr else self.DigitiserBinWidth), None, lfac, rfac)

    def get_combined_cut(self, n=1000):
        """returns: cut of [n] pulser and [n] signal events."""
        return self.Ana.make_event_cut(concatenate([self.get_events()[:n], self.Ana.get_events()[:n]]))

    def get_signal_indices(self):  # for Waveform
        """ :returns: indices that are smaller than the smallest signals """
        return where(array(self.get_values()) > self.Ana.get_min_signal() - 5)[0]

    def get_min_signal(self, name=None):
        h = self.draw_distribution(name, show=False, save=False)
        return h.GetBinCenter(h.FindFirstBinAbove(h.GetMaximum() * .01))

    @save_pickle('Rate')
    def rate(self, _redo=False):
        return calc_eff(values=self.Run.get_tree_vec('pulser', dtype=bool, cut=self.Ana.Cut.include('beam stops', 'event range')))

    def get_rate_stability(self, bin_size=10, b=None, show=False):
        return self.Draw.pull(self.draw_rate(bin_size, show=False, cut=self.Ana.Cut['beam stops']), b, show=show)

    def get_t_bins(self, bin_size=None):
        return bins.make(*ax_range(self.SignalRegion, 0, .5, .5), choose(bin_size, default=self.Waveform.BinWidth))

    def get_pulse_height(self, corr=True, beam_on=True, bw=None, redo=False):
        return self.get_distribution_fit(corr, beam_on, bw, _redo=redo)[1]

    def get_sigma(self, corr=True, beam_on=True, bw=None, redo=False):
        return self.get_distribution_fit(corr, beam_on, bw, _redo=redo)[2]

    def get_pedestal(self, par=1, beam_on=True, redo=False):
        pickle_path = self.make_simple_pickle_path('Pedestal', str(int(beam_on)))
        fit = do_pickle(pickle_path, partial(self.draw_pedestal_fit, show=False, prnt=False, redo=redo, beam_on=beam_on), redo=redo)
        return fit[par]

    def get_pedestal_mean(self, redo=False):
        return self.get_pedestal(redo=redo)

    def get_pedestal_sigma(self, redo=False):
        return self.get_pedestal(par=2, redo=redo)

    def get_values(self, cut=None):
        return self.get_tree_vec(var=self.get_signal_var(cut=self.Cut(cut)), cut=self.Cut(cut), dtype='f4')

    def get_cft(self, ind=None):
        cft = self.Peaks.get_cft(ind=ind)
        tmin, tmax = self.SignalRegion
        cft = array([[t for t in lst if tmin <= t <= tmax] for lst in cft])
        return array([lst[0] if len(lst) else 0 for lst in cft])

    def get_peak_times(self, cut=None):
        t = self.Peaks.get_times(flat=True, fit=True, cut=choose(cut, self.get_event_cut()))
        return t if self.Type == 'extern' else t[::2]

    def get_peak_time(self, sigma=False):
        return mean_sigma(self.get_peak_times())[sigma]

    @save_pickle('Meansigma', suf_args=0)
    def get_mean_sigma(self, cut=None, _redo=False):
        return fit_fwhm(self.Draw.distribution(self.get_values(cut), show=False))[1:]

    def get_fraction(self, show=False, prnt=True):
        """ :returns the fitted value of the fraction of pulser events with event range and beam interruptions cuts and its fit error. """
        h = self.draw_rate(show=show, cut=self.Ana.Cut.get('beam stops'), prnt=prnt)
        format_statbox(h, fit=True)
        format_histo(h, markersize=None)
        fit = h.Fit('pol0', 'qs')
        return fit2u(FitRes(fit), par=0)

    def get_real_fraction(self):
        """ :returns the real fraction of pulser events compared to signal events. (pulser in rate / flux / area)  """
        in_rate = 40 if self.Ana.get_flux().n < 20 else 100
        return in_rate / (self.Ana.get_flux() * self.DUT.ActiveArea)

    def get_signal_name(self, peak_int=None):
        return self.Ana.get_signal_name('', peak_int, 'pulser')

    def get_all_signal_names(self):
        return self.Ana.get_all_signal_names('pulser')

    def get_short_name(self, signal=None):
        return self.get_all_signal_names()[choose(signal, self.SignalName)]

    def pedestal(self, cut=None):
        return self.Pedestal.get_mean(cut=self.Cut(cut))

    def get_signal_var(self, name=None, ped_corr=True, cut=None):
        return f'{self.Polarity} * {choose(name, self.SignalName)}{f"- {self.pedestal(cut).n}" if ped_corr else ""}'  # sign of pedestal is already corrected
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_rate(self, bw=10, cut='', rel_t=True, hist=False, **dkw):
        """ Shows the fraction of pulser events as a function of time. Peaks appearing in this graph are most likely beam interruptions. """
        x, y = self.get_tree_vec(var=[self.get_t_var(), 'pulser'], cut=self.Cut(choose(cut, self.Ana.Cut.get('beam stops'))))
        b = self.Bins.get_raw_time(bw)
        f = partial(self.Draw.profile, draw_opt='hist', **Draw.mode(3)) if hist else self.Draw.efficiency
        return f(x, y * 100 if hist else 1, b, 'Pulser Rate', **prep_kw(dkw, y_tit='Pulser Fraction [%]', gridy=True, y_range=[0, 105], **self.get_t_args(rel_t), file_name='PulRate', stats=False))

    def draw_beam_stop_cut(self):
        t, s0 = self.Ana.Cut.interruption_times(), self.Run.ev2t(self.Ana.Cut.get_event_range()[0])
        t = concatenate([[[t[0][0], s0]], t[all(t > s0, axis=1)]])  # remove interruptions that are covered by event range
        b = [Draw.box(t0, -5, t1, 200, line_color=2, width=2, fillcolor=2, style=7, opacity=.1) for t0, t1 in t]
        return b + [Draw.legend(b[:1], ['excluded events'], 'lf', margin=.45, w=.3, scale=2)]

    @save_pickle('PulserPulseHeight', suf_args='all')
    def get_pulse_height_trend(self, bin_size=None, cut=None, _redo=False):
        x, y = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var()], cut=self.Cut(cut))
        m, s = self.get_mean_sigma(cut)
        ph_cut = (m - 4 * s < y) & (y < m + 4 * s)
        return self.Draw.profile(x[ph_cut], y[ph_cut], self.Bins.get_time(bin_size, cut), 'Pulser Pulse Height', show=False, graph=True, markersize=.7)

    def draw_pulse_height(self, bin_size=None, cut=None, redo=False, **kwargs):
        """ Shows the average pulse height of the pulser as a function of time """
        g = self.get_pulse_height_trend(bin_size, cut, _redo=redo)
        fit = FitRes(g.Fit('pol0', 'qs'))
        kwargs = prep_kw(kwargs, gridy=True, lm=.14, y_off=1.7, stats=set_statbox(fit=True, entries=True), y_range=ax_range(graph_y(g, err=False), 0, .6, 1))
        g = self.Draw(g, **self.get_t_args(), **kwargs, file_name='PulserPulseHeight')
        return g, fit

    @save_pickle('Disto', suf_args='all')
    def get_distribution(self, name=None, corr=True, beam_on=True, bw=None, _redo=False):
        cut = self.Ana.Cut.get_pulser(beam_on=beam_on)()
        x, (m, s) = self.get_tree_vec(self.get_signal_var(name, corr, cut), cut), self.Pedestal()
        if mean(x) > m + 3 * s:
            x = x[x > m + 2 * s]  # filter out very low signals
        m, s = mean_sigma(x[x < mean(x) + 10], err=False)
        return self.Draw.distribution(x, bins.make(m - 7 * s, m + 10 * s, choose(bw, max(.2, self.Bins.find_width(x)))), 'Pulser Pulse Height', x_tit='Pulse Height [mV]', show=False)

    def draw_distribution(self, name=None, corr=True, beam_on=True, bw=None, redo=False, **kwargs):
        """ Shows the distribution of the pulser integrals. """
        h = self.get_distribution(name, corr, beam_on, bw, _redo=redo)
        rx = ax_range(h=h, thresh=h.GetMaximum() * .02, fl=.3, fh=.6)
        return self.Draw.distribution(h, file_name='PulserDistribution', **prep_kw(kwargs, x_range=rx))

    def get_fit_range(self, h, lsig=2.5, rsig=1):
        f0, same_polarity = fit_fwhm(h), self.Polarity == self.Ana.Polarity
        if f0.get_chi2() > 50:
            return get_fwhm(h, ret_edges=True, err=False)
        return uarr2n([f0[1] + i * f0[2] for i in ([-lsig, rsig] if same_polarity else [-rsig, lsig])])  # fit left tail if same pol and right tail otherwise

    @save_pickle('Fit', suf_args='all')
    def get_distribution_fit(self, corr=True, beam_on=True, bw=None, _redo=False):
        h = self.get_distribution(corr=corr, beam_on=beam_on, bw=bw, _redo=_redo)
        fit = FitRes(h.Fit('gaus', 'qs0', '', *self.get_fit_range(h)))
        return fit

    def draw_distribution_fit(self, corr=True, beam_on=True, bw=None, redo=False, **kwargs):
        fit = self.get_distribution_fit(corr, beam_on, bw, _redo=redo)
        h = self.draw_distribution(corr=corr, beam_on=beam_on, bw=bw, redo=redo, **kwargs)
        f = self.Draw.function(Draw.make_f('gp0', 'gaus', 0, 500, pars=fit.Pars, npx=100, line_style=7), draw_opt='same')
        h.GetListOfFunctions().Add(deepcopy(f))
        self.Draw.function(Draw.make_f('gp1', 'gaus', *self.get_fit_range(h), pars=fit.Pars, npx=100), draw_opt='same')
        format_statbox(h, fit=True, all_stat=True, form='.2f')
        self.Draw.save_plots('PulserDistributionFit', **kwargs)
        return fit

    def draw_pedestal(self, redo=False, **kwargs):
        return self.Ana.Pedestal.draw_distribution(name=self.PedestalName, cut=self.Cut(), redo=redo, **kwargs, prefix='Pulser')

    def draw_pedestal_fit(self, beam_on=True, redo=False, **kwargs):
        kwargs = prep_kw(kwargs, prefix='Pulser')
        return self.Ana.Pedestal.draw_dist_fit(name=self.PedestalName, cut=self.Ana.Cut.get_pulser(beam_on=beam_on)(), redo=redo, **kwargs)

    def compare_pedestal(self, show=True):
        histos = [self.Ana.Pedestal.draw_distribution(show=False), self.draw_pedestal(show=False)]
        self.Draw.stack(histos, 'Comparison of Pulser and Signal Pedestal', ['Signal', 'Pulser'], scale=True, show=show, lw=2, rebin=3)

    def draw_hit_efficiency(self, bin_size=200, show=True):
        from pad.alignment import PadAlignment
        return PadAlignment(self.Run.Converter).draw(bin_size=bin_size, show=show)

    def draw_signal_vs_peaktime(self, bin_size=None, x=None, y=None, show=True):
        return self.Ana.draw_ph_peaktime(x=choose(x, self.Peaks.get_from_tree()), y=choose(y, self.get_values()), xbins=self.get_t_bins(bin_size), show=show)

    def draw_peak_times(self, bin_size=None, x_range=None, y_range=None, draw_ph=False):
        ind = self.get_signal_indices()
        return self.Peaks.draw_signal(bin_size, x_range=x_range, y_range=y_range, x=array(self.Peaks.get_from_tree())[ind], y=array(self.get_values())[ind], draw_ph=draw_ph)

    def draw_cft(self, bin_size=None, show=True):
        h = self.Draw.distribution(self.get_cft(self.get_signal_indices()), self.get_t_bins(bin_size), 'Pulser Constant Fraction Times', x_tit='Constant Fraction Time [ns]', show=show)
        h.Fit('gaus')
        format_statbox(h, fit=True, entries=True)
    # endregion DRAW
    # ----------------------------------------
