# --------------------------------------------------------
#       PULSER ANALYSIS
# created on June 16th 2016 by M. Reichmann
# --------------------------------------------------------

from ROOT import TProfile, gROOT, THStack, TCut
from InfoLegend import InfoLegend
from analysis import *
from binning import make_bins
from numpy import argmax
from waveform import Waveform
from peaks import PeakAnalysis


class PulserAnalysis(Analysis):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.Channel
        self.Tree = self.Ana.Tree
        self.AnaCut = self.Ana.Cut
        self.set_save_directory(self.Ana.SubDir)
        self.Cut = self.AnaCut.get_pulser()
        self.Polarity = self.Ana.PulserPolarity
        self.SignalName = self.load_signal_name()
        self.SignalRegion = array(self.Run.IntegralRegions[self.Ana.DUT.Number - 1]['pulser'])
        self.PedestalName = self.load_pedestal_name()
        self.Type = self.load_type()

        self.Bins = self.Ana.Bins
        self.DUT = self.Ana.DUT
        self.StartEvent = self.Ana.StartEvent
        self.DigitiserBinWidth = self.Ana.DigitiserBinWidth
        self.Waveform = Waveform(self)
        self.Pedestal = self.Ana.Pedestal
        self.PeakName = self.Ana.get_peak_name(type_='pulser')
        self.Peaks = PeakAnalysis(self)
        self.Peaks.Threshold = 14
        self.InfoLegend = InfoLegend(pad_analysis)
        self.set_pickle_sub_dir('Pulser')

    # ----------------------------------------
    # region INIT
    def load_signal_name(self, peak_int=None):
        num = self.Ana.get_signal_number('', peak_int, 'pulser')
        return self.Ana.SignalDefinition.format(pol=self.Polarity, num=num)

    def load_pedestal_name(self, peak_int=None):
        region = self.Config.get('BASIC', 'pulser pedestal') if self.Config.has_option('BASIC', 'pulser pedestal') else 'ac'
        return self.Ana.get_pedestal_name(region, self.Ana.PeakIntegral if None else peak_int)

    def load_type(self):
        return str(self.Ana.Run.RunInfo['pulser']) if 'pulser' in self.Ana.Run.RunInfo else None
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_events(self, cut=None):
        return self.Ana.get_events(cut)

    def get_signal_indices(self):
        """ :returns: indices that are smaller then the smallest signals """
        return where(array(self.get_values()) < self.Ana.get_min_signal() - 5)[0]

    def get_rate(self):
        values = self.Run.get_root_vec(dtype=bool, var='pulser', cut=self.AnaCut.CutStrings.get('beam_interruptions'))
        rate = calc_eff(values=values)
        self.info('pulser rate: {:1.1f} ({:1.1f}) %'.format(rate.n, rate.s))
        return rate

    def get_t_bins(self, bin_size):
        xmin, xmax = self.SignalRegion * self.DigitiserBinWidth
        return make_bins(xmin, xmax, choose(bin_size, default=self.DigitiserBinWidth))

    def get_min_signal(self, name=None):
        h = self.draw_distribution(name, show=False)
        return h.GetBinCenter(h.FindFirstBinAbove(h.GetMaximum() * .01))

    def get_pulse_height(self, corr=True, bin_width=.2, redo=False):
        pickle_path = self.make_simple_pickle_path('HistoFit', '{}_{}'.format(int(corr), 'BeamOn'))
        fit = do_pickle(pickle_path, partial(self.draw_distribution_fit, show=False, prnt=False, corr=corr, redo=redo, bin_width=bin_width), redo=redo)
        return make_ufloat(fit, par=1)

    def get_pedestal(self, par=1, redo=False):
        pickle_path = self.make_simple_pickle_path('Pedestal')
        fit = do_pickle(pickle_path, partial(self.draw_pedestal, show=False, prnt=False, redo=redo))
        return make_ufloat(fit, par=par)

    def get_pedestal_mean(self, redo=False):
        return self.get_pedestal(par=1, redo=redo)

    def get_pedestal_sigma(self, redo=False):
        return self.get_pedestal(par=2, redo=redo)

    def get_values(self, cut=None, redo=False):
        cut = self.Cut(cut)
        return do_hdf5(self.make_simple_hdf5_path('V', cut.GetName()), self.Run.get_root_vec, redo, var=self.generate_signal_name(cut=cut), cut=cut, dtype='f2')

    def get_cft(self, ind=None):
        cft = self.Peaks.get_cft(ind=ind)
        tmin, tmax = self.SignalRegion * self.DigitiserBinWidth
        cft = array([[t for t in lst if tmin <= t <= tmax] for lst in cft])
        return array([lst[0] if len(lst) else 0 for lst in cft])

    def get_ms(self, cut=None):
        values = array(self.get_values(cut))
        h = self.draw_disto(values, '', self.Ana.Bins.get_pad_ph(bin_width=.5), show=False)
        mpv = get_hist_args(h)[argmax(get_hist_vec(h))].n  # find MPV
        fit = h.Fit('gaus', 'qs0', '', mpv - 10, mpv + 10)
        return fit.Parameter(1), fit.Parameter(2)

    def get_fraction(self, show=False, prnt=True):
        """ :returns the fitted value of the fraction of pulser events with event range and beam interruptions cuts and its fit error. """
        cut = self.AnaCut.generate_custom(include=['beam_interruptions'], prnt=prnt)
        self.format_statbox(only_fit=True, x=.9, w=.2)
        h = self.draw_rate(show=show, cut=cut, prnt=prnt)
        format_histo(h, 'Fit Result', markersize=None)
        fit = h.Fit('pol0', 'qs')
        self.info('The fraction of pulser events is: {0:5.2f} +- {1:4.2f} %'.format(fit.Parameter(0), fit.ParError(0)))
        return FitRes(fit)

    def get_real_fraction(self):
        """ :returns the real fraction of pulser events compared to signal events. (pulser in rate / flux / area)  """
        in_rate = 40 if self.Ana.get_flux().n < 10 else 100
        diamond_size = make_ufloat((.4, .05)) ** 2
        particle_rate = self.Ana.get_flux() * diamond_size
        return in_rate / particle_rate

    def get_signal_name(self, peak_integral=None):
        return self.load_signal_name(peak_int=peak_integral)

    def get_all_signal_names(self):
        return self.Ana.get_all_signal_names('pulser')

    def generate_signal_name(self, name=None, event_corr=False, off_corr=True, cut=None):
        return self.Ana.generate_signal_name(choose(name, self.SignalName), evnt_corr=event_corr, off_corr=off_corr, cut=self.Cut(cut))

    def get_attenuator(self):
        return self.Ana.get_attenuator()

    def get_irradiation(self):
        return self.Ana.get_irradiation()
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_rate(self, evts_per_bin=1000, cut=None, vs_time=True, rel_t=True, show=True, prnt=True):
        """ Shows the fraction of pulser events as a function of the event number. Peaks appearing in this graph are most likely beam interruptions. """
        cut = '' if cut is None else TCut(cut)
        set_root_output(False)
        h = TProfile('hpr', 'Pulser Rate', *self.Ana.Bins.get_raw(evts_per_bin, vs_time=vs_time, t_from_event=True))
        self.Tree.Draw('pulser*100:{v}>>hpr'.format(v=self.Ana.get_t_var() if vs_time else 'Entry$'), cut, 'goff')
        format_histo(h, x_tit='Time [hh:mm]' if vs_time else 'Event Number', y_tit='Pulser Fraction [%]', y_off=.8, fill_color=self.FillColor, y_range=[0, 105], markersize=.7, stats=0,
                     t_ax_off=self.Ana.Run.StartTime if rel_t else 0)
        self.save_histo(h, 'PulserRate', show, lm=.08, draw_opt='bare', x=1.5, y=.75, prnt=prnt)
        return h

    def draw_pulse_height(self, bin_size=None, cut=None, y_range=None, show=True, redo=False):
        """ Shows the average pulse height of the pulser as a function of time """
        def f():
            values, times = array(self.get_values(cut)), self.Run.get_root_vec(var=self.Ana.get_t_var(), cut=self.Cut(cut))
            m, s = self.get_ms(cut)
            p0 = TProfile('ppph', 'Pulser Pulse Height', *self.Ana.Bins.get_time(bin_size))
            cut_ind = (m - 4 * s < values) & (values < m + 4 * s)
            fill_hist(p0, x=times[cut_ind], y=values[cut_ind])
            return p0
        p = do_pickle(self.make_simple_pickle_path('PHT', self.Cut(cut).GetName()), f, redo=redo)
        fit = p.Fit('pol0', 'qs')
        ph = get_hist_vec(p, err=False)
        y_range = choose(y_range, increased_range([min(ph[ph > 0]), max(ph)], .3, .7))
        format_histo(p, x_tit='Time [hh:mm]', y_tit='Pulse Height [mV]', y_off=1.7, fill_color=self.FillColor, t_ax_off=self.Run.StartTime, y_range=y_range)
        self.save_histo(p, 'PulserPulserHeight{}'.format(bin_size), show, gridy=True, lm=.14)
        return p, FitRes(fit)

    def draw_distribution(self, name=None, corr=True, beam_on=True, bin_width=.2, show=True, redo=False):
        """ Shows the distribution of the pulser integrals. """
        def f():
            cut = self.AnaCut.get_pulser(beam_on=beam_on)()
            var = self.generate_signal_name(name, event_corr=False, off_corr=corr, cut=cut)
            values = self.Run.get_root_vec(var=var, cut=cut)
            m, s = mean_sigma(values[values < mean(values) + 10])
            return self.draw_disto(values, 'Pulser Pulse Height', make_bins(m - 3 * s, m + 5 * s, bin_width), x_tit='Pulse Height [mV]', show=False)
        suf = '{corr}_{beam}_{}'.format(self.get_all_signal_names()[choose(name, self.SignalName)], corr='ped_corr' if corr else '', beam='BeamOff' if not beam_on else 'BeamOn')
        h = do_pickle(self.make_simple_pickle_path('Disto', suf), f, redo=redo)
        self.draw_histo(h, show, lm=.12)
        return h

    def draw_distribution_fit(self, show=True, redo=False, corr=True, beam_on=True, bin_width=.2, prnt=True):
        suffix = '{corr}_{beam}'.format(corr='ped_corr' if corr else '', beam='BeamOff' if not beam_on else 'BeamOn')
        pickle_path = self.make_pickle_path('Pulser', 'HistoFit', self.Run.Number, self.DUT.Number, suf=suffix)
        self.format_statbox(only_fit=True, w=.3, h=.2)
        h = self.draw_distribution(show=show, corr=corr, beam_on=beam_on, bin_width=bin_width, redo=redo)
        h.SetName('Fit Result')
        same_pols = self.Polarity == self.Ana.Polarity

        def f():
            xmax = h.GetBinCenter(h.GetMaximumBin())
            full_fit = h.Fit('gaus', 'qs0', '', xmax - 5, xmax + 5)
            xmin, xmax = [full_fit.Parameter(1) + i * full_fit.Parameter(2) for i in ([-2, .5] if same_pols else [-.5, 2])]
            fit_func = h.Fit('gaus', 'qs{0}'.format('' if show else '0'), '', xmin, xmax)
            return FitRes(fit_func)

        fit = do_pickle(pickle_path, f, redo=redo)
        f2 = deepcopy(gROOT.GetFunction('gaus'))
        f2.SetLineStyle(7)
        f2.SetRange(0, 500)
        h.GetListOfFunctions().Add(f2)
        set_drawing_range(h)
        self.save_plots('PulserDistributionFit', show=show, prnt=prnt)
        self.Ana.server_pickle(pickle_path, fit)
        return fit

    def draw_pedestal(self, show=True, save=True, prnt=True, redo=False):
        return self.Ana.Pedestal.draw_disto(name=self.PedestalName, cut=self.Cut(), show=show, save=save, prnt=prnt, redo=redo)

    def draw_pedestal_fit(self, show=True, save=True, prnt=True, redo=False):
        return self.Ana.Pedestal.draw_disto_fit(name=self.PedestalName, cut=self.Cut(), show=show, save=save, prnt=prnt, redo=redo)

    def compare_pedestal(self, show=True):
        h1, h2 = self.Ana.Pedestal.draw_disto(show=False), self.draw_pedestal(show=False)
        legend = self.make_legend(.7)
        names = ['Signal', 'Pulser']
        stack = THStack('spc', 'Comparison of Pulser and Signal Pedestal')
        for i, h in enumerate([h1, h2]):
            h.Scale(1 / h.GetMaximum())
            color = get_color_gradient(2)[i]
            format_histo(h, color=color, fill_color=color, lw=2, stats=0, y_range=[0, 1.1], opacity=.5)
            stack.Add(h)
            legend.AddEntry(h, names[i], 'l')
        format_histo(stack, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.3, draw_first=True)
        self.save_histo(stack, 'PulserPedestalComparison', show, lm=.12, leg=legend, draw_opt='nostack')
        self.reset_colors()

    def draw_hit_efficiency(self, xbins=200, show=True):
        xbins = self.Ana.Bins.get_pulser(xbins) if type(xbins) is int else xbins
        p = TProfile('pa{}'.format(self.Run.Number), 'Hit Efficiency at Pulser Events', *xbins)
        self.Ana.Tree.Draw('(@col.size()>1)*100:Entry$>>pa{}'.format(self.Run.Number), 'pulser', 'goff')
        format_histo(p, x_tit='Event Number', y_tit='Hit Efficiency [%]', y_off=1.3, stats=0, y_range=[0, 105], fill_color=self.FillColor)
        self.save_histo(p, 'PulserHitEfficiency', show, self.Ana.TelSaveDir, draw_opt='hist', prnt=show, rm=.08)
        return p

    def draw_signal_vs_peaktime(self, bin_size=None, x=None, y=None, show=True):
        return self.Ana.draw_signal_vs_peaktime(x=choose(x, self.Peaks.get_from_tree()), y=choose(y, self.get_values()), xbins=self.get_t_bins(bin_size), show=show)

    def draw_peak_times(self, bin_size=None, x_range=None, y_range=None, draw_ph=False):
        ind = self.get_signal_indices()
        return self.Peaks.draw_signal(bin_size, x_range=x_range, y_range=y_range, x=array(self.Peaks.get_from_tree())[ind], y=array(self.get_values())[ind], draw_ph=draw_ph)

    def draw_cft(self, bin_size=None, show=True):
        h = self.draw_disto(self.get_cft(self.get_signal_indices()), 'Pulser Constant Fraction Times', self.get_t_bins(bin_size), x_tit='Constant Fraction Time [ns]', show=show)
        self.format_statbox(fit=True, entries=True, h=.2)
        h.Fit('gaus')
        update_canvas()
    # endregion DRAW
    # ----------------------------------------
