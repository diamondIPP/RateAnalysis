# --------------------------------------------------------
#       PULSER ANALYSIS
# created on June 16th 2016 by M. Reichmann
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import TProfile, gROOT, THStack, TCut, TH2F
from Utils import *
from copy import deepcopy
from numpy import mean
from InfoLegend import InfoLegend


class PulserAnalysis(Elementary):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.channel
        self.Tree = self.Ana.tree
        self.Cut = self.Ana.Cut
        self.save_dir = self.Ana.save_dir
        self.PulserCut = self.Cut.generate_pulser_cut()
        self.Polarity = self.Ana.PulserPolarity
        self.SignalName = self.get_signal_name()
        self.PedestalName = self.get_pedestal_name()
        self.Type = self.load_type()

        self.DiamondName = self.Ana.DiamondName
        self.DiamondNumber = self.Ana.DiamondNumber
        self.RunNumber = self.Ana.RunNumber
        self.InfoLegend = InfoLegend(pad_analysis)

    def get_signal_name(self, peak_int=None):
        num = self.Ana.get_signal_number('', peak_int, 'pulser')
        return self.Ana.SignalDefinition.format(pol=self.Polarity, num=num)

    def get_pedestal_name(self, peak_int=None):
        region = self.ana_config_parser.get('BASIC', 'pulser_pedestal') if self.ana_config_parser.has_option('BASIC', 'pulser_pedestal') else 'ac'
        return self.Ana.get_pedestal_name(region, self.Ana.PeakIntegral if None else peak_int)

    def load_type(self):
        return str(self.Ana.RunInfo['pulser']) if 'pulser' in self.Ana.RunInfo else None

    def draw_rate(self, evts_per_bin=1000, cut=None, vs_time=True, rel_t=True, show=True, prnt=True):
        """ Shows the fraction of pulser events as a function of the event number. Peaks appearing in this graph are most likely beam interruptions. """
        cut = '' if cut is None else TCut(cut)
        set_root_output(False)
        h = TProfile('hpr', 'Pulser Rate', *self.Ana.Plots.get_binning(evts_per_bin, time_bins=vs_time))
        self.Tree.Draw('pulser*100:{v}>>hpr'.format(v='time / 1000.' if vs_time else 'Entry$'), cut, 'goff')
        self.format_histo(h, x_tit='Time [hh:mm]' if vs_time else 'Event Number', y_tit='Pulser Fraction [%]', y_off=.8, fill_color=self.FillColor, y_range=[0, 105], markersize=.7, stats=0,
                          t_ax_off=self.Ana.Run.StartTime if rel_t else 0)
        self.save_histo(h, 'PulserRate', show, lm=.08, draw_opt='bare', x=1.5, y=.75, prnt=prnt)
        return h

    def calc_fraction(self, show=False, prnt=True):
        """ :returns the fitted value of the fraction of pulser events with event range and beam interruptions cuts and its fit error. """
        cut = self.Cut.generate_special_cut(included=['beam_interruptions'], prnt=prnt)
        self.set_statbox(only_fit=True, n_entries=2, x=.9, w=.2)
        h = self.draw_rate(show=show, cut=cut, prnt=prnt)
        self.format_histo(h, 'Fit Result', markersize=None)
        fit = h.Fit('pol0', 'qs')
        self.log_info('The fraction of pulser events is: {0:5.2f} +- {1:4.2f} %'.format(fit.Parameter(0), fit.ParError(0)), prnt=prnt)
        return FitRes(fit)

    def calc_real_fraction(self):
        in_rate = 40 if self.Ana.get_flux().n < 10 else 100
        diamond_size = make_ufloat((.4, .05)) ** 2
        particle_rate = self.Ana.get_flux() * diamond_size
        return in_rate / particle_rate

    def draw_pulse_height(self, bin_size=10000, y_range=None, show=True, redo=False):
        """ Shows the average pulse height of the pulser as a function of time """
        pickle_path = self.make_pickle_path('Pulser', 'PH', self.RunNumber, self.Ana.DiamondNumber, bin_size)
        print self.RunNumber

        def f():
            gr = self.make_tgrapherrors('gpph', 'Pulser Pulse Height')
            h = TH2F('temp', '', *[v for info in [self.Ana.get_time_bins(bin_size), self.Ana.Plots.get_ph_bins(bin_width=1)] for v in info])
            signal = self.Ana.generate_signal_name(self.SignalName, evnt_corr=False, off_corr=True, cut=self.PulserCut)
            self.Tree.Draw('{sig}:time/1000.>>temp'.format(sig=signal), self.PulserCut, 'goff')
            for xbin in xrange(2, h.GetNbinsX() + 1):  # first bin is always empty
                py = h.ProjectionY('_py{}'.format(xbin), xbin, xbin)
                m = py.GetBinCenter(py.GetMaximumBin())
                f1 = py.Fit('gaus', 'qs0', '', m - 10, m + 10)
                m, s = (f1.Parameter(j) for j in [1, 2])
                f2 = py.Fit('gaus', 'qs0', '', m - 4 * s, m + 1.5 * s)
                gr.SetPoint(xbin - 2, h.GetXaxis().GetBinCenter(xbin), f2.Parameter(1))
                gr.SetPointError(xbin - 2, h.GetXaxis().GetBinWidth(xbin) / 2., f2.ParError(1))
            return gr

        g = do_pickle(pickle_path, f, redo=redo)
        self.set_statbox(n_entries=2, only_fit=True, form='1.2f')
        fit_res = g.Fit('pol0', 'qs')
        values = [g.GetY()[i] for i in xrange(g.GetN()) if g.GetY()[i]]
        y_range = increased_range([min(values), max(values)], .7, .7) if y_range is None else y_range
        self.format_histo(g, x_tit='Time [hh:mm]', y_tit='Pulse Height [au]', y_off=1.7, fill_color=self.FillColor, y_range=y_range, t_ax_off=self.Ana.Run.StartTime)
        self.save_histo(g, 'PulserPulserHeight{}'.format(bin_size), show, gridy=True, draw_opt='ap', lm=.14, save=not file_exists(pickle_path))
        return g, FitRes(fit_res)

    def find_range(self, corr):
        n, i = 0, 0
        while not n:
            n = self.Ana.tree.Draw(self.Ana.generate_signal_name(self.SignalName, off_corr=corr, evnt_corr=False, cut=self.PulserCut), self.PulserCut, 'goff', 10000, self.Ana.StartEvent + i)
            i += 10000
        values = sorted([self.Ana.tree.GetV1()[i] for i in xrange(n)])
        ran = mean(values[:5]), mean(values[-5:])
        return increased_range(ran, .1, .3)

    def draw_distribution(self, show=True, corr=True, beam_on=True, bin_width=10., events=None, start=None, stats=False, redo=False, prnt=True):
        """ Shows the distribution of the pulser integrals. """
        cut = self.Cut.generate_pulser_cut(beam_on=beam_on)
        x = self.find_range(corr)
        h = self.Ana.draw_signal_distribution(cut=cut, sig=self.SignalName, show=False, off_corr=corr, evnt_corr=False, bin_width=bin_width, events=events,
                                              start=start, redo=redo, x_range=x, prnt=prnt, save=False)
        self.format_histo(h, name='p_hd', stats=stats, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.save_histo(h, 'PulserDistribution', show, lm=.12, prnt=prnt)
        return h

    def draw_distribution_fit(self, show=True, redo=False, corr=True, beam_on=True, events=None, start=None, bin_width=.1, prnt=True):
        start_string = '_{0}'.format(start) if start is not None else ''
        events_string = '_{0}'.format(events) if events is not None else ''
        suffix = '{corr}_{beam}{st}{ev}'.format(corr='ped_corr' if corr else '', beam='BeamOff' if not beam_on else 'BeamOn', st=start_string, ev=events_string)
        pickle_path = self.make_pickle_path('Pulser', 'HistoFit', self.RunNumber, self.DiamondNumber, suf=suffix)
        set_statbox(.95, .88, entries=4, only_fit=True, w=.5)
        h = self.draw_distribution(show=show, corr=corr, beam_on=beam_on, bin_width=bin_width, events=events, start=start, stats=True, redo=redo, prnt=prnt)
        h.SetName('Fit Result')
        same_pols = self.Polarity == self.Ana.Polarity
        xmax = h.GetBinCenter(h.GetMaximumBin())
        full_fit = h.Fit('gaus', 'qs0', '', xmax - 5, xmax + 5)
        xmin, xmax = [full_fit.Parameter(1) + i * full_fit.Parameter(2) for i in ([-2, .5] if same_pols else [-.5, 2])]
        fit_func = h.Fit('gaus', 'qs{0}'.format('' if show else '0'), '', xmin, xmax)
        f = deepcopy(gROOT.GetFunction('gaus'))
        f.SetLineStyle(7)
        f.SetRange(0, 500)
        h.GetListOfFunctions().Add(f)
        set_drawing_range(h)
        self.save_plots('PulserDistributionFit', show=show, prnt=prnt)
        fit = FitRes(fit_func)
        server_pickle(pickle_path, fit)
        return fit

    def draw_peak_timing(self, show=True, corr=False):
        self.Ana.draw_peak_timing('', 'pulser', cut=self.PulserCut, show=show, draw_cut=False, corr=corr)

    def draw_pedestal(self, show=True, save=True, prnt=True, redo=False):
        return self.Ana.Pedestal.draw_disto_fit(name=self.PedestalName, cut=self.PulserCut, show=show, save=save, prnt=prnt, redo=redo)

    def compare_pedestal(self, show=True):
        self.Ana.draw_pedestal_disto_fit(show=False)
        h1 = deepcopy(self.Ana.PedestalHisto)
        self.draw_pedestal(show=False)
        h2 = self.Ana.PedestalHisto
        legend = self.make_legend(.7)
        names = ['Signal', 'Pulser']
        stack = THStack('spc', 'Comparison of Pulser and Signal Pedestal')
        means = []
        for i, h in enumerate([h1, h2]):
            h.Scale(1 / h.GetMaximum())
            self.format_histo(h, color=self.get_color(), lw=2, stats=0, y_range=[0, 1.1])
            stack.Add(h)
            legend.AddEntry(h, names[i], 'l')
            fit = h.Fit('gaus', 'qs0')
            means.append(fit.Parameter(1))
        self.format_histo(stack, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.3, draw_first=True)
        self.save_histo(stack, 'PulserPedestalComparison', show, lm=.12, l=legend, draw_opt='nostack')
        self.reset_colors()

    def draw_waveforms(self, n=1, ch=None, start_event=None, show=True, t_corr=True):
        self.Ana.Waveform.draw(n=n, start_event=start_event, cut=self.PulserCut, show=show, t_corr=t_corr, channel=ch)

    def save_pulser_shapes(self, n_pics=5, show=True):
        events_spacing = (self.Ana.EndEvent - self.Ana.StartEvent) / n_pics
        start_events = [self.Ana.StartEvent + events_spacing * i for i in xrange(n_pics)]
        for i, start in enumerate(start_events):
            self.count = 0
            self.draw_waveforms(n=1000, start_event=start, show=show)
            self.save_plots('WaveForms{0}'.format(i), sub_dir='{0}/WaveForms'.format(self.save_dir))

    def save_felix(self):
        self.save_dir = self.Ana.save_dir
        self.set_save_directory('PlotsFelix')

    def draw_hit_efficiency(self, xbins=200, show=True):
        xbins = self.Ana.Plots.get_pulser_bins(xbins) if type(xbins) is int else xbins
        p = TProfile('pa{}'.format(self.RunNumber), 'Hit Efficiency at Pulser Events', *xbins)
        self.Ana.tree.Draw('(@col.size()>1)*100:Entry$>>pa{}'.format(self.RunNumber), 'pulser', 'goff')
        self.format_histo(p, x_tit='Event Number', y_tit='Hit Efficiency [%]', y_off=1.3, stats=0, y_range=[0, 105], fill_color=self.FillColor)
        self.save_histo(p, 'PulserHitEfficiency', show, self.Ana.TelSaveDir, draw_opt='hist', prnt=show, rm=.08)
        return p
