# --------------------------------------------------------
#       PULSER ANALYSIS
# created on June 16th 2016 by M. Reichmann
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import TProfile, gROOT, THStack
from Utils import *
from copy import deepcopy


class PulserAnalysis(Elementary):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        # self.Ana = PadAnalysis(5, 5)
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.Run = self.Ana.run
        self.Channel = self.Ana.channel
        self.Tree = self.Ana.tree
        self.Cut = self.Ana.Cut
        self.save_dir = self.Ana.save_dir
        self.PulserCut = self.Cut.generate_pulser_cut()
        self.Polarity = self.Ana.PulserPolarity

        self.ROOTObjects = []

    def draw_rate(self, evts_per_bin=1000, cut=None, show=True):
        """ Shows the fraction of pulser events as a function of the event number. Peaks appearing in this graph are most likely beam interruptions. """
        cut = '' if cut is None else cut
        entries = self.Run.n_entries
        nbins = entries / evts_per_bin
        h = TProfile('hpr', 'Pulser Rate', nbins, 0, entries)
        self.Tree.Draw('pulser*100:Entry$>>hpr', cut, 'goff')
        self.format_histo(h, x_tit='Event Number', y_tit='Pulser Fraction [%]', y_off=1.3, fill_color=self.FillColor, stats=0)
        self.save_histo(h, 'PulserRate', show, lm=.13, draw_opt='hist')
        return h

    def calc_fraction(self):
        """ :returns the fitted value of the fraction of pulser events with event range and beam interruptions cuts and its fit error. """
        cut = self.Cut.generate_special_cut(included_cuts=['beam_interruptions', 'event_range'])
        h = self.draw_rate(show=False, cut=cut)
        fit = h.Fit('pol0', 'qs0')
        print 'The fraction of pulser events is: {0:5.2f} +- {1:4.2f}%'.format(fit.Parameter(0), fit.ParError(0))
        return fit.Parameter(0), fit.ParError(0)

    def draw_pulseheight(self, binning=20000, draw_opt='histe', show=True):
        """ Shows the average pulse height of the pulser as a function of event numbers. """
        entries = self.Run.n_entries
        nbins = entries / binning
        h = TProfile('hpph', 'Pulser Pulse Height', nbins, 0, entries)
        signal = self.Ana.generate_signal_name(self.Ana.PulserName, evnt_corr=False, off_corr=True, cut=self.PulserCut)
        self.Tree.Draw('{sig}:Entry$>>hpph'.format(sig=signal), self.PulserCut, 'goff')
        values = [h.GetBinContent(i) for i in xrange(h.GetNbinsX()) if h.GetBinContent(i)]
        diff = max(values) - min(values)
        self.format_histo(h, x_tit='Event Number', y_tit='Pulse Height [au]', y_off=1.7, stats=0, y_range=[min(values) - diff, max(values) + diff * 2], fill_color=self.FillColor)
        self.draw_histo(h, '', show, gridy=True, draw_opt=draw_opt, lm=.15)
        h.Draw('same')
        self.save_plots('PulserPulserHeight', self.save_dir)
        return h

    def draw_pulseheight_fit(self, show=True):
        """ Shows the pulse height fit for the pulser. """
        set_statbox(.88, .88, only_fit=True)
        h = self.draw_pulseheight(show=show)
        h.SetName('Fit Result')
        h.SetStats(1)
        fit = h.Fit('pol0', 'qs')
        self.save_plots('PulserPulserHeightFit')
        return fit.Parameter(0), fit.ParError(0)

    def draw_distribution(self, show=True, corr=True, beam_on=True, binning=700, events=None, start=None):
        """ Shows the distribution of the pulser integrals. """
        cut = self.Cut.generate_pulser_cut(beam_on=beam_on)
        h = self.Ana.show_signal_histo(cut=cut, sig=self.Ana.PulserName, show=False, off_corr=corr, evnt_corr=False, binning=binning, events=events, start=start)
        self.format_histo(h, stats=0, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.3)
        self.save_histo(h, 'PulserDistribution', show, logy=True, lm=.12)
        return h

    def draw_distribution_fit(self, show=True, corr=True, beam_on=True, events=None, start=None, binning=350):
        start_string = '_{0}'.format(start) if start is not None else ''
        events_string = '_{0}'.format(events) if events is not None else ''
        suffix = '{corr}_{beam}{st}{ev}'.format(corr='_ped_corr' if corr else '', beam='BeamOff' if not beam_on else 'BeamOn', st=start_string, ev=events_string)
        pickle_path = self.Ana.PickleDir + 'Pulser/HistoFit_{tc}_{run}_{dia}{suf}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.Ana.run_number, dia=self.Ana.diamond_name, suf=suffix)

        def func():
            set_statbox(.88, .88, only_fit=True)
            h = self.draw_distribution(show=show, corr=corr, beam_on=beam_on, binning=binning, events=events, start=start)
            h.SetStats(1)
            same_pols = self.Polarity == self.Ana.Polarity
            h.GetXaxis().SetRangeUser(20, h.GetXaxis().GetXmax())
            x_min = 10 if same_pols else h.GetBinCenter(h.GetMaximumBin() - 2)
            h.GetXaxis().UnZoom()
            x_max = h.GetBinCenter(h.GetMaximumBin() + 2) if same_pols else h.GetBinCenter(h.GetNbinsX() - 1)
            fit_func = h.Fit('gaus', 'qs{0}'.format('' if show else '0'), '', x_min, x_max)
            f = deepcopy(gROOT.GetFunction('gaus'))
            f.SetLineStyle(7)
            f.SetRange(0, 500)
            f.Draw('same')
            h.SetName('Fit Results')
            self.Ana.RootObjects.append(f)
            self.save_plots('PulserDistributionFit')
            return fit_func

        fit = func() if show else None
        return self.do_pickle(pickle_path, func, fit)

    def draw_peak_timing(self, show=True):
        self.Ana.draw_peak_timing('', 'pulser', ucut=self.PulserCut, show=show, draw_cut=False)

    def draw_pedestal(self, show=True, fit=True, x_range=None):
        return self.Ana.show_pedestal_histo(cut=self.PulserCut, show=show, fit=fit, x_range=x_range)

    def compare_pedestal(self, show=True):
        self.Ana.show_pedestal_histo(show=False, fit=False, x_range=[-20, 20])
        h1 = deepcopy(self.Ana.PedestalHisto)
        self.draw_pedestal(show=False, fit=False, x_range=[-20, 20])
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

    def draw_waveform(self, n=1, start_event=None, add_buckets=False, y_range=None, show=True):
        fit = self.draw_distribution_fit(show=False)
        if y_range is None:
            y_min = fit.Parameter(1) / 50 * 50 + 50 if self.Polarity < 0 else -150
            y_max = +150 if self.Polarity < 0 else fit.Parameter(1) * 50 / 50 + 50
            print y_min, y_max
        else:
            y_min, y_max = y_range[0], y_range[1]
        self.Ana.draw_waveforms(n=n, start_event=start_event, add_buckets=add_buckets, cut_string=self.PulserCut, fixed_range=[y_min, y_max], show=show)

    def save_pulser_shapes(self, n_pics=5, fixed_range=None, show=True):
        frange = [-100, 50] if fixed_range is None else fixed_range
        events_spacing = (self.Ana.EndEvent - self.Ana.StartEvent) / n_pics
        start_events = [self.Ana.StartEvent + events_spacing * i for i in xrange(n_pics)]
        for i, start in enumerate(start_events):
            self.count = 0
            self.draw_waveform(n=1000, start_event=start, y_range=frange, show=show)
            self.save_plots('WaveForms{0}'.format(i), sub_dir='{0}/WaveForms'.format(self.save_dir))

    def draw_pulser_vs_time(self, n_points=5, _mean=True, show=True, corr=True, events=5000):
        events_spacing = (self.EndEvent - self.StartEvent) / n_points
        start_events = [self.StartEvent + events_spacing * i for i in xrange(n_points)]
        mode = 'Mean' if _mean else 'Sigma'
        gr = self.make_tgrapherrors('gr', '{0} of Pulser vs Time'.format(mode))
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        for i, start in enumerate(start_events):
            fit = self.calc_pulser_fit(show=False, start=start, events=events, binning=200, corr=corr)
            par = 1 if _mean else 2
            gr.SetPoint(i, (self.run.get_time_at_event(start) - self.run.startTime) / 60e3, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
        gROOT.SetBatch(0) if show else gROOT.SetBatch(1)
        c = TCanvas('c', '{0} of Pulser vs Time'.format(mode), 1000, 1000)
        gr.Draw('alp')
        self.save_plots('Pulser{0}VsTime'.format(mode), sub_dir=self.save_dir)
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.histos.append([gr, c])
