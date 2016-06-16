# --------------------------------------------------------
#       PULSER ANALYSIS
# created on June 16th 2016 by M. Reichmann
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import TProfile
from Utils import *


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
        # todo beamon nicht vergessen!
        """ Shows the distribution of the pulser integrals. """
        h = self.Ana.show_signal_histo(cut=self.PulserCut, sig=self.Ana.PulserName, show=False, off_corr=corr, evnt_corr=False, binning=binning, events=events, start=start)
        self.format_histo(h, stats=0, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.3)
        self.save_histo(h, 'PulserDistribution', show, logy=True, lm=.12)
        return h
