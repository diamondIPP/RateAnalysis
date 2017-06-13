#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import TH1F, TCut
from Utils import set_statbox, fit_poissoni
from math import pi


class PeakAnalysis(Elementary):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.Run = self.Ana.run
        self.Channel = self.Ana.channel
        self.Tree = self.Ana.tree
        self.AllCut = self.Ana.AllCuts
        self.save_dir = self.Ana.save_dir

    def draw_n_peaks(self, spec=False, show=True, fit=False):
        h = TH1F('h_pn', 'Number of Peaks', 10, 0, 10)
        draw_var = '@peaks{ch}_x.size()>>h_pn' if spec and self.Run.has_branch('peaks{ch}_x'.format(ch=self.Channel)) else 'n_peaks[{ch}] - 1>>h_pn'
        self.Tree.Draw(draw_var.format(ch=self.Channel), self.AllCut, 'goff')
        set_statbox(only_fit=True, entries=4, w=.3) if fit else set_statbox(only_entries=True)
        self.format_histo(h, x_tit='Number of Peaks', y_tit='Number of Entries', y_off=1.4, fill_color=836, lw=2, fill_style=3004)
        self.save_histo(h, 'PeakNumbers', show, logy=True, lm=.11)
        return h

    def draw_n_peaks_fit(self, show=True, spec=False):
        h = self.draw_n_peaks(spec, show=False, fit=True)
        fit = fit_poissoni(h, show=show)
        self.format_histo(h, 'Fit Result')
        self.save_histo(h, 'PeakNumbersFit', show, logy=True, lm=.11)
        self.get_flux(fit=fit)
        return fit

