#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import TH1F


class PeakAnalysis(Elementary):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.Run = self.Ana.run
        self.Channel = self.Ana.channel
        self.Tree = self.Ana.tree
        self.AllCut = self.Ana.AllCuts
        self.save_dir = self.Ana.save_dir

    def draw_n_peaks(self, show=True, p1=0.7, p2=1):
        h = TH1F('h_pn', 'Number of Peaks', 21, -.5, 20.5)
        h1 = TH1F('h_pn1', 'Number of Peaks', 21, -.5, 20.5)
        self.Tree.Draw('@peaks{ch}_x.size()>>h_pn'.format(ch=self.Channel), self.AllCut, 'goff')
        self.format_histo(h, x_tit='number of peaks', y_tit='number of entries', y_off=1.5, fill_color=836, lw=2, fill_style=3004)
        self.save_histo(h, 'PeakNumbers', show, logy=True)
        # while h1.GetBinContent(2) != h.GetBinContent(2):
        #     # h1.Fill(gRandom.Poisson(24 * self.get_flux() / 5e4 * .5 * .5 * p2) + gRandom.Binomial(1, p1))
        #     # h1.Fill(gRandom.Poisson(24 * self.get_flux() / 5e4 * .5 * .5 * p2) + 1)
        #     h1.Fill(gRandom.Poisson(24 * f / 5e4 * .5 * .5 * p2) + 1)
        # self.format_histo(h1, x_tit='number of peaks', y_tit='Number of Entries', y_off=1.5, fill_color=896, lw=2)
        # h1.SetFillStyle(3005)
        # h1.Draw('same')
        # self.histos.append(h1)

