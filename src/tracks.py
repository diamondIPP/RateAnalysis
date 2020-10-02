#!/usr/bin/env python
# --------------------------------------------------------
#       analysis class for the telescope tracks
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TProfile
from analysis import Analysis, Draw, format_histo, ufloat, calc_eff
from helpers.draw import format_statbox


class Tracks(Analysis):

    def __init__(self, analysis):

        super().__init__(analysis.TCString, results_dir=analysis.Run.Number, pickle_dir='Tracks')
        self.Ana = analysis
        self.Tree = self.Ana.Tree
        self.Bins = self.Ana.Bins

    def get_efficiency(self):
        return calc_eff(values=self.Ana.get_root_vec(var='n_tracks'))

    def draw_n(self, show=True):
        p = TProfile('ptt', 'Number of Tracks vs. Time', *self.Bins.get_raw_time(10))
        self.Tree.Draw('n_tracks:{}>>ptt'.format(self.Ana.get_t_var()), '', 'goff')
        format_statbox(all_stat=True)
        format_histo(p, x_tit='Time [hh:mm}', y_tit='Number of Tracks', y_off=1.3, t_ax_off=0, fill_color=Draw.FillColor)
        self.Draw(p, 'NTracksTime', draw_opt='hist', show=show)
