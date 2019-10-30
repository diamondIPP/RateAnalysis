#!/usr/bin/env python
# --------------------------------------------------------
#       parent class for the analysis of a single device under test
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from telescope_analysis import *
from CurrentInfo import Currents


class DUTAnalyis(TelecopeAnalysis):
    def __init__(self, run_number, diamond_nr, test_campaign=None, tree=None, verbose=False):

        TelecopeAnalysis.__init__(self, run_number, test_campaign, tree=tree, verbose=verbose)

        self.DiamondNumber = self.load_diamond_nr(diamond_nr)
        self.DiamondName = self.Run.DiamondNames[diamond_nr - 1]
        self.Bias = self.Run.Bias[diamond_nr - 1]

        self.update_config()
        self.set_save_directory(join(self.DiamondName, str(self.RunNumber).zfill(3)))

        self.Currents = Currents(self)

    def update_config(self):
        pass

    def load_diamond_nr(self, diamond_nr):
        if diamond_nr not in self.Run.DiamondNumbers:
            critical('wrong diamond number "{}". The following diamond numbers are valid: {}'.format(diamond_nr, self.Run.DiamondNumbers))
        return diamond_nr

    def draw_current(self, relative_time=False, averaging=1, show=True, v_range=None):
        self.Currents.draw_indep_graphs(rel_time=relative_time, averaging=averaging, show=show, v_range=v_range)

    def get_current(self):
        return self.Currents.get_current()

    def get_irradiation(self):
        return self.Run.get_irradiations()[self.DiamondNumber - 1]

    def get_attenuator(self):
        return False

    def draw_fid_cut(self, scale=1):
        self.Cut.draw_fid_cut(scale)

