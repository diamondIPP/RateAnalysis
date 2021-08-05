#!/usr/bin/env python
# --------------------------------------------------------
#       Class do overwrite methods in case of voltage scans
# created on August 14th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from helpers.draw import array, ufloat


def make_volage_scan(base_class):

    class VoltageScan(base_class):

        def __init__(self, run_plan, dut_nr, testcampaign=None, load_tree=True, verbose=False):
            super().__init__(run_plan, dut_nr, testcampaign, load_tree, verbose)

        def get_voltages(self, runs=None):
            return array([ufloat(ana.DUT.Bias, abs(ana.DUT.Bias * 1e-3)) for ana in self.get_analyses(runs)])

        def get_x_var(self, vs_time=False, *args, **kwargs):
            return self.get_times() if vs_time else self.get_voltages()

        @staticmethod
        def get_x_tit(vs_time):
            return 'Time [hh:mm]' if vs_time else 'Voltage [V]'

        @staticmethod
        def get_x_draw(vs_time=False):
            return {'grid': vs_time}

        def get_x_args(self, vs_time=False, rel_time=False, x_range=None, draw=False):
            hist_kwargs = {'x_tit': self.get_x_tit(vs_time), 't_ax_off': self.get_tax_off(vs_time, rel_time), 'x_off': None if vs_time else 1.2}
            return {**hist_kwargs, **self.get_x_draw(vs_time)} if draw else hist_kwargs

        def get_repr_error(self, *args, **kwargs):
            return 0

        def draw_legend(self, graphs, **kwargs):
            return super(VoltageScan, self).draw_legend(graphs, x=.75)

    return VoltageScan
