#!/usr/bin/env python
# --------------------------------------------------------
#       Class do overwrite methods in case of voltage scans
# created on August 14th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from plotting.draw import array, ufloat, prep_kw
from pixel.collection import PixCollection
from pad.collection import PadCollection, AnalysisCollection


def v_scan(cls):
    class VScan(cls):
        def __init__(self, run_plan, dut_nr, testcampaign=None, load_tree=True, verbose=False):
            super().__init__(run_plan, dut_nr, testcampaign, load_tree, verbose)

        def get_voltages(self, runs=None):
            return array([ufloat(ana.DUT.Bias, abs(ana.DUT.Bias * 1e-3)) for ana in self.get_analyses(runs)])

        def get_x_var(self, vs_time=False, *args, **kwargs):  # noqa
            return self.get_times() if vs_time else self.get_voltages()

        @staticmethod
        def get_x_tit(vs_time=False, vs_irrad=False, e_field=False):  # noqa
            return 'Time [hh:mm]' if vs_time else 'Electric Field [V/#mum]' if e_field else 'Voltage [V]'

        @staticmethod
        def get_x_draw(vs_time=False):
            return {'grid': vs_time}

        @staticmethod
        def get_x_args(vs_time=False, rel_time=False, vs_irrad=False, draw=False, e_field=False, **kwargs):
            kw = prep_kw(kwargs, x_tit=VScan.get_x_tit(vs_time, 0, e_field), t_ax_off=AnalysisCollection.get_tax_off(vs_time, rel_time), x_off=None if vs_time or vs_irrad else 1.1)
            return {**kw, **VScan.get_x_draw(vs_time)} if draw else kw

        def get_repr_error(self, *args, **kwargs):  # noqa
            return 0

        def draw_legend(self, graphs, **kwargs):  # noqa
            return super(VScan, self).make_legend(graphs, x=.75)

    return VScan


PixVScan = v_scan(PixCollection)
PadVScan = v_scan(PadCollection)
VoltageScan = v_scan(PadCollection)   # statically the same as pixel
