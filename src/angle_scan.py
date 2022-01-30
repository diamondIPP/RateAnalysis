#!/usr/bin/env python
# --------------------------------------------------------
#       Class do overwrite methods in case of angle scans
# created on January 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import insert

from pad.collection import PadCollection, AnalysisCollection
from pixel.collection import PixCollection
from plotting.draw import array, ufloat, prep_kw


def a_scan(cls):
    class AScan(cls):
        def __init__(self, run_plan, dut_nr, testcampaign=None, load_tree=True, verbose=False):
            super().__init__(run_plan, dut_nr, testcampaign, load_tree, verbose)

        @property
        def info_header(self):
            return insert(super().info_header, 3, 'Angle [deg]').tolist()

        @property
        def info_rows(self):
            return insert(super().info_rows, 3, [f'{a.n:>10}' for a in self.get_angles()], axis=1)

        def get_angles(self):
            return array([ufloat(run.Info['angle'], .5) for run in self.Ensemble.Runs])

        def get_x_var(self, vs_time=False, *args, **kwargs):  # noqa
            return self.get_times() if vs_time else self.get_angles()

        @staticmethod
        def get_x_tit(vs_time=False, vs_irrad=False):  # noqa
            return 'Time [hh:mm]' if vs_time else 'Rotation Angle [deg]'

        @staticmethod
        def get_x_draw(vs_time=False):
            return {'grid': vs_time}

        @staticmethod
        def get_x_args(vs_time=False, rel_time=False, vs_irrad=False, draw=False, **kwargs):
            kw = prep_kw(kwargs, x_tit=AScan.get_x_tit(vs_time), t_ax_off=AnalysisCollection.get_tax_off(vs_time, rel_time), x_off=None if vs_time or vs_irrad else 1.1)
            return {**kw, **AScan.get_x_draw(vs_time)} if draw else kw

        def get_repr_error(self, *args, **kwargs):  # noqa
            return 0

        def draw_legend(self, graphs, **kwargs):  # noqa
            return super(AScan, self).draw_legend(graphs, x=.75)

    return AScan


PixAScan = a_scan(PixCollection)
PadAScan = a_scan(PadCollection)
AngleScan = a_scan(PadCollection)   # statically the same as pixel
