#!/usr/bin/env python
# --------------------------------------------------------
#       Class to draw the info legend for an analysis class
# created on Jan 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from helpers.utils import make_tc_str, timedelta, flux2str, make_irr_string, mean, plural, uarr2n
from numpy import sort
from plotting.info import Info
from plotting.draw import get_window_ratio


class AnaInfo(Info):
    def __init__(self, draw):
        super().__init__(draw)
        self.Ana = draw.Analysis
        Info.ShowLegend = self.Draw.Config.get_value('SAVE', 'info legend', default=False)

    @property
    def has_dut(self):
        return hasattr(self.Ana, 'DUT')

    @property
    def is_collection(self):
        return hasattr(self.Ana, 'Runs') or hasattr(self.Ana, 'IsCollection') and self.Ana.IsCollection

    def draw_legend(self):
        run_str, info_str = self.get_run_string(), self.get_info_string()
        width = max(len(run_str), len(info_str)) * get_window_ratio() / 80
        show = self.ShowLegend and self.is_active
        return self.Draw.legend([0] * 3, [run_str, self.get_dia_string(), info_str], '', x1=.005, y1=.003, w=width, nentries=3, clean=False, scale=.75, margin=.05, fix=True, show=show)

    def is_active(self):
        if self.is_collection:
            return self.Ana.Ana.LoadTree if hasattr(self.Ana, 'Ana') and hasattr(self.Ana.Ana, 'LoadTree') else self.Ana.LoadTree if hasattr(self.Ana, 'LoadTree') else True
        return hasattr(self.Ana, 'Run') and self.Ana.Run.Number is not None and bool(self.Ana.Tree.Hash() if hasattr(self.Ana, 'Tree') else self.Ana.Ana.Tree.Hash())

    def get_run_string(self):
        dur = sum([ana.Run.Duration for ana in self.Ana.get_analyses()], timedelta()) if self.is_collection else self.Ana.Run.Duration
        dur -= timedelta(microseconds=dur.microseconds)
        nr = f's {"{}-{}".format(*self.Ana.Runs[[0, -1]])}' if self.is_collection else f' {self.Ana.Run.Number}'
        run_str = f'Run{nr}: {self.get_rate_string()}, {dur}'
        return run_str + ('' if self.is_collection else f' ({self.Ana.Run.NEvents} evts)')

    def get_dia_string(self):
        dia_str = self.Ana.DUT.Name if self.has_dut else ', '.join(dut.Name for dut in self.Ana.Run.DUTs)
        return f'{plural("Detector", self.has_dut)}: {dia_str} ({make_tc_str(self.Ana.TCString)})'

    def get_rate_string(self):
        if self.is_collection:
            fluxes = uarr2n(self.Ana.get_fluxes(pbar=False))
            return flux2str(mean(fluxes), prec=0) if 'voltage' in self.Ana.Type else ' - '.join(flux2str(i) for i in sort(fluxes)[[0, -1]])
        return flux2str(self.Ana.Run.Flux.n)

    def get_info_string(self):
        voltage = '{0:+4.0f}V'.format(self.Ana.DUT.Bias) if self.has_dut else '/'.join('{0:+4.0f}V'.format(dut.Bias) for dut in self.Ana.Run.DUTs)
        irradiation = make_irr_string(self.Ana.DUT.Irradiation[self.Ana.TCString]) if self.has_dut else '/'.join(make_irr_string(dut.get_irradiation(self.Ana.TCString)) for dut in self.Ana.Run.DUTs)
        attenuator = 'Att: {}'.format(str(self.Ana.DUT.Attenuator)) if self.has_dut and self.Ana.DUT.Attenuator and 'Pad' in self.Ana.__class__.__name__ else ''
        return 'Info: {v}, {i}{a}'.format(v=voltage, i=irradiation, a=', {}'.format(attenuator) if attenuator else '')
