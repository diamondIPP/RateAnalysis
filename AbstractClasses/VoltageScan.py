#!/usr/bin/env python
# --------------------------------------------------------
#       Class do analysis of the voltage scans as subclass of AnalysisCollection
# created on August 14th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import gStyle
from InfoLegend import InfoLegend
from Utils import make_ufloat


class VoltageScan(Elementary):

    def __init__(self, ana_collection):
        self.Ana = ana_collection
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.save_dir = self.Ana.save_dir
        self.collection = self.Ana.collection
        self.RunPlan = self.Ana.RunPlan
        self.DiamondName = self.Ana.DiamondName
        self.DiamondNumber = self.Ana.DiamondNumber
        self.InfoLegend = InfoLegend(ana_collection)

    def get_voltages(self):
        return [make_ufloat((ana.Bias, abs(ana.Bias * 1e-4))) for ana in self.collection.itervalues()]

    def draw_all(self, redo=False):
        self.draw_pulse_height(show=False, redo=redo)
        self.draw_pulser_pulse_height(show=False, redo=redo)
        self.draw_pedestals(show=False, redo=redo)
        self.draw_pedestals(show=False, sigma=True, redo=redo)
        self.draw_pedestals(show=False, pulser=True, redo=redo)

    def draw_pedestals(self, cut=None, pulser=False, sigma=False, show=True, redo=False):

        mode = 'Sigma' if sigma else 'Mean'
        pedestals = self.Ana.get_pedestals(cut, pulser, redo)
        y_values = [dic['sigma' if sigma else 'ph'] for dic in pedestals.itervalues()]
        gr = self.make_tgrapherrors('g_vpd', '{p}Pedestal {m} vs. Bias Voltage'.format(p='Pulser' if pulser else '', m=mode), x=self.get_voltages(), y=y_values, color=self.get_color())
        self.format_histo(gr, x_tit='Voltage [V]', y_tit='Pulse Height [au]', y_off=1.4)
        self.save_histo(gr, '{s}Pedestal{m}Voltage'.format(s='Pulser' if pulser else '', m=mode), show, lm=.13, draw_opt='ap')
        self.reset_colors()

    def draw_efficiency(self, show=True):
        x = self.get_voltages()
        efficiencies = [ana.get_hit_efficiency() for ana in self.collection.itervalues()]
        g = self.make_tgrapherrors('gev', 'Efficiency vs. Voltage', x=x, y=efficiencies, color=self.colors[2])
        self.format_histo(g, x_tit='Voltage [V]', y_tit='Hit Efficiency [%]', y_off=1.3, y_range=[0, 108], x_range=[min(x).n - 10, max(x).n + 10])
        self.draw_histo(g, show=show, gridy=True)
        self.draw_preliminary()
        self.save_plots('EfficiencyVoltage')

    def draw_pulser_pulse_height(self, binning=10000, redo=False, show=True):
        return self.draw_pulse_height(binning, pulser=True, redo=redo, show=show)

    def draw_pulse_height(self, binning=None, pulser=False, redo=False, show=True):

        marker_size = 1
        gStyle.SetEndErrorSize(4)
        ph = self.Ana.Pulser.get_pulse_heights(redo=redo) if pulser else self.Ana.get_pulse_heights(binning, redo)
        y_values = [dic['ph'] for dic in ph.itervalues()]
        x = self.get_voltages()
        rel_sys_error = 0  # TODO think of something else here...
        fac = 1 if self.Ana.DUTType == 'pad' else 1e-3
        y_values = [make_ufloat((v.n, v.s + rel_sys_error * v.n)) * fac for v in y_values]
        g = self.make_tgrapherrors('gStatError', 'stat. error', self.colors[2], x=x, y=y_values, marker_size=marker_size)
        ytit = 'Pulse Height [{}]'.format('mV' if self.Ana.DUTType == 'pad' else 'ke')
        self.format_histo(g, x_tit='Voltage [V]', y_tit=ytit, y_off=1.4, draw_first=True, y_range=[0, max(y_values).n * 1.2], x_range=[min(x).n - 10, max(x).n + 10])
        self.draw_histo(g, draw_opt='apl', lm=.14, show=show)
        self.draw_preliminary()
        self.save_plots('{s}VoltageScan'.format(s='Signal' if not pulser else 'Pulser'))
        return g
