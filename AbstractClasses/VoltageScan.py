#!/usr/bin/env python
# --------------------------------------------------------
#       Class do analysis of the voltage scans as subclass of AnalysisCollection
# created on August 14th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import gStyle, TMultiGraph
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
        g = self.make_tgrapherrors('gev', 'Efficiency vs. Voltage')
        for i, (key, ana) in enumerate(self.collection.iteritems()):
            fit = ana.draw_hit_efficiency(show=False)
            x = ana.Run.RunInfo['dia{nr}hv'.format(nr=self.DiamondNumber)]
            s, e = (fit.Parameter(0), fit.ParError(0))
            g.SetPoint(i, x, s)
            g.SetPointError(i, 0, e)
        self.format_histo(g, x_tit='Voltage [V]', y_tit='Hit Efficiency [%]', y_off=1.3)
        self.save_histo(g, 'EfficiencyVoltage', show=show)

    def draw_pulser_pulse_height(self, binning=10000, redo=False, show=True):
        return self.draw_pulse_height(binning, pulser=True, redo=redo, show=show)

    def draw_pulse_height(self, binning=None, pulser=False, redo=False, first_last=True, show=True):

        marker_size = 1
        gStyle.SetEndErrorSize(4)
        ph = self.Ana.Pulser.get_pulse_heights(redo=redo) if pulser else self.Ana.get_pulse_heights(binning, redo)
        y_values = [dic['ph'] for dic in ph.itervalues()]
        x_values = self.get_voltages()
        g = self.make_tgrapherrors('gStatError', 'stat. error', self.colors[2], x=self.get_voltages(), y=y_values, marker_size=marker_size)
        rel_sys_error = 0  # TODO think of something else here...
        y_values = [make_ufloat((v.n, v.s + rel_sys_error * v.n)) for v in y_values]
        g_errors = self.make_tgrapherrors('gerr', 'full error', marker=0, color=602, marker_size=0, x=x_values, y=y_values)
        g_first = self.make_tgrapherrors('g1', 'first run', marker=22, color=2, marker_size=marker_size, x=[x_values[0].n], y=[y_values[0].n])
        g_last = self.make_tgrapherrors('g2', 'last run', marker=23, color=2, marker_size=marker_size, x=[x_values[-1].n], y=[y_values[-1].n])
        graphs = [g, g_errors]
        graphs += [g_first, g_last] if first_last else []
        l = self.make_legend(.75, .37, nentries=len(graphs))
        mg = TMultiGraph('mg_ph', 'Pulse Height - {dia}'.format(dia=self.DiamondName))
        for gr in graphs:
            l.AddEntry(gr, gr.GetTitle(), 'l' if gr.GetName() == 'gerr' else 'p')
            mg.Add(gr, 'pl')
        mg.GetListOfFunctions().Add(l)
        self.reset_colors()
        max_y = max(graphs[0].GetY()[i] for i in xrange(graphs[0].GetN()))
        self.format_histo(mg, x_tit='Voltage [V]', y_tit='Pulse Height [au]', y_off=1.4, draw_first=True, y_range=[0, max_y * 1.2])
        self.save_histo(mg, '{s}VoltageScan'.format(s='Signal' if not pulser else 'Pulser'), draw_opt='a', lm=.14, show=show)
        return mg
