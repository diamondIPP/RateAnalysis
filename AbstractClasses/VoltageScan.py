#!/usr/bin/env python
# --------------------------------------------------------
#       Class do analysis of the voltage scans as subclass of AnalysisCollection
# created on August 14th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import gStyle, TMultiGraph
from InfoLegend import InfoLegend


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

    def draw_all(self):
        self.draw_pulse_height(show=False)
        self.draw_pulser_pulse_height(show=False)
        self.draw_pedestals(show=False)
        self.draw_pedestals(show=False, sigma=True)
        self.draw_pedestals(show=False, pulser=True)

    def draw_pedestals(self, cut=None, pulser=False, sigma=False, show=True):

        mode = 'Sigma' if sigma else 'Mean'
        gr = self.make_tgrapherrors('g_vpd', '{p}Pedestal {m} vs. Bias Voltage'.format(p='Pulser' if pulser else '', m=mode), self.get_color())

        for i, (key, ana) in enumerate(self.collection.iteritems()):
            fit_par = ana.Pedestal.draw_disto_fit(cut=cut, show=False) if not pulser else ana.Pulser.draw_pedestal(show=False)
            x = ana.Run.RunInfo['dia{nr}hv'.format(nr=self.DiamondNumber)]
            par = 2 if sigma else 1
            gr.SetPoint(i, x, fit_par.Parameter(par))
            gr.SetPointError(i, 0, fit_par.ParError(par))

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

    def draw_pulse_height(self, binning=30000, pulser=False, redo=False, show=True, cut=None):
        gr1 = self.make_tgrapherrors('gStatError', 'stat. error', self.get_color())
        gStyle.SetEndErrorSize(4)
        gr_first = self.make_tgrapherrors('gFirst', 'first run', marker=22, color=2, marker_size=2)
        gr_last = self.make_tgrapherrors('gLast', 'last run', marker=23, color=2, marker_size=2)
        gr_errors = self.make_tgrapherrors('gFullError', 'stat. + repr. error', marker=0, color=602, marker_size=0)

        # flux_errors = self.draw_ph_distributions_below_flux(flux=80, show=False, save_plot=False)
        # rel_sys_error = flux_errors[1] / flux_errors[0]
        rel_sys_error = 0
        i, j = 0, 0
        for key, ana in self.collection.iteritems():
            try:
                fit1 = ana.draw_pulse_height(binning=binning, corr=True, save=redo, show=False, redo=redo)[1] if not pulser else ana.Pulser.draw_distribution_fit(show=False, prnt=False)
            except TypeError:
                fit1 = ana.draw_pulse_height(bin_size=binning, cut=cut, show=False)
            x = ana.Run.RunInfo['dia{nr}hv'.format(nr=self.DiamondNumber)]
            s, e = (fit1.Parameter(0), fit1.ParError(0)) if not pulser else (fit1.Parameter(1), fit1.ParError(1))
            gr1.SetPoint(i, x, s)
            self.log_info('{x}\t{s:5.2f} {e:3.2f}'.format(x=x, s=s, e=e))
            gr1.SetPointError(i, 0, e)
            gr_errors.SetPoint(i, x, s)
            gr_errors.SetPointError(i, 0, e + rel_sys_error * s)
            # set special markers for the first and last run
            if i == 0:
                gr_first.SetPoint(0, x, s)
            if j == len(self.collection) - 1:
                gr_last.SetPoint(0, x, s)
            i += 1
            j += 1
        graphs = [gr_errors, gr1]
        gr_line = gr1.Clone()
        self.format_histo(gr_line, name='gLine', color=920)
        graphs += [gr_first, gr_last]
        legend = self.make_legend(.65, .35, nentries=len(graphs))
        # gr1.SetName('data') if len(graphs) < 5 else self.do_nothing()

        mg = TMultiGraph('mg_ph', '' + self.DiamondName)
        # mg.Add(gr_line, 'l')
        for gr in graphs:
            if gr.GetName().startswith('gFull'):
                legend.AddEntry(gr, gr.GetTitle(), 'l')
            else:
                legend.AddEntry(gr, gr.GetTitle(), 'p')
            mg.Add(gr, 'p')
        self.format_histo(mg, x_tit='Voltage [V]', y_tit='Pulse Height [au]', y_off=1.3, draw_first=True)
        self.save_histo(mg, '{s}VoltageScan'.format(s='Signal' if not pulser else 'Pulser'), draw_opt='a', lm=.12, show=show)
        self.reset_colors()
