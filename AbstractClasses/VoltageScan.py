#!/usr/bin/env python
# --------------------------------------------------------
#       Class do analysis of the voltage scans as subclass of AnalysisCollection
# created on August 14th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary


class VoltageScan(Elementary):

    def __init__(self, ana_collection):
        self.Ana = ana_collection
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.save_dir = self.Ana.save_dir
        self.collection = self.Ana.collection
        self.RunPlan = self.Ana.RunPlan
        self.DiamondName = self.Ana.DiamondName
        self.DiamondNumber = self.Ana.DiamondNumber

    def draw_pedestals(self, cut=None, pulser=False, sigma=False, show=True):

        mode = 'Sigma' if sigma else 'Mean'
        gr = self.make_tgrapherrors('g_vpd', '{p}Pedestal {m} vs. Bias Voltage'.format(p='Pulser' if pulser else '', m=mode), self.get_color())

        for i, (key, ana) in enumerate(self.collection.iteritems()):
            fit_par = ana.Pedestal.draw_disto_fit(cut=cut, show=False) if not pulser else ana.Pulser.draw_pedestal(show=False)
            x = ana.run.RunInfo['dia{nr}hv'.format(nr=self.DiamondNumber)]
            par = 2 if sigma else 1
            gr.SetPoint(i, x, fit_par.Parameter(par))
            gr.SetPointError(i, 0, fit_par.ParError(par))

        self.format_histo(gr, x_tit='Voltage [V]', y_tit='Pulse Height [au]', y_off=1.4)
        self.save_histo(gr, '{s}Pedestal{m}Voltage'.format(s='Pulser' if pulser else '', m=mode), show, lm=.13, draw_opt='ap')
