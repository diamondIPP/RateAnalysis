#!/usr/bin/env python
# --------------------------------------------------------
#       Pulser sub class for AnalysisCollection
# created by Michael Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TMultiGraph
from Elementary import Elementary
from InfoLegend import InfoLegend
from Utils import *
from functools import partial


class PulserCollection(Elementary):

    def __init__(self, ana_collection):
        Elementary.__init__(self, verbose=ana_collection.verbose)
        self.Analysis = ana_collection

        self.Collection = self.Analysis.collection
        self.DiamondName = self.Analysis.DiamondName
        self.RunPlan = self.Analysis.RunPlan
        self.InfoLegend = InfoLegend(ana_collection)
        self.save_dir = self.Analysis.save_dir

    def get_pulse_height_graph(self, sigma=False, vs_time=False, corr=True, beam_on=True, redo=False):

        self.log_info('Getting pulser pulse heights{}'.format(' vs time' if vs_time else ''))
        marker_size = 2
        par = 2 if sigma else 1
        y_values = []
        self.start_pbar(self.Analysis.NRuns)
        for i, ana in enumerate(self.Collection.itervalues()):
            y_values.append(make_ufloat(ana.Pulser.draw_distribution_fit(corr=corr, beam_on=beam_on, redo=redo, show=False, prnt=False), par))
            # add pedestal sigma as error
            cut = ana.Cut.generate_pulser_cut(beam_on)
            pedestal_fit = ana.Pedestal.draw_disto_fit(cut=cut, save=False)
            y_values[i] += make_ufloat((0, pedestal_fit.ParError(par)))
            self.ProgressBar.update(i + 1)
        self.ProgressBar.finish()
        x_values = [make_ufloat(ana.run.get_time() if vs_time else ana.get_flux()) for ana in self.Collection.itervalues()]
        g = self.make_tgrapherrors('g_pph', 'data', self.get_color(), marker_size=marker_size, x=x_values, y=y_values)
        g_first = self.make_tgrapherrors('g1', 'first run', marker=22, color=2, marker_size=marker_size, x=[x_values[0].n], y=[y_values[0].n])
        g_last = self.make_tgrapherrors('g2', 'last run', marker=23, color=2, marker_size=marker_size, x=[x_values[-1].n], y=[y_values[-1].n])
        graphs = [g, g_first, g_last]
        legend = self.make_legend(.17, .35, nentries=3, x2=.4)
        mg = TMultiGraph('mg_pph', 'Pulser Pulse Height vs {mod} - {dia}'.format(mod='Time' if vs_time else 'Flux', dia=self.DiamondName))
        for gr in graphs:
            legend.AddEntry(gr, gr.GetTitle(), 'l' if gr.GetName() == 'gerr' else 'p')
            mg.Add(gr, 'p')
        mg.GetListOfFunctions().Add(legend)
        self.reset_colors()
        if vs_time:
            g = mg.GetListOfGraphs()[1]
            for i, (ana, x) in enumerate(zip(self.Collection.itervalues(), x_values)):
                y, ey = g.GetY()[i], g.GetErrorY(i)
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(self.draw_tlatex(x.n, y + ey * 1.2, '{:1.0f}'.format(ana.get_flux()[0]), color=1, align=21, size=.02))
        return mg

    def draw_pulse_heights(self, sigma=False, corr=True, beam_on=True, vs_time=False, do_fit=False, save_comb=True, show=True, redo=False):

        mode = 'Time' if vs_time else 'Flux'
        pickle_path = self.make_pickle_path('Pulser', 'PulseHeights', self.RunPlan, self.DiamondName, '{}_{}'.format(mode, sigma))

        f = partial(self.get_pulse_height_graph, sigma=sigma, vs_time=vs_time, corr=corr, beam_on=beam_on, redo=redo)
        mg = do_pickle(pickle_path, f, redo=redo)

        y_values = [mg.GetListOfGraphs()[0].GetY()[i] for i in xrange(mg.GetListOfGraphs()[0].GetN())]
        y_range = increased_range([min(y_values), max(y_values)], .3, .3)
        self.format_histo(mg, x_tit=self.Analysis.make_x_tit(vs_time), y_tit='{} [au]'.format('Sigma' if sigma else 'Pulser Pulse Height'), draw_first=True, y_range=y_range, y_off=1.75, x_off=1.3)
        mg.GetXaxis().SetLimits(1, 40000) if not vs_time else do_nothing()
        self.save_histo(mg, 'Pulser{mean}{a}{b}'.format(mean='Sigma' if sigma else 'Mean', a=corr, b=beam_on), lm=.14, draw_opt='A', logx=False if vs_time else True, show=False)
        mg1 = mg.Clone()
        mg1.GetListOfGraphs()[0].SetLineColor(self.colors[0])
        self.draw_legend()
        if do_fit:
            set_statbox(only_fit=True)
            mg.GetListOfGraphs()[0].Fit('pol0', 'q')
        if save_comb:
            ymin = increased_range([min(y_values), max(y_values)], .5)[0]
            self.save_combined_pulse_heights(mg, mg1, ymin, show, name='CombinedPulserPulseHeights', pulser_leg=self.draw_legend)
            self.ROOTObjects.append(mg1)
        return mg

    def draw_scaled_pulse_heights(self, sigma=False, vs_time=False, show=True, redo=False, y_range=None):

        mode = 'Time' if vs_time else 'Flux'
        pickle_path = self.make_pickle_path('Pulser', 'PulseHeights', self.RunPlan, self.DiamondName, '{}_{}'.format(mode, sigma))
        f = partial(self.get_pulse_height_graph, sigma, vs_time, redo=redo)
        mg = do_pickle(pickle_path, f, redo=redo)
        scale_multigraph(mg)
        xtit = 'Time [hh:mm]' if vs_time else 'Flux [kHz/cm^{2}]'
        y_range = [.95, 1.05] if y_range is None else y_range
        self.format_histo(mg, x_tit=xtit, y_tit='Scaled Pulser Pulse Height', y_off=1.75, x_off=1.3, draw_first=True, t_ax_off=0 if vs_time else None, y_range=y_range, ndivx=503, center_y=1)
        mg.GetXaxis().SetLimits(1, 40000) if not vs_time else do_nothing()
        move_legend(mg.GetListOfFunctions()[0], .16, .20)
        self.save_histo(mg, 'ScaledPulseHeights', show, lm=.14, draw_opt='a', logx=not vs_time, grid=vs_time, gridy=True, bm=.18)
        self.draw_irradiation(make_irr_string(self.Analysis.selection.get_irradiation()))

    def draw_legend(self):
        try:
            typ = self.Analysis.FirstAnalysis.RunInfo['pulser']
            pol = 'positive' if self.Analysis.FirstAnalysis.PulserPolarity > 0 else 'negative'
            sig = 'positive' if self.Analysis.FirstAnalysis.Polarity > 0 else 'negative'
            l1 = self.make_legend(.17, .88, nentries=3, margin=.05, clean=True, x2=.5)
            l1.AddEntry(0, 'Pulser Type:', '')
            l1.AddEntry(0, typ, '').SetTextAlign(12)
            l1.AddEntry(0, 'Pulser Polarity:', '')
            l1.AddEntry(0, pol, '').SetTextAlign(12)
            l1.AddEntry(0, 'Signal Polarity:', '')
            l1.AddEntry(0, sig, '').SetTextAlign(12)
            l1.AddEntry(0, 'Pulser Ped. Substr.:', '')
            l1.AddEntry(0, 'yes', '').SetTextAlign(12)
            l1.SetNColumns(2)
            l1.Draw()
            self.ROOTObjects.append(l1)
        except KeyError:
            pass
