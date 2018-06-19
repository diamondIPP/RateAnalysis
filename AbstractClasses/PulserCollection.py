#!/usr/bin/env python
# --------------------------------------------------------
#       Pulser sub class for AnalysisCollection
# created by Michael Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TMultiGraph, TH1F
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

    def get_pulse_height_graph(self, sigma=False, vs_time=False, corr=True, beam_on=True, redo=False, legend=True, show_flux=True):

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
        graphs = [g] if vs_time else [g, g_first, g_last]
        l = self.make_legend(.17, .35, nentries=3, x2=.4)
        mg = TMultiGraph('mg_pph', 'Pulser Pulse Height vs {mod} - {dia}'.format(mod='Time' if vs_time else 'Flux', dia=self.DiamondName))
        for gr in graphs:
            l.AddEntry(gr, gr.GetTitle(), 'l' if gr.GetName() == 'gerr' else 'p')
            mg.Add(gr, 'p')
        if legend:
            mg.GetListOfFunctions().Add(l)
        self.reset_colors()
        if vs_time and show_flux:
            g = mg.GetListOfGraphs()[0]
            for i, (ana, x) in enumerate(zip(self.Collection.itervalues(), x_values)):
                y, ey = g.GetY()[i], g.GetErrorY(i)
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(self.draw_tlatex(x.n, y + ey * 1.2, '{:1.0f}'.format(ana.get_flux()[0]), color=1, align=21, size=.02))
        return mg

    def draw_pulse_height(self, evts_per_bin=10000, show=True, rel_t=True):
        """Draw time evolution of the pulse height"""
        histos = [ana.Pulser.draw_pulse_height(evts_per_bin, show=False, fit=False) for ana in self.Collection.itervalues()]
        h1 = TH1F('hpra', 'Pulser Rate for Run Plan {n}'.format(n=self.RunPlan), *self.Analysis.get_binning(evts_per_bin, t_bins=True))
        i_bin = 0
        for h in histos:
            for i in xrange(1, h.GetNbinsX() + 1):
                h1.SetBinContent(i_bin + 1, h.GetBinContent(i))
                h1.SetBinError(i_bin + 1, h.GetBinError(i))
                i_bin += 1
            i_bin += 1
        self.format_histo(h1, x_tit='Time [hh:mm]', y_tit='Pulser Pulse Height [au]', y_off=.8, fill_color=self.FillColor, stats=0, y_range=[0, 105])
        set_time_axis(h1, off=self.Analysis.FirstAnalysis.run.StartTime if rel_t else 0)
        self.save_histo(h1, 'PulserPulseHeight{}'.format(evts_per_bin), show=show, draw_opt='hist', x_fac=1.5, y_fac=.75, lm=.065)

    def draw_pulse_heights(self, sigma=False, corr=True, beam_on=True, vs_time=False, do_fit=False, save_comb=True, show=True, redo=False):

        mode = 'Time' if vs_time else 'Flux'
        pickle_path = self.make_pickle_path('Pulser', 'PulseHeights', self.RunPlan, self.DiamondName, '{}_{}_{}{}'.format(mode, sigma, corr, beam_on))

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

    def draw_scaled_pulse_heights(self, sigma=False, vs_time=False, redo=False, y_range=None, scale=1, show=True):

        mode = 'Time' if vs_time else 'Flux'
        pickle_path = self.make_pickle_path('Pulser', 'PulseHeights', self.RunPlan, self.DiamondName, '{}_{}'.format(mode, sigma))
        f = partial(self.get_pulse_height_graph, sigma, vs_time, redo=redo)
        mg = do_pickle(pickle_path, f, redo=redo)
        scale_multigraph(mg, scale)
        xtit = 'Time [hh:mm]' if vs_time else 'Flux [kHz/cm^{2}]'
        y_range = [.95, 1.05] if y_range is None and scale == 1 else y_range
        self.format_histo(mg, x_tit=xtit, y_tit='Scaled Pulser Pulse Height', y_off=1.75, x_off=1.3, draw_first=True, t_ax_off=0 if vs_time else None, y_range=y_range, ndivx=503, center_y=1)
        mg.GetXaxis().SetLimits(1, 40000) if not vs_time else do_nothing()
        move_legend(mg.GetListOfFunctions()[0], .16, .20)
        self.save_histo(mg, 'ScaledPulseHeights', show, lm=.14, draw_opt='a', logx=not vs_time, grid=vs_time, gridy=True, bm=.18)
        self.draw_irradiation(make_irr_string(self.Analysis.selection.get_irradiation()))
        return mg

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

    def draw_distributions(self, show=True, corr=True):

        set_root_output(show)
        legend = self.make_legend(nentries=self.Analysis.NRuns)
        c = self.make_canvas('c_ph', 'Pulser Histos', 1, 1)
        c.SetRightMargin(.03)
        for i, ana in enumerate(self.Collection.itervalues()):
            h = ana.Pulser.draw_distribution(show=False, corr=corr, prnt=False, binning=5)
            c.cd()
            h.SetTitle('Pulser Distributions {0}Corrected'.format('Pedestal' if corr else 'Un'))
            x_range = increased_range([h.GetBinCenter(h.FindFirstBinAbove(2) * 10 / 10 - 20), h.GetBinCenter(h.FindLastBinAbove(2) * 10 / 10 + 10)], 0, .3)
            self.format_histo(h, fill_color=0, stats=0, x_range=x_range, color=self.get_color(), lw=2)
            h.Scale(1 / h.GetMaximum())
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm^{{2}}'.format(ana.get_flux()[0]), 'l')
            self.ROOTObjects.append(h)
        legend.Draw()
        self.save_plots('AllPulserHistos{0}'.format('Uncorrected' if not corr else ''))
        self.reset_colors()

    def draw_all_pulse_heights(self, sigma=False):
        graphs = [self.get_pulse_height_graph(sigma=sigma, corr=bool(x), beam_on=bool(y)) for x, y in zip([1, 1, 0, 0], [1, 0, 1, 0])]
        y_range = increased_range(find_graph_margins(graphs), .1, .1)
        c = self.make_canvas('c_apph', 'Pulser Info', 1, 1)
        c.Divide(2, 2)
        for i, gr in enumerate(graphs, 1):
            pad = c.cd(i)
            self.format_histo(gr, y_off=1.3, draw_first=True, y_range=y_range)
            gr.GetXaxis().SetLimits(1, 40000)
            pad.SetLogx()
            pad.SetBottomMargin(.15)
            gr.Draw('a')
        self.ROOTObjects.append([graphs, c])
        self.save_plots('AllPulserOverview{0}'.format('Sigma' if sigma else 'Mean'))

    def draw_rate(self, evts_per_bin=1000, cut=None, rel_t=True, show=True):
        histos = [ana.Pulser.draw_rate(evts_per_bin, show=False, cut=cut, vs_time=True) for ana in self.Collection.itervalues()]
        h1 = TH1F('hpra', 'Pulser Rate for Run Plan {n}'.format(n=self.RunPlan), *self.Analysis.get_binning(evts_per_bin, t_bins=True))
        i_bin = 0
        for h in histos:
            for i in xrange(1, h.GetNbinsX() + 1):
                h1.SetBinContent(i_bin + 1, h.GetBinContent(i))
                h1.SetBinError(i_bin + 1, h.GetBinError(i))
                i_bin += 1
            i_bin += 1
        self.format_histo(h1, x_tit='Time [hh:mm]', y_tit='Pulser Rate [%]', y_off=.8, fill_color=self.FillColor, stats=0, y_range=[0, 105])
        set_time_axis(h1, off=self.Analysis.FirstAnalysis.run.StartTime if rel_t else 0)
        self.save_histo(h1, 'AllPulserRate', show=show, draw_opt='hist', x_fac=1.5, y_fac=.75, lm=.065)

    def draw_rates(self, show=True, vs_time=False, real=False):
        mode = self.Analysis.get_mode(vs_time)
        x_values, y_values = [], []
        self.start_pbar(self.Analysis.NRuns)
        for i, ana in enumerate(self.Collection.itervalues(), 1):
            x_values += [make_ufloat(ana.run.get_time() if vs_time else ana.get_flux())]
            y_values += [ana.Pulser.calc_real_fraction() if real else make_ufloat(ana.Pulser.calc_fraction(prnt=False))]
            self.ProgressBar.update(i)
        g = self.make_tgrapherrors('g_pr', 'Pulser Fraction vs {mod} '.format(mod=mode), x=x_values, y=y_values, color=self.colors[0], marker_size=2)
        self.format_histo(g, x_tit=self.Analysis.make_x_tit(vs_time), y_tit='Pulser Fraction [%]')
        self.save_histo(g, 'PulserFraction{0}'.format(mode), show, logx=not vs_time, draw_opt='ap')
        return g

    def calc_error(self, fs11, fsh13):
        runs = self.Analysis.get_runs_by_collimator(fs11=fs11, fsh13=fsh13)
        pulse_heights = [make_ufloat(self.Collection[run].Pulser.draw_distribution_fit(show=False), par=1) for run in runs]
        return mean_sigma(pulse_heights)

    def calc_all_errors(self):
        collimator_settings = set([(ana.run.RunInfo['fs11'], ana.run.RunInfo['fs13']) for key, ana in self.Collection.iteritems()])
        fits = {}
        for s in collimator_settings:
            retval = self.calc_error(fs11=s[0], fsh13=s[1])
            fits[s] = retval
        return fits
