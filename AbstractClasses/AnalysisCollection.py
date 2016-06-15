# ==============================================
# IMPORTS
# ==============================================
import json
import os
from ConfigParser import ConfigParser
from argparse import ArgumentParser
from collections import OrderedDict
from numpy import sqrt, log, array, zeros
from time import time

from ROOT import gROOT, TCanvas, TLegend, TExec, gStyle, TMultiGraph, THStack

from CurrentInfo import Currents
from Elementary import Elementary
from Extrema import Extrema2D
from PadAnalysis import SignalAnalysis
from RunSelection import RunSelection
from TelescopeAnalysis import Analysis
from Utils import *


# ==============================================
# MAIN CLASS
# ==============================================
class AnalysisCollection(Elementary):
    """
    An object of this class contains several analysis of runs.
    It gives the ability to compare the data from different runs.
    """
    current_run_number = -1

    def __init__(self, list_of_runs, diamonds=None, verbose=False):
        Elementary.__init__(self, verbose=verbose)

        # dict where all analysis objects are saved
        self.collection = OrderedDict()
        self.selection = list_of_runs if isinstance(list_of_runs, RunSelection) else self.make_runselection(list_of_runs)

        self.runs = self.load_runs(list_of_runs)
        self.diamonds = self.load_diamonds(diamonds, list_of_runs)
        self.min_max_rate_runs = self.get_high_low_rate_runs()

        self.generate_slope_pickle()
        self.generate_threshold_pickle()

        self.add_analyses()

        self.signalValues = None

        # important information
        self.diamond_name = self.get_first_analysis().diamond_name
        self.bias = self.get_first_analysis().bias

        # root stuff
        self.run_plan = list_of_runs.selected_runplan if isinstance(list_of_runs, RunSelection) else '-'
        self.save_dir = '{dia}/runplan{plan}'.format(tc=self.TESTCAMPAIGN[2:], plan=self.run_plan, dia=self.diamond_name)
        self.RootObjects = []
        # important plots
        self.FWHM = None
        self.PulseHeight = None
        self.Pedestal = None
        self.PeakDistribution = None

        # current information
        self.StartTime = float(self.get_first_analysis().run.log_start.strftime('%s'))
        self.channel = self.get_first_analysis().channel
        self.Currents = Currents(self)

    def __del__(self):
        print 'deleting AnalysisCollection...'
        for nr, ana in self.collection.iteritems():
            ana.__del__()
            print '  deleted Analysis of Run {nr}'.format(nr=nr)
        print 'AnalyisCollection deleted'

    def close_files(self):
        for ana in self.collection.itervalues():
            ana.run.RootFile.Close()

    # ============================================
    # region INIT
    def load_run_config(self):
        return self.load_run_configs(0)

    def add_analyses(self):
        """
        Creates and adds Analysis objects with run numbers in runs.
        """
        for run, dia in sorted(zip(self.runs, self.diamonds)):
            ch = 0 if dia == 1 or dia == 3 else 3
            analysis = SignalAnalysis(run, ch, self.min_max_rate_runs)
            self.collection[analysis.run.run_number] = analysis
            self.current_run_number = analysis.run.run_number

    @staticmethod
    def load_runs(run_list):
        if type(run_list) is list:
            return run_list
        elif isinstance(run_list, RunSelection):
            return run_list.get_selected_runs()
        else:
            raise ValueError('listOfRuns has to be of type list or instance of RunSelection')

    def load_diamonds(self, diamonds, run_list):
        dias = diamonds
        assert type(dias) is list or dias in [1, 2, 3], '"diamonds" has to be 1, 2, 3, or None (0x1: diamond1, 0x2: diamond2)'
        if dias is not None:
            if type(dias) is not list:
                dias = [dias] * len(self.runs)
        else:
            dias = [3] * len(run_list) if type(run_list) is list else run_list.get_selected_diamonds()
        return dias

    def get_high_low_rate_runs(self):
        keydict = ConfigParser()
        keydict.read('Configuration/KeyDict_{tc}.cfg'.format(tc=self.TESTCAMPAIGN))
        path = self.run_config_parser.get('BASIC', 'runinfofile')
        flux_name = keydict.get('KEYNAMES', 'measured flux')
        f = open(path, 'r')
        run_log = json.load(f)
        fluxes = {}
        print 'RUN   FLUX'
        for run in self.runs:
            flux = run_log[str(run)][flux_name]
            print '{run:3d} {flux:6.1f}'.format(run=run, flux=flux)
            fluxes[flux] = run
        return {'min': fluxes[min(fluxes)], 'max': fluxes[max(fluxes)]}

    def generate_slope_pickle(self):
        picklepath = 'Configuration/Individual_Configs/Slope/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.min_max_rate_runs['min'])
        if os.path.exists(picklepath):
            return
        Analysis(self.min_max_rate_runs['min'])

    def generate_threshold_pickle(self):
        picklepath = 'Configuration/Individual_Configs/Cuts/SignalThreshold_{tc}_{run}_{ch}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.min_max_rate_runs['max'], ch=0)
        if os.path.exists(picklepath):
            return
        Analysis(self.min_max_rate_runs['max'])

    @staticmethod
    def make_runselection(run_list):
        assert type(run_list) is list, 'run list argument has to be a list!'
        selection = RunSelection()
        selection.select_runs(run_list, do_assert=True)
        return selection
    # endregion

    def create_all_single_run_plots(self):
        for key, ana in self.collection.iteritems():
            print 'Create Plots for Run ', key
            ana.compare_consecutive_cuts(scale=False)
            ana.compare_consecutive_cuts(scale=True)
            ana.show_cut_contributions(show=False)
            ana.draw_bucket_pedestal(show=False)

    # ============================================
    # region SIGNAL/PEDESTAL
    def draw_ph_with_currents(self, show=True):
        ph = self.draw_pulse_heights(show=False, vs_time=True, fl=False)
        self.Currents.set_graphs()
        cur = self.Currents.CurrentGraph.Clone()
        cur.SetLineColor(899)
        scale = ph.GetListOfGraphs()[0].GetY()[0]
        pul = self.draw_pulser_info(show=False, do_fit=False, vs_time=True, scale=scale)
        pul.SetLineColor(859)
        pul.SetMarkerColor(859)

        entries = [1, 3, 1]
        legends = [self.make_legend(.8, .7, nentries=entries[i], scale=1.7, w=.2) for i in xrange(3)]
        for gr in ph.GetListOfGraphs():
            if gr.GetName() != 'flux':
                legends[1].AddEntry(gr, gr.GetName(), 'p')
        dummy_gr = self.make_tgrapherrors('g', 'g', width=2)
        legends[1].AddEntry(0, 'flux in kHz/cm^{2}', 'p')
        legends[1].AddEntry(dummy_gr, 'duration', 'l')
        legends[0].AddEntry(pul, 'pulser', 'p')
        legends[2].AddEntry(cur, 'current', 'l')

        gROOT.SetBatch(1) if not show else do_nothing()
        c = TCanvas('c', 'c', 1500, 1000)
        margins = [[.075, .05, 0, .1], [.075, .05, 0, 0], [.075, .05, .15, 0]]
        pads = [self.Currents.make_tpad('p{0}'.format(i + 1), 'p{0}'.format(i + 1), pos=[0, i / 3., 1, (i + 1) / 3.], gridx=True, margins=margins[2 - i]) for i in xrange(3)]
        for pad in pads:
            pad.Draw()
        # c.Divide(1, 3)
        draw_opts = ['pl', '', 'l']
        # y_tits = ['Pulse Height [au] ', 'Pulse Height [au] ', 'Current [nA] ']
        y_tits = ['pulse height [au] ', 'pulse height [au] ', 'current [nA] ']

        pad = None
        for i, gr in enumerate([pul, ph, cur], 1):
            pad = pads[3 - i]
            pad.cd()
            pad.SetMargin(*margins[i - 1])
            ymin, ymax = gr.GetYaxis().GetXmin(), gr.GetYaxis().GetXmax()
            if i == 1:
                ymin, ymax = scale_margins(pul, ph)
            if i == 2:
                ymin, ymax = scale_margins(ph, pul)
            self.Currents.draw_frame(pad, ymin, ymax, y_tits[i - 1])
            if i == 2:
                self.Currents.draw_time_axis(ymin, opt='t')
            if i == 3:
                self.Currents.draw_time_axis(ymax, opt='t')
            gr.Draw(draw_opts[i - 1])
            # pad.SetLogx()
            legends[i - 1].Draw()

        run_info = self.collection.values()[0].run.get_runinfo(self.channel, pad=pad)
        run_info[0].SetX2NDC(.25)
        run_info[0].SetY2NDC(.145)
        run_info[0].Draw()
        run_info[1].Draw()
        c.Update()
        gROOT.SetBatch(0)
        self.save_canvas(c, self.save_dir, 'PhPulserCurrent')
        self.RootObjects.append([ph, cur, pul, c, legends, pads])

    def draw_pulse_heights(self, binning=20000, flux=True, raw=False, all_corr=False, show=True, vs_time=False, fl=True):
        flux = False if vs_time else flux
        mode = 'Run'
        if flux:
            mode = 'Flux'
        if vs_time:
            mode = 'Time'
        prefix = 'Pulse Height vs {mod} - '.format(mod=mode)
        gr1 = self.make_tgrapherrors('eventwise', prefix + 'eventwise correction', self.get_color())
        gr2 = self.make_tgrapherrors('binwise', prefix + 'binwise correction', self.get_color())
        gr3 = self.make_tgrapherrors('mean ped', prefix + 'mean correction', self.get_color())
        gr4 = self.make_tgrapherrors('raw', prefix + 'raw', self.get_color())
        gr5 = self.make_tgrapherrors('flux', 'bla', 1, width=1, marker_size=0)
        gStyle.SetEndErrorSize(4)
        gr_first = self.make_tgrapherrors('first run', prefix + 'first', marker=22, color=2)
        gr_first.SetMarkerSize(2)
        gr_last = self.make_tgrapherrors('last run', prefix + 'last', marker=23, color=2)
        gr_last.SetMarkerSize(2)
        not_found_for = False in [coll.run.FoundForRate for coll in self.collection.itervalues()]

        i, j = 0, 0
        for key, ana in self.collection.iteritems():
            fit1 = ana.draw_pulse_height(binning, evnt_corr=True, show=False)
            fit2 = ana.draw_pulse_height(binning, bin_corr=True, show=False)
            fit3 = ana.draw_pulse_height(binning, off_corr=True, show=False, evnt_corr=False)
            fit4 = ana.draw_pulse_height(binning, evnt_corr=False, show=False)
            # print '\033[1A', '{0:5.1f}'.format(fit1.Parameter(0))
            x = key
            if flux:
                x = ana.run.RunInfo['measured flux'] if not_found_for else ana.run.flux
            if vs_time:
                x_err = ana.run.duration.seconds / 2.
                x = int(ana.run.log_start.strftime('%s')) + x_err - self.StartTime
                gr5.SetPoint(i, x, fit1.Parameter(0))
                gr5.SetPointError(i, x_err, 0)
                l1 = self.make_tlatex(gr5.GetX()[i] - x_err, gr5.GetY()[i] + .03, '{0:5.0f}'.format(ana.run.flux), color=1, align=10, size=.04)
                gr1.GetListOfFunctions().Add(l1)
            if fit1.Parameter(0) > 10:
                gr1.SetPoint(i, x, fit1.Parameter(0))
                gr2.SetPoint(i, x, fit2.Parameter(0))
                gr3.SetPoint(i, x, fit3.Parameter(0))
                gr4.SetPoint(i, x, fit4.Parameter(0))
                gr1.SetPointError(i, 0, fit1.ParError(0))
                gr2.SetPointError(i, 0, fit2.ParError(0))
                gr3.SetPointError(i, 0, fit3.ParError(0))
                gr4.SetPointError(i, 0, fit4.ParError(0))
                # set special markers for the first and last run
                if i == 0:
                    gr_first.SetPoint(0, x, fit1.Parameter(0))
                if j == len(self.collection) - 1:
                    gr_last.SetPoint(0, x, fit1.Parameter(0))
                i += 1
            j += 1
        gROOT.SetBatch(1)
        graphs = [gr1]
        if fl:
            graphs += [gr_first, gr_last]
        if all_corr:
            graphs += [gr2, gr3]
        if raw:
            graphs.append(gr4)
        legend = self.make_legend(.65, .35, nentries=len(graphs))
        gr1.SetName('data') if len(graphs) < 4 else self.do_nothing()

        mg = TMultiGraph('mg_ph', prefix + self.diamond_name,)
        for gr in graphs:
            legend.AddEntry(gr, gr.GetName(), 'p')
            mg.Add(gr, 'lp')
        mg.Draw('a')

        gROOT.SetBatch(0)
        ymin = mg.GetYaxis().GetXmin()
        ymax = mg.GetYaxis().GetXmax()

        # small range
        self.format_histo(mg, color=None, x_tit=mode + ' [kHz/cm^{2}]' if flux else '', y_tit='Pulse Height [au]', y_off=1.75, x_off=1.3)
        mg.GetYaxis().SetRangeUser(ymin - (ymax - ymin) * .3, ymax)
        if vs_time:
            mg.Add(gr5, '[]')
            mg.Add(gr5, 'p')
            mg.GetXaxis().SetTimeDisplay(1)
            mg.GetXaxis().SetTimeFormat('%H:%M%F2000-02-28 23:00:00')
            mg.GetXaxis().SetLabelSize(.03)
        x_vals = sorted([gr1.GetX()[i] for i in xrange(gr1.GetN())])
        mg.GetXaxis().SetLimits(x_vals[0] * 0.8, x_vals[-1] * 1.2) if flux else self.do_nothing()
        self.RootObjects.append(self.save_histo(mg, 'PulseHeight{mod}'.format(mod=mode.title()), False, self.save_dir, lm=.14, draw_opt='A', l=legend, logx=True if flux else 0,
                                                grid=1 if vs_time else 0))

        # no zero suppression
        mg1 = mg.Clone()
        mg1.SetName('mg1_ph')
        mg1.GetYaxis().SetRangeUser(0, ymax * 1.1)
        self.RootObjects.append(self.save_histo(mg1, 'PulseHeightZero{mod}'.format(mod=mode.title()), False, self.save_dir, lm=.14, draw_opt='A', l=legend, logx=True if flux else 0))

        gROOT.SetBatch(1) if not show else self.do_nothing()
        c = TCanvas('c_phall', 'Rate Scan', 2000, 1000)
        c.Divide(2)
        for i, gr in enumerate([mg, mg1], 1):
            pad = c.cd(i)
            pad.SetLogx() if flux else self.do_nothing()
            pad.SetMargin(.13, .1, .15, .1)
            gr.Draw('a')
            legend.Draw()
        self.RootObjects.append(c)
        self.save_plots('PHOverview{mod}'.format(mod=mode), self.save_dir)
        gROOT.SetBatch(0)
        self.reset_colors()

        self.PulseHeight = gr1
        self.save_combined_pulse_heights(mg, mg1, legend)
        return mg, mg1

    def save_combined_pulse_heights(self, mg, mg1, legend):
        gROOT.SetBatch(0)
        c = TCanvas('c', 'c', 2000, 2000)
        margins = [.13, .13, .15, .1]
        p0 = self.Currents.make_tpad('p0', 'p0', margins=margins, logx=True)
        p1 = self.Currents.make_tpad('p1', 'p1', margins=margins, logx=True, transparent=True)
        for pad in [p0, p1]:
            pad.Draw()

        x = [mg.GetXaxis().GetXmin(), mg.GetXaxis().GetXmax()]
        y = [mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()]
        # first graph
        y0 = [y[0] - (y[1] - y[0]) * .3, y[1]]
        draw_frame(p0, x, y0, base=True, x_tit='flux [kHz/cm^{2}]', y_tit='pulse height [au]', y_off=1.7, x_off=1.2)
        mg.Draw()
        legend.Draw()
        # second graph
        y1 = [0, y[1] * 1.1]
        draw_frame(p1, x, y1)
        self.draw_x_axis(x[1], y1[0], y1[1], 'pulse height [au]', off=1.5)
        mg1.Draw()
        self.save_plots('CombinedPulseHeights', self.save_dir)
        self.RootObjects.append([c])

    def draw_pedestals(self, region='ab', peak_int='2', flux=True, all_regions=False, sigma=False, show=True, cut=None, beam_on=True):
        legend = TLegend(0.7, 0.3, 0.98, .7)
        legend.SetName('l1')
        mode = 'Flux' if flux else 'Run'
        y_val = 'Sigma' if sigma else 'Mean'
        gr1 = self.make_tgrapherrors('pedestal', 'Pedestal {y} in {reg}'.format(y=y_val, reg=region + peak_int))
        regions = self.get_first_analysis().run.pedestal_regions
        graphs = [self.make_tgrapherrors('pedestal', 'Pedestal {y} in {reg}'.format(y=y_val, reg=reg + peak_int), color=self.get_color()) for reg in regions]
        gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        i = 0
        par = 2 if sigma else 1
        cut_string = None
        for key, ana in self.collection.iteritems():
            print 'getting pedestal for run {n}...'.format(n=key)
            cut_string = ana.Cut.generate_pulser_cut(beam_on=beam_on) if cut == 'pulser' else cut
            fit_par = ana.show_pedestal_histo(region, peak_int, cut=cut_string, draw=False)
            x = ana.run.flux if flux else key
            gr1.SetPoint(i, x, fit_par.Parameter(par))
            gr1.SetPointError(i, 0, fit_par.ParError(par))
            if all_regions:
                for reg, gr in zip(regions, graphs):
                    fit_par = ana.show_pedestal_histo(reg, peak_int, draw=False)
                    gr.SetPoint(i, x, fit_par.Parameter(par))
                    gr.SetPointError(i, 0, fit_par.ParError(par))
            i += 1
        if show:
            gROOT.SetBatch(0)
            gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        c = TCanvas('c', 'Pedestal vs Run', 1000, 1000)
        if flux:
            c.SetLogx()
        self.format_histo(gr1, color=None, x_tit=self.make_x_tit(mode, flux), y_tit='Mean Pedestal [au]')
        gr1.Draw('alp')
        if all_regions:
            for i, gr in enumerate(graphs):
                legend.AddEntry(gr, str(regions.values()[i]), 'p')
                gr.Draw('alp') if not i else gr.Draw('lp')
            legend.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.Pedestal = gr1
        self.RootObjects.append([graphs, legend, c])
        save_name = 'Pedestal_{mod}{cut}'.format(mod=mode, cut='' if cut is None else cut_string.GetName())
        self.save_plots(save_name, sub_dir=self.save_dir)
        self.reset_colors()
        return gr1

    def draw_signal_distributions(self, show=True, off=3):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        self.reset_colors()
        stack = THStack('hsd', 'Pulse Height Distributions')
        legend = self.make_legend(nentries=self.get_number_of_analyses())
        histos = [ana.show_signal_histo(show=False) for ana in self.collection.itervalues()]
        for i, h in enumerate(histos):
            self.format_histo(h, lw=2, color=self.get_color())
            h.Scale(1 / h.GetMaximum())
            stack.Add(h)
            legend.AddEntry(h, '{0:06.1f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux()), 'l')
        self.format_histo(stack, y_off=1.55, draw_first=True, x_tit='Pulse Height [au]', y_tit='Number of Entries')
        self.RootObjects.append(self.save_histo(stack, 'SignalDistributions', False, self.save_dir, lm=.13, draw_opt='nostack', l=legend))
        log_stack = stack.Clone()
        log_stack.SetMaximum(off)
        log_stack.SetNameTitle('hsdl', 'Signal Distribution LogY')
        self.RootObjects.append(self.save_histo(log_stack, 'SignalDistributionsLogY', False, self.save_dir, lm=.13, draw_opt='nostack', logy=True, l=legend))
        gROOT.SetBatch(1) if not show else self.do_nothing()
        c = TCanvas('c_sd1', 'Signal Distributions', 1500, 750)
        c.Divide(2)
        legends = [legend, legend.Clone()]
        for i, s in enumerate([stack, log_stack], 1):
            pad = c.cd(i)
            pad.SetLogy() if i == 2 else self.do_nothing()
            s.Draw('nostack')
            legends[i - 1].Draw()
        self.RootObjects.append([c, legends])
        gROOT.SetBatch(0)

    def draw_snrs(self, flux=True, draw=True):
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        gr = self.make_tgrapherrors('gr', 'SNR vs {mode}'.format(mode=mode))
        i = 0
        for key, ana in self.collection.iteritems():
            snr = ana.calc_snr()
            x = ana.run.flux if flux else key
            gr.SetPoint(i, x, snr[0])
            gr.SetPointError(i, 0, snr[1])
            i += 1
        if draw:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        if flux:
            c.SetLogx()
        gr.Draw('ap')
        self.save_plots('AllSNRs', canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        gROOT.SetBatch(0)

    # endregion

    # ============================================
    # region PULSER
    def draw_pulser_info(self, flux=True, show=True, mean=True, corr=True, beam_on=True, vs_time=False, do_fit=True, scale=1):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        flux = False if vs_time else flux
        mode = 'Flux' if flux else 'Run'
        mode = 'Time' if vs_time else mode
        title = '{mean} of Pulser vs {mod} ({ped}, {beam})'.format(mean='Mean' if mean else 'Sigma', mod=mode,
                                                                   ped='pedcorrected' if corr else 'uncorrected', beam='BeamOff' if not beam_on else 'BeamOn')
        gr = self.make_tgrapherrors('gr', title)
        i = 0
        y0 = None
        for key, ana in self.collection.iteritems():
            x = ana.run.flux if flux else key
            fit = ana.calc_pulser_fit(show=False, corr=corr, beam_on=beam_on)
            par = 1 if mean else 2
            cut = ana.Cut.generate_pulser_cut(beam_on)
            ped_fit = ana.show_pedestal_histo(cut=cut, draw=False)
            ped_err = ped_fit.ParError(par)
            if vs_time:
                xerr = ana.run.duration.seconds / 2.
                x = int(ana.run.log_start.strftime('%s')) + xerr - self.StartTime
            y = fit.Parameter(par)
            y0 = y if y0 is None else y0
            yerr = sqrt(pow(fit.ParError(par), 2) + pow(ped_err, 2))
            if scale != 1:
                y *= scale / y0
                yerr *= scale / y0
            if ana.IsAligned:
                gr.SetPoint(i, x, y)
                gr.SetPointError(i, 0, yerr)
                i += 1
        if vs_time:
            gr.GetXaxis().SetTimeDisplay(1)
            gr.GetXaxis().SetTimeFormat('%H:%M%F2000-02-28 23:00:00')
            gr.GetXaxis().SetLabelSize(.03)
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Overview', 1000, 1000)
        if do_fit:
            gStyle.SetOptFit(1)
            gr.Fit('pol0', 'q')
        c.SetLeftMargin(.125)
        c.SetLogx() if flux else self.do_nothing()
        self.format_histo(gr, x_tit=self.make_x_tit(mode, flux), y_tit='{mean} [au]'.format(mean='Mean' if mean else 'Sigma'), y_off=1.8)
        gr.Draw('alp')
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('Pulser{mean}{a}{b}'.format(mean='Mean' if mean else 'Sigma', a=corr, b=beam_on), sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        return gr

    def draw_pulser_histos(self, show=True, corr=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = {i: ana.show_pulser_histo(show=False, corr=corr) for i, ana in enumerate(self.collection.itervalues()) if ana.IsAligned}
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Histos', 1000, 1000)
        legend = TLegend(.13, .88 - self.get_number_of_analyses() * 0.03, .33, .88)
        histos[0].SetTitle('Pulser Distributions {0}Corrected'.format('Pedestal' if corr else 'Un'))
        for i, h in histos.iteritems():
            h.SetStats(0)
            h.GetXaxis().SetRangeUser(h.GetBinCenter(h.FindFirstBinAbove(2) * 10 / 10 - 20), h.GetBinCenter(h.FindLastBinAbove(2) * 10 / 10 + 10))
            h.Scale(1 / h.GetMaximum())
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux()) + '^{2}', 'l')
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.save_plots('AllPulserHistos{0}'.format('Uncorrected' if not corr else ''), sub_dir=self.save_dir)
        self.RootObjects.append([c, legend] + histos.values())
        z.reset_colors()

    def draw_all_pulser_info(self, mean=True):
        graphs = [self.draw_pulser_info(show=False, mean=mean, corr=x, beam_on=y) for x, y in zip([1, 1, 0, 0], [1, 0, 1, 0])]
        margins = self.find_graph_margins(graphs)
        c = TCanvas('c', 'Pulser Info', 1500, 1500)
        c.Divide(2, 2)
        for i, gr in enumerate(graphs, 1):
            gr.GetYaxis().SetRangeUser(*margins)
            self.format_histo(gr, y_off=1.3)
            pad = c.cd(i)
            pad.SetLogx()
            pad.SetBottomMargin(.15)
            gr.Draw('alp')
        gROOT.SetBatch(0)
        self.RootObjects.append([graphs, c])
        self.save_plots('AllPulserOverview{0}'.format('Mean' if mean else 'Sigma'), sub_dir=self.save_dir)

    def compare_pedestals(self, show=True):
        gr1 = self.draw_pedestals(show=False)
        gr2 = self.draw_pedestals(cut='pulser', show=False)
        gr3 = self.draw_pedestals(cut='pulser', show=False, beam_on=False)
        graphs = [gr1, gr2, gr3]
        margins = self.find_graph_margins(graphs)
        gROOT.SetBatch(0) if show else gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Pedestal Comparison', 1000, 1000)
        c.SetLogx()
        legend = TLegend(.7, .78, .88, .88)
        names = ['Signal', 'Pulser', 'BeamOff']
        for i, gr in enumerate(graphs):
            gr.GetYaxis().SetRangeUser(*margins)
            self.format_histo(gr, color=self.get_color())
            gr.Draw('lp') if i else gr.Draw('alp')
            legend.AddEntry(gr, names[i], 'pl')
        legend.Draw()
        gROOT.SetBatch(0)
        self.save_plots('PulserPedestalComparison', sub_dir=self.save_dir)
        self.RootObjects.append([c, graphs, legend])

    def draw_pulser_rates(self, show=True, flux=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        title = 'Pulser Rate vs {mod} '.format(mod=mode)
        gr = self.make_tgrapherrors('gr', title)
        for i, (key, ana) in enumerate(self.collection.iteritems()):
            x = ana.run.flux if flux else key
            fit = ana.fit_pulser_rate(show=False)
            gr.SetPoint(i, x, fit.Parameter(0))
            gr.SetPointError(i, 0, fit.ParError(0))
        gROOT.SetBatch(0) if show else self.do_nothing()
        c = TCanvas('c', 'Pulser Rates', 1000, 1000)
        c.SetLogx() if flux else self.do_nothing()
        self.format_histo(gr, x_tit=self.make_x_tit(mode, flux), y_tit='Pulser Rate [Hz]')
        gr.Draw('alp')
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('PulserRate{0}'.format(mode), sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        return gr

    # endregion

    # ============================================
    # region CUTS

    def draw_bucket_info(self, flux=True, show=True, mean=True):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        gr1 = self.make_tgrapherrors('gr1', '', color=self.get_color())
        prefix = 'Number of Bucket Cut Events' if not mean else 'Mean Pulse Height with Different Bucket Cuts'
        gr2 = self.make_tgrapherrors('gr2', '{pref} vs {mod}'.format(pref=prefix, mod=mode), color=self.get_color())
        gr3 = self.make_tgrapherrors('gr3', '', color=self.get_color())
        i = 0
        for key, ana in self.collection.iteritems():
            x = ana.run.flux if flux else key
            if not mean:
                n = ana.show_bucket_numbers(show=False)
                gr1.SetPoint(i, x, n['new'] / n['all'] * 100)
                gr2.SetPoint(i, x, n['old'] / n['all'] * 100)
            else:
                info = ana.show_bucket_means(show=False, plot_histos=False)
                gr1.SetPoint(i, x, info['new'][0])
                gr2.SetPoint(i, x, info['old'][0])
                gr3.SetPoint(i, x, info['no'][0])
                gr1.SetPointError(i, 0, info['new'][1])
                gr2.SetPointError(i, 0, info['old'][1])
                gr3.SetPointError(i, 0, info['no'][1])
            i += 1
        c = TCanvas('c', 'Bucket Numbers', 1000, 1000)
        c.SetLeftMargin(.13)
        if flux:
            c.SetLogx()
        self.format_histo(gr2, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='Events [%]' if not mean else 'Mean [au]', y_off=1.7, color=None)
        gr2.Draw('apl')
        gr1.Draw('pl')
        if mean:
            gr3.Draw('pl')
        leg = TLegend(.2, .8, .35, .9)
        leg.AddEntry(gr2, 'old cut', 'pl')
        leg.AddEntry(gr1, 'new cut', 'pl')
        if mean:
            leg.AddEntry(gr3, 'no bucket', 'pl')
        leg.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('{mode}_' + mode, sub_dir=self.save_dir)
        self.RootObjects.append([c, gr1, gr2, gr3, leg])

    # endregion

    # ============================================
    # region PEAK VALUES
    def draw_signal_peaks(self, flux=True, draw=True, pulser=False):
        """
        Shows the means of the signal peak distribution.
        :param flux:
        :param draw:
        :param pulser:
        """
        mode = 'Flux' if flux else 'Run'
        signal = 'Pulser' if pulser else 'Signal'
        prefix = 'Mean of {pul} Peaks: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.bias, pul='Pulser' if pulser else 'Signal')
        gr = self.make_tgrapherrors('gr', prefix)
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.fit_peak_values(draw=False, pulser=pulser)
            x = ana.run.flux if flux else key
            gr.SetPoint(i, x, fit.Parameter(1))
            gr.SetPointError(i, 0, fit.ParError(1))
            i += 1
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''))
        c = TCanvas('c', 'Mean of {0} Peaks'.format(signal), 1000, 1000)
        c.SetLogx()
        gr.Draw('alp')
        if not draw:
            c.Close()
        self.RootObjects.append([gr, c])

    def draw_signal_fwhm(self, flux=True, draw=True):
        """
        Shows the FWHM of the signal peak distribution.
        :param flux:
        :param draw:
        """
        mode = 'Flux' if flux else 'Run'
        prefix = 'FWHM of Signal Peaks: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.bias)
        gr = self.make_tgrapherrors('gr1', prefix)
        i = 0
        for key, ana in self.collection.iteritems():
            fwhm = ana.calc_peak_value_fwhm()
            x = ana.run.flux if flux else key
            gr.SetPoint(i, x, fwhm)
            i += 1
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''))
        c = TCanvas('c', 'FWHM of Signal Peaks', 1000, 1000)
        gr.Draw('alp')
        if not draw:
            c.Close()
        self.RootObjects.append([gr, c])

    # endregion

    # ============================================
    # region 2D SIGNAL MAP
    def draw_mean_fwhm(self, saveplots=True, flux=True, draw=True):
        """
        Creates the FWHM Distribution of all selected MeanSignalHistograms
        :param saveplots: if True saves the plot
        :param flux: draw vs flux if True else vs run
        :param draw:
        """
        if not draw:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        prefix = 'FWHM of Mean Signal Histogram: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.bias)
        gr = self.make_tgrapherrors('pedestal', prefix)
        conversion_factor = 2 * sqrt(2 * log(2))  # sigma to FWHM
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.fit_mean_signal_distribution()
            x = ana.run.flux if flux else key
            gr.SetPoint(i, x, fit.Parameter(2) * conversion_factor)
            gr.SetPointError(i, 0, fit.ParError(2) * conversion_factor)
            i += 1
        c = TCanvas('c', 'FWHM', 1000, 1000)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='Flux [kHz/cm2]', y_tit='FWHM [au]', y_off=1.1)
        gr.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        if saveplots:
            self.save_plots('Mean_FWHM_' + mode, canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        self.FWHM = gr

    def save_signal_maps(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        graphs = []
        for i, ana in enumerate(self.collection.values()):
            h = ana.draw_signal_map(show=False)
            h.SetContour(50)
            graphs.append(h)
        # find min/max
        glob_min = int(min([gr.GetMinimum() for gr in graphs])) / 5 * 5
        glob_max = (int(max([gr.GetMaximum() for gr in graphs])) + 5) / 5 * 5
        print glob_min, glob_max

        c = TCanvas('sig_map', 'Signal Maps', 1000, 1000)
        c.SetTheta(55)
        c.SetPhi(20)
        for i, gr in enumerate(graphs):
            gr.GetZaxis().SetRangeUser(glob_min, glob_max)
            gr.Draw('surf2')
            self.save_plots('map{}'.format(i), canvas=c, sub_dir=self.save_dir, ind=i)
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        # return graphs

    def draw_signal_spreads(self, flux=True, draw=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        gr = self.make_tgrapherrors('gr', 'Relative Spread vs {mode}'.format(mode=mode))
        i = 0
        for key, ana in self.collection.iteritems():
            rel_spread = ana.calc_signal_spread()
            x = ana.run.flux if flux else key
            gr.SetPoint(i, x, rel_spread[0])
            gr.SetPointError(i, 0, rel_spread[1])
            i += 1
        if draw:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='Relative Spread [%]')
        gr.Draw('ap')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('RelativeSpread', canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        gROOT.SetBatch(0)

    def show_peak_distribution(self, show=True):
        """
        Shows the positions of the peaks of the 2D map.
        :param show:
        """
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        # create an overall VotingHistogram
        ana = self.get_first_analysis()
        ana.draw_mean_signal_distribution(show=False)
        extrema = Extrema2D(ana.SignalMapHisto, ana.MeanSignalHisto)
        self.PeakDistribution = extrema.create_voting_histo()
        for run, ana in self.collection.iteritems():
            if ana.IsAligned:
                ana.find_2d_extrema(histo=self.PeakDistribution, show=False)
            else:
                print 'Run {run} is not aligned...'.format(run=run)
        c = TCanvas('c', 'Voting Histos', 1600, 800)
        c.Divide(2, 1)
        # new_pal = ar.array('i', [kYellow, kYellow, kOrange, kOrange - 3, kOrange + 7, kRed])
        ex = [TExec('ex1', 'gStyle->SetPalette(56);'), TExec('ex2', 'gStyle->SetPalette(51)')]
        if show:
            gROOT.SetBatch(0)
        for i, histo in enumerate(self.PeakDistribution.itervalues(), 1):
            c.cd(i)
            histo.Draw('col')
            ex[i - 1].Draw()
            histo.Draw('colz same')
        self.save_plots('PeakDistribution', sub_dir=self.save_dir)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.RootObjects.append([c, ex])

    # endregion

    # ====================================================================================
    # region BEAM PROFILE
    def draw_beam_info(self, mean=True, flux=True, show=True, direction='x', fit_margin=.6):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        title = 'Mean' if mean else 'Sigma'
        gr = self.make_tgrapherrors('gr', '{tit} of the Beam Profile in {dir}'.format(tit=title, dir=direction.title()))
        i = 0
        for key, ana in self.collection.iteritems():
            x = ana.run.flux if flux else key
            fit = ana.fit_beam_profile(direction, show=False, fit_margin=fit_margin)
            par = 1 if mean else 2
            gr.SetPoint(i, x, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
            i += 1
        c = TCanvas('c', 'Beam Profile', 1000, 1000)
        c.SetLeftMargin(.125)
        c.SetLogx() if flux else self.do_nothing()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='{tit} [cm]'.format(tit=title), y_off=1.8)
        gr.Draw('alp')
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('BeamProfile{tit}{dir}{mar}'.format(tit=title, dir=direction.title(), mar=fit_margin * 100), sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        return gr

    def draw_xy_profiles(self, flux=True, show=True, fitx=.4, fity=.7):
        gr1 = self.draw_beam_info(mean=True, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr2 = self.draw_beam_info(mean=True, flux=flux, show=False, direction='y', fit_margin=fity)
        gr3 = self.draw_beam_info(mean=False, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr4 = self.draw_beam_info(mean=False, flux=flux, show=False, direction='y', fit_margin=fity)
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulse Height Distribution', 1500, 1500)
        c.Divide(2, 2)
        for i, gr in enumerate([gr1, gr2, gr3, gr4], 1):
            self.format_histo(gr, y_off=1.3)
            pad = c.cd(i)
            pad.SetLogx() if flux else self.do_nothing()
            pad.SetBottomMargin(.15)
            gr.Draw('alp')
        gROOT.SetBatch(0)
        self.save_plots('BeamProfileOverview', sub_dir=self.save_dir)
        self.RootObjects.append([c, gr1, gr2, gr3, gr4])

    def draw_beam_profiles(self, show=True, direction='x'):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [ana.draw_beam_profile(mode=direction, show=False, fit=False) for ana in self.collection.itervalues()]
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        legend = TLegend(.4, .6 - self.get_number_of_analyses() * 0.03, .6, .6)
        for i, h in enumerate(histos):
            h.SetStats(0)
            self.normalise_histo(h)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux()) + '^{2}', 'l')
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.save_plots('AllBeamProfiles{mod}'.format(mod=direction.title()), sub_dir=self.save_dir)
        self.RootObjects.append([c, legend] + histos)

    # endregion

    # ====================================================================================
    # region TRACKS
    def show_chi2s(self, mode=None, show=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        self.reset_colors()
        histos = [ana.show_chi2(mode=mode, show=False) for ana in self.collection.itervalues()]
        yq = zeros(1)
        histos[0].GetQuantiles(1, yq, array([.9]))
        legend = self.make_legend(nentries=self.get_number_of_analyses())
        stack = THStack('hx2', '#chi^{{2}}{mode}'.format(mode=' in ' + mode if mode is not None else ''))
        for i, h in enumerate(histos):
            self.format_histo(h, stats=0, color=self.get_color(), lw=2)
            h.GetXaxis().SetRangeUser(0, yq[0])
            self.normalise_histo(h, to100=True)
            stack.Add(h)
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux()), 'l')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.format_histo(stack, x_tit='#chi^{2}', y_tit='Probability [%]', y_off=1.5, draw_first=True)
        stack.GetXaxis().SetRangeUser(0, yq[0])
        mode = '' if mode is None else mode
        self.RootObjects.append(self.save_histo(stack, 'AllChi2{mod}'.format(mod=mode.upper()), show, self.save_dir, lm=.15, draw_opt='nostack', l=legend))

    def draw_all_chi2s(self, show=True):
        self.show_chi2s(show=show)
        self.show_chi2s('x', show)
        self.show_chi2s('y', show)

    def show_angles(self, mode='x', show=True):
        self.reset_colors()
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [ana.show_angle(mode=mode, show=False) for ana in self.collection.itervalues()]
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')

        legend = self.make_legend(nentries=self.get_number_of_analyses())
        stack = THStack('has', 'Track Angles in {mode}'.format(mode=mode.title()))
        for i, h in enumerate(histos):
            self.format_histo(h, stats=0, color=self.get_color())
            self.normalise_histo(h, to100=True)
            stack.Add(h)
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux()), 'l')
        self.format_histo(stack, x_tit='Angle [deg]', y_tit='Probability [%]', y_off=1.5, draw_first=True)
        stack.GetXaxis().SetRangeUser(-3, 4)
        self.RootObjects.append(self.save_histo(stack, 'AllTrackAngles{mod}'.format(mod=mode.title()), show, self.save_dir, lm=.15, draw_opt='nostack', l=legend))

    def draw_both_angles(self, show=True):
        self.show_angles('x', show)
        self.show_angles('y', show)

    def show_angle_peaks(self, mode='x', sigma=False, flux=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        tit = '{mod} of Track Angle Distribution in {dir}'.format(mod='Sigma' if sigma else 'Mean', dir=mode.title())
        gr = self.make_tgrapherrors('gr', tit)
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.calc_angle_fit(mode, show=False)
            x = ana.run.flux if flux else key
            par = 2 if sigma else 1
            gr.SetPoint(i, x, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
            i += 1
        c = TCanvas('c', 'Angle Peaks', 1000, 1000)
        c.SetLeftMargin(.12)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod='Flux' if flux else 'Run', unit=' [kHz/cm2]' if flux else ''), y_tit='Sigma [deg]' if sigma else 'Mean [deg]', y_off=1.5)
        gr.Draw('ap')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('TrackAngle{0}s{mod}'.format('Sigma' if sigma else 'Mean', mod=mode.title()), canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])

    # endregion

    def make_signal_analysis(self, saveplots=True):
        """
        Run all available signal analyises together and plot them in an overview.
        :param saveplots:
        """
        start_time = time()
        self.draw_pulse_heights(show=False)
        self.draw_pedestals(show=False)
        self.draw_mean_fwhm(draw=False)
        c = TCanvas('c', 'overview', 800, 1000)
        c.Divide(1, 3)
        plots = [self.PulseHeight, self.Pedestal, self.FWHM]
        for i, plot in enumerate(plots, 1):
            pad = c.cd(i)
            pad.SetLogx()
            plot.Draw()
        if saveplots:
            self.save_plots('Overview Plot', sub_dir=self.save_dir, canvas=c)
        self.RootObjects.append(c)

        print '\nThe preanalysis for this selection took', self.print_elapsed_time(start_time)

    def select_runs_in_range(self, start, stop):
        new_collection = OrderedDict()
        for key, ana in self.collection.iteritems():
            if start <= key <= stop:
                new_collection[key] = ana
        if not new_collection:
            print 'You did not select any run! No changes were made!'
        else:
            self.collection = new_collection

    def set_channel(self, ch):
        """
        Sets the channels to be analysed by the SignalAnalysisobjects.
        :param ch: int (0 || 3)
        """
        for ana in self.collection.values():
            ana.set_channel(ch)

    def get_fluxes(self):
        flux = OrderedDict()
        for key, ana in self.collection.iteritems():
            flux[key] = ana.run.get_flux()
        return flux

    def get_first_analysis(self):
        return self.collection.values()[0]

    def get_last_analysis(self):
        return self.collection.values()[-1]

    def get_run_numbers(self):
        """ :return: sorted list of run numbers in AnalysisCollection instance """
        return sorted(self.collection.keys())

    def get_number_of_analyses(self):
        """ :return: number of analyses that the analysis collection object contains """
        return len(self.collection)

    def show_information(self):
        print "ANALYSIS COLLECTION INFO:"
        self.get_first_analysis().print_info_header()
        for ana in self.collection.itervalues():
            ana.print_information(header=False)

    @staticmethod
    def make_x_tit(mode, flux):
        return '{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else '')


if __name__ == "__main__":
    st = time()
    main_parser = ArgumentParser()
    main_parser.add_argument('runplan', nargs='?', default=3)
    main_parser.add_argument('dia', nargs='?', default=1, type=int)
    main_parser.add_argument('-tc', '--testcampaign', nargs='?', default='')
    args = main_parser.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    run_plan = args.runplan
    diamond = args.dia
    a = Elementary(tc)
    sel = RunSelection(testcampaign=tc)
    sel.select_runs_from_runplan(run_plan)
    message = 'STARTING PAD-ANALYSIS COLLECTION OF RUNPLAN {0}'.format(run_plan)
    print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(message)) * '=', msg=message)
    a.print_testcampaign()

    z = AnalysisCollection(sel, diamond)
    z.print_elapsed_time(st, 'Instantiation')
