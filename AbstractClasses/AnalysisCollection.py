# ==============================================
# IMPORTS
# ==============================================
import ROOT
from ROOT import gROOT, TCanvas, TLegend, TExec
import os
import json
import numpy as np
from numpy import sqrt, log, array, zeros
from time import time
from collections import OrderedDict
from PadAnalysis import SignalAnalysis
from ConfigParser import ConfigParser
from argparse import ArgumentParser
from Elementary import Elementary
from Extrema import Extrema2D
from TelescopeAnalysis import Analysis
from RunSelection import RunSelection


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
        self.save_dir = '{tc}_Runplan{plan}_{dia}'.format(tc=self.TESTCAMPAIGN[2:], plan=self.run_plan, dia=self.diamond_name)
        self.canvases = {}
        self.histos = {}
        # important plots
        self.FWHM = None
        self.PulseHeight = None
        self.Pedestal = None
        self.PeakDistribution = None

    def __del__(self):
        print "deleting AnalysisCollection..."
        for runnumber in self.collection.keys():
            print "in AnalysisCollection.__del__ : deleting Analysis of Run ", runnumber
            self.collection[runnumber].__del__()
        for obj in [self.PulseHeight, self.Pedestal, self.FWHM, self.PeakDistribution]:
            self.del_rootobj(obj)
        print "AnalyisCollection deleted"

    # ============================================
    # region INIT
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
        parser = ConfigParser()
        parser.read('Configuration/RunConfig_' + self.TESTCAMPAIGN + '.cfg')
        keydict = ConfigParser()
        keydict.read('Configuration/KeyDict_{tc}.cfg'.format(tc=self.TESTCAMPAIGN))
        path = parser.get('BASIC', 'runinfofile')
        flux_name = keydict.get('KEYNAMES', 'measured flux')
        f = open(path, 'r')
        run_log = json.load(f)
        fluxes = {}
        for run in self.runs:
            flux = run_log[str(run)][flux_name]
            print run, flux
            fluxes[flux] = run
        min_flux = min(fluxes)
        max_flux = max(fluxes)
        return {'min': fluxes[min_flux], 'max': fluxes[max_flux]}

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

    # ============================================
    # region SIGNAL/PEDESTAL
    def draw_pulse_heights(self, binning=20000, flux=True, raw=False, all_corr=False, draw=True):
        legend = TLegend(0.79, 0.13, 0.98, .34)
        legend.SetName('l1')
        mode = 'Flux' if flux else 'Run'
        prefix = 'Pulse Height {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.bias)
        gr1 = self.make_tgrapherrors('eventwise', prefix + 'eventwise correction', self.get_color())
        gr2 = self.make_tgrapherrors('binwise', prefix + 'binwise correction', self.get_color())
        gr3 = self.make_tgrapherrors('mean ped', prefix + 'mean correction', self.get_color())
        gr4 = self.make_tgrapherrors('raw', prefix + 'raw', self.get_color())

        gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        i = 0
        for key, ana in self.collection.iteritems():
            print 'getting ph for run', key
            fit = ana.draw_pulse_height(binning, show=False)
            fit1 = ana.draw_pulse_height(binning, show=False)
            fit2 = ana.draw_pulse_height(binning, eventwise_corr=True, show=False)
            ped = ana.show_pedestal_histo(draw=False)
            x = ana.run.flux if flux else key
            gr1.SetPoint(i, x, fit1.Parameter(0))
            gr2.SetPoint(i, x, fit2.Parameter(0))
            gr3.SetPoint(i, x, fit.Parameter(0) - ana.Polarity * ped.Parameter(1))
            gr4.SetPoint(i, x, fit.Parameter(0))
            gr1.SetPointError(i, 0, fit1.ParError(0))
            gr2.SetPointError(i, 0, fit2.ParError(0))
            gr3.SetPointError(i, 0, fit.ParError(0) + ped.ParError(1))
            gr4.SetPointError(i, 0, fit.ParError(0))
            i += 1
        if draw:
            gROOT.SetBatch(0)
            gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        graphs = [gr1]
        if all_corr:
            graphs += [gr2, gr3]
        if raw:
            graphs.append(gr4)
        c = TCanvas('c1', 'ph', 1000, 1000)
        c.SetLeftMargin(0.14)
        if flux:
            c.SetLogx()
        for i, gr in enumerate(graphs):
            self.histos[i] = gr
            legend.AddEntry(gr, gr.GetName(), 'lp')
            if not i:
                self.format_histo(gr, title=prefix, color=None, x_tit='Flux [kHz/cm2]', y_tit='Pulse Height [au]', y_off=2)
                gr.Draw('alp')
            else:
                gr.Draw('lp')
        self.histos['legend'] = legend
        if all_corr or raw:
            legend.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        self.save_plots('PulseHeight_' + mode, 'png', canvas=c, sub_dir=self.save_dir)
        self.save_plots('PulseHeight_' + mode, 'root', canvas=c, sub_dir=self.save_dir)
        self.canvases[0] = c
        self.PulseHeight = gr1

    def draw_pedestals(self, region='ab', peak_int='2', flux=True, all_regions=False, sigma=False, draw=True):
        legend = TLegend(0.7, 0.3, 0.98, .7)
        legend.SetName('l1')
        mode = 'Flux' if flux else 'Run'
        y_val = 'Sigma' if sigma else 'Mean'
        prefix = '{y} of Pedestal {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.diamond_name, bias=self.bias, y=y_val)
        gr1 = self.make_tgrapherrors('pedestal', prefix + 'in {reg}'.format(reg=region + peak_int))
        graphs = []
        regions = self.get_first_analysis().run.pedestal_regions
        for reg in regions:
            graphs.append(self.make_tgrapherrors('pedestal', prefix + 'in {reg}'.format(reg=reg + peak_int), color=self.get_color()))
        gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        i = 0
        par = 2 if sigma else 1
        for key, ana in self.collection.iteritems():
            print 'getting pedestal for run {n}...'.format(n=key)
            fit_par = ana.show_pedestal_histo(region, peak_int)
            flux = ana.run.flux
            x = ana.run.flux if flux else key
            gr1.SetPoint(i, x, fit_par.Parameter(par))
            gr1.SetPointError(i, 0, fit_par.ParError(par))
            if all_regions:
                for reg, gr in zip(regions, graphs):
                    fit_par = ana.show_pedestal_histo(reg, peak_int)
                    gr.SetPoint(i, x, fit_par.Parameter(par))
                    gr.SetPointError(i, 0, fit_par.ParError(par))
            i += 1
        if draw:
            gROOT.SetBatch(0)
            gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        c = TCanvas('c', 'Pedestal vs Run', 1000, 1000)
        c.SetLeftMargin(0.14)
        if flux:
            c.SetLogx()
        gr1.Draw('ap')
        if all_regions:
            for i, gr in enumerate(graphs):
                legend.AddEntry(gr, str(regions.values()[i]), 'p')
                self.format_histo(gr, title=prefix, color=None, x_tit='Flux [kHz/cm2]', y_tit='Pulse Height [au]', y_off=2)
                self.histos[i + 1] = gr
                if not i:
                    gr.Draw('ap')
                else:
                    gr.Draw('p')
            legend.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.Pedestal = gr1
        self.histos['legend'] = legend
        self.canvases[0] = c
        self.save_plots('Pedestal_' + mode, 'png', canvas=self.canvases[0], sub_dir=self.save_dir)
        self.save_plots('Pedestal_' + mode, 'root', canvas=self.canvases[0], sub_dir=self.save_dir)

    def draw_signal_distributions(self, show=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1) if not show else self.do_nothing()
        histos = [ana.show_signal_histo(show=False) for ana in self.collection.itervalues()]
        c = TCanvas('c', 'Signal Distributions', 1000, 1000)
        c.SetLeftMargin(.13)
        legend = TLegend(.7, .8 - self.get_number_of_analyses() * 0.03, .9, .9)
        for i, h in enumerate(histos):
            h.SetStats(0)
            h.Scale(1 / h.GetMaximum())
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux()) + '^{2}', 'l')
        legend.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('SignalDistributions', sub_dir=self.save_dir)
        self.histos[0] = [c, histos, legend]

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
        self.canvases[0] = c
        self.histos[0] = gr
        gROOT.SetBatch(0)

    # endregion

    # ============================================
    # region PULSER
    def draw_pulser_info(self, flux=True, show=True, mean=True, corr=True, at_jumps=False):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        title = '{mean} of Pulser vs {mod} ({ped}, {beam})'.format(mean='Mean' if mean else 'Sigma', mod=mode, ped='pedcorrected' if corr else 'uncorrected', beam='BeamOff' if at_jumps else 'BeamOn')
        gr = self.make_tgrapherrors('gr', title)
        i = 0
        for key, ana in self.collection.iteritems():
            x = ana.run.flux if flux else key
            fit = ana.calc_pulser_fit(show=False, corr=corr, at_jumps=at_jumps)
            par = 1 if mean else 2
            gr.SetPoint(i, x, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
            i += 1
        c = TCanvas('c', 'Pulser Overview', 1000, 1000)
        c.SetLeftMargin(.125)
        c.SetLogx() if flux else self.do_nothing()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='{mean} [au]'.format(mean='Mean' if mean else 'Sigma'), y_off=1.8)
        gr.Draw('alp')
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('Pulser{mean}{a}{b}'.format(mean='Mean' if mean else 'Sigma', a=corr, b=at_jumps), sub_dir=self.save_dir)
        self.histos[0] = [c, gr]
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
        gr2 = self.make_tgrapherrors('gr2', '{pref} vs {mod}'.format(pref=prefix,  mod=mode), color=self.get_color())
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
        self.histos[0] = [c, gr1, gr2, gr3, leg]
    # endregion

    # ============================================
    # region PEAK VALUES
    def draw_signal_peaks(self, flux=True, draw=True):
        """
        Shows the means of the signal peak distribution.
        :param flux:
        :param draw:
        """
        mode = 'Flux' if flux else 'Run'
        prefix = 'Mean of Signal Peaks: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.bias)
        gr = self.make_tgrapherrors('gr', prefix)
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.fit_peak_values(draw=False)
            x = ana.run.flux if flux else key
            gr.SetPoint(i, x, fit.Parameter(1))
            gr.SetPointError(i, 0, fit.ParError(1))
            i += 1
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''))
        c = TCanvas('c', 'Mean of Signal Peaks', 1000, 1000)
        gr.Draw('alp')
        self.canvases = c
        if not draw:
            c.Close()
        self.histos[0] = gr

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
        self.canvases = c
        if not draw:
            c.Close()
        self.histos[0] = gr

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
        self.canvases[0] = c
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
        self.canvases[0] = c
        self.histos[0] = gr
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
        self.histos[0] = ex
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.canvases[0] = c
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
        self.histos[0] = [c, gr]
        return gr

    def draw_xy_profiles(self, flux=True, show=True, fitx=.4, fity=.7):
        gr1 = self.draw_beam_info(mean=True, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr2 = self.draw_beam_info(mean=True, flux=flux, show=False, direction='y', fit_margin=fity)
        gr3 = self.draw_beam_info(mean=False, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr4 = self.draw_beam_info(mean=False, flux=flux, show=False, direction='y', fit_margin=fity)
        c = TCanvas('c', 'Pulse Height Distribution', 1500, 1500)
        c.Divide(2, 2)
        for i, gr in enumerate([gr1, gr2, gr3, gr4], 1):
            self.format_histo(gr, y_off=1.3)
            pad = c.cd(i)
            pad.SetLogx() if flux else self.do_nothing()
            pad.SetBottomMargin(.15)
            gr.Draw('alp')
        self.save_plots('BeamProfileOverview', sub_dir=self.save_dir)
        self.histos[1] = [c, gr1, gr2, gr3, gr4]

    def draw_beam_profiles(self, show=True, direction='x'):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [ana.draw_beam_profile(mode=direction, show=False, fit=False) for ana in self.collection.itervalues()]
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
        self.save_plots('AllBeamProfiles{mod}'.format(mod=direction.title()), sub_dir=self.save_dir)
        self.histos[0] = [c, legend] + histos

    # endregion

    # ====================================================================================
    # region TRACKS
    def show_chi2s(self, mode=None):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [ana.show_chi2(mode=mode, show=False) for ana in self.collection.itervalues()]
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        histos[0].SetStats(0)
        yq = zeros(1)
        histos[0].GetQuantiles(1, yq, array([.9]))
        legend = TLegend(.7, .8 - self.get_number_of_analyses() * 0.03, .9, .9)
        for i, h in enumerate(histos):
            h.GetXaxis().SetRangeUser(0, yq[0])
            self.normalise_histo(h)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux()) + '^{2}', 'l')
            self.histos[i] = h
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        mode = '' if mode is None else mode
        self.save_plots('AllChi2{mod}'.format(mod=mode.upper()), canvas=c, sub_dir=self.save_dir)
        self.canvases[0] = c
        self.histos['legend'] = legend

    def show_angles(self, mode='x'):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [ana.show_angle(mode=mode, show=False) for ana in self.collection.itervalues()]
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        legend = TLegend(.7, .8 - self.get_number_of_analyses() * 0.03, .9, .9)
        for i, h in enumerate(histos):
            h.SetStats(0)
            self.normalise_histo(h)
            h.SetLineColor(self.get_color())
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux()) + '^{2}', 'l')
            self.histos[i] = h
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('AllTrackAngles{mod}'.format(mod=mode.upper()), sub_dir=self.save_dir)
        self.canvases[0] = c
        self.histos['legend'] = legend

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
        self.histos[0] = [gr, c]

    # endregion

    def make_signal_analysis(self, saveplots=True):
        """
        Run all available signal analyises together and plot them in an overview.
        :param saveplots:
        """
        start_time = time()
        self.draw_pulse_heights(draw=False)
        self.draw_pedestals(draw=False)
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
        self.canvases[0] = c

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

    def PeakSignalEvolution(self, channel, NMax=3, NMin=3, OnThisCanvas=None, BinRateEvolution=False):
        '''
        Shows a rate scan of individual bins. For the plot NMax maxima
        and NMin minima are chosen and its mean signal evolution
        is shown as a function of run number (i.e. rate).
        :param channel:
        :param NMax:
        :param NMin:
        :param OnThisCanvas:
        :param BinRateEvolution:
        :return:
        '''
        if self.get_number_of_analyses() == 0: return 0

        if OnThisCanvas != None:
            assert (isinstance(OnThisCanvas, ROOT.TCanvas)), "OnThisCanvas has to be a TCanvas object"
        print "Signal Evolution start"
        if not hasattr(self, "PeakPadMax"):
            self.PeakComparison(channel=channel, show=False)

        def BinsAreNearby(x1, y1, x2, y2, R):
            d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2
            if d2 <= R ** 2:
                return True
            else:
                return False

        # Find the separated peaks / lows (binnumbers) to consider in Signal Evolution
        def FindPeakBins(PeakPad, N, maximum=False):
            self.PeakPadMaxPad.CalculateMeanSignalDistribution()
            peakbins = [-1] * N  # container to store the binnumbers of the separated maximas found
            peakPad = PeakPad  # copy.deepcopy(PeakPad)

            PeakPad2 = self.PeakPadMaxPad.meansignaldistribution  # peaksearch due to mean signal content in bins
            peakPad2 = PeakPad2  # copy.deepcopy(PeakPad2)
            i = 0
            while i < int(N):
                if i == 3 and maximum:
                    peakPad = peakPad2

                maxcount = peakPad.GetMaximum()  # counts of maximum

                if maxcount < 1:
                    break
                peakbins[i] = peakPad.GetMaximumBin()  # binnumber with hightest counts
                coordinates = peakPad.GetBinCenter(peakbins[i])
                peakPad.Fill(coordinates[0], coordinates[1], -maxcount)  # remove content of maximum bin

                # if the binnumber is already in a neighborhood of a found peak, don't use it:
                IsInNBHD = False
                for k in xrange(i):
                    # IsInNBHD |= peakbins[k] in peakPad.GetBinsInNbhd(peakbins[i], include_center=True, extended=True)
                    other_coordinates = peakPad.GetBinCenter(peakbins[k])
                    IsInNBHD |= BinsAreNearby(coordinates[0], coordinates[1], other_coordinates[0], other_coordinates[1], 0.05)
                if IsInNBHD:
                    pass
                else:
                    i += 1
            self.verbose_print("Number of separated extrema to look at: {0:0.0f}".format(i))
            peakbins = peakbins[:i]
            return peakbins

        peakbins = FindPeakBins(self.PeakPadMax, NMax, maximum=True)
        lowbins = FindPeakBins(self.PeakPadMin, NMin)

        # Rate Time Evolution of Bin:
        if BinRateEvolution:
            tmpSaveDir = self.save_directory

            high1 = peakbins[0]
            high2 = peakbins[1]
            high3 = peakbins[2]

            low1 = lowbins[0]
            low2 = lowbins[1]
            low3 = lowbins[2]

            for run_number in self.collection.keys():
                self.set_save_directory(tmpSaveDir + str(run_number) + "/")
                self.collection[run_number].RateTimeEvolution(time_spacing=5, save=True, binnumber=high1, nameExtension="High1")
                self.collection[run_number].RateTimeEvolution(time_spacing=5, save=True, binnumber=high2, nameExtension="High2")
                self.collection[run_number].RateTimeEvolution(time_spacing=5, save=True, binnumber=high3, nameExtension="High3")

                self.collection[run_number].RateTimeEvolution(time_spacing=5, save=True, binnumber=low1, nameExtension="Low1")
                self.collection[run_number].RateTimeEvolution(time_spacing=5, save=True, binnumber=low2, nameExtension="Low2")
                self.collection[run_number].RateTimeEvolution(time_spacing=5, save=True, binnumber=low3, nameExtension="Low3")

            self.set_save_directory(tmpSaveDir)

        # Fill all graphs of all separated peaks / lows
        def FillGraphDict(self, GraphDict, ListOfBins):
            runnumbers = self.collection.keys()
            runnumbers.sort()
            for peakbin in ListOfBins:
                GraphDict[peakbin] = ROOT.TGraphErrors()
                GraphDict[peakbin].SetNameTitle("MaxGraph_" + str(peakbin), "Evolution of Signal Response during Rate Scan")
                # signals = []
                i = 0
                for runnumber in runnumbers:
                    self.collection[runnumber].Pads[channel].listOfBins[peakbin].CreateBinSignalHisto(saveplot=True, savedir=self.SaveDirectory + str(runnumber) + "/", show_fit=False)
                    mean = self.collection[runnumber].Pads[channel].listOfBins[peakbin].BinSignalHisto.GetMean()
                    error = self.collection[runnumber].Pads[channel].listOfBins[peakbin].BinSignalHisto.GetRMS() / np.sqrt(
                        self.collection[runnumber].Pads[channel].listOfBins[peakbin].BinSignalHisto.GetEntries())
                    # mpv = self.collection[runnumber].Pads[channel].listOfBins[peakbin].Fit['MPV']
                    # signals.append(mpv)
                    GraphDict[peakbin].SetPoint(i, runnumber, mean)
                    GraphDict[peakbin].SetPointError(i, 0, error)
                    i += 1

        MaxGraphs = {}
        FillGraphDict(self, MaxGraphs, peakbins)
        MinGraphs = {}
        FillGraphDict(self, MinGraphs, lowbins)
        theseMaximas = []
        theseMinimas = []

        # Prepare for drawing: Settings, create Canvas, create Legend
        if len(MaxGraphs) > 0:
            marker = 20
            npeaks = len(MaxGraphs)
            PeakSignalEvolutionCanvas = ROOT.gROOT.GetListOfCanvases().FindObject("PeakSignalEvolutionCanvas")
            if not PeakSignalEvolutionCanvas:
                PeakSignalEvolutionCanvas = ROOT.TCanvas("PeakSignalEvolutionCanvas", "Signal Evolution Canvas")
            PeakSignalEvolutionCanvas.cd()
            legend = ROOT.TLegend(0.1, 0.1, 0.3, 0.35)

            # determine the signal range for y axis:
            MaxSignals = []
            MinSignals = []
            for peakbin in peakbins:
                MaxSignals.append(MaxGraphs[peakbin].GetYaxis().GetXmax())
                MinSignals.append(MaxGraphs[peakbin].GetYaxis().GetXmin())
            MaxRange_peak = 1.1 * np.array(MaxSignals).max()
            MinRange_peak = 0.9 * np.array(MinSignals).min()

            for peaknr in xrange(npeaks):
                MaxGraphs[peakbins[peaknr]].SetMarkerStyle(marker)
                MaxGraphs[peakbins[peaknr]].SetMarkerColor(ROOT.kRed)
                MaxGraphs[peakbins[peaknr]].SetLineColor(ROOT.kRed)
                MaxGraphs[peakbins[peaknr]].Draw("SAME LP")
                legend.AddEntry(MaxGraphs[peakbins[peaknr]], "high" + str(peaknr + 1), "lp")

                theseMaximas += [self.PeakPadMax.GetBinCenter(peakbins[peaknr])]

                marker += 1
        else:
            PeakSignalEvolutionCanvas = ROOT.gROOT.GetListOfCanvases().FindObject("PeakSignalEvolutionCanvas")
            if not PeakSignalEvolutionCanvas:
                PeakSignalEvolutionCanvas = ROOT.TCanvas("PeakSignalEvolutionCanvas", "Signal Evolution Canvas")
            PeakSignalEvolutionCanvas.cd()
            legend = ROOT.TLegend(0.1, 0.1, 0.2, 0.25)
            MaxRange_peak = None
            MinRange_peak = None

        if len(MinGraphs) > 0:
            marker = 20
            nlows = len(MinGraphs)

            # determine the signal range for y axis:
            MaxSignals = []
            MinSignals = []
            for lowbin in lowbins:
                MaxSignals.append(MinGraphs[lowbin].GetYaxis().GetXmax())
                MinSignals.append(MinGraphs[lowbin].GetYaxis().GetXmin())
            MaxRange_low = 1.1 * np.array(MaxSignals).max()
            MinRange_low = 0.9 * np.array(MinSignals).min()

            for lownr in xrange(nlows):
                MinGraphs[lowbins[lownr]].SetMarkerStyle(marker)
                MinGraphs[lowbins[lownr]].SetMarkerColor(ROOT.kBlue)
                MinGraphs[lowbins[lownr]].SetLineColor(ROOT.kBlue)
                # MinGraphs[lowbins[lownr+1]].Draw("SAME LP")
                legend.AddEntry(MinGraphs[lowbins[lownr]], "low" + str(lownr + 1), "lp")

                theseMinimas += [self.PeakPadMin.GetBinCenter(lowbins[lownr])]

                marker += 1
        else:
            MaxRange_low = None
            MinRange_low = None

        # Prepare for drawing: Evaluate the Range in y direction:
        MaxRange = np.array([i for i in [MaxRange_low, MaxRange_peak] if i != None])
        MinRange = np.array([i for i in [MinRange_low, MinRange_peak] if i != None])
        if len(MaxRange) > 0:
            MaxRange = MaxRange.max()
        else:
            MaxRange = 200
        if len(MinRange) > 0:
            MinRange = 0.8 * MinRange.min()
        else:
            MinRange = 0

        # Prepare for drawing: Individual Print options:
        NumbersOfGraphs = len(MaxGraphs) + len(MinGraphs)
        DrawOptions = ["SAME LP"] * NumbersOfGraphs
        try:
            DrawOptions[0] = "ALP"
        except IndexError:  # if neither maxima nor minima found
            pass

        # Prepare for drawing rate:
        runnumbers = self.collection.keys()
        runnumbers.sort()
        first = runnumbers[0]
        last = runnumbers[-1]
        ratebins = last - first + 1
        RateHisto = ROOT.TH1D("RateHisto", "Rate Histogram", ratebins, first - 0.5, last + 0.5)
        for runnumber in runnumbers:
            rate_kHz = self.collection[runnumber].get_flux()
            print "runnumber: ", self.collection[runnumber].run.run_number, " == ", runnumber, " rate_kHz: ", rate_kHz
            assert (self.collection[runnumber].run.run_number == runnumber)
            RateHisto.Fill(runnumber, rate_kHz)
        RateHisto.GetXaxis().SetTitle("Run Number")
        RateHisto.GetYaxis().SetTitle("Rate / kHz")

        # Draw everything:
        i = 0  # i-th draw option
        if len(MaxGraphs) > 0:
            MaxGraphs[peakbins[0]].GetYaxis().SetRangeUser(MinRange, MaxRange)
            MaxGraphs[peakbins[0]].GetXaxis().SetTitle("Run Number")
            MaxGraphs[peakbins[0]].GetYaxis().SetTitle("Mean Signal Response")
            MaxGraphs[peakbins[0]].Draw(DrawOptions[i])
            i += 1
            for peaknr in xrange(npeaks - 1):
                MaxGraphs[peakbins[peaknr + 1]].Draw(DrawOptions[i])
                i += 1
        if len(MinGraphs) > 0:
            if i == 0:
                MinGraphs[lowbins[0]].GetYaxis().SetRangeUser(MinRange, MaxRange)
                MinGraphs[lowbins[0]].GetXaxis().SetTitle("Run Number")
                MinGraphs[lowbins[0]].GetYaxis().SetTitle("Mean Signal Response")
            MinGraphs[lowbins[0]].Draw(DrawOptions[i])
            i += 1
            for lownr in xrange(nlows - 1):
                MinGraphs[lowbins[lownr + 1]].Draw(DrawOptions[i])
                i += 1
        legend.Draw()
        self.save_plots("PeakSignalEvolution.png")
        self.save_plots("PeakSignalEvolution.root")
        raw_input("waiting in AnalysisCollection->Line 566")
        pad = PeakSignalEvolutionCanvas.GetPad(0)
        RateHisto.SetStats(0)
        RateHisto.Draw("SAME HIST Y+")  # include in plot instead of second plot
        pad.SetLogy()
        self.save_plots("PeakSignalEvolution_Rate.png")

        # show the selected bins in another canvas:
        if OnThisCanvas:
            ROOT.gStyle.SetPalette(53)  # Dark Body Radiator palette
            OnThisCanvas.cd(1)

            if len(MaxGraphs) > 0:
                for peaknr in xrange(npeaks):
                    maxima = self.PeakPadMax.GetBinCenter(peakbins[peaknr])
                    text = ROOT.TText()
                    text.SetTextColor(ROOT.kRed)
                    text.DrawText(maxima[0] - 0.02, maxima[1] - 0.005, 'high' + str(peaknr + 1))

            if len(MinGraphs) > 0:
                for lownr in xrange(nlows):
                    minima = self.PeakPadMin.GetBinCenter(lowbins[lownr])
                    text = ROOT.TText()
                    text.SetTextColor(ROOT.kBlue)
                    text.DrawText(minima[0] - 0.01, minima[1] - 0.005, 'low' + str(lownr + 1))

            OnThisCanvas.Update()
            self.save_plots("IIa-2_neutron_SignalDistribution_MAXSearch.png")
            raw_input("wait")
        print "Highs: ", theseMaximas
        print "Lows: ", theseMinimas

    # In [4]: a = coll.collection[445]
    #
    # In [5]: a.ShowSignalMaps(False)
    # In [7]: c1 = ROOT.gROOT.FindObject("signal_canvas{run}")
    # In [8]: pad = c1.cd(1)
    # In [10]: a._DrawMinMax(pad, channel, theseMaximas, theseMinimas)
    # --> ADD number to high low labels..

if __name__ == "__main__":
    main_parser = ArgumentParser()
    main_parser.add_argument('runplan', nargs='?', default=3, type=int)
    main_parser.add_argument('dia', nargs='?', default=1, type=int)
    args = main_parser.parse_args()
    run_plan = args.runplan
    diamond = args.dia
    sel = RunSelection()
    sel.select_runs_from_runplan(run_plan)
    z = AnalysisCollection(sel, diamond)
