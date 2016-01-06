# ==============================================
# IMPORTS
# ==============================================
import ROOT
from ROOT import gROOT, TCanvas, TGraphErrors, TLegend
from AbstractClasses.ATH2D import ATH2D
from AbstractClasses.BinCollection import BinCollection
from AbstractClasses.newAnalysis import Analysis
from AbstractClasses.RunSelection import RunSelection
import os
import numpy as np
from numpy import sqrt, log
from Elementary import Elementary
from time import time
from collections import OrderedDict
from signal_analysis import SignalAnalysis
from ConfigParser import ConfigParser
import json
from argparse import ArgumentParser


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
        self.lowest_rate_run = self.get_lowest_rate_run()

        self.generate_slope_pickle()

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

    def __del__(self):
        print "deleting AnalysisCollection.."
        for runnumber in self.collection.keys():
            print "in AnalysisCollection.__del__ : deleting Analysis of Run ", runnumber
            self.collection[runnumber].__del__()
        if hasattr(self, "FWHMcanvas"):
            canvas = ROOT.gROOT.FindObject("FWHMcanvas")
            if canvas:
                canvas.Close()
                del canvas
        if hasattr(self, "fwhm_histo"):
            ROOT.gROOT.Delete("fwhm_histo")
        print "AnalyisCollection deleted"

    # ============================================
    # region INIT
    def add_analyses(self):
        """
        Creates and adds Analysis objects with run numbers in runs.
        """
        for run, dia in sorted(zip(self.runs, self.diamonds)):
            ch = 0 if dia == 1 or dia == 3 else 3
            analysis = SignalAnalysis(run, ch, self.lowest_rate_run)
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

    def get_lowest_rate_run(self):
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
        return fluxes[min_flux]

    def generate_slope_pickle(self):
        picklepath = 'Configuration/Individual_Configs/Slope/{tc}_{run}_{ch}_Slope.pickle'.format(tc=self.TESTCAMPAIGN, run=self.lowest_rate_run, ch=0)
        if os.path.exists(picklepath):
            return
        Analysis(self.lowest_rate_run)

    @staticmethod
    def make_runselection(run_list):
        assert type(run_list) is list, 'run list argument has to be a list!'
        selection = RunSelection()
        selection.select_runs(run_list, do_assert=True)
        return selection
    # endregion

    # ============================================
    # region ANALYSIS
    def draw_pulse_heights(self, binning=10000, flux=True, raw=False, all_corr=False, draw=True):
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
            fit = ana.draw_pulse_height(binning)
            fit1 = ana.draw_pulse_height(binning, ped_corr=True)
            fit2 = ana.draw_pulse_height(binning, eventwise_corr=True)
            ped = ana.show_pedestal_histo()
            x = ana.run.flux if flux else key
            gr1.SetPoint(i, x, fit1.GetParameters()[0])
            gr2.SetPoint(i, x, fit2.GetParameters()[0])
            gr3.SetPoint(i, x, fit.GetParameters()[0] - ana.polarity * ped.Parameter(1))
            gr4.SetPoint(i, x, fit.GetParameters()[0])
            gr1.SetPointError(i, 0, fit1.GetParError(0))
            gr2.SetPointError(i, 0, fit2.GetParError(0))
            gr3.SetPointError(i, 0, fit.GetParError(0) + ped.ParError(1))
            gr4.SetPointError(i, 0, fit.GetParError(0))
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

    def draw_mean_fwhm(self, saveplots=True, flux=True, draw=True):
        """
        Creates the FWHM Distribution of all selected MeanSignalHistograms
        :param saveplots: if True saves the plot
        :param flux: draw vs flux if True else vs run
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

        print '\nThe preanalysis for this selection took', self.elapsed_time(start_time)
    # endregion

    def select_runs_in_range(self, start, stop):
        new_collection = OrderedDict()
        for key, ana in self.collection.iteritems():
            if start <= key <= stop:
                new_collection[key] = ana
        if not new_collection:
            print 'You did not select any run! No changes were made!'
        else:
            self.collection = new_collection

    def get_first_analysis(self):
        return self.collection.values()[0]

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

    def ShowPulserRates(self):
        '''
        Execute the ShowPulserRate method for all runs (i.e. Analysis
        objects) in AnalysisCollection.
        :return:
        '''
        runnumbers = self.get_run_numbers()
        for run in runnumbers:
            self.collection[run].ShowPulserRate()

    def CreateSigmaMPVPlot(self):
        '''
        Analysis.FindMaxima() has to be executed for all contributing analysis
        :param saveplots:
        :param savename:
        :param ending:
        :return:
        '''
        if self.get_number_of_analyses() == 0: return 0

        # self.AnalysisCollection.collection[run_number].MaximaAnalysis
        canvas = ROOT.TCanvas('MPVSigmaCanvas', 'MPVSigmaCanvas')
        canvas.cd()
        graph = TGraphErrors()
        graph.SetNameTitle("graph", "MPV vs Sigma of underlying Landau")

        count = 0
        for run_number in self.collection:
            MPVs, Sigmas, MPVErrs, SigmaErrs = self.collection[run_number].GetMPVSigmas(minimum_counts=300)
            for i in xrange(len(MPVs)):
                graph.SetPoint(count, MPVs[i], Sigmas[i])
                graph.SetPointError(count, MPVErrs[i], SigmaErrs[i])
                count += 1
        graph.SaveAs("Total_MPV_Sigma_graph.root")
        graph.Draw('AP')
        ROOT.gPad.Print("Results/MPV_Sigma_graph.png")
        ROOT.gPad.Print("Results/MPV_Sigma_graph.root")
        self.if_wait("MPV vs Sigma shown...")

    def SignalHeightScan(self, channel):  # improve!
        if self.get_number_of_analyses() == 0: return 0

        # tmp = self.ShowAndWait
        # self.ShowAndWait = True
        SignalHeightScanCanvas = ROOT.TCanvas("SignalHeightScanCanvas", "SignalHeightScan Canvas")
        SignalHeightScanCanvas.cd()

        SignalHeightScanGraph = ROOT.TGraph()
        SignalHeightScanGraph.SetNameTitle("SignalHeightScanGraph", "Signal Height Scan")

        runnumbers = self.collection.keys()
        runnumbers.sort()

        count = 0
        for runnumber in runnumbers:
            if not self.collection[runnumber].TimingAlignmentFailed:
                SignalHeightScanGraph.SetPoint(count, runnumber, self.collection[runnumber].extremaResults[channel]['SignalHeight'])
                count += 1
            else:
                print "INFO: Run number {0} excluded in SignalHeightScan plot due to bad timing alignment !"
        SignalHeightScanGraph.SaveAs(self.save_directory + "SignalHeightGraph.root")
        SignalHeightScanGraph.GetXaxis().SetTitle("Run Number")
        SignalHeightScanGraph.GetYaxis().SetTitle("Reconstructed Signal Height")
        SignalHeightScanGraph.Draw("AP*")
        self.save_plots("SignalHeightGraph.png")
        self.if_wait("SignalHeightScan shown...")
        # self.ShowAndWait = tmp

    def PeakComparison(self, channel, show=True):
        if self.get_number_of_analyses() == 0: return 0

        print "PeakComparision start"
        if show:
            self.peakComparisonCanvasMax = ROOT.TCanvas("peakComparisonCanvasMax", "PeakComparisonCanvas")
            self.peakComparisonCanvasMin = ROOT.TCanvas("peakComparisonCanvasMin", "PeakComparisonCanvas")

        runnumbers = self.collection.keys()
        if not hasattr(self.collection[runnumbers[0]], "Pads"):
            self.collection[runnumbers[0]].LoadTrackData()
        pad_attributes = self.collection[runnumbers[0]].Pads[channel].Get2DAttributes()

        self.PeakPadMax = ATH2D("PeakPadMax", "Peak distribution over all selected runs", *pad_attributes)
        self.PeakPadMaxPad = BinCollection(self.collection[runnumbers[0]], channel, *pad_attributes)  # CHANGE NAME !
        self.PeakPadMin = ATH2D("PeakPadMin", "Low distribution over all selected runs", *pad_attributes)

        for runnumber in runnumbers:
            analysis = self.collection[runnumber]
            if analysis.extremaResults[channel]["FoundMaxima"] == None: analysis.FindMaxima(channel=channel, show=False)
            if analysis.extremaResults[channel]["FoundMinima"] == None: analysis.FindMinima(channel=channel, show=False)
            maxima = analysis.extremaResults[channel]["FoundMaxima"]
            minima = analysis.extremaResults[channel]["FoundMinima"]
            if maxima != None:
                for peak in maxima:
                    self.PeakPadMax.Fill(*peak)
                    if not hasattr(analysis.Pads[channel], "meansignaldistribution"):
                        analysis.Pads[channel].CalculateMeanSignalDistribution()
                    signal_ = analysis.Pads[channel].meansignaldistribution.GetBinContent(analysis.Pads[channel].GetBinNumber(*peak))
                    self.PeakPadMaxPad.Fill(peak[0], peak[1], signal_)
            else:
                print "WARNING: No Maxima results found in run ", runnumber, ". PeakComparisonMax will be incomplete."
            if minima != None:
                for peak in minima:
                    self.PeakPadMin.Fill(*peak)
            else:
                print "WARNING: No Minima results found in run ", runnumber, ". PeakComparisonMin will be incomplete."
        if show:
            ROOT.gStyle.SetPalette(55)  # Rainbow palette
            self.peakComparisonCanvasMax.cd()
            self.PeakPadMax.Draw("COLZ")
            self.save_plots("PeakPadMax.png")
            self.peakComparisonCanvasMin.cd()
            self.PeakPadMin.Draw("COLZ")
            self.save_plots("PeakPadMin.png")

            # raw_input("peakpad")

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

    def get_run_numbers(self):
        """ :return: sorted list of run numbers in AnalysisCollection instance """
        return sorted(self.collection.keys())

    def get_number_of_analyses(self):
        """ :return: number of analyses that the analysis collection object contains """
        return len(self.collection)

    def ShowInfo(self):
        print "ANALYSIS COLLECTION INFO:"
        print "\tRuns: \tDiamond1 \tBias1 \tSelected \tDiamond2 \tBias2 \tSelected \tType"
        contentstring = ""
        for run in self.get_run_numbers():
            contentstring += "\t{run} \t{Diamond1} \t{Bias1} \t{Selected1} \t\t{Diamond2} \t{Bias2} \t{Selected2} \t\t{Type}\n".format(
                run=str(run).zfill(3), Diamond1=self.collection[run].run.get_diamond_name(0).ljust(8), Bias1=str(self.collection[run].run.bias[0]).zfill(5),
                Selected1=str(self.collection[run].run.analyzeCh[0]).ljust(5),
                Diamond2=self.collection[run].run.get_diamond_name(3).ljust(8), Bias2=str(self.collection[run].run.bias[3]).zfill(5), Selected2=str(self.collection[run].run.analyzeCh[3]).ljust(5),
                Type=self.collection[run].run.RunInfo["type"])
        print contentstring

    def MakeGlobalPedestalCorrections(self, channel=None):
        for run in self.collection.keys():
            self.collection[run].MakeGlobalPedestalCorrection(channel=channel)

    def SetIndividualCuts(self, showOverview=True, savePlot=False):
        for run in self.collection.keys():
            self.collection[run].SetIndividualCuts(showOverview=showOverview, savePlot=savePlot)

    def AnalyzePedestalContribution(self, channel, normalize=False, refactor=5):
        '''
        Example:
            sel = RunSelection()
            sel.SelectRunsFromRunPlan(12)
            sel.UnSelectUnlessInRange(439, 446)
            coll = AnalysisCollection(sel)
            coll.AnalyzePedestalContribution(3)
        :param channel:
        :param refactor:
        :return:
        '''
        self.set_save_directory("Results/Pedestal_Analysis/")
        self.pedestalresults = {
            "full": ROOT.TGraph(),
            "no_tail": ROOT.TGraph(),
            "no_ped": ROOT.TGraph()
        }
        colors = {
            "full": ROOT.kRed,
            "no_tail": ROOT.kBlue,
            "no_ped": ROOT.kGreen
        }
        for key in self.pedestalresults.keys():
            self.pedestalresults[key].SetNameTitle("Scan_Mean_" + key, "Scan_Mean_" + key)
            self.pedestalresults[key].SetLineColor(colors[key])

        i = 0

        means = {
            "full": [],
            "no_tail": [],
            "no_ped": []
        }
        for run in self.collection.keys():
            self.collection[run].CalculateSNR(channel=channel, savePlots=False)
            fullsignalhisto = self.collection[run].snr_canvas.GetPrimitive(
                "{dia}_SNRSignalHisto{run}".format(dia=self.collection[run].run.diamondname[channel], run=self.collection[run].run.run_number))
            self.save_plots(savename="SignalHisto_full_{run}{ch}.png".format(run=run, ch=channel), canvas=self.collection[run].snr_canvas)
            self.save_plots(savename="SignalHisto_full_{run}{ch}.root".format(run=run, ch=channel), sub_dir="root", canvas=self.collection[run].snr_canvas)
            mean_full = fullsignalhisto.GetMean()
            mean, mean_nopedestal = self.collection[run].AnalyzePedestalContribution(channel=channel, refactor=refactor)
            self.save_plots(savename="SignalHisto_fit_{run}{ch}.png".format(run=run, ch=channel), canvas=self.collection[run].signalpedestalcanvas)
            self.save_plots(savename="SignalHisto_fit_{run}{ch}.root".format(run=run, ch=channel), sub_dir="root", canvas=self.collection[run].signalpedestalcanvas)

            means["full"] += [mean_full]
            means["no_tail"] += [mean]
            means["no_ped"] += [mean_nopedestal]

            if normalize:
                factor = mean_nopedestal
            else:
                factor = 1.

            self.pedestalresults["full"].SetPoint(i, self.collection[run].get_flux(), mean_full / factor)
            self.pedestalresults["no_tail"].SetPoint(i, self.collection[run].get_flux(), mean / factor)
            self.pedestalresults["no_ped"].SetPoint(i, self.collection[run].get_flux(), mean_nopedestal / factor)
            i += 1

        self.pedestal_analysis_canvas = ROOT.TCanvas("pedestal_analysis_canvas", "pedestal_analysis_canvas")
        self.pedestalresults["full"].Draw("ALP")
        self.pedestalresults["no_tail"].Draw("LP")
        self.pedestalresults["no_ped"].Draw("LP")
        self.pedestal_analysis_canvas.Update()

        for key in means.keys():
            print key, "  -  ", means[key]


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
