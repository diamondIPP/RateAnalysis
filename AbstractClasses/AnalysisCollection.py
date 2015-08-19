import ROOT
from AbstractClasses.ATH2D import ATH2D
from AbstractClasses.BinCollection import BinCollection
from AbstractClasses.RunClass import Run
from AbstractClasses.newAnalysis import Analysis
import types as t
import os
import numpy as np
from Elementary import Elementary
from ROOT import TGraphErrors

class AnalysisCollection(Elementary):
    '''
    An object of this class contains several analysis of runs.
    It gives the ability to compare the data from different runs.
    '''
    current_run_number = -1

    def __init__(self, listOfRuns=None, diamonds=None, verbose = False):
        Elementary.__init__(self, verbose=verbose)
        self.collection = {} # dict where all analysis objects are saved
        if listOfRuns != None:
            assert(type(listOfRuns) is t.ListType), "listOfRuns has to be of type list"
            self.AddRuns(listOfRuns, diamonds=diamonds)


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

    def AddAnalysis(self,analysis_obj):
        '''
        Adds an Analysis object to the analysis collection object
        :param analysis_obj: Analysis Object of type "Analysis"
        :return: -
        '''

        self.collection[analysis_obj.run.run_number] = analysis_obj
        self.current_run_number = analysis_obj.run.run_number

    def AddRuns(self, list, diamonds=None):
        assert(type(list) is t.ListType), "argument has to be a list of run numbers"
        if diamonds == None: diamonds=3
        assert(diamonds in [1,2,3]), "'diamonds' has to be 1, 2, 3, or None (0x1: diamond1, 0x2: diamond2)"
        for runnr in list:
            self.AddAnalysis(Analysis(Run(runnr, diamonds)))

    def SetChannels(self, channels):
        runnumbers = self.GetRunNumbers()
        for runnumber in runnumbers:
            self.collection[runnumber].run.SetChannels(channels=channels)

    def RemoveBeamInterruptions(self):
        runnumbers = self.GetRunNumbers()
        for runnumber in runnumbers:
            self.collection[runnumber].RemoveBeamInterruptions()

    def CreateFWHMPlot(self, saveplots = True, savename = 'FWHM_Histo', ending = 'png'):
        '''
        Creates the FWHM Distribution of all the MeanSignalHistogram histograms from all
        Analysis object inside the analysis collection
        :param saveplots: if True saves the plot
        :param savename:  filename if saveplots = True
        :param ending:  file typ if saveplots = True
        :return: -
        '''

        self.FWHMcanvas = ROOT.TCanvas("FWHMcanvas", "FWHM")
        self.fwhm_histo = ROOT.TH1D("fwhm_histo", "FWHM Distribution of "+str(self.GetNumberOfAnalyses())+" runs",50,0,100)

        for run in AnalysisCollection.collection:
            self.fwhm_histo.Fill(self.CalculateFWHM(print_result=False,run_number=run))
        self.FWHMcanvas.cd()
        self.fwhm_histo.GetXaxis().SetTitle('FWHM')
        self.fwhm_histo.Draw()
        self.FWHMcanvas.Update()

        if saveplots:
            # Results directories:
            resultsdir = 'Results/' # eg. 'Results/run_364/'
            if not os.path.exists(resultsdir): # if directory doesn't exist, create it!
                os.makedirs(resultsdir)

            ROOT.gPad.Print(resultsdir+savename+'.'+ending)
            ROOT.gPad.Print(resultsdir+savename+'.'+'root')

        #raw_input("wait")

    #
    def CalculateFWHM(self, print_result = True, run_number = None):
        '''
        Calculates the FWHM of the Mean Signal Histogram (Histogram of
        mean signal response of 2D Signal response distribution)
        :param print_result: if True prints FWHM in console
        :param run_number:  run number to analyze - if None it takes
                            current run number from AnalysisCollection object
        :return: FWHM
        '''
        channel = 0
        if run_number == None:
            run_number = AnalysisCollection.current_run_number
        assert(type(run_number) == t.IntType and 0 < run_number < 1000), "Invalid run number"

        analysis_obj = AnalysisCollection.collection[run_number]

        if not hasattr(analysis_obj, "MeanSignalHisto"):
            analysis_obj.CreateMeanSignalHistogram(channel=channel)

        maximum = analysis_obj.MeanSignalHisto[channel].GetMaximum()
        low_bin = analysis_obj.MeanSignalHisto[channel].FindFirstBinAbove(maximum/2.)
        high_bin = analysis_obj.MeanSignalHisto[channel].FindLastBinAbove(maximum/2.)

        fwhm = analysis_obj.MeanSignalHisto[channel].GetBinCenter(high_bin) - analysis_obj.MeanSignalHisto[channel].GetBinCenter(low_bin)

        if print_result:
            print "FWHM of run ",run_number," is: ",fwhm

        return fwhm

    def MakePreAnalysis(self, channel=None, mode="mean"):
        assert(channel in [0,3, None]), "invalid channel: channel has to be either 0, 3 or None"
        runnumbers = self.GetRunNumbers()

        for run in runnumbers:
            if channel == None:
                self.collection[run].MakePreAnalysis(mode=mode)
            else:
                self.collection[run].MakePreAnalysis(channel=channel, mode=mode)

    def ShowSignalVSRate(self, canvas=None, diamonds=None): #, method="mean"
        '''
        Draws the signal vs rate scan into the canvas. If no canvas is given, it will create
        a new canvas attached to the intance as self.ratecanvas .
        If no diamonds are selected in particular, then the active diamonds of the first run will
        be selected for the rate scan of all runs
        :param canvas: optional. A canvas to draw the ratescan into
        :param diamonds: 0x1: diamond1 0x2: diamond2
        :return:
        '''
        assert(diamonds in [1,2,3,None]), "wrong diamonds selection: 0x1: diamond1, 0x2: diamond2"
        if canvas==None:
            self.ratecanvas = ROOT.TCanvas("signalvsratecanvas", "signalvsratecanvas")
            ROOT.SetOwnership(self.ratecanvas, False)
            self.ratelegend = ROOT.TLegend(0.1, 0.1, 0.4, 0.4)
            axisoption = "A"
        else:
            self.ratecanvas = canvas
            self.ratelegend = canvas.FindObject("TPave")
            axisoption = ""

        tmpcanvas = ROOT.TCanvas("tmpcanvas", "tmpcanvas")
        tmpsignalhisto = ROOT.TH1D("tmpsignalhisto", "tmpsignalhisto", 600, -100, 500)
        runnumbers = self.GetRunNumbers()

        if diamonds == None:
            channels = self.collection[runnumbers[0]].run.GetChannels() # get channels from first run
        elif diamonds == 1:
            channels = [0]
        elif diamonds == 2:
            channels = [3]
        else:
            channels = [0,3]

        self.graphs = {}
        results = {}
        for channel in channels:
            tmpcanvas.cd()
            color = self.GetNewColor()
            self.graphs[channel] = ROOT.TGraphErrors()
            self.graphs[channel].SetNameTitle("graphCh0"+self.collection[runnumbers[0]].run.diamondname[channel], "Signal Rate Scan")
            ROOT.SetOwnership(self.graphs[channel], False)
            i = -1

            for runnumber in runnumbers:
                i += 1
                #runnumber = self.collection[runnumber].run.run_number
                results[runnumber] = {}
                print "Signal VS Rate: Processing Run {run} (Rate: {rate}) - Channel {channel}".format(run=runnumber, channel=channel, rate=self.collection[runnumber].run.RunInfo["measured flux"])
                results[runnumber][channel] = {}
                self.collection[runnumber].run.tree.Draw((self.collection[runnumber].signaldefinition+">>tmpsignalhisto").format(channel=channel), self.collection[runnumber].cut.format(channel=channel), "", 2000000, self.collection[runnumber].excludefirst)
                results[runnumber][channel]["mean"] = tmpsignalhisto.GetMean()
                results[runnumber][channel]["error"] = tmpsignalhisto.GetRMS()/np.sqrt(tmpsignalhisto.GetEntries())
                self.graphs[channel].SetPoint(i, self.collection[runnumber].run.RunInfo["measured flux"], results[runnumber][channel]["mean"])
                self.graphs[channel].SetPointError(i, 0, results[runnumber][channel]["error"])

            #save graph:
            self.graphs[channel].SaveAs(self.graphs[channel].GetName()+".root")

            self.ratecanvas.cd()
            if axisoption == "A":
                self.ratecanvas.SetLogx()
                self.ratecanvas.SetGridx()
                self.ratecanvas.SetGridy()
                self.graphs[channel].GetYaxis().SetRangeUser(0, 200)
                self.graphs[channel].GetYaxis().SetTitle("Mean Signal ({signal})".format(signal=self.collection[runnumbers[0]].signaldefinition.format(channel="")))
                self.graphs[channel].GetXaxis().SetTitle("Rate / kHz")
                self.graphs[channel].GetXaxis().SetLimits(1, 7000)
            self.graphs[channel].SetLineColor(color)
            self.graphs[channel].Draw(axisoption+"LP")
            self.ratelegend.AddEntry(self.graphs[channel], self.collection[runnumbers[0]].run.diamondname[channel]+" "+str(self.collection[runnumbers[0]].run.bias[channel])+"V ({startrun}-{endrun})".format(startrun=runnumbers[0], endrun=runnumbers[-1]), "lep")
            axisoption = ""

        self.ratelegend.Draw("SAME")
        #self.ratecanvas.Modified()
        self.ratecanvas.Update()
        tmpcanvas.Close()
        self.ShowAndWait = True
        self.IfWait("Signal VS Rate shown")


    def CreateSigmaMPVPlot(self):
        '''
        Analysis.FindMaxima() has to be executed for all contributing analysis
        :param saveplots:
        :param savename:
        :param ending:
        :return:
        '''
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
        self.IfWait("MPV vs Sigma shown...")

    def SignalHeightScan(self): # improve!
        #tmp = self.ShowAndWait
        #self.ShowAndWait = True
        SignalHeightScanCanvas = ROOT.TCanvas("SignalHeightScanCanvas", "SignalHeightScan Canvas")
        SignalHeightScanCanvas.cd()

        SignalHeightScanGraph = ROOT.TGraph()
        SignalHeightScanGraph.SetNameTitle("SignalHeightScanGraph", "Signal Height Scan")

        runnumbers = self.collection.keys()
        runnumbers.sort()

        count = 0
        for runnumber in runnumbers:
            if not self.collection[runnumber].TimingAlignmentFailed:
                SignalHeightScanGraph.SetPoint(count, runnumber, self.collection[runnumber].ExtremaResults['SignalHeight'])
                count += 1
            else:
                print "INFO: Run number {0} excluded in SignalHeightScan plot due to bad timing alignment !"
        SignalHeightScanGraph.SaveAs(self.SaveDirectory+"SignalHeightGraph.root")
        SignalHeightScanGraph.GetXaxis().SetTitle("Run Number")
        SignalHeightScanGraph.GetYaxis().SetTitle("Reconstructed Signal Height")
        SignalHeightScanGraph.Draw("AP*")
        self.SavePlots("SignalHeightGraph.png")
        self.IfWait("SignalHeightScan shown...")
        #self.ShowAndWait = tmp

    def PeakComparison(self, show = True):
        print "PeakComparision start"
        if show:
            PeakComparisonCanvasMax = ROOT.TCanvas("PeakComparisonCanvasMax", "PeakComparisonCanvas")
            PeakComparisonCanvasMin = ROOT.TCanvas("PeakComparisonCanvasMin", "PeakComparisonCanvas")

        runnumbers = self.collection.keys()
        pad_attributes = self.collection[runnumbers[0]].ExtremeAnalysis.Pad.Get2DAttributes()

        self.PeakPadMax = ATH2D("PeakPadMax", "Peak distribution over all selected runs", *pad_attributes)
        self.PeakPadMaxPad = BinCollection(*pad_attributes) # CHANGE NAME !
        self.PeakPadMin = ATH2D("PeakPadMin", "Low distribution over all selected runs", *pad_attributes)

        for runnumber in runnumbers:
            analysis = self.collection[runnumber]
            maxima = analysis.ExtremeAnalysis.ExtremaResults["FoundMaxima"]
            minima = analysis.ExtremeAnalysis.ExtremaResults["FoundMinima"]
            if maxima != None:
                for peak in maxima:
                    self.PeakPadMax.Fill(*peak)
                    if not hasattr(analysis.ExtremeAnalysis.Pad, "meansignaldistribution"):
                        analysis.ExtremeAnalysis.Pad.CalculateMeanSignalDistribution()
                    signal_ = analysis.ExtremeAnalysis.Pad.meansignaldistribution.GetBinContent(analysis.ExtremeAnalysis.Pad.GetBinNumber(*peak))
                    self.PeakPadMaxPad.Fill(peak[0], peak[1], signal_)
            else:
                print "WARNING: No Maxima results found in run ", runnumber, ". PeakComparisonMax will be incomplete."
            if minima != None:
                for peak in minima:
                    self.PeakPadMin.Fill(*peak)
            else:
                print "WARNING: No Minima results found in run ", runnumber, ". PeakComparisonMin will be incomplete."
        if show:
            ROOT.gStyle.SetPalette(55) # Rainbow palette
            PeakComparisonCanvasMax.cd()
            self.PeakPadMax.Draw("COLZ")
            self.SavePlots("PeakPadMax.png")
            PeakComparisonCanvasMin.cd()
            self.PeakPadMin.Draw("COLZ")
            self.SavePlots("PeakPadMin.png")

        # raw_input("peakpad")

    def PeakSignalEvolution(self, NMax = 3, NMin = 3, OnThisCanvas = None, BinRateEvolution = False):
        if OnThisCanvas != None:
            assert(isinstance(OnThisCanvas, ROOT.TCanvas)), "OnThisCanvas has to be a TCanvas object"
        print "Signal Evolution start"
        if not hasattr(self, "PeakPadMax"):
            self.PeakComparison(show = False)

        def BinsAreNearby(x1,y1,x2,y2, R):
            d2 = (x1-x2)**2 + (y1-y2)**2
            if d2 <= R**2:
                return True
            else:
                return False

        # Find the separated peaks / lows (binnumbers) to consider in Signal Evolution
        def FindPeakBins(PeakPad, N, maximum=False):
            self.PeakPadMaxPad.CalculateMeanSignalDistribution()
            peakbins = [-1]*N # container to store the binnumbers of the separated maximas found


            PeakPad2 = self.PeakPadMaxPad.meansignaldistribution # peaksearch due to mean signal content in bins
            i = 0
            while i<int(N):
                if i == 3 and maximum:
                    PeakPad = PeakPad2

                maxcount = PeakPad.GetMaximum() # counts of maximum

                if maxcount < 1:
                    break
                peakbins[i] = PeakPad.GetMaximumBin() # binnumber with hightest counts
                coordinates = PeakPad.GetBinCenter(peakbins[i])
                PeakPad.Fill(coordinates[0], coordinates[1], -maxcount) # remove content of maximum bin

                # if the binnumber is already in a neighborhood of a found peak, don't use it:
                IsInNBHD = False
                for k in xrange(i):
                    # IsInNBHD |= peakbins[k] in PeakPad.GetBinsInNbhd(peakbins[i], include_center=True, extended=True)
                    other_coordinates = PeakPad.GetBinCenter(peakbins[k])
                    IsInNBHD |= BinsAreNearby(coordinates[0], coordinates[1], other_coordinates[0], other_coordinates[1], 0.05)
                if IsInNBHD:
                    pass
                else:
                    i += 1
            self.VerbosePrint("Number of separated extrema to look at: {0:0.0f}".format(i))
            peakbins = peakbins[:i]
            return peakbins

        peakbins = FindPeakBins(self.PeakPadMax, NMax, maximum=True)
        lowbins = FindPeakBins(self.PeakPadMin, NMin)

        # Rate Time Evolution of Bin:
        if BinRateEvolution:
            tmpSaveDir = self.SaveDirectory

            high1 = peakbins[0]
            high2 = peakbins[1]
            high3 = peakbins[2]

            low1 = lowbins[0]
            low2 = lowbins[1]
            low3 = lowbins[2]

            for run_number in self.collection.keys():
                self.SetSaveDirectory(tmpSaveDir+str(run_number)+"/")
                self.collection[run_number].RateTimeEvolution(time_spacing = 5, save = True, binnumber = high1, nameExtension = "High1")
                self.collection[run_number].RateTimeEvolution(time_spacing = 5, save = True, binnumber = high2, nameExtension = "High2")
                self.collection[run_number].RateTimeEvolution(time_spacing = 5, save = True, binnumber = high3, nameExtension = "High3")

                self.collection[run_number].RateTimeEvolution(time_spacing = 5, save = True, binnumber = low1, nameExtension = "Low1")
                self.collection[run_number].RateTimeEvolution(time_spacing = 5, save = True, binnumber = low2, nameExtension = "Low2")
                self.collection[run_number].RateTimeEvolution(time_spacing = 5, save = True, binnumber = low3, nameExtension = "Low3")

            self.SetSaveDirectory(tmpSaveDir)

        # Fill all graphs of all separated peaks / lows
        def FillGraphDict(self, GraphDict, ListOfBins):
            runnumbers = self.collection.keys()
            runnumbers.sort()
            for peakbin in ListOfBins:
                GraphDict[peakbin] = ROOT.TGraphErrors()
                GraphDict[peakbin].SetNameTitle("MaxGraph_"+str(peakbin), "Evolution of Signal Response during Rate Scan")
                # signals = []
                i = 0
                for runnumber in runnumbers:
                    if not self.collection[runnumber].TimingAlignmentFailed:
                        self.collection[runnumber].ExtremeAnalysis.Pad.ListOfBins[peakbin].CreateBinSignalHisto(saveplot = True, savedir=self.SaveDirectory+str(runnumber)+"/",show_fit = False)
                        mean = self.collection[runnumber].ExtremeAnalysis.Pad.ListOfBins[peakbin].BinSignalHisto.GetMean()
                        error = self.collection[runnumber].ExtremeAnalysis.Pad.ListOfBins[peakbin].BinSignalHisto.GetRMS()/np.sqrt(self.collection[runnumber].ExtremeAnalysis.Pad.ListOfBins[peakbin].BinSignalHisto.GetEntries())
                        #mpv = self.collection[runnumber].ExtremeAnalysis.Pad.ListOfBins[peakbin].Fit['MPV']
                        # signals.append(mpv)
                        GraphDict[peakbin].SetPoint(i, runnumber, mean)
                        GraphDict[peakbin].SetPointError(i, 0, error)
                        i += 1
        MaxGraphs = {}
        FillGraphDict(self, MaxGraphs, peakbins)
        MinGraphs = {}
        FillGraphDict(self, MinGraphs, lowbins)

        # Prepare for drawing: Settings, create Canvas, create Legend
        if len(MaxGraphs)>0:
            marker = 20
            npeaks = len(MaxGraphs)
            PeakSignalEvolutionCanvas = ROOT.gROOT.GetListOfCanvases().FindObject("PeakSignalEvolutionCanvas")
            if not PeakSignalEvolutionCanvas:
                PeakSignalEvolutionCanvas = ROOT.TCanvas("PeakSignalEvolutionCanvas", "Signal Evolution Canvas")
            PeakSignalEvolutionCanvas.cd()
            legend = ROOT.TLegend(0.1,0.1,0.3,0.35)

            # determine the signal range for y axis:
            MaxSignals = []
            MinSignals = []
            for peakbin in peakbins:
                MaxSignals.append(MaxGraphs[peakbin].GetYaxis().GetXmax())
                MinSignals.append(MaxGraphs[peakbin].GetYaxis().GetXmin())
            MaxRange_peak = 1.1*np.array(MaxSignals).max()
            MinRange_peak = 0.9*np.array(MinSignals).min()

            for peaknr in xrange(npeaks):
                MaxGraphs[peakbins[peaknr]].SetMarkerStyle(marker)
                MaxGraphs[peakbins[peaknr]].SetMarkerColor(ROOT.kRed)
                MaxGraphs[peakbins[peaknr]].SetLineColor(ROOT.kRed)
                MaxGraphs[peakbins[peaknr]].Draw("SAME LP")
                legend.AddEntry(MaxGraphs[peakbins[peaknr]], "high"+str(peaknr+1), "lp")
                marker += 1
        else:
            PeakSignalEvolutionCanvas = ROOT.gROOT.GetListOfCanvases().FindObject("PeakSignalEvolutionCanvas")
            if not PeakSignalEvolutionCanvas:
                PeakSignalEvolutionCanvas = ROOT.TCanvas("PeakSignalEvolutionCanvas", "Signal Evolution Canvas")
            PeakSignalEvolutionCanvas.cd()
            legend = ROOT.TLegend(0.1,0.1,0.2,0.25)
            MaxRange_peak = None
            MinRange_peak = None

        if len(MinGraphs)>0:
            marker = 20
            nlows = len(MinGraphs)

            # determine the signal range for y axis:
            MaxSignals = []
            MinSignals = []
            for lowbin in lowbins:
                MaxSignals.append(MinGraphs[lowbin].GetYaxis().GetXmax())
                MinSignals.append(MinGraphs[lowbin].GetYaxis().GetXmin())
            MaxRange_low = 1.1*np.array(MaxSignals).max()
            MinRange_low = 0.9*np.array(MinSignals).min()

            for lownr in xrange(nlows):
                MinGraphs[lowbins[lownr]].SetMarkerStyle(marker)
                MinGraphs[lowbins[lownr]].SetMarkerColor(ROOT.kBlue)
                MinGraphs[lowbins[lownr]].SetLineColor(ROOT.kBlue)
                # MinGraphs[lowbins[lownr+1]].Draw("SAME LP")
                legend.AddEntry(MinGraphs[lowbins[lownr]], "low"+str(lownr+1), "lp")
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
            MinRange = 0.8*MinRange.min()
        else:
            MinRange = 0


        # Prepare for drawing: Individual Print options:
        NumbersOfGraphs = len(MaxGraphs) + len(MinGraphs)
        DrawOptions = ["SAME LP"]*NumbersOfGraphs
        try:
            DrawOptions[0] = "ALP"
        except IndexError: # if neither maxima nor minima found
            pass

        # Prepare for drawing rate:
        runnumbers = self.collection.keys()
        runnumbers.sort()
        first = runnumbers[0]
        last = runnumbers[-1]
        ratebins = last - first + 1
        RateHisto = ROOT.TH1D("RateHisto", "Rate Histogram", ratebins, first - 0.5, last + 0.5)
        for runnumber in runnumbers:
            rate_raw = self.collection[runnumber].RunInfo["rate_raw"]
            rate_kHz = 1.*rate_raw/100
            print "runnumber: ", self.collection[runnumber].run.run_number," == ",  runnumber, " rate_kHz: ", rate_kHz
            RateHisto.Fill(runnumber, rate_kHz)
        RateHisto.GetXaxis().SetTitle("Run Number")
        RateHisto.GetYaxis().SetTitle("Rate / kHz")


        # Draw everything:
        i = 0 # i-th draw option
        if len(MaxGraphs) > 0:
            MaxGraphs[peakbins[0]].GetYaxis().SetRangeUser(MinRange, MaxRange)
            MaxGraphs[peakbins[0]].GetXaxis().SetTitle("Run Number")
            MaxGraphs[peakbins[0]].GetYaxis().SetTitle("Mean Signal Response")
            MaxGraphs[peakbins[0]].Draw(DrawOptions[i])
            i += 1
            for peaknr in xrange(npeaks-1):
                MaxGraphs[peakbins[peaknr+1]].Draw(DrawOptions[i])
                i += 1
        if len(MinGraphs) > 0:
            if i == 0:
                MinGraphs[lowbins[0]].GetYaxis().SetRangeUser(MinRange, MaxRange)
                MinGraphs[lowbins[0]].GetXaxis().SetTitle("Run Number")
                MinGraphs[lowbins[0]].GetYaxis().SetTitle("Mean Signal Response")
            MinGraphs[lowbins[0]].Draw(DrawOptions[i])
            i += 1
            for lownr in xrange(nlows-1):
                MinGraphs[lowbins[lownr+1]].Draw(DrawOptions[i])
                i += 1
        legend.Draw()
        self.SavePlots("PeakSignalEvolution.png")
        self.SavePlots("PeakSignalEvolution.root")
        raw_input("waiting in AnalysisCollection->Line 375")
        pad = PeakSignalEvolutionCanvas.GetPad(0)
        RateHisto.SetStats(0)
        RateHisto.Draw("HIST Y+") # include in plot instead of second plot
        pad.SetLogy()
        self.SavePlots("PeakSignalEvolution_Rate.png")

        # show the selected bins in another canvas:
        if OnThisCanvas:
            ROOT.gStyle.SetPalette(53) # Dark Body Radiator palette
            OnThisCanvas.cd(1)

            if len(MaxGraphs) > 0:
                for peaknr in xrange(npeaks):
                    maxima = self.PeakPadMax.GetBinCenter(peakbins[peaknr])
                    text = ROOT.TText()
                    text.SetTextColor(ROOT.kRed)
                    text.DrawText(maxima[0]-0.02, maxima[1]-0.005, 'high'+str(peaknr+1))

            if len(MinGraphs) > 0:
                for lownr in xrange(nlows):
                    minima = self.PeakPadMin.GetBinCenter(lowbins[lownr])
                    text = ROOT.TText()
                    text.SetTextColor(ROOT.kBlue)
                    text.DrawText(minima[0]-0.01, minima[1]-0.005, 'low'+str(lownr+1))

            OnThisCanvas.Update()
            self.SavePlots("IIa-2_neutron_SignalDistribution_MAXSearch.png")
            raw_input("wait")

    def GetRunNumbers(self):
        '''
        :return: sorted list of run numbers in AnalysisCollection instance
        '''
        runnumbers = self.collection.keys()
        runnumbers.sort()
        return runnumbers

    def GetNumberOfAnalyses(self):
        '''
        :return: number of analyses that the analysis collection object contains
        '''
        return len(self.collection.keys())

    def ShowInfo(self):
        print "ANALYSIS COLLECTION INFO:"
        print "\tRuns: \tDiamond1 \tBias1 \tSelected \tDiamond2 \tBias2 \tSelected \tType"
        contentstring = ""
        for run in self.GetRunNumbers():
            contentstring += "\t{run} \t{Diamond1} \t{Bias1} \t{Selected1} \t\t{Diamond2} \t{Bias2} \t{Selected2} \t\t{Type}\n".format(
                run=str(run).zfill(3), Diamond1=self.collection[run].run.GetDiamondName(0).ljust(8), Bias1=str(self.collection[run].run.bias[0]).zfill(5), Selected1=str(self.collection[run].run.analyzeCh[0]).ljust(5),
                Diamond2=self.collection[run].run.GetDiamondName(3).ljust(8), Bias2=str(self.collection[run].run.bias[3]).zfill(5), Selected2=str(self.collection[run].run.analyzeCh[3]).ljust(5),
                Type=self.collection[run].run.RunInfo["type"])
        print contentstring
