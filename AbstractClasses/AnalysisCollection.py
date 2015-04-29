import ROOT
from AbstractClasses.ATH2D import ATH2D
import types as t
import os
from Elementary import Elementary
from ROOT import TGraphErrors

class AnalysisCollection(Elementary):
    '''
    An object of this class contains several analysis of runs.
    It gives the ability to compare the data from different runs.
    '''
    collection = {}
    current_run_number = -1

    def __init__(self, verbose = False):
        Elementary.__init__(self, verbose=verbose)

    def __del__(self):
        for runnumber in self.collection.keys():
            self.collection[runnumber].__del__()
        if hasattr(self, "FWHMcanvas"):
            ROOT.gROOT.Delete("FWHMcanvas")
        if hasattr(self, "fwhm_histo"):
            ROOT.gROOT.Delete("fwhm_histo")

    def AddAnalysis(self,analysis_obj):
        '''
        Adds an Analysis object to the analysis collection object
        :param analysis_obj: Analysis Object of type "Analysis"
        :return: -
        '''

        AnalysisCollection.collection[analysis_obj.run_object.run_number] = analysis_obj
        AnalysisCollection.current_run_number = analysis_obj.run_object.run_number

    #
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
        if run_number is None:
            run_number = AnalysisCollection.current_run_number
        assert(type(run_number) is t.IntType and 0 < run_number < 1000), "Invalid run number"

        analysis_obj = AnalysisCollection.collection[run_number]

        assert(analysis_obj.MeanSignalHistoIsCreated), "Histogram not created yet or not found"

        maximum = analysis_obj.MeanSignalHisto.GetMaximum()
        low_bin = analysis_obj.MeanSignalHisto.FindFirstBinAbove(maximum/2.)
        high_bin = analysis_obj.MeanSignalHisto.FindLastBinAbove(maximum/2.)

        fwhm = analysis_obj.MeanSignalHisto.GetBinCenter(high_bin) - analysis_obj.MeanSignalHisto.GetBinCenter(low_bin)

        if print_result:
            print "FWHM of run ",run_number," is: ",fwhm

        return fwhm

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
        tmp = self.ShowAndWait
        self.ShowAndWait = True
        SignalHeightScanCanvas = ROOT.TCanvas("SignalHeightScanCanvas", "SignalHeightScan Canvas")
        SignalHeightScanCanvas.cd()

        SignalHeightScanGraph = ROOT.TGraph()
        SignalHeightScanGraph.SetNameTitle("SignalHeightScanGraph", "Signal Height Scan")

        runnumbers = self.collection.keys()
        runnumbers.sort()

        count = 0
        for runnumber in runnumbers:

            SignalHeightScanGraph.SetPoint(count, runnumber, self.collection[runnumber].ExtremaResults['SignalHeight'])
            count += 1
        SignalHeightScanGraph.SaveAs("SignalHeightGraph.root")
        SignalHeightScanGraph.GetXaxis().SetTitle("Run Number")
        SignalHeightScanGraph.GetYaxis().SetTitle("Reconstructed Signal Height")
        SignalHeightScanGraph.Draw("AP*")
        self.IfWait("SignalHeightScan shown...")
        self.ShowAndWait = tmp

    def PeakComparison(self, show = True):
        print "PeakComparision start"
        if show:
            PeakComparisonCanvasMax = ROOT.TCanvas("PeakComparisonCanvasMax", "PeakComparisonCanvas")
            PeakComparisonCanvasMin = ROOT.TCanvas("PeakComparisonCanvasMin", "PeakComparisonCanvas")

        runnumbers = self.collection.keys()
        pad_attributes = self.collection[runnumbers[0]].ExtremeAnalysis.Pad.Get2DAttributes()

        self.PeakPadMax = ATH2D("PeakPadMax", "Peak distribution over all selected runs", *pad_attributes)
        self.PeakPadMin = ATH2D("PeakPadMin", "Low distribution over all selected runs", *pad_attributes)

        for runnumber in runnumbers:
            analysis = self.collection[runnumber]
            maxima = analysis.ExtremeAnalysis.ExtremaResults["FoundMaxima"]
            minima = analysis.ExtremeAnalysis.ExtremaResults["FoundMinima"]
            if maxima is not None:
                for peak in maxima:
                    self.PeakPadMax.Fill(*peak)
            else:
                print "WARNING: No Maxima results found in run ", runnumber, ". PeakComparisonMax will be incomplete."
            if minima is not None:
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

        raw_input("peakpad")

    def PeakSignalEvolution(self, NMax = 3, NMin = 3, OnThisCanvas = None):
        if OnThisCanvas is not None:
            assert(isinstance(OnThisCanvas, ROOT.TCanvas)), "OnThisCanvas has to be a TCanvas object"
        print "Signal Evolution start"
        if not hasattr(self, "PeakPadMax"):
            self.PeakComparison(show = False)

        # Find the separated peaks (binnumbers) to consider in Signal Evolution
        peakbins = [-1]*NMax # container to store the binnumbers of the separated maximas found
        i = 0
        while i<int(NMax):
            maxcount = self.PeakPadMax.GetMaximum() # counts of maximum
            if maxcount < 1:
                break
            peakbins[i] = self.PeakPadMax.GetMaximumBin() # binnumber with hightest counts
            coordinates = self.PeakPadMax.GetBinCenter(peakbins[i])
            self.PeakPadMax.Fill(coordinates[0], coordinates[1], -maxcount) # remove content of maximum bin
            # if the binnumber is already in a neighborhood of a found peak, don't use it:
            IsInNBHD = False
            for k in xrange(i):
                IsInNBHD |= peakbins[k] in self.PeakPadMax.GetBinsInNbhd(peakbins[i], include_center=True)
            if IsInNBHD:
                pass
            else:
                i += 1
        self.VerbosePrint("Number of separated peaks to look at: {0:0.0f}".format(i))
        peakbins = peakbins[:i]

        # Fill all graphs of all separated peaks
        MaxGraphs = {}
        for peakbin in peakbins:
            MaxGraphs[peakbin] = ROOT.TGraph()
            MaxGraphs[peakbin].SetNameTitle("MaxGraph_"+str(peakbin), "Evolution of Signal Response during Rate Scan")
            signals = []
            runnumbers = self.collection.keys()
            runnumbers.sort()
            i = 0
            for runnumber in runnumbers:
                self.collection[runnumber].ExtremeAnalysis.Pad.ListOfBins[peakbin].CreateBinSignalHisto(saveplot = True, savedir=self.SaveDirectory+str(runnumber)+"/",show_fit = True)
                mpv = self.collection[runnumber].ExtremeAnalysis.Pad.ListOfBins[peakbin].Fit['MPV']
                signals.append(mpv)
                MaxGraphs[peakbin].SetPoint(i, runnumber, mpv)
                i += 1

        # Print all Graphs of all peaks into the same canvas
        if len(MaxGraphs)>0:
            marker = 20
            npeaks = len(MaxGraphs)
            PeakSignalEvolutionCanvas = ROOT.TCanvas("PeakSignalEvolutionCanvas", "Signal Evolution Canvas")
            PeakSignalEvolutionCanvas.cd()
            legend = ROOT.TLegend(0.1,0.1,0.4,0.3)
            MaxGraphs[peakbins[0]].SetMarkerStyle(marker)
            MaxGraphs[peakbins[0]].SetMarkerColor(ROOT.kRed)
            MaxGraphs[peakbins[0]].SetLineColor(ROOT.kRed)
            MaxGraphs[peakbins[0]].Draw("ALP")
            i = 1
            legend.AddEntry(MaxGraphs[peakbins[0]], "high"+str(i), "lp")
            marker += 1
            i+=1
            for peaknr in xrange(npeaks-1):
                MaxGraphs[peakbins[peaknr+1]].SetMarkerStyle(marker)
                MaxGraphs[peakbins[peaknr+1]].SetMarkerColor(ROOT.kRed)
                MaxGraphs[peakbins[peaknr+1]].SetLineColor(ROOT.kRed)
                MaxGraphs[peakbins[peaknr+1]].Draw("SAME LP")
                legend.AddEntry(MaxGraphs[peakbins[peaknr+1]], "high"+str(i), "lp")
                marker += 1
                i += 1

            legend.Draw()
            self.SavePlots("PeakSignalEvolution.png")
            self.SavePlots("PeakSignalEvolution.root")

        # show the selected bins in another canvas:
        OnThisCanvas.cd(1)
        i = 0
        for peakbin in peakbins:
            maxima = self.PeakPadMax.GetBinCenter(peakbin)
            text = ROOT.TText()
            text.SetTextColor(ROOT.kRed)
            text.DrawText(maxima[i][0]-0.02, maxima[i][1]-0.005, 'high'+str(i))
            i += 1

        else:
            print "No Maxima found over all considered Runs"


    def GetNumberOfAnalyses(self):
        '''
        :return: number of analyses that the analysis collection object contains
        '''
        return len(AnalysisCollection.collection.keys())
