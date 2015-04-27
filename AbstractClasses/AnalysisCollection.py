import ROOT
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

    def PeakComparison(self):
        PeakComparisonCanvas = ROOT.TCanvas("PeakComparisonCanvas", "PeakComparisonCanvas")
        PeakComparisonCanvas.cd()

        runnumbers = self.collection.keys()
        pad_attributes = self.collection[runnumbers[0]].ExtremeAnalysis.Pad.Get2DAttributes()

        PeakPad = ROOT.TH2D("PeakPad", "Peak distribution over all selected runs", *pad_attributes)

        for runnumber in runnumbers:
            analysis = self.collection[runnumber]
            maxima = analysis.ExtremeAnalysis.ExtremaResults["FoundMaxima"]
            if maxima is not None:
                for peak in maxima:
                    PeakPad.Fill(*peak)
            else:
                print "WARNING: No Maxima results found in run ", runnumber, ". PeakComparison will be incomplete."

        PeakPad.Draw("COLZ")
        self.SavePlots("PeakPad.png")

        raw_input("peakpad")

    def GetNumberOfAnalyses(self):
        '''
        :return: number of analyses that the analysis collection object contains
        '''
        return len(AnalysisCollection.collection.keys())
