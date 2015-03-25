import ROOT
from RunClass import Run
import os
import types as t
from BinCollection import BinCollection
from ConfigClass import *

#from Configuration.initialize_ROOT import initialize_ROOT

class Analysis(object):

    Signal2DDistribution = ROOT.TH2D()

    def __init__(self, run_object, config_object = Pad2DHistConfig(50)):
        #initialize_ROOT()
        assert(run_object.run_number > 0), "No run selected, choose run.SetRun(run_nr) before you pass the run object"
        self.run_object = run_object
        self.config_object = config_object
        self.TrackingPadAnalysisROOTFile = run_object.TrackingPadAnalysis['ROOTFile']
        self.Signal2DDistribution = ROOT.TH2D()
        self.Signal2DDistribution.SetDirectory(0) # is needed because of garbage collection

        self.signal_canvas = ROOT.TCanvas()
        ROOT.SetOwnership(self.signal_canvas, False)
        self.MeanSignalHisto = ROOT.TH1D()

        # loading data file
        assert (os.path.exists(self.TrackingPadAnalysisROOTFile)), 'cannot find '+self.TrackingPadAnalysisROOTFile
        self.rootfile = ROOT.TFile(self.TrackingPadAnalysisROOTFile)

        self.MeanSignalHistoIsCreated = False


    def DoAnalysis(self,minimum_bincontent = 1):
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer" # bins with less hits are ignored

        self.track_info = self.rootfile.Get('track_info') # Get TTree called "track_info"

        bins = self.config_object.bins_x
        xmin = self.run_object.diamond.Position['xmin']
        xmax = self.run_object.diamond.Position['xmax']
        ymin = self.run_object.diamond.Position['ymin']
        ymax = self.run_object.diamond.Position['ymax']

        self.Pad = BinCollection(bins,xmin,xmax,bins,ymin,ymax)

        # fill two 2-dim histograms to collect the hits and signal strength
        for i in xrange(self.track_info.GetEntries()):
            self.track_info.GetEntry(i)

            x_ = self.track_info.track_x
            y_ = self.track_info.track_y
            signal_ = abs(self.track_info.integral50)

            self.Pad.Fill(x_, y_, signal_)

        self.Signal2DDistribution = self.Pad.GetMeanSignalDistribution(minimum_bincontent)

    def CreatePlots(self,saveplots = False,savename = '2DSignalDistribution',ending='png',saveDir = 'Results/'):
        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        self.signal_canvas.Clear()

        self.signal_canvas.SetName("signal_canvas")
        self.signal_canvas.SetTitle("signal distribution")
        # Plot the Signal2D TH2D histogram
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        self.Signal2DDistribution.SetStats(False)

        self.Signal2DDistribution.Draw('colz')

        if saveplots:
            self.SavePlots(savename, ending, saveDir)

    def CreateMeanSignalHistogram(self,saveplots = False,savename = 'MeanSignalDistribution',ending='png',saveDir = 'Results/'):
        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        self.signal_canvas.Clear()

        minimum = self.Signal2DDistribution.GetMinimum()
        maximum = self.Signal2DDistribution.GetMaximum()
        self.MeanSignalHisto = ROOT.TH1D("MeanSignalHisto","Mean Signal Histogram",100,minimum,maximum)

        self.signal_canvas.SetName("signal_canvas")
        self.signal_canvas.SetTitle("Mean Signal Distribution")

        nbins = (self.Signal2DDistribution.GetNbinsX()+2)*(self.Signal2DDistribution.GetNbinsY()+2)
        for i in xrange(nbins):
            bincontent = self.Signal2DDistribution.GetBinContent(i)
            if bincontent != 0.:
                self.MeanSignalHisto.Fill(bincontent)

        self.MeanSignalHisto.Draw()

        if saveplots:
            self.SavePlots(savename, ending, saveDir)

        self.MeanSignalHistoIsCreated = True

    def CreateBoth(self,saveplots = False,savename = 'SignalDistribution',ending='png',saveDir = 'Results/'):

        self.combined_canvas = ROOT.TCanvas("combined_canvas","Combined Canvas",1000,500)
        self.combined_canvas.Divide(2,1)

        self.CreatePlots(False)
        self.CreateMeanSignalHistogram(False)

        self.combined_canvas.cd(1)
        self.Signal2DDistribution.Draw('colz')
        self.combined_canvas.cd(2)
        # ROOT.gStyle.SetPalette(53)
        # ROOT.gStyle.SetNumberContours(999)
        ROOT.gStyle.SetHistFillColor(6)#7)
        ROOT.gStyle.SetHistFillStyle(1001)#3003)
        self.MeanSignalHisto.UseCurrentStyle()
        self.MeanSignalHisto.Draw()
        self.combined_canvas.cd()

        savename = self.run_object.diamond.Specifications['Name']+'_'+self.run_object.diamond.Specifications['Irradiation']+'_'+savename+'_'+str(self.run_object.run_number) # diamond_irradiation_savename_runnr
        if saveplots:
            self.SavePlots(savename, ending, saveDir)
            self.SavePlots(savename, 'root', saveDir)
        ROOT.gStyle.SetHistFillColor(7)
        ROOT.gStyle.SetHistFillStyle(3003)

    def CreateHitsDistribution(self):
        #ROOT.gStyle.SetNumberContours(10)
        canvas = ROOT.TCanvas('canvas', 'Hits',500,500)
        canvas.cd()
        self.Pad.counthisto.SetStats(False)
        self.Pad.counthisto.Draw('colz')#'surf2')
        #self.Pad.counthisto.Draw('CONT1 SAME')
        raw_input('hits distribution')

    def SavePlots(self, savename, ending, saveDir):
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        resultsdir = saveDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)

        ROOT.gPad.Print(resultsdir+savename+'.'+ending)

class AnalysisCollectionClass(object):

    collection = {}
    current_run_number = -1

    def __init__(self):
        pass

    def AddAnalysis(self,analysis_obj):

        AnalysisCollectionClass.collection[analysis_obj.run_object.run_number] = analysis_obj
        AnalysisCollectionClass.current_run_number = analysis_obj.run_object.run_number


    def CreateFWHMPlot(self, saveplots = True, savename = 'FWHM_Histo', ending = 'png'):

        self.canvas = ROOT.TCanvas("canvas", "FWHM")
        self.fwhm_histo = ROOT.TH1D("fwhm_histo", "FWHM Distribution of "+str(self.GetNumberOfAnalyses())+" runs",50,0,100)

        for run in AnalysisCollectionClass.collection:
            self.fwhm_histo.Fill(self.CalculateFWHM(False,run))
        self.canvas.cd()
        self.fwhm_histo.Draw()

        if saveplots:
            # Results directories:
            resultsdir = 'Results/' # eg. 'Results/run_364/'
            if not os.path.exists(resultsdir):
                os.makedirs(resultsdir)

            ROOT.gPad.Print(resultsdir+savename+'.'+ending)
            ROOT.gPad.Print(resultsdir+savename+'.'+'root')

        #raw_input("wait")

    def CalculateFWHM(self, print_result = True, run_number = None):
        if run_number is None:
            run_number = AnalysisCollectionClass.current_run_number
        assert(type(run_number) is t.IntType and 0 < run_number < 1000), "Invalid run number"

        analysis_obj = AnalysisCollectionClass.collection[run_number]

        assert(analysis_obj.MeanSignalHistoIsCreated), "Histogram not created yet or not found"

        maximum = analysis_obj.MeanSignalHisto.GetMaximum()
        low_bin = analysis_obj.MeanSignalHisto.FindFirstBinAbove(maximum/2.)
        high_bin = analysis_obj.MeanSignalHisto.FindLastBinAbove(maximum/2.)

        fwhm = analysis_obj.MeanSignalHisto.GetBinCenter(high_bin) - analysis_obj.MeanSignalHisto.GetBinCenter(low_bin)

        if print_result:
            print "FWHM of run ",run_number," is: ",fwhm

        return fwhm

    def GetNumberOfAnalyses(self):
        return len(AnalysisCollectionClass.collection.keys())
