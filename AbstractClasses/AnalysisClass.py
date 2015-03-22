import ROOT
from RunClass import Run
import os
import types as t
from ConfigClass import *

#from Configuration.initialize_ROOT import initialize_ROOT

class Analysis(object):

    Signal2DDistribution = {
        'Histogram':''

    }

    def __init__(self, run_object, config_object = Pad2DHistConfig(50,-1,1)):
        #initialize_ROOT()
        self.run_object = run_object
        self.config_object = config_object
        self.TrackingPadAnalysisROOTFile = run_object.TrackingPadAnalysis['ROOTFile']
        self.Signal2DDistribution['Histogram'] = ROOT.TH2D("Signal2D",
                                                        "2D Signal distribution",
                                                        config_object.bins_x,
                                                        config_object.min_x,
                                                        config_object.max_x,
                                                        config_object.bins_y,
                                                        config_object.min_y,
                                                        config_object.max_y
                                                        )
        self.Signal2DDistribution['Histogram'].SetDirectory(0) # is needed because of garbage collection

        self.signal_canvas = ROOT.TCanvas()
        ROOT.SetOwnership(self.signal_canvas, False)
        self.MeanSignalHisto = ROOT.TH1D()
        self.signal_sum = ROOT.TH2D("signal_sum",
                                    "Sum of signal distribution",
                                    config_object.bins_x,
                                    config_object.min_x,
                                    config_object.max_x,
                                    config_object.bins_y,
                                    config_object.min_y,
                                    config_object.max_y
                                    )
        self.signal_sum.SetDirectory(0)

        self.signal_counts = self.signal_sum.Clone("Signal Counts")
        self.signal_counts.SetDirectory(0)

        # loading data file
        assert (os.path.exists(self.TrackingPadAnalysisROOTFile)), 'cannot find '+self.TrackingPadAnalysisROOTFile
        self.rootfile = ROOT.TFile(self.TrackingPadAnalysisROOTFile)

        self.MeanSignalHistoIsCreated = False


    def DoAnalysis(self,minimum_bincontent = 1):
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer"
        minimum_statistics = minimum_bincontent # bins with less hits are ignored

        self.track_info = self.rootfile.Get('track_info') # Get TTree called "track_info"

        bins = self.config_object.bins_x
        xmin = self.config_object.min_x
        xmax = self.config_object.max_x
        ymin = self.config_object.min_y
        ymax = self.config_object.max_y

        # fill two 2-dim histograms to collect the hits and signal strength
        for i in xrange(self.track_info.GetEntries()):

            self.track_info.GetEntry(i)
            self.signal_sum.Fill(self.track_info.track_x, self.track_info.track_y, self.track_info.integral50)
            self.signal_counts.Fill(self.track_info.track_x, self.track_info.track_y, 1)

        # go through every bin, calculate the average signal strength and fill the main 2D hist
        binwidth_x = 1.*(xmax-xmin)/bins
        binwidth_y = 1.*(ymax-ymin)/bins
        current_pos_x = xmin + 1.*binwidth_x/2.
        for bin_x in xrange(1,bins+1):

            current_pos_y = ymin + 1.*binwidth_y/2.

            for bin_y in xrange(1,bins+1):

                binsignalsum = abs(self.signal_sum.GetBinContent(bin_x, bin_y))
                binsignalcount = self.signal_counts.GetBinContent(bin_x, bin_y)

                if binsignalcount >= minimum_statistics :
                    self.Signal2DDistribution['Histogram'].Fill(current_pos_x, current_pos_y, abs(binsignalsum/binsignalcount))

                current_pos_y += binwidth_y

            current_pos_x += binwidth_x

    def CreatePlots(self,saveplots = False,savename = '2DSignalDistribution',ending='png',saveDir = 'Results/'):
        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        self.signal_canvas.Clear()

        self.signal_canvas.SetName("signal_canvas")
        self.signal_canvas.SetTitle("signal distribution")
        # Plot the Signal2D TH2D histogram
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        self.Signal2DDistribution['Histogram'].SetStats(False)

        self.Signal2DDistribution['Histogram'].Draw('colz')

        if saveplots:
            self.SavePlots(savename, ending, saveDir)

    def CreateMeanSignalHistogram(self,saveplots = False,savename = 'MeanSignalDistribution',ending='png',saveDir = 'Results/'):
        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        self.signal_canvas.Clear()

        minimum = self.Signal2DDistribution['Histogram'].GetMinimum()
        maximum = self.Signal2DDistribution['Histogram'].GetMaximum()
        self.MeanSignalHisto = ROOT.TH1D("MeanSignalHisto","Mean Signal Histogram",100,minimum,maximum)

        self.signal_canvas.SetName("signal_canvas")
        self.signal_canvas.SetTitle("Mean Signal Distribution")

        nbins = (self.Signal2DDistribution['Histogram'].GetNbinsX()+2)*(self.Signal2DDistribution['Histogram'].GetNbinsY()+2)
        for i in xrange(nbins):
            bincontent = self.Signal2DDistribution["Histogram"].GetBinContent(i)
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
        self.Signal2DDistribution['Histogram'].Draw('colz')
        self.combined_canvas.cd(2)
        # ROOT.gStyle.SetPalette(53)
        # ROOT.gStyle.SetNumberContours(999)
        ROOT.gStyle.SetHistFillColor(7)
        ROOT.gStyle.SetHistFillStyle(3003)
        self.MeanSignalHisto.UseCurrentStyle()
        self.MeanSignalHisto.Draw()

        savename = self.run_object.diamond.Specifications['Name']+'_'+self.run_object.diamond.Specifications['Irradiation']+'_'+savename+'_'+str(self.run_object.run_number) # diamond_irradiation_savename_runnr
        if saveplots:
            self.SavePlots(savename, ending, saveDir)
            self.SavePlots(savename, 'root', saveDir)

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
