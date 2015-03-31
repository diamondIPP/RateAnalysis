import ROOT
from RunClass import Run
import os
import types as t
from BinCollection import BinCollection
from ConfigClass import *

#from Configuration.initialize_ROOT import initialize_ROOT
class blah(object):
    def __init__(self):
        pass


class Analysis(object):
    '''
    An Analysis Object contains all Data and Results of a single run.
    '''

    Signal2DDistribution = ROOT.TH2D()

    def __init__(self, run_object, config_object = Config()):
        '''
        Initializes the Analysis object.
        :param run_object: run object of type "Run"
        :param config_object: config object of type "Pad2DHistConfig"
        :return: -
        '''

        #initialize_ROOT()
        assert(run_object.run_number > 0), "No run selected, choose run.SetRun(run_nr) before you pass the run object"
        self.run_object = run_object
        self.config_object = config_object
        self.config_object.SetWindowFromDiamond(self.run_object.diamond)
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
        '''
        Create a bin collection object as self.Pad and load data from ROOT TTree
        into the Pad object. Then get the 2-dim signal distribution from self.Pad
        :param minimum_bincontent: Bins with less hits are ignored
        :return: -
        '''

        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer" # bins with less hits are ignored

        self.track_info = self.rootfile.Get('track_info') # Get TTree called "track_info"

        # create a bin collection object:
        self.Pad = BinCollection(*self.config_object.Get2DAttributes())

        # fill two 2-dim histograms to collect the hits and signal strength
        for i in xrange(self.track_info.GetEntries()):
            # read the ROOT TTree
            self.track_info.GetEntry(i)
            x_ = self.track_info.track_x
            y_ = self.track_info.track_y
            signal_ = abs(self.track_info.integral50)

            self.Pad.Fill(x_, y_, signal_)

        self.Signal2DDistribution = self.Pad.GetMeanSignalDistribution(minimum_bincontent)

    def CreatePlots(self,saveplots = False,savename = '2DSignalDistribution',ending='png',saveDir = 'Results/'):
        '''
        Creates 2D Signal Distribution plot
        :param saveplots: if True, save the plot
        :param savename: filename if saveplots = True
        :param ending: datatype of file if saveplots = True
        :param saveDir: directory to save the plot - has to end with '/'
        :return: -
        '''

        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        self.signal_canvas.Clear()

        self.signal_canvas.SetName("signal_canvas")
        self.signal_canvas.SetTitle("signal distribution")
        # Plot the Signal2D TH2D histogram
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        self.Signal2DDistribution.SetStats(False)

        self.Signal2DDistribution.Draw("SPEC dm(2,10) pa(1,1,1) ci(1,1,1) a(15,45,0) s(1,1)")
        raw_input('2d drawn')
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

    def CreateBoth(self,saveplots = False,savename = 'SignalDistribution',ending='png',saveDir = 'Results/',PS=False):
        self.combined_canvas = ROOT.TCanvas("combined_canvas","Combined Canvas",1000,500)
        self.combined_canvas.Divide(2,1)

        self.CreatePlots(False)
        self.CreateMeanSignalHistogram(False)

        self.combined_canvas.cd(1)
        self.Signal2DDistribution.GetXaxis().SetTitle('pos x / cm')
        self.Signal2DDistribution.GetYaxis().SetTitle('pos y / cm')
        self.Signal2DDistribution.GetYaxis().SetTitleOffset(1.4)
        #ROOT.gStyle.SetPalette(55)
        # ROOT.gStyle.SetNumberContours(999)
        self.Signal2DDistribution.Draw('colz')#"CONT1Z")#)'colz')
        self.combined_canvas.cd(2)

        if PS: #if photoshop mode, fill histogram pink
            ROOT.gStyle.SetHistFillColor(6)
            ROOT.gStyle.SetHistFillStyle(1001)
        else:
            ROOT.gStyle.SetHistFillColor(7)
            ROOT.gStyle.SetHistFillStyle(3003)

        self.MeanSignalHisto.UseCurrentStyle()
        self.MeanSignalHisto.GetXaxis().SetTitle('Signal response')
        self.MeanSignalHisto.Draw()
        self.combined_canvas.cd()

        savename = self.run_object.diamond.Specifications['Name']+'_'+self.run_object.diamond.Specifications['Irradiation']+'_'+savename+'_'+str(self.run_object.run_number) # diamond_irradiation_savename_runnr
        if saveplots:
            self.SavePlots(savename, ending, saveDir)
            self.SavePlots(savename, 'root', saveDir)
        if PS:
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
    '''
    An object of this class contains several analysis of runs.
    It gives the ability to compare the data from different runs.
    '''
    collection = {}
    current_run_number = -1

    def __init__(self):
        pass

    def AddAnalysis(self,analysis_obj):
        '''
        Adds an Analysis object to the analysis collection object
        :param analysis_obj: Analysis Object of type "Analysis"
        :return: -
        '''

        AnalysisCollectionClass.collection[analysis_obj.run_object.run_number] = analysis_obj
        AnalysisCollectionClass.current_run_number = analysis_obj.run_object.run_number

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

        self.canvas = ROOT.TCanvas("canvas", "FWHM")
        self.fwhm_histo = ROOT.TH1D("fwhm_histo", "FWHM Distribution of "+str(self.GetNumberOfAnalyses())+" runs",50,0,100)

        for run in AnalysisCollectionClass.collection:
            self.fwhm_histo.Fill(self.CalculateFWHM(False,run))
        self.canvas.cd()
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
        '''
        :return: number of analyses that the analysis collection object contains
        '''
        return len(AnalysisCollectionClass.collection.keys())
