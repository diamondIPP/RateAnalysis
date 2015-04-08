import ROOT
from RunClass import Run
import os
from BinCollection import BinCollection
from ConfigClass import *


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

        self.MeanSignalHisto = ROOT.TH1D()

        # loading data file
        assert (os.path.exists(self.TrackingPadAnalysisROOTFile)), 'cannot find '+self.TrackingPadAnalysisROOTFile
        self.rootfile = ROOT.TFile(self.TrackingPadAnalysisROOTFile)

        self.Checklist = { # True if Plot was created
            'DoAnalysis': False,
            'MeanSignalHisto': False,
            'HitsDistribution': False,
        }

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
            calibflag = self.track_info.calibflag
            if calibflag == 0:
                self.Pad.Fill(x_, y_, signal_)

        self.Pad.MakeFits()

        self.Signal2DDistribution = self.Pad.GetMeanSignalDistribution(minimum_bincontent)
        self.Signal2DDistribution.SetStats(False)
        self.Signal2DDistribution.GetXaxis().SetTitle('pos x / cm')
        self.Signal2DDistribution.GetYaxis().SetTitle('pos y / cm')
        self.Signal2DDistribution.GetYaxis().SetTitleOffset(1.4)

        self.Checklist['DoAnalysis'] = True

    def CreatePlots(self,saveplots = False,savename = '2DSignalDistribution',ending='png',saveDir = 'Results/', show3d = False):
        '''
        Creates 2D Signal Distribution plot
        :param saveplots: if True, save the plot
        :param savename: filename if saveplots = True
        :param ending: datatype of file if saveplots = True
        :param saveDir: directory to save the plot - has to end with '/'
        :return: -
        '''
        self.signal_canvas = ROOT.TCanvas()
        ROOT.SetOwnership(self.signal_canvas, False)
        self.signal_canvas.Clear()

        self.signal_canvas.SetName("signal_canvas")
        self.signal_canvas.SetTitle("signal distribution")
        # Plot the Signal2D TH2D histogram
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        if show3d:
            self.Signal2DDistribution.Draw("SPEC dm(2,10) pa(1,1,1) ci(1,1,1) a(15,45,0) s(1,1)")
            savename += '3D'
        else:
            self.Signal2DDistribution.Draw('colz')
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

        self.Checklist['MeanSignalHisto'] = True

    def CreateBoth(self,saveplots = False,savename = 'SignalDistribution',ending='png',saveDir = 'Results/',PS=False):
        self.combined_canvas = ROOT.TCanvas("combined_canvas","Combined Canvas",1000,500)
        self.combined_canvas.Divide(2,1)

        self.CreatePlots(False)
        self.CreateMeanSignalHistogram(False)

        self.combined_canvas.cd(1)
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

    def CreateHitsDistribution(self,saveplot = False):
        #ROOT.gStyle.SetNumberContours(10)
        canvas = ROOT.TCanvas('canvas', 'Hits',500,500)
        canvas.cd()
        self.Pad.counthisto.SetStats(False)
        self.Pad.counthisto.Draw('colz')#'surf2')
        #self.Pad.counthisto.Draw('CONT1 SAME')
        if saveplot:
            self.SavePlots('Hits_Distribution', 'png', 'Results/')
        raw_input('hits distribution')

    def SavePlots(self, savename, ending, saveDir):
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        resultsdir = saveDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)

        ROOT.gPad.Print(resultsdir+savename+'.'+ending)

    def FindMaxima(self,show=False):
        minimum_bincontent = 10
        self.MaximaAnalysis = Analysis(self.run_object,Config(200))
        self.MaximaAnalysis.DoAnalysis(minimum_bincontent)
        self.MaximaAnalysis.CreatePlots(saveplots=True, show3d=True)
        self.MaximaAnalysis.Pad.FindMaxima(minimum_bincount=minimum_bincontent,show=show)

    def GetMPVSigmas(self,show = False, minimum_counts = 10):
        '''
        Returns mpvs and signals of current run
        :return:
        '''
        MPVs = []
        MPVErrs = []
        Sigmas = []
        SigmaErrs = []
        for bin in self.Pad.ListOfBins:
            entries = bin.GetEntries()
            if entries >= minimum_counts:
                MPVs.append(bin.Fit['MPV'])
                MPVErrs.append(bin.Fit['MPVErr'])
                Sigmas.append(bin.Fit['Sigma'])
                SigmaErrs.append(bin.Fit['SigmaErr'])

        if show:
            canvas = ROOT.TCanvas('MPVSigmaCanvas', 'MPVSigmaCanvas')
            canvas.cd()
            graph = ROOT.TGraphErrors()
            count = 0
            for i in xrange(len(MPVs)):
                graph.SetPoint(count, MPVs[i], Sigmas[i])
                graph.SetPointError(count, MPVErrs[i], SigmaErrs[i])
                count += 1
            graph.Draw('AP')
            raw_input("MPV vs Sigma shown...")

        return MPVs, Sigmas, MPVErrs, SigmaErrs

    def ExportMC(self, MCDir = 'MCInputs/'):
        '''
        Hit distribution
        :return:
        '''
        if not os.path.exists(MCDir):
            os.makedirs(MCDir)
        self.Pad.counthisto.SaveAs(MCDir+str(self.run_object.run_number)+'counthisto.root')
