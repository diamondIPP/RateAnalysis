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

    def __init__(self, run_object, config_object = Config(), verbose = False):
        '''
        Initializes the Analysis object.
        :param run_object: run object of type "Run"
        :param config_object: config object of type "Pad2DHistConfig"
        :return: -
        '''

        #initialize_ROOT()
        assert(run_object.run_number > 0), "No run selected, choose run.SetRun(run_nr) before you pass the run object"
        self.verbose = verbose
        self.check_offset = True
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

        self.MCResults ={
            'TrueNPeaks': None,
            'FoundNMaxima': None,
            'Ninjas': None,
            'Ghosts': None
        }

    def VerbosePrint(self, *args):
        if self.verbose:
            # Print each argument separately so caller doesn't need to
            # stuff everything to be printed into a single string
            for arg in args:
               print arg,
            print
        else:
            pass

    def DoAnalysis(self,minimum_bincontent = 1):
        '''
        Create a bin collection object as self.Pad and load data from ROOT TTree
        into the Pad object. Then get the 2-dim signal distribution from self.Pad
        :param minimum_bincontent: Bins with less hits are ignored
        :return: -
        '''
        print 'Check 1 in DoAnalysis'
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer" # bins with less hits are ignored

        if not self.run_object.IsMonteCarlo:
            self.track_info = self.rootfile.Get('track_info') # Get TTree called "track_info"
        print 'Check 2 in DoAnalysis'
        # create a bin collection object:
        # HERE IT CRASHES:
        self.Pad = BinCollection(self, *self.config_object.Get2DAttributes())

        # Check for misalignment in track vs signal using calibflag:
        if self.check_offset and not self.run_object.IsMonteCarlo:
            offset_canvas = ROOT.TCanvas("offset_canvas", "offset_canvas")
            th1 = ROOT.TH1D("th1", "calibration offset", 21, -10.5, 10.5)
            self.track_info.Draw("calib_offset>>th1", "calibflag==1 && calib_offset < 50")
            pad = offset_canvas.GetPad(0)
            pad.SetLogy()
            offset_canvas.Update()
            offset = int(th1.GetBinCenter(th1.GetMaximumBin()))
            if abs(th1.GetMean())>0.3:
                print "\nINFO: BAD TIMING ALIGNMENT !\n"
            print "MOST COMMON OFFSET: {0:0.0F}\n".format(offset)
        else:
            offset = 0
        print 'Check 3 in DoAnalysis'
        # fill two 2-dim histograms to collect the hits and signal strength
        # if we have an offset, pick values from different events in tree
        if not self.run_object.IsMonteCarlo:
            if offset == 0:
                for i in xrange(self.track_info.GetEntries()):
                    # read the ROOT TTree
                    self.track_info.GetEntry(i)
                    x_ = self.track_info.track_x
                    y_ = self.track_info.track_y
                    signal_ = abs(self.track_info.integral50)
                    calibflag = self.track_info.calibflag
                    if calibflag == 0:
                        self.Pad.Fill(x_, y_, signal_)
            else:
                if offset > 0:
                    pluscorrection = offset
                    minuscorrection = 0
                else:
                    pluscorrection = 0
                    minuscorrection = abs(offset)

                for i in xrange(self.track_info.GetEntries()-abs(offset)):
                    # read the ROOT TTree
                    self.track_info.GetEntry(i+pluscorrection)
                    x_ = self.track_info.track_x
                    y_ = self.track_info.track_y
                    self.track_info.GetEntry(i+minuscorrection)
                    signal_ = abs(self.track_info.integral50)
                    calibflag = self.track_info.calibflag
                    if calibflag == 0:
                        self.Pad.Fill(x_, y_, signal_)
        else: # run is Monte Carlo:
            GoOn = True
            d = 0
            while GoOn and d<5: # loop for different MC Signal Distributions
                try:
                    print 'Check 4 in DoAnalysis'
                    self.VerbosePrint('try a Signal Distribution')
                    if not self.run_object.DataIsMade:
                        print 'Check 5 in DoAnalysis'
                        self.run_object.Simulate(save=True, draw=True) # if draw=False the first distribution will be taken
                    for i in xrange(self.run_object.NumberOfHits):
                        print 'Check 6 in DoAnalysis'
                        x_ = self.run_object.Data['track_x'][i]
                        y_ = self.run_object.Data['track_y'][i]
                        signal_ = self.run_object.Data['integral50'][i]
                        self.Pad.Fill(x_, y_, signal_)

                except IndexError:
                    pass
                else:
                    GoOn = False
                d += 1


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

    def CreateHitsDistribution(self,saveplot = False, drawoption = 'colz'): # add palette!
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        canvas = ROOT.TCanvas('canvas', 'Hits',500,500)
        canvas.cd()
        self.Pad.counthisto.SetStats(False)
        self.Pad.counthisto.Draw(drawoption)#'surf2')
        #self.Pad.counthisto.Draw('CONT1 SAME')
        if saveplot:
            self.SavePlots('Hits_Distribution', 'png', 'Results/')
        raw_input('hits distribution')
        self.Checklist['HitsDistribution'] = True

    def SavePlots(self, savename, ending, saveDir):
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        resultsdir = saveDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)

        ROOT.gPad.Print(resultsdir+savename+'.'+ending)

    def FindMaxima(self,show=False):
        minimum_bincontent = 10
        print 'find maxima started'
        try:
            print self.MaximaAnalysis
        except AttributeError:
            print 'self.MaximaAnalysis not yet defined'
        else:
            del self.MaximaAnalysis
            print '\nself.MaximaAnalysis deleted\n'
        self.MaximaAnalysis = Analysis(self.run_object,Config(200))
        print 'self.MaximaAnalysis defined'
        self.MaximaAnalysis.DoAnalysis(minimum_bincontent) # HERE IT CRASHES
        print 'DoAnalysis executed'
        if show:
            self.MaximaAnalysis.CreatePlots(saveplots=True, show3d=False)
        self.MaximaAnalysis.Pad.FindMaxima(minimum_bincount=minimum_bincontent,show=show)
        print 'end of findmaxima'

    def WaldWolfowitzRunsTest(self):
        pass

    def Chi2Test(self):
        pass

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
            raw_input("MPV vs Sigma shown. Press Enter to go on.")

        return MPVs, Sigmas, MPVErrs, SigmaErrs

    def ExportMC(self, MCDir = 'MCInputs/'):
        '''
        Hit distribution
        :return:
        '''
        if not self.run_object.IsMonteCarlo:
            if not os.path.exists(MCDir):
                os.makedirs(MCDir)
            self.Pad.counthisto.SaveAs(MCDir+str(self.run_object.run_number)+'counthisto.root')
            self.VerbosePrint("CountHisto exported..")
        else:
            print "INFO: Monte Carlo run can not be exported as MC input"
