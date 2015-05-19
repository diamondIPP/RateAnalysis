import ROOT
from ROOT import gROOT
from RunClass import Run
from AbstractClasses.Elementary import Elementary
from AbstractClasses.Langau import Langau
import os
from copy import deepcopy
import numpy as np
from BinCollection import BinCollection
from ConfigClass import *
from array import array
import types as t
import ConfigParser
from ConfigParser import NoSectionError


class Analysis(Elementary):
    '''
    An Analysis Object contains all Data and Results of a SINGLE run.
    '''


    def __init__(self, run_object, config_object = Config(), verbose = False):
        '''
        Initializes the Analysis object.
        :param run_object: run object of type "Run"
        :param config_object: config object of type "Pad2DHistConfig"
        :return: -
        '''
        Elementary.__init__(self, verbose=verbose)
        #initialize_ROOT()
        assert(run_object.run_number > 0), "No run selected, choose run.SetRun(run_nr) before you pass the run object"
        self.run_object = run_object
        self.RunInfo = deepcopy(run_object.RunInfo)
        self.config_object = config_object
        self.config_object.SetWindowFromDiamond(self.run_object.diamond)
        self.TrackingPadAnalysisROOTFile = run_object.TrackingPadAnalysis["ROOTFile"]
        # self.Signal2DDistribution = ROOT.TH2D()
        # self.Signal2DDistribution.SetDirectory(0) # is needed because of garbage collection


        # loading data file
        assert (os.path.exists(self.TrackingPadAnalysisROOTFile)), "cannot find "+self.TrackingPadAnalysisROOTFile
        self.rootfile = ROOT.TFile(self.TrackingPadAnalysisROOTFile)

        self.Checklist = { # True if Plot was created # -> USED?!
            "DoAnalysis": False,
            "MeanSignalHisto": False,
            "HitsDistribution": False,
        }

        self.ExtremaResults ={
            "TrueNPeaks": None, # MC true number of peaks
            "FoundNMaxima": None, # number of Maxima found
            "FoundMaxima": None, # Maxima found as list [(x1, y1), (x2, y2), ...]
            "FoundNMinima": None, # number of Minima found
            "FoundMinima": None, # Minima found as list [(x1, y1), (x2, y2), ...]
            "Ninjas": None, # number of true peaks not found (only available when MC Run)
            "Ghosts": None, # nunmber of wrong Maxima found (only available when MC Run)
            "SignalHeight": 0.
        }

        self.SignalHistoFitResults = {
            "FitFunction": None,
            "Peak": None,
            "FWHM": None,
            "Chi2": None,
            "NDF": None,

        }
        self.TimingAlignmentFailed = False

    def __del__(self):
        self.VerbosePrint("Deleting Analysis..")
        if hasattr(self, "Pad"):
            self.Pad.__del__()
            del self.Pad
            print "Pad of Analysis deleted"
            # delattr(self, "Pad")
        # gROOT.GetListOfFiles().ls()
        rootfile = gROOT.GetListOfFiles().FindObject(self.TrackingPadAnalysisROOTFile)
        if rootfile:
            rootfile.Close()
        if hasattr(self, "ExtremeAnalysis"):
            self.ExtremeAnalysis.__del__()
            print "ExtremeAnalysis deleted"
        if hasattr(self, "combined_canvas"):
            canvas = ROOT.gROOT.GetListOfCanvases().FindObject(self.combined_canvas_name)
            canvas.Close()
            del canvas
        if hasattr(self, "MeanSignalHisto"):
            ROOT.gROOT.Delete(self.MeanSignalHisto_name)
        self.VerbosePrint("Analysis deleted")




    def LoadConfig(self): # BUG: CRASH WHEN Config loaded after many runs (CWD lost)
        # GoOn = True
        # configfile = "Configuration/AnalysisConfig.cfg"
        # while GoOn:
        #     try:
        #         parser = ConfigParser.ConfigParser()
        #         output = parser.read(configfile)
        #         print " ---- Config Parser Read ---- \n - ", output, " -\n"
        #         ShowAndWait = parser.getboolean("DISPLAY","ShowAndWait")
        #         SaveMCData = parser.getboolean("SAVE","SaveMCData")
        #         check_offset = parser.getboolean("DO-ANALYSIS","check_offset")
        #         self.ShowAndWait = ShowAndWait
        #         self.SaveMCData = SaveMCData
        #         self.check_offset = check_offset
        #     except NoSectionError:
        #         raw_input("\n ---- FAILED TO LOAD THE CONFIG FILE ---- \n\n")
        #         # printanswer = raw_input("print CWD? ")
        #         # if printanswer in ["yes", "1", '1']:
        #         #     print "CWD: ", os.getcwd()
        #         configfile = raw_input("Type the path of the config file again to try one more time: ")
        #         answer = raw_input("Try one more time?")
        #         GoOn = answer in ["yes", "True", "true", "1", '1']
        #         if not GoOn:
        #             self.ShowAndWait = bool(int(raw_input("Decided to not go on. Set `ShowAndWait` manually (1 or 0): ")))
        #             self.SaveMCData = bool(int(raw_input("Set `SaveMCData` manually (1 or 0): ")))
        #             self.check_offset = bool(int(raw_input("Set `check_offset` manually (1 or 0): ")))
        #     else:
        #         GoOn = False
        self.ShowAndWait = False
        self.SaveMCData = True
        self.check_offset = True

    def CheckOffset(self):
        if self.check_offset and not self.run_object.IsMonteCarlo:
            offset_canvas = ROOT.TCanvas("offset_canvas", "offset_canvas")
            th1 = ROOT.TH1D("th1", "calibration offset", 21, -10.5, 10.5)
            self.track_info.Draw("calib_offset>>th1", "calibflag==1 && calib_offset < 50")
            pad = offset_canvas.GetPad(0)
            pad.SetLogy()
            offset = int(th1.GetBinCenter(th1.GetMaximumBin()))
            if th1.GetEntries() == 0:
                print "Offset analysis failed."
                offset_canvas.Update()
                self.SavePlots("Offset.png")
                return 0
            good = 1.*th1.GetMaximum()/th1.GetEntries()
            th1.SetTitle(th1.GetTitle()+" ({0:0.2f}% in peak)".format(100.*good))
            offset_canvas.Update()
            self.SavePlots("Offset.png")
            if good < 0.3: # less than 30% aligned hits
                print "\nINFO: BAD TIMING ALIGNMENT RUN!\n"
                self.TimingAlignmentFailed = True
            print "MOST COMMON OFFSET: {0:0.0F}\n".format(offset)
            return offset
        else:
            return 0

    def DoAnalysis(self,minimum_bincontent = 1):
        '''
        Create a bin collection object as self.Pad and load data from ROOT TTree
        into the Pad object. Then get the 2-dim signal distribution from self.Pad
        :param minimum_bincontent: Bins with less hits are ignored
        :return: -
        '''
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer" # bins with less hits are ignored
        self.minimum_bincontent = minimum_bincontent
        if not self.run_object.IsMonteCarlo:
            self.track_info = self.rootfile.Get("track_info") # Get TTree called "track_info"
        # create a bin collection object:
        self.Pad = BinCollection(*self.config_object.Get2DAttributes(), parent_analysis_obj=self)

        # Check for misalignment in track vs signal using calibflag:
        offset = self.CheckOffset()

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
                    self.VerbosePrint("try a Signal Distribution")
                    if not self.run_object.DataIsMade:
                        self.run_object.Simulate() # if draw=False the first distribution will be taken
                    for i in xrange(self.run_object.NumberOfHits):
                        x_ = self.run_object.Data["track_x"][i]
                        y_ = self.run_object.Data["track_y"][i]
                        signal_ = self.run_object.Data["integral50"][i]
                        self.Pad.Fill(x_, y_, signal_)

                except IndexError:
                    pass
                else:
                    GoOn = False
                d += 1


        self.Pad.MakeFits()

        self.Signal2DDistribution = self.Pad.GetMeanSignalDistribution(self.minimum_bincontent)
        self.Signal2DDistribution.SetDirectory(0)
        self.Signal2DDistribution.SetStats(False)
        self.Signal2DDistribution.GetXaxis().SetTitle("pos x / cm")
        self.Signal2DDistribution.GetYaxis().SetTitle("pos y / cm")
        self.Signal2DDistribution.GetYaxis().SetTitleOffset(1.4)

        self.Checklist["DoAnalysis"] = True

    def CreatePlots(self,saveplots = False,savename = "2DSignalDistribution",ending="png",saveDir = "Results/", show3d = False):
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
        ROOT.gStyle.SetPalette(53) # Dark Body Radiator palette
        ROOT.gStyle.SetNumberContours(999)
        if show3d:
            self.Signal2DDistribution.Draw("SPEC dm(2,10) pa(1,1,1) ci(1,1,1) a(15,45,0) s(1,1)")
            savename += "3D"
        else:
            self.Signal2DDistribution.Draw("colz")
        self.IfWait("2d drawn")
        if saveplots:
            self.SavePlots(savename, ending, saveDir)

    def CreateMeanSignalHistogram(self, saveplots = False, savename = "MeanSignalDistribution",ending="png",saveDir = "Results/", show = False):
        if not hasattr(self, "Pad"):
            print "Performing auto analysis with minimum binhits 100"
            self.DoAnalysis(minimum_bincontent=100)

        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        if show:
            if hasattr(self, "signal_canvas"):
                self.signal_canvas.Clear()
            else:
                self.signal_canvas = ROOT.TCanvas()
                ROOT.SetOwnership(self.signal_canvas, False)
                ROOT.gStyle.SetPalette(53) # Dark Body Radiator palette
                ROOT.gStyle.SetNumberContours(999)
        # print self.Signal2DDistribution
        minimum = self.Signal2DDistribution.GetMinimum()
        maximum = self.Signal2DDistribution.GetMaximum()
        self.MeanSignalHisto_name = "MeanSignalHisto"+str(self.GLOBAL_COUNT)
        self.MeanSignalHisto = ROOT.TH1D(self.MeanSignalHisto_name, "Mean Signal Histogram", 100, minimum,maximum)
        self.GLOBAL_COUNT += 1

        if show:
            self.signal_canvas.SetName("signal_canvas")
            self.signal_canvas.SetTitle("Mean Signal Distribution")

        nbins = (self.Signal2DDistribution.GetNbinsX()+2)*(self.Signal2DDistribution.GetNbinsY()+2)
        for i in xrange(nbins):
            bincontent = self.Signal2DDistribution.GetBinContent(i)
            binhits = self.Pad.counthisto.GetBinContent(i)
            if binhits >= self.minimum_bincontent:
                self.MeanSignalHisto.Fill(bincontent)
        if show:
            self.MeanSignalHisto.Draw()

        if saveplots:
            self.SavePlots(savename, ending, saveDir)

        self.Checklist["MeanSignalHisto"] = True

    def CreateBoth(self,saveplots = False,savename = "SignalDistribution",ending="png",saveDir = None, PS=False, test=""):
        self.combined_canvas_name = "combined_canvas"+test
        self.combined_canvas = ROOT.gROOT.GetListOfCanvases().FindObject(self.combined_canvas_name)
        if not self.combined_canvas:
            self.combined_canvas = ROOT.TCanvas(self.combined_canvas_name,"Combined Canvas",1000,500)
            ROOT.SetOwnership(self.combined_canvas, False)
            self.combined_canvas.Divide(2,1)

        #self.CreatePlots(False)
        self.CreateMeanSignalHistogram(saveplots=False, show=False)

        self.combined_canvas.cd(1)
        #ROOT.gStyle.SetPalette(55)
        # ROOT.gStyle.SetNumberContours(999)
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        self.Signal2DDistribution.Draw("colz")#"CONT1Z")#)'colz')
        self.combined_canvas.cd(2)

        if PS: #if photoshop mode, fill histogram pink
            ROOT.gStyle.SetHistFillColor(6)
            ROOT.gStyle.SetHistFillStyle(1001)
        else:
            ROOT.gStyle.SetHistFillColor(7)
            ROOT.gStyle.SetHistFillStyle(3003)

        self.MeanSignalHisto.UseCurrentStyle()
        self.MeanSignalHisto.GetXaxis().SetTitle("Signal response")
        self.MeanSignalHisto.Draw()
        self.combined_canvas.cd()
        self.combined_canvas.Update()

        savename = self.run_object.diamond.Specifications["Name"]+"_"+self.run_object.diamond.Specifications["Irradiation"]+"_"+savename+"_"+str(self.run_object.run_number) # diamond_irradiation_savename_runnr
        if saveplots:
            self.SavePlots(savename, ending, saveDir)
            self.SavePlots(savename, "root", saveDir)
        if PS:
            ROOT.gStyle.SetHistFillColor(7)
            ROOT.gStyle.SetHistFillStyle(3003)
        self.IfWait("Combined 2D Signal DistributionsShown")

    def CreateHitsDistribution(self,saveplot = False, drawoption = "colz", RemoveLowStatBins = 0): # add palette!
        if RemoveLowStatBins > 0:
            if type(RemoveLowStatBins) == t.BooleanType:
                RemoveLowStatBins = 10
            bins = (self.Pad.counthisto.GetNbinsX() + 2)*(self.Pad.counthisto.GetNbinsY() + 2)
            for bin in xrange(bins):
                if self.Pad.counthisto.GetBinContent(bin) < RemoveLowStatBins:
                    coordinates = self.Pad.GetBinCenter(bin)
                    content = self.Pad.counthisto.GetBinContent(bin)
                    self.Pad.counthisto.Fill(coordinates[0], coordinates[1], -content)
            extension = "_"+str(RemoveLowStatBins)
        else:
            extension = ""

        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        canvas = ROOT.TCanvas("canvas", "Hits", 500, 500) # adjust the width slightly
        canvas.cd()
        self.Pad.counthisto.SetStats(False)
        self.Pad.counthisto.Draw(drawoption)#"surf2")
        #self.Pad.counthisto.Draw("CONT1 SAME")
        if saveplot:
            self.SavePlots("Hits_Distribution"+extension, "png")
        self.IfWait("Hits Distribution shown")
        self.Checklist["HitsDistribution"] = True

    def FindMaxima(self,show=False, binning = 200, minimum_bincontent = 30):
        # add a config file for minimum_bincontent
        if not hasattr(self, "ExtremeAnalysis"):
            self.ExtremeAnalysis = Analysis(self.run_object,Config(binning),verbose=self.verbose)
            self.ExtremeAnalysis.DoAnalysis(minimum_bincontent)
        if show:
            # self.ExtremeAnalysis.CreatePlots()
            self.ExtremeAnalysis.CreateBoth(saveplots=False, savename="SignalDistribution_MAXSearch", test="maxima") # destroys combined_canvas
        self.ExtremeAnalysis.Pad.FindMaxima(minimum_bincount=minimum_bincontent,show=show)

    def FindMinima(self,show=False, binning = 200, minimum_bincontent = 30):
        if not hasattr(self, "ExtremeAnalysis"):
            self.ExtremeAnalysis = Analysis(self.run_object,Config(binningsize=binning),verbose=self.verbose)
            self.ExtremeAnalysis.DoAnalysis(minimum_bincontent)
        if show:
            self.ExtremeAnalysis.CreatePlots()
            # self.ExtremeAnalysis.CreateBoth(saveplots=True, savename="SignalDistribution_MINSearch") # distroys combined_canvas
        self.ExtremeAnalysis.Pad.FindMinima(minimum_bincount=minimum_bincontent,show=show)

    def WaldWolfowitzRunsTest(self):
        pass

    def Chi2Test(self):
        pass

    def GetMPVSigmas(self,show = False, minimum_counts = 10):
        '''
        Returns MPVs and signals form all Bins of current run as Lists
        :return:
        '''
        MPVs = []
        MPVErrs = []
        Sigmas = []
        SigmaErrs = []
        for bin in self.Pad.ListOfBins:
            entries = bin.GetEntries()
            if entries >= minimum_counts:
                MPVs.append(bin.Fit["MPV"])
                MPVErrs.append(bin.Fit["MPVErr"])
                Sigmas.append(bin.Fit["Sigma"])
                SigmaErrs.append(bin.Fit["SigmaErr"])

        if show:
            canvas = ROOT.TCanvas("MPVSigmaCanvas", "MPVSigmaCanvas")
            canvas.cd()
            graph = ROOT.TGraphErrors()
            count = 0
            for i in xrange(len(MPVs)):
                graph.SetPoint(count, MPVs[i], Sigmas[i])
                graph.SetPointError(count, MPVErrs[i], SigmaErrs[i])
                count += 1
            graph.Draw("AP")
            self.IfWait("MPV vs Sigma shown")

        return MPVs, Sigmas, MPVErrs, SigmaErrs

    def GetSignalHeight(self, min_percent = 5, max_percent = 99):
        if not hasattr(self, "ExtremeAnalysis"):
            self.FindMaxima()
        if not hasattr(self.ExtremeAnalysis, "MeanSignalHisto"):
            self.ExtremeAnalysis.CreateMeanSignalHistogram()
        q = array('d', [1.*min_percent/100., 1.*max_percent/100.])
        y = array('d', [0,0])
        self.ExtremeAnalysis.MeanSignalHisto.GetQuantiles(2, y, q)
        SignalHeight = y[1]/y[0]-1.
        self.VerbosePrint('\nApproximated Signal Amplitude: {0:0.0f}% - ({1:0.0f}%/{2:0.0f}% Quantiles approximation)\n'.format(100.*(SignalHeight), max_percent, min_percent))
        self.ExtremaResults['SignalHeight'] = SignalHeight
        self.ExtremeAnalysis.ExtremaResults['SignalHeight'] = SignalHeight
        return SignalHeight

    def ShowTotalSignalHistogram(self, save = False, scale = False, showfit=False):
        if  hasattr(self, "Pad"):
            self.Pad.CreateTotalSignalHistogram(saveplot=save, scale=scale, showfit=showfit)
        elif hasattr(self, "ExtremeAnalysis") and hasattr(self.ExtremeAnalysis, "Pad"):
            self.ExtremeAnalysis.Pad.CreateTotalSignalHistogram(saveplot=save, scale=scale, showfit=showfit)
            if showfit:
                self.GetSignalHistoFitResults() # Get the Results from self.Extremeanalysis.SignalHistoFitResults to self.SignalHistoFitResults
        else:
            self.DoAnalysis()
            self.Pad.CreateTotalSignalHistogram(saveplot=save, scale=scale, showfit=showfit)

    def GetSignalHistoFitResults(self):
        if self.SignalHistoFitResults["Peak"] != None:
            FitFunction = self.SignalHistoFitResults["FitFunction"]
            Peak = self.SignalHistoFitResults["Peak"]
            FWHM = self.SignalHistoFitResults["FWHM"]
            Chi2 = self.SignalHistoFitResults["Chi2"]
            NDF = self.SignalHistoFitResults["NDF"]
            return FitFunction, Peak, FWHM, Chi2, NDF
        elif hasattr(self, "ExtremeAnalysis") and self.ExtremeAnalysis.SignalHistoFitResults["Peak"] != None:
            FitFunction = self.ExtremeAnalysis.SignalHistoFitResults["FitFunction"]
            Peak = self.ExtremeAnalysis.SignalHistoFitResults["Peak"]
            FWHM = self.ExtremeAnalysis.SignalHistoFitResults["FWHM"]
            Chi2 = self.ExtremeAnalysis.SignalHistoFitResults["Chi2"]
            NDF = self.ExtremeAnalysis.SignalHistoFitResults["NDF"]
            self.SignalHistoFitResults["FitFunction"] = FitFunction
            self.SignalHistoFitResults["Peak"] = Peak
            self.SignalHistoFitResults["FWHM"] = FWHM
            self.SignalHistoFitResults["Chi2"] = Chi2
            self.SignalHistoFitResults["NDF"] = NDF
            return FitFunction, Peak, FWHM, Chi2, NDF
        else:
            self.ShowTotalSignalHistogram(save=False, showfit=True)
            if (self.SignalHistoFitResults["Peak"] != None) or (hasattr(self, "ExtremeAnalysis") and self.ExtremeAnalysis.SignalHistoFitResults["Peak"] != None):
                self.GetSignalHistoFitResults()
            else:
                assert(False), "BAD SignalHistogram Fit, Stop program due to possible infinity loop"

    def SignalTimeEvolution(self, Mode="Mean", show=True, time_spacing = 3, save = True): # not show: save evolution data / Comment more
        '''
        Creates Signal vs time plot. The time is bunched into time buckets of width time_spaceing,
        from all Signals inside one time bucket the Mean or the MPV is evaluated and plotted in a TGraph.
        The Mean is taken from the histogramm of the particular time bucket. The MPV is simply the
        Position (Center) of the Bin containing the most entries in the time bucket histogram.
        :param Mode: What should be plotted as y: either "Mean" or "MPV"
        :param show:
        :param time_spacing: time bucket width in minutes
        :param save:
        :return:
        '''
        assert(Mode in ["Mean", "MPV"]), "Wrong Mode, Mode has to be `Mean` or `MPV`"
        if not self.run_object.IsMonteCarlo:
            results = {}
            SignalEvolution = ROOT.TGraphErrors()
            SignalEvolution.SetNameTitle("SignalEvolution", "Signal Response Evolution")

            self.track_info.GetEntry(1)
            starttime = self.track_info.time_stamp
            self.track_info.GetEntry(self.track_info.GetEntries())
            # endtime = self.track_info.time_stamp
            ok_deltatime = 0
            TimeERROR = False
            for i in xrange(self.track_info.GetEntries()-1):
                i += 1
                # read the ROOT TTree
                self.track_info.GetEntry(i)
                signal_ = abs(self.track_info.integral50)
                time_ = self.track_info.time_stamp
                deltatime_min = (time_-starttime)/60.
                if deltatime_min < -1: # possible time_stamp reset
                    deltatime_min = ok_deltatime + time_/60.
                    TimeERROR = True
                else:
                    ok_deltatime = deltatime_min
                time_bucket = int(deltatime_min)/int(time_spacing)*int(time_spacing)
                calibflag = self.track_info.calibflag
                if calibflag == 0:
                    try:
                        results[time_bucket].append(signal_)
                    except KeyError:
                        results[time_bucket] = [signal_]
            if TimeERROR:
                print "\n\n\n\n\nWARNING: Error occured in time flow of Run "+str(self.run_object.run_number)+". Possible reset of timestamp during data taking!\n\n\n\n\n"

            time_buckets = results.keys()
            time_buckets.sort()
            count = 0
            for t in time_buckets:
                histo = ROOT.TH1D("SignalEvolution_time_"+str(t), "Signal Response Histogram for t = {0:0.0f}-{1:0.0f}min".format(t, t+int(time_spacing)), 500, 0, 500)
                for i in xrange(len(results[t])):
                    histo.Fill(results[t][i])
                if Mode == "Mean":
                    SignalEvolution.SetPoint(count, t, histo.GetMean())
                    SignalEvolution.SetPointError(count, 0, histo.GetRMS()/np.sqrt(histo.GetEntries()))
                    signalname = "Mean of Signal Response"
                elif Mode == "MPV":
                    SignalEvolution.SetPoint(count, t, histo.GetBinCenter(histo.GetMaximumBin()))
                    SignalEvolution.SetPointError(count, 0, 0)
                    signalname = "MPV of Signal Respose"
                count += 1
                del histo
                ROOT.gROOT.Delete("SignalEvolution_time_"+str(t))

            SignalEvolution.GetXaxis().SetTitle("Time / min")
            SignalEvolution.GetYaxis().SetTitle(signalname+" / ADC Units")
            SignalEvolution.Draw("ALP*")
            if save:
                self.SavePlots("SignalTimeEvolution"+Mode+".png")
            if TimeERROR:
                SignalEvolution.SaveAs(self.SaveDirectory+"ERROR_SignalTimeEvolution"+Mode+".root")
            self.IfWait("Showing Signal Time Evolution..")


    def ExportMC(self, MCDir = "MCInputs/"):
        '''
        Hit distribution
        :return:
        '''
        if not self.run_object.IsMonteCarlo:
            if not os.path.exists(MCDir):
                os.makedirs(MCDir)
            self.Pad.counthisto.SaveAs(MCDir+str(self.run_object.run_number)+"counthisto.root")
            self.VerbosePrint("CountHisto exported..")
        else:
            print "INFO: Monte Carlo run can not be exported as MC input"
