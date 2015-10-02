import ROOT
from ROOT import gROOT
from RunClass import Run
import numpy as np
from AbstractClasses.Elementary import Elementary
from AbstractClasses.PreAnalysisPlot import PreAnalysisPlot
from AbstractClasses.ConfigClass import *
from AbstractClasses.RunClass import Run
from AbstractClasses.Cut import Cut
from AbstractClasses.Langau import Langau
import os
import copy, collections, numpy
import sys
import json
import numpy as np
from BinCollection import BinCollection
from array import array
import pickle
import types as t
import ConfigParser
from ConfigParser import NoSectionError

class Analysis(Elementary):
    '''
    An Analysis Object contains all Data and Results of a SINGLE run.
    '''

    def __init__(self, run, diamonds=None, verbose = False, maskfilename=""):
        '''
        Initializes the Analysis object.
        An Analysis Object collects all the information and data for the analysis of
        one single run.
        The most important fields are:
            .run        - instance of Run class, containing run info and data
            .cut        - dict containing two instances of Cut class
            .RunInfo    - dict containing run infos
            .Pads       - dict containing two intances of
        The most important methods are:
            .CreateMeanSignalHistogram(channel)
            .Draw(varexp)
            .DrawRunInfo()
            .GetCut(channel)
            .GetEventAtTime(dt)
            .MakePreAnalysis()
            .ShowPulserRate()

        :param run:
            run object of type "Run"
        :param diamonds:
            An integer number defining the diamonds activated for analysis:
            0x1=ch0 (diamond 1)
            0x2=ch3 (diamond 2)
            e.g. 1 -> diamond 1, 2 -> diamond 2, 3 -> diamonds 1&2

        :param verbose:
            if True, verbose printing is activated
        :return: -
        '''
        Elementary.__init__(self, verbose=verbose)
        if not isinstance(run, Run):
            assert (type(run) is t.IntType), "run has to be either an instance of Run or run number (int)"
            if diamonds != None:
                run = Run(run, diamonds=diamonds, maskfilename=maskfilename)
            else:
                run = Run(run, maskfilename=maskfilename)
        else:
            assert(run.run_number > 0), "No run selected, choose run.SetRun(run_nr) before you pass the run object"
        self.run = run
        self.run.analysis = self
        self.config_object = {
            0: BinCollectionConfig(run=self.run, channel=0),
            3: BinCollectionConfig(run=self.run, channel=3)
        }
        self.RunInfo = copy.deepcopy(run.RunInfo)

        self.Checklist = { # True if Plot was created # -> USED?!
            "RemoveBeamInterruptions": False,
            "LoadTrackData": False,
            "MeanSignalHisto": {
                0: False,
                3: False
            },
            "HitsDistribution": False,
            "FindExtrema": {
                "Maxima": {
                    0: False,
                    3: False
                },
                "Minima": {
                    0: False,
                    3: False
                }
            }
        }

        extremaResults ={
           "TrueNPeaks": None, # MC true number of peaks
           "FoundNMaxima": None, # number of Maxima found
           "FoundMaxima": None, # Maxima found as list [(x1, y1), (x2, y2), ...]
           "FoundNMinima": None, # number of Minima found
           "FoundMinima": None, # Minima found as list [(x1, y1), (x2, y2), ...]
           "Ninjas": None, # number of true peaks not found (only available when MC Run)
           "Ghosts": None, # nunmber of wrong Maxima found (only available when MC Run)
           "SignalHeight": 0.
        }
        self.extremaResults = {
            0: copy.deepcopy(extremaResults),
            3: copy.deepcopy(extremaResults)
        }

        signalHistoFitResults = {
           "FitFunction": None,
           "Peak": None,
           "FWHM": None,
           "Chi2": None,
           "NDF": None,
        }
        self.signalHistoFitResults = {
            0: copy.deepcopy(signalHistoFitResults),
            3: copy.deepcopy(signalHistoFitResults)
        }
        self.cut = {
            0: Cut(self, 0),
            3: Cut(self, 3)
        }

    def LoadConfig(self):
        configfile = "Configuration/AnalysisConfig_"+self.TESTCAMPAIGN+".cfg"
        parser = ConfigParser.ConfigParser()
        output = parser.read(configfile)
        print " ---- Config Parser Read ---- \n - ", output, " -\n"


        self.showAndWait = parser.getboolean("DISPLAY","ShowAndWait")
        self.saveMCData = parser.getboolean("SAVE","SaveMCData")
        self.signalname = parser.get("BASIC", "signalname")
        self.pedestal_correction = parser.getboolean("BASIC", "pedestal_correction")
        self.pedestalname = parser.get("BASIC", "pedestalname")
        # self.cut = parser.get("BASIC", "cut")
        # self.excludefirst = parser.getint("BASIC", "excludefirst")
        # self.excludeBeforeJump = parser.getint("BASIC", "excludeBeforeJump")
        # self.excludeAfterJump = parser.getint("BASIC", "excludeAfterJump")

        self.loadMaxEvent = parser.getint("TRACKING", "loadMaxEvent")
        self.minimum_bincontent = parser.getint("TRACKING", "min_bincontent")
        self.minimum_bincontent = parser.getint("TRACKING", "min_bincontent")
        self.PadsBinning = parser.getint("TRACKING", "padBinning")


        if not self.pedestal_correction:
            self.signaldefinition = self.signalname+"[{channel}]"
        else:
            self.signaldefinition = self.signalname+"[{channel}]-"+self.pedestalname+"[{channel}]"

    def GetCut(self, channel, gen_PulserCut=True, gen_EventRange=True, gen_ExcludeFirst=True):
        '''
        Returns the full Cut string, which corresponds to the given
        channel (i.e. diamond). If some of the other arguments are
        False, these cuts will be ignored.
        :param channel:
        :param gen_PulserCut:
        :param gen_EventRange:
        :param gen_ExcludeFirst:
        :return:
        '''
        return self.cut[channel].GetCut(gen_PulserCut=gen_PulserCut, gen_EventRange=gen_EventRange, gen_ExcludeFirst=gen_ExcludeFirst)

    def MakePreAnalysis(self, channel=None, mode="mean", binning=5000, savePlot=True, canvas=None, setyscale_sig=None, setyscale_ped=None):
        '''
        Creates an overview plot, containing the signal time-evolution
        as an average graph and as a 2-dimensional histogram.
        An additional third plot shows the pedestal time-evolution.
        For these plots is no tracking information required.
        :return:
        '''
        if channel != None:
            channels = [channel]
        else:
            channels = self.run.GetChannels()
        #c1 = ROOT.TCanvas("c1", "c1")
        if not hasattr(self, "preAnalysis"):
            self.preAnalysis = {}

        for ch in channels:
            self.preAnalysis[ch] = PreAnalysisPlot(analysis=self, channel=ch, canvas=canvas, binning=binning)
            self.preAnalysis[ch].Draw(mode=mode, savePlot=savePlot, setyscale_sig=setyscale_sig, setyscale_ped=setyscale_ped)

    def ShowPulserRate(self, binning=2000, canvas=None, savePlot=True):
        '''
        Shows the fraction of accepted pulser events as a function of
        event numbers. Peaks appearing in this graph are most likely
        beam interruptions.
        :param binning:
        :param canvas:
        :return:
        '''
        assert(binning>=100), "binning too low"
        binning = int(binning)

        if canvas == None:
            self.pulserRateCanvas = ROOT.TCanvas("pulserratecanvas{run}".format(run=self.run.run_number), "Pulser Rate Canvas")
        else:
            self.pulserRateCanvas = canvas
        self.pulserRateCanvas.cd()

        self.pulserRateGraph = ROOT.TGraph()
        self.pulserRateGraph.SetNameTitle("pulserrategraph{run}".format(run=self.run.run_number), "Pulser Rate")
        nbins = int(self.run.tree.GetEntries())/binning

        for i in xrange(nbins):
            pulserevents = self.run.tree.Draw("1", "pulser", "", binning, i*binning)
            pulserrate = 1.*pulserevents/binning
            self.pulserRateGraph.SetPoint(i, (i+0.5)*binning, pulserrate)
        self.pulserRateGraph.Draw("AL")
        self.pulserRateGraph.GetXaxis().SetTitle("Event Number")
        self.pulserRateGraph.GetYaxis().SetTitle("Fraction of Pulser Events")
        self.pulserRateGraph.GetYaxis().SetTitleOffset(1.2)
        self.DrawRunInfo(canvas=self.pulserRateCanvas)
        self.pulserRateCanvas.Update()
        if savePlot: self.SavePlots("Run{run}_PulserRate.png".format(run=self.run.run_number), canvas=self.pulserRateCanvas)

    def DrawRunInfo(self, channel=None, canvas=None, diamondinfo=True, showcut=False, comment=None, infoid="", userHeight=None, userWidth=None):
        '''
        Draws the run infos inside the canvas. If no canvas is given, it
        will be drawn into the active Pad. If The channel number will be
        given, channel number and diamond name will be drawn.
        :param channel:
        :param canvas:
        :param diamondinfo:
        :param showcut:
        :param comment:
        :param infoid:
        :param userHeight:
        :param userWidth:
        :return:
        '''
        self.run.DrawRunInfo(channel=channel, canvas=canvas, diamondinfo=diamondinfo, showcut=showcut, comment=comment, infoid=infoid, userHeight=userHeight, userWidth=userWidth)

    def ShowFFT(self, drawoption="", cut=None, channel=None, startevent=0, endevent=10000000, savePlots=True, canvas=None):
        '''
        Draws the fft_mean VS 1/fft_max scatter plot.
        For the drawoption ROOT TH2D Drawoption can be used, as well as
        "mc" for a special multicolor plot, which differentiates between
        different waveform types.
        The cut can be None, for using a default cut, True for using the
        cut from the Analysis instance or can be set customly.
        :param drawoption:
            "" draws all events with the given cut
            "mc" for multicolor: draw different colors for different
                waveform types
            "col" gives colored TH2D plot (also possible: "colz")
        :param cut:
            if True: apply cut from Analysis instance
            if None: "!pulser" set as default (except for "mc"
                drawoption)
            custom: cut="!pulser&&!is_saturated[{channel}]"
                (--> syntax "[{channel}]" important for vectors!)
        :param channel:
            if None: takes the activated channels (diamonds) from
                Analysis instance
            else: 0 or 3 to select one particular channel
        :param startevent:
            first event number to draw
        :param endevent:
            last event number to draw
        :param savePlots:
            True to save the generated scatter plots
        :return:
        '''
        if channel != None:
            channels = [channel]
        else:
            channels = self.run.GetChannels()

        if drawoption in ["MC", "mc", "Mc", "mC"]:
            drawoption = ""
            multicolor = True
            if cut == None:
                cut = ""
        else:
            multicolor = False

        if cut == None:
            cut = "!pulser"

        events = endevent-startevent
        if events<10000000:
            if endevent >= 10000000: endevent = self.run.tree.GetEntries()
            comment = "Events: {start}-{end}".format(start=startevent, end=endevent)

        else:
            comment = None
        if canvas==None:
            self.fftcanvas = ROOT.TCanvas("fftcanvas", "fftcanvas", len(channels)*500, 500)
        else:
            #assert(isinstance(canvas, ROOT.TCanvas))
            self.fftcanvas = canvas

        self.fftcanvas.Divide(len(channels), 1)

        if len(channels) ==2:
            namesuffix = ""
        else:
            namesuffix = "_Ch{ch}".format(ch=channels[0])

        self.fftcanvas.cd()
        i = 1
        self.fftHistos = {}
        self.fftlegend = {}
        ROOT.gStyle.SetOptStat(0)
        for ch in channels:
            if cut == True:
                thiscut = self.GetCut(channel=ch)
                thisusercut = self.GetUserCutString(channel=ch)
                events = self.GetNEventsCut(channel=ch)
                startevent = self.GetMinEventCut(channel=ch)
            else:
                thiscut = cut.format(channel=ch)
                thisusercut = thiscut
            self.fftHistos[ch] = ROOT.TH2D("fft_ch{channel}".format(channel=ch), "FFT {"+thisusercut+"}", 5000, 2e-6, 0.0025, 5000, 1e1, 1e4)
            self.fftcanvas.cd(i)
            ROOT.gPad.SetLogy()
            ROOT.gPad.SetLogx()
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}".format(channel=ch), thiscut, "", events, startevent)
            self.fftHistos[ch].SetTitle("{diamond} ".format(diamond=self.run.diamondname[ch])+self.fftHistos[ch].GetTitle())
            self.fftHistos[ch].GetXaxis().SetTitle("1/fft_max")
            self.fftHistos[ch].GetYaxis().SetTitle("fft_mean")
            #self.fftHistos[ch].Draw(drawoption)
            if multicolor:
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_sat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), "is_saturated[{channel}]".format(channel=ch), "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_sat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), (thiscut+"&&is_saturated[{channel}]").format(channel=ch), "", events, startevent)
                saturated_histo = ROOT.gROOT.FindObject("fft_ch{channel}_sat".format(channel=ch))
                saturated_histo.SetMarkerStyle(1)
                saturated_histo.SetMarkerColor(6)
                saturated_histo.SetFillColor(6)
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_med(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), "abs(median[{channel}])>8".format(channel=ch), "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_med(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), (thiscut+"&&abs(median[{channel}])>8").format(channel=ch), "", events, startevent)
                median_histo = ROOT.gROOT.FindObject("fft_ch{channel}_med".format(channel=ch))
                median_histo.SetMarkerStyle(1)
                median_histo.SetMarkerColor(8)
                median_histo.SetFillColor(8)
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_flat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), "sig_spread[{channel}]<10".format(channel=ch), "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_flat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), (thiscut+"&&sig_spread[{channel}]<10").format(channel=ch), "", events, startevent)
                flat_histo = ROOT.gROOT.FindObject("fft_ch{channel}_flat".format(channel=ch))
                flat_histo.SetMarkerStyle(1)
                flat_histo.SetMarkerColor(4)
                flat_histo.SetFillColor(4)
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_pulser(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), "pulser", "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_pulser(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), (thiscut+"&&pulser").format(channel=ch), "", events, startevent)
                pulser_histo = ROOT.gROOT.FindObject("fft_ch{channel}_pulser".format(channel=ch))
                pulser_histo.SetMarkerStyle(1)
                pulser_histo.SetMarkerColor(2)
                pulser_histo.SetFillColor(2)
                self.fftHistos[ch].Draw("")
                saturated_histo.Draw("same")
                median_histo.Draw("same")
                pulser_histo.Draw("same")
                flat_histo.Draw("same")
                self.fftlegend[ch] = ROOT.TLegend(0.1,0.1,0.3,0.3)
                self.fftlegend[ch].AddEntry(pulser_histo, "pulser", "f")
                self.fftlegend[ch].AddEntry(saturated_histo, "saturated", "f")
                self.fftlegend[ch].AddEntry(median_histo, "|median|>8", "f")
                self.fftlegend[ch].AddEntry(flat_histo, "sig_spread<10", "f")
                self.fftlegend[ch].Draw("same")
                self.fftcanvas.Update()
                self.DrawRunInfo(channel=ch, comment=comment, canvas=self.fftcanvas)
            else:
                self.fftHistos[ch].Draw(drawoption)
                self.fftcanvas.Update()
                self.DrawRunInfo(channel=ch, comment=comment, canvas=self.fftcanvas)
            i+=1
        self.fftcanvas.Update()
        if savePlots:
            savename = "Run{run}_FFT"+namesuffix
            self.SavePlots(savename.format(run=self.run.run_number), ending="png", canvas=self.fftcanvas, subDir="png/")
            #self.SavePlots(savename.format(run=self.run.run_number), ending="root", canvas=self.fftcanvas, subDir="root/")

        self.IfWait("FFT shown..")

    def GetIncludedEvents(self, maxevent, channel=None):
        '''
        Returns list of all event numbers, which are neither excluded by
        beaminerruptions, by the excludeFirst cut, or are outside of a
        event range window.
        :return: list of included event numbers
        '''
        if channel == None:
            minevent0 = self.cut[0].GetMinEvent()
            minevent3 = self.cut[3].GetMinEvent()
            minevent = min(minevent0, minevent3)
            if maxevent == None:
                maxevent0 = self.cut[0].GetMaxEvent()
                maxevent3 = self.cut[3].GetMaxEvent()
                maxevent = max(maxevent0, maxevent3)
            excluded = [i for i in np.arange(0, minevent)] # first events
            if self.cut[0]._cutTypes["noBeamInter"] and self.cut[3]._cutTypes["noBeamInter"]:
                self.cut[0].GetBeamInterruptions()
                for i in xrange(len(self.cut[0].jumpsRanges["start"])):
                    excluded += [i for i in np.arange(self.cut[0].jumpsRanges["start"][i], self.cut[0].jumpsRanges["stop"][i]+1)] # events around jumps
            excluded.sort()
            all_events = np.arange(0, maxevent)
            included = np.delete(all_events, excluded)
            return included
            # cut0 = self.cut[0].GetIncludedEvents(maxevent=maxevent)
            # cut3 = self.cut[3].GetIncludedEvents(maxevent=maxevent)
            # assert(cut0 == cut3), "Included Events not The same for both channels, select a particular channel"
            # return cut0
        else:
            assert(channel in [0,3])
            return self.cut[channel].GetIncludedEvents(maxevent=maxevent)

    def GetMinEventCut(self, channel=None):
        if channel == None:
            opt0 = self.cut[0].GetMinEvent()
            opt3 = self.cut[3].GetMinEvent()
            assert(opt0==opt3), "GetMinEvent Not the same for both cut channels. Choose a specific channel."
            return opt0
        else:
            assert(channel in [0,3])
            return self.cut[channel].GetMinEvent()

    def GetMaxEventCut(self, channel=None):
        '''
        Returns the highest event number, which satisfies the cut
        conditions. It is recommended to handle the channel ID. If no
        channel will be passed and the event numbers differ for both
        channels, an assertion error will be raised.
        :param channel:
        :return:
        '''
        if channel == None:
            opt0 = self.cut[0].GetMaxEvent()
            opt3 = self.cut[3].GetMaxEvent()
            assert(opt0==opt3), "GetMaxEvent Not the same for both cut channels. Choose a specific channel."
            return opt0
        else:
            assert(channel in [0,3])
            return self.cut[channel].GetMaxEvent()

    def GetNEventsCut(self, channel=None):
        '''
        Returns the number of events, which satisfies the cut
        conditions. It is recommended to handle the channel ID. If no
        channel will be passed and the event numbers differ for both
        channels, an assertion error will be raised.
        :param channel:
        :return:
        '''
        if channel == None:
            opt0 = self.cut[0].GetNEvents()
            opt3 = self.cut[3].GetNEvents()
            assert(opt0==opt3), "GetNEvents Not the same for both cut channels. Choose a specific channel."
            return opt0
        else:
            assert(channel in [0,3])
            return self.cut[channel].GetNEvents()

    def ShowSignalHisto(self, channel=None, canvas=None, drawoption="", cut="", color=None, normalized=True, drawruninfo=True, binning=600, xmin=None, xmax=None, savePlots=False, logy=False, gridx=False):
        '''
        Shows a histogram of the signal responses for the selected
        channels.
        :param channel:
        :param canvas:
        :param drawoption:
        :param cut:
        :param color:
        :param normalized:
        :param drawruninfo:
        :param binning:
        :param xmin:
        :param xmax:
        :param savePlots:
        :param logy:
        :param gridx:
        :return:
        '''
        self._ShowHisto(self.signaldefinition, channel=channel, canvas=canvas, drawoption=drawoption, cut=cut, color=color, normalized=normalized, infoid="SignalHisto", drawruninfo=drawruninfo, binning=binning, xmin=xmin, xmax=xmax,savePlots=savePlots, logy=logy, gridx=gridx)

    def ShowPedestalHisto(self, channel=None, canvas=None, drawoption="", color=None, normalized=True, drawruninfo=False, binning=600, xmin=None, xmax=None, savePlots=False, logy=False, cut=""):
        '''
        Shows a histogram of the pedestal signal responses for the
        selected channels.
        :param channel:
        :param canvas:
        :param drawoption:
        :param color:
        :param normalized:
        :param drawruninfo:
        :param binning:
        :param xmin:
        :param xmax:
        :param savePlots:
        :param logy:
        :return:
        '''
        self._ShowHisto(self.pedestalname+"[{channel}]", channel=channel, canvas=canvas, drawoption=drawoption, color=color, cut=cut, normalized=normalized, infoid="PedestalHisto", drawruninfo=drawruninfo, binning=binning, xmin=xmin, xmax=xmax,savePlots=savePlots, logy=logy)

    def ShowMedianHisto(self, channel=None, canvas=None, drawoption="", color=None, normalized=True, drawruninfo=True, binning=100, xmin=-30, xmax=30, savePlots=False, logy=True):
        '''
        Shows a histogram of the median for the selected
        channels.
        :param channel:
        :param canvas:
        :param drawoption:
        :param color:
        :param normalized:
        :param drawruninfo:
        :param binning:
        :param xmin:
        :param xmax:
        :param savePlots:
        :param logy:
        :return:
        '''
        self.medianhisto = self._ShowHisto("median[{channel}]", logy=logy, infoid="median", drawruninfo=drawruninfo, xmin=xmin, xmax=xmax, binning=binning, channel=channel, canvas=canvas, drawoption=drawoption, color=color, normalized=normalized, savePlots=False)
        if canvas==None:
            canvas = ROOT.gROOT.FindObject("mediancanvas")

        canvas.SetGridx()

        if self.medianhisto:
            self.medianhisto.GetXaxis().SetNdivisions(20)

        canvas.Update()
        if savePlots:
            if channel != None:
                dia = "_"+self.run.diamondname[channel]
            else:
                dia = ""
            self.SavePlots("Run{run}_MedianHisto{dia}.png".format(run=self.run.run_number, dia=dia))

    def ShowSignalPedestalHisto(self, channel, canvas=None, savePlots=True, cut="", normalized=True, drawruninfo=True, binning=600, xmin=None, xmax=None, logy=False, gridx=True):
        if canvas == None:
            self.signalpedestalcanvas = ROOT.TCanvas("signalpedestalcanvas{run}{channel}".format(run=self.run.run_number, channel=channel), "signalpedestalcanvas")
        else:
            self.signalpedestalcanvas = canvas

        self.ResetColorPalette()
        self.ShowPedestalHisto(channel=channel, canvas=self.signalpedestalcanvas, savePlots=False, cut=cut, normalized=normalized, drawruninfo=False, binning=binning, xmin=xmin, xmax=xmax, logy=logy)
        self.ShowSignalHisto(channel=channel, canvas=self.signalpedestalcanvas, drawoption="sames", savePlots=False, cut=cut, normalized=normalized, drawruninfo=False, binning=binning, xmin=xmin, xmax=xmax, logy=logy, gridx=gridx)

        if drawruninfo:
            self.DrawRunInfo(channel=channel, canvas=self.signalpedestalcanvas, infoid="signalpedestal"+str(self.run.run_number)+str(channel))

        if savePlots:
            for ending in ["png", "root"]:
                self.SavePlots(savename="Run{run}_SignalPedestal_Ch{channel}.{end}".format(run=self.run.run_number, channel=channel, end=ending), subDir=ending, canvas=self.signalpedestalcanvas)

        self.IfWait("ShowSignalPedeslHisto")


    def _ShowHisto(self, signaldef, channel=None, canvas=None, drawoption="", cut="", color=None, normalized=True, infoid="histo", drawruninfo=False, binning=600, xmin=None, xmax=None, savePlots=False, logy=False, gridx=False):
        '''

        :param channel:
        :param canvas:
        :param drawoption:
        :param normalized:
        :return:
        '''
        if channel != None:
            channels = [channel]
        else:
            channels = self.run.GetChannels()
        if canvas == None:
            canvas = ROOT.TCanvas(infoid+"canvas")
            self.ResetColorPalette()
        else:
            pass
            #drawoption = "sames"
        canvas.cd()

        if xmin==None:
            xmin = -100
        if xmax==None:
            xmax = 500

        ROOT.gStyle.SetOptStat(1)

        if color == None: color = self.GetNewColor()
        for ch in channels:
            if len(channels)>1 and drawoption=="" and ch==3:
                drawoption = "sames"
                color = self.GetNewColor()
            if cut == "":
                thiscut = self.GetCut(ch)
                thisusercut = self.GetUserCutString(channel=ch)
            else:
                thiscut = cut.format(channel=ch)
                thisusercut = thiscut
            print "making "+infoid+" using\nSignal def:\n\t{signal}\nCut:\n\t({usercut})\n\t{cut}".format(signal=signaldef, usercut=thisusercut, cut=thiscut)
            self.run.tree.Draw((signaldef+">>{infoid}{run}({binning}, {min}, {max})").format(infoid=(self.run.diamondname[ch]+"_"+infoid), channel=ch, run=self.run.run_number, binning=binning, min=xmin, max=xmax), thiscut, drawoption, self.GetNEventsCut(channel=ch), self.GetMinEventCut(channel=ch))
            canvas.Update()
            if logy: canvas.SetLogy()
            if gridx: canvas.SetGridx()
            histoname = "{infoid}{run}".format(infoid=(self.run.diamondname[ch]+"_"+infoid), run=self.run.run_number)
            histo = ROOT.gROOT.FindObject(histoname)

            if histo:
                histo.GetXaxis().SetTitle(signaldef.format(channel=""))
                histo.SetLineColor(color)
                histo.Draw(drawoption)
                histo.SetTitle("{signal} {cut}".format(signal=infoid, cut="{"+thisusercut+"}"))
                stats = histo.FindObject("stats")
                if stats:
                    stats.SetTextColor(color)
                    if channels.index(ch) == 1: # shift second statbox down
                        point = stats.GetBBoxCenter()
                        point.SetY(150)
                        stats.SetBBoxCenter(point)
            canvas.Modified()

            if normalized:
                histo.Scale(1./histo.GetMaximum())
                histo.Draw(drawoption)

        if drawruninfo:
            if len(channels) == 2:
                self.DrawRunInfo(canvas=canvas)
            else:
                self.DrawRunInfo(channel=channels[0], canvas=canvas)
        canvas.Update()
        self.IfWait(infoid+" shown")
        if savePlots:
            self.SavePlots("Run{run}_{signal}.png".format(run=self.run.run_number, signal=infoid), canvas=canvas, subDir=infoid+"/png/")
            self.SavePlots("Run{run}_{signal}.root".format(run=self.run.run_number, signal=infoid), canvas=canvas, subDir=infoid+"/root/")
        return histo

    def ShowDiamondCurrents(self):
        pass

    def LoadTrackData(self, minimum_bincontent=None): # min_bincontent in config file
        '''
        Create a bin collection object as self.Pads and load data from ROOT TTree
        into the Pad object. Then get the 2-dim signal distribution from self.Pads
        :param minimum_bincontent: Bins with less hits are ignored
        :return: -
        '''
        print "Loading Track information with \n\tmin_bincontent: {mbc}\n\tfirst: {first}\n\tmaxevent: {me}".format(mbc=self.minimum_bincontent, first=self.GetMinEventCut(channel=0), me=self.loadMaxEvent)
        if minimum_bincontent != None: self.minimum_bincontent = minimum_bincontent
        assert (self.minimum_bincontent > 0), "minimum_bincontent has to be a positive integer" # bins with less hits are ignored

        # create a bin collection object:
        self.Pads = {}
        for ch in [0,3]: # self.run.GetChannels():
            self.Pads[ch] = BinCollection(self, ch, *self.config_object[ch].Get2DAttributes())

        # fill two 2-dim histograms to collect the hits and signal strength
        x_ = {}
        y_ = {}
        channels = [0,3] # self.run.GetChannels()
        includedCh = {}
        if self.loadMaxEvent > 0:
            included = self.GetIncludedEvents(maxevent=self.loadMaxEvent) # all event numbers without jump events and initial cut
            includedCh[0] = self.GetIncludedEvents(maxevent=self.loadMaxEvent, channel=0)
            includedCh[3] = self.GetIncludedEvents(maxevent=self.loadMaxEvent, channel=3)
        else:
            included = self.GetIncludedEvents() # all event numbers without jump events and initial cut
            includedCh[0] = self.GetIncludedEvents(channel=0)
            includedCh[3] = self.GetIncludedEvents(channel=3)

        if not self.pedestal_correction:
            signaldef = "self.run.tree.{signal}[{channel}]"
        else:
            signaldef = "self.run.tree.{signal}[{channel}] - self.run.tree.{pedestal}[{channel}]"

        print "Loading tracking informations..."
        print "Signal def:\n\t{signal}".format(signal=signaldef)
        try:
            self.run.tree.diam1_track_x
        except AttributeError:
            print "\n\n--------------- ERROR --------------"
            print "NO TRACKING INFORMATION IN ROOT FILE"
            print "------------------------------------\n\n"

        cutfunctions = {}
        #cutfunctions = lambda pulser, is_saturated, n_tracks, fft_mean, INVfft_max: 1
        exec("cutfunctions[0] = {cutf}".format(cutf=self.cut[0].GetCutFunctionDef()))
        exec("cutfunctions[3] = {cutf}".format(cutf=self.cut[3].GetCutFunctionDef()))

        for i in included:
            if i%10000==0:print "processing Event ", i
            # read the ROOT TTree
            self.run.tree.GetEntry(i)
            x_[0] = self.run.tree.diam1_track_x
            y_[0] = self.run.tree.diam1_track_y
            x_[3] = self.run.tree.diam2_track_x
            y_[3] = self.run.tree.diam2_track_y
            for ch in channels:
                signal_ = eval(signaldef.format(signal=self.signalname, channel=ch, pedestal=self.pedestalname))
                pulser = self.run.tree.pulser
                is_saturated = self.run.tree.is_saturated[ch]
                fft_mean = self.run.tree.fft_mean[ch]
                INVfft_max = 1./self.run.tree.fft_max[ch]
                n_tracks = self.run.tree.n_tracks
                sig_time = self.run.tree.sig_time[ch]
                sig_spread = self.run.tree.sig_spread[ch]
                median = self.run.tree.median[ch]
                if cutfunctions[ch](pulser, is_saturated, n_tracks, fft_mean, INVfft_max, sig_time, sig_spread, median): #(not pulser and not is_saturated and fft_mean>50 and fft_mean<500 and INVfft_max>1e-4):
                    self.Pads[ch].Fill(x_[ch], y_[ch], signal_)
        print "Tracking information of ",len(included), " events loaded"


        #self.Pads[channel].MakeFits()
        self.Signal2DDistribution = {}
        for ch in channels:
            self.Signal2DDistribution[ch] = self.Pads[ch].GetMeanSignalDistribution(self.minimum_bincontent)
            self.Signal2DDistribution[ch].SetDirectory(0)
            self.Signal2DDistribution[ch].SetStats(False)
            self.Signal2DDistribution[ch].GetXaxis().SetTitle("pos x / cm")
            self.Signal2DDistribution[ch].GetYaxis().SetTitle("pos y / cm")
            self.Signal2DDistribution[ch].GetYaxis().SetTitleOffset(1.4)

        self.Checklist["LoadTrackData"] = True

    def ShowHitMap(self, channel=None, saveplot = False, drawoption = "colz", RemoveLowStatBins = 0):
        '''
        Shows a 2-dimensional (TH2D) histogram of the hit distributions
        on the pad.
        :param channel:
        :param saveplot:
        :param drawoption:
        :param RemoveLowStatBins:
        :return:
        '''
        if channel != None:
            assert (channel in self.run.GetChannels()), "channel not selected"
            channels = [channel]
        else:
            channels = self.run.GetChannels()

        if not self.Checklist["LoadTrackData"]:
            self.LoadTrackData()

        if RemoveLowStatBins > 0:
            if type(RemoveLowStatBins) == t.BooleanType:
                RemoveLowStatBins = 10
            bins = (self.Pads[channel].counthisto.GetNbinsX() + 2)*(self.Pads[channel].counthisto.GetNbinsY() + 2)
            for bin in xrange(bins):
                if self.Pads[channel].counthisto.GetBinContent(bin) < RemoveLowStatBins:
                    coordinates = self.Pads[channel].GetBinCenter(bin)
                    content = self.Pads[channel].counthisto.GetBinContent(bin)
                    self.Pads[channel].counthisto.Fill(coordinates[0], coordinates[1], -content)
            extension = "_min"+str(RemoveLowStatBins)
        else:
            extension = ""

        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        self.hitmapcanvas = ROOT.TCanvas("hitmapcanvas", "Hits", len(channels)*500, 500) # adjust the width slightly
        self.hitmapcanvas.Divide(len(channels), 1)

        for ch in channels:
            self.hitmapcanvas.cd(channels.index(ch)+1)
            self.Pads[ch].counthisto.SetStats(False)
            self.Pads[ch].counthisto.Draw(drawoption)#"surf2")
            #self.Pads[ch].counthisto.Draw("CONT1 SAME")
            if saveplot:
                self.SavePlots(("Run{run}_HitMap"+extension).format(run=self.run.run_number), "png", canvas=self.hitmapcanvas)
            self.hitmapcanvas.Update()
        self.IfWait("Hits Distribution shown")
        self.Checklist["HitsDistribution"] = True

    def SetDiamondPosition(self, diamonds=3):
        '''
        With this method the diamond window is set. The diamond mask
        info is stored in DiamondPositions/ and will be loaded the next
        time an Analysis object with the same mask file (silicon pixel
        mask) is created.
        The margins has to be set in the terminal window and can be
        found by the HitMaps, which will pop up.
        To cancel the process, just pass a non-number character.
        :param diamonds:
        :return:
        '''
        tmp = copy.deepcopy(self.run.analyzeCh)
        self.SetChannels(diamonds)
        self.ShowHitMap()
        maskname = self.run.RunInfo["mask"][:-4]
        for ch in [0,3]:
            try:
                xmin = float(raw_input("Set x_min of Diamond {dia}: ".format(dia=self.run.diamondname[ch])))
                xmax = float(raw_input("Set x_max of Diamond {dia}: ".format(dia=self.run.diamondname[ch])))
                ymin = float(raw_input("Set y_min of Diamond {dia}: ".format(dia=self.run.diamondname[ch])))
                ymax = float(raw_input("Set y_max of Diamond {dia}: ".format(dia=self.run.diamondname[ch])))
                window = [xmin, xmax, ymin, ymax]
                # save list to file
                windowfile = open("DiamondPositions/{testcampaign}_{mask}_{dia}.pickle".format(testcampaign=self.TESTCAMPAIGN, mask=maskname, dia=ch), "wb")
                pickle.dump(window, windowfile)
                windowfile.close()
            except ValueError:
                pass
        self.run.analyzeCh = copy.deepcopy(tmp)

    def SetChannels(self, diamonds):
        '''
        Sets the channels (i.e. diamonds) to be analyzed. This is just a
        short cut for: self.run.SetChannels(diamonds)
        :param diamonds:
        :return:
        '''
        self.run.SetChannels(diamonds=diamonds)

    def GetRate(self):
        '''

        :return:
        '''
        return self.run.GetRate()

    def ShowSignalMaps(self, draw_minmax=True, saveplots = False, savename = "Run{run}_SignalMaps",ending="png",saveDir = "Results/", show3d = False):
        '''
        Shows a 2-dimensional mean signal response map for each channel
        which is activated for analysis.
        :param saveplots: if True, save the plot
        :param savename: filename if saveplots = True
        :param ending: datatype of file if saveplots = True
        :param saveDir: directory to save the plot - has to end with '/'
        :return: -
        '''
        if not self.Checklist["LoadTrackData"]:
            self.LoadTrackData()

        channels = self.run.GetChannels()
        self.signal_canvas = ROOT.TCanvas("signal_canvas{run}", "Mean Signal Maps", len(channels)*500, 500)
        self.signal_canvas.Divide(len(channels), 1)
        if len(channels) ==2:
            namesuffix = ""
        else:
            namesuffix = "_Ch{ch}".format(ch=channels[0])
        self.signal_canvas.cd(1)
        ROOT.SetOwnership(self.signal_canvas, False)

        ROOT.gStyle.SetPalette(53) # Dark Body Radiator palette
        ROOT.gStyle.SetNumberContours(999)

        for channel in self.run.GetChannels():
            pad = self.signal_canvas.cd(channels.index(channel)+1)
            # Plot the Signal2D TH2D histogram

            if show3d:
                self.Signal2DDistribution[channel].Draw("SPEC dm(2,10) pa(1,1,1) ci(1,1,1) a(15,45,0) s(1,1)")
                savename += "3D"
            else:
                self.Signal2DDistribution[channel].Draw("colz")

            if draw_minmax and not show3d:
                self._DrawMinMax(pad = self.signal_canvas.cd(channels.index(channel)+1), channel=channel)

            if not show3d: self.DrawRunInfo(canvas=pad, infoid="signalmap{run}{ch}".format(run=self.run.run_number, ch=channel))

        self.signal_canvas.Update()
        self.IfWait("2d drawn")
        if saveplots:
            savename = savename.format(run=self.run.run_number)+namesuffix
            self.SavePlots(savename, ending, canvas=self.signal_canvas, subDir=ending)
            self.SavePlots(savename, "root", canvas=self.signal_canvas, subDir="root")

    def _DrawMinMax(self, pad, channel, theseMaximas=None, theseMinimas=None):
        pad.cd()

        if theseMaximas==None:
            if self.extremaResults[channel]['FoundMaxima'] == None: self.FindMaxima(channel=channel)
        else:
            assert(type(theseMaximas) is t.ListType)
        if theseMinimas==None:
            if self.extremaResults[channel]['FoundMinima'] == None: self.FindMinima(channel=channel)
        else:
            assert(type(theseMinimas) is t.ListType)

        if self.run.IsMonteCarlo:
                print "Run is MONTE CARLO"
                if self.run.SignalParameters[0] > 0:
                    height = self.run.SignalParameters[4]
                    self.Pads[channel].GetSignalInRow(height, show=True)
                #self.combined_canvas.cd(1)
                self.Pads[channel].MaximaSearch.real_peaks.SetMarkerColor(ROOT.kBlue)
                self.Pads[channel].MaximaSearch.real_peaks.SetLineColor(ROOT.kBlue)
                self.Pads[channel].MaximaSearch.real_peaks.Draw('SAME P0')
        # self.Pads[channel].MaximaSearch.found_extrema.SetMarkerColor(ROOT.kGreen+2)
        # self.Pads[channel].MaximaSearch.found_extrema.Draw('SAME P0')
        if theseMinimas==None:
            minima = self.extremaResults[channel]['FoundMinima']
        else:
            minima = theseMinimas

        for i in xrange(len(minima)):
            text = ROOT.TText()
            text.SetTextColor(ROOT.kBlue-4)
            text.DrawText(minima[i][0]-0.01, minima[i][1]-0.005, 'low')

        if theseMaximas==None:
            maxima = self.extremaResults[channel]['FoundMaxima']
        else:
            maxima = theseMaximas

        for i in xrange(len(maxima)):
            text = ROOT.TText()
            text.SetTextColor(ROOT.kRed)
            text.DrawText(maxima[i][0]-0.02, maxima[i][1]-0.005, 'high')
        if len(maxima)*len(minima)>0:
            maxbin = self.Pads[channel].GetBinByCoordinates(*(maxima[0]))
            maxbin.FitLandau()
            minbin = self.Pads[channel].GetBinByCoordinates(*(minima[0]))
            minbin.FitLandau()
            print '\nApproximated Signal Amplitude: {0:0.0f}% - (high/low approximation)\n'.format(100.*(maxbin.Fit['MPV']/minbin.Fit['MPV']-1.))

    def CreateMeanSignalHistogram(self, channel, saveplots = False, savename = "MeanSignalDistribution",ending="png",saveDir = "Results/", show = False):
        '''
        Creates a histogram containing all the mean signal responses of
        the two-dimensional mean signal map. (i.e. the mean signal value
        of each bin of the 2D distribution is one entry in the
        histogram)
        :param channel:
        :param saveplots:
        :param savename:
        :param ending:
        :param saveDir:
        :param show:
        :return:
        '''
        if not hasattr(self, "Pads"):
            print "Performing auto analysis"
            self.LoadTrackData(minimum_bincontent=self.minimum_bincontent)
        if saveplots:
            show = True

        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        if show:
            if hasattr(self, "signal_canvas") and bool(self.signal_canvas):
                self.signal_canvas.Clear()
            else:
                self.signal_canvas = ROOT.TCanvas("signal_canvas{run}", "Mean Signal Maps")
                ROOT.SetOwnership(self.signal_canvas, False)
                ROOT.gStyle.SetPalette(53) # Dark Body Radiator palette
                ROOT.gStyle.SetNumberContours(999)
        # print self.Signal2DDistribution
        minimum = self.Signal2DDistribution[channel].GetMinimum()
        maximum = self.Signal2DDistribution[channel].GetMaximum()
        self.MeanSignalHisto_name = "MeanSignalHisto"+str(self.run.run_number)+str(channel)
        if not hasattr(self, "MeanSignalHisto"):
            self.MeanSignalHisto = {}
        self.MeanSignalHisto[channel] = ROOT.TH1D(self.MeanSignalHisto_name, "Run{run}: {diamond} Mean Signal Histogram".format(run=self.run.run_number, diamond=self.run.diamondname[channel]), 100, minimum,maximum)
        ROOT.SetOwnership(self.MeanSignalHisto[channel], False)

        nbins = (self.Signal2DDistribution[channel].GetNbinsX()+2)*(self.Signal2DDistribution[channel].GetNbinsY()+2)
        for i in xrange(nbins):
            bincontent = self.Signal2DDistribution[channel].GetBinContent(i)
            binhits = self.Pads[channel].counthisto.GetBinContent(i)
            if binhits >= self.minimum_bincontent and bincontent != 0: # bincontent != 0: bins in overflow frame give 0 content by ROOT.TH2D
                self.MeanSignalHisto[channel].Fill(bincontent)
        if show:
            self.MeanSignalHisto[channel].Draw()
            self.signal_canvas.Update()
        else:
            print "INFO: MeanSignalHisto created as attribute: self.MeanSignalHisto (ROOT.TH1D)"

        if saveplots:
            self.SavePlots(savename, ending, saveDir, canvas=self.signal_canvas)

        self.Checklist["MeanSignalHisto"][channel] = True

    def ShowExtendedSignalMaps(self, draw_minmax=True, saveplots = False, savename = "SignalDistribution",ending="png",saveDir = None, PS=False, test=""):
        '''
        Shows the 2-dimensional mean signal map, as well as the
        corresponding histogram of the mean signals inside this map.
        (i.e. the mean of each bin corresponds to one entry in the
        histogram)
        :param draw_minmax:
        :param saveplots:
        :param savename:
        :param ending:
        :param saveDir:
        :param PS:
        :param test:
        :return:
        '''
        self.combined_canvas_name = "combined_canvas"+test
        self.combined_canvas = ROOT.gROOT.GetListOfCanvases().FindObject(self.combined_canvas_name)

        channels = self.run.GetChannels()

        if not self.combined_canvas:
            self.combined_canvas = ROOT.TCanvas(self.combined_canvas_name,"Combined Canvas",1000,len(channels)*500)
            ROOT.SetOwnership(self.combined_canvas, False)
            self.combined_canvas.Divide(2,len(channels))


        for channel in channels:
            #self.CreatePlots(False)
            self.CreateMeanSignalHistogram(channel=channel, saveplots=False, show=False)

            self.combined_canvas.cd(1+channels.index(channel)*2) # 2D Signal Map Pad
            #ROOT.gStyle.SetPalette(55)
            # ROOT.gStyle.SetNumberContours(999)
            ROOT.gStyle.SetPalette(53)
            ROOT.gStyle.SetNumberContours(999)
            self.Signal2DDistribution[channel].Draw("colz")#"CONT1Z")#)'colz')
            if draw_minmax:
                self._DrawMinMax(pad = self.combined_canvas.cd(1+channels.index(channel)*2), channel=channel)

            self.combined_canvas.cd(2+channels.index(channel)*2) # Histo Pad

            if PS: #if photoshop mode, fill histogram pink
                ROOT.gStyle.SetHistFillColor(6)
                ROOT.gStyle.SetHistFillStyle(1001)
            else:
                ROOT.gStyle.SetHistFillColor(7)
                ROOT.gStyle.SetHistFillStyle(3003)

            self.MeanSignalHisto[channel].UseCurrentStyle()
            self.MeanSignalHisto[channel].GetXaxis().SetTitle("Signal response")
            self.MeanSignalHisto[channel].Draw()
        self.combined_canvas.cd()
        self.combined_canvas.Update()

        savename = self.run.diamondname[channel]+"_"+savename+"_"+str(self.run.run_number) # diamond_irradiation_savename_runnr
        if saveplots:
            self.SavePlots(savename, ending, saveDir, canvas=self.combined_canvas)
            self.SavePlots(savename, "root", saveDir, canvas=self.combined_canvas)
        if PS:
            ROOT.gStyle.SetHistFillColor(7)
            ROOT.gStyle.SetHistFillStyle(3003)
        self.IfWait("Combined 2D Signal DistributionsShown")

    def FindMaxima(self, channel, show=False):
        '''
        Conducts a maxima search using the class FindMaxima. The Results
        are stored in the dict:
            self.extremaResults[channel]
        :param channel:
        :param show:
        :return:
        '''
        if not hasattr(self, "Pads"):
            self.LoadTrackData()
        self.Pads[channel].FindMaxima(minimum_bincount=self.minimum_bincontent, show=show)

    def FindMinima(self, channel, show=False):
        if not hasattr(self, "Pads"):
            self.LoadTrackData()
        self.Pads[channel].FindMinima(minimum_bincount=self.minimum_bincontent, show=show)

    def GetMPVSigmas(self, channel, show = False):
        '''
        Returns MPVs and sigmas form all Bins of current run as Lists.
        If shown=True then a scatter plot is shown
        :return:
        '''
        if not hasattr(self, "Pads"):
            self.LoadTrackData()

        MPVs = []
        MPVErrs = []
        Sigmas = []
        SigmaErrs = []
        for bin in self.Pads[channel].listOfBins:
            entries = bin.GetEntries()
            if entries >= self.minimum_bincontent:
                if bin.Fit["MPV"] == None:
                    bin.FitLandau()
                MPVs.append(bin.Fit["MPV"])
                MPVErrs.append(bin.Fit["MPVErr"])
                Sigmas.append(bin.Fit["Sigma"])
                SigmaErrs.append(bin.Fit["SigmaErr"])

        if show:
            canvas = ROOT.TCanvas("MPVSigmaCanvas", "MPVSigmaCanvas")
            canvas.cd()
            graph = ROOT.TGraphErrors()
            graph.SetNameTitle("mpvsigmascatterplot", "Landau Fit Parameters of all Bins")
            count = 0
            for i in xrange(len(MPVs)):
                graph.SetPoint(count, MPVs[i], Sigmas[i])
                graph.SetPointError(count, MPVErrs[i], SigmaErrs[i])
                count += 1
            graph.Draw("AP")
            graph.GetYaxis().SetTitle("Sigma of Landau fit")
            graph.GetXaxis().SetTitle("MPV of Landau fit")
            self.DrawRunInfo(channel=channel)
            canvas.Update()
            self.IfWait("MPV vs Sigma shown")

        return MPVs, Sigmas, MPVErrs, SigmaErrs

    def GetSignalHeight(self, channel, min_percent = 5, max_percent = 99):
        '''
        Calculates the spread of mean signal response from the 2D signal response map.
        The spread is taken from min_percent quantile to max_percent quantile.
        Definition: SignalHeight := max_percent_quantile / min_percent_quantile - 1
        :param channel:
        :param min_percent:
        :param max_percent:
        :return:
        '''
        if self.Checklist["FindExtrema"]["Maxima"][channel]:
            self.FindMaxima(channel=channel)

        if not self.Checklist["MeanSignalHisto"][channel]:
            self.CreateMeanSignalHistogram(channel=channel)
        q = array('d', [1.*min_percent/100., 1.*max_percent/100.])
        y = array('d', [0,0])
        self.MeanSignalHisto[channel].GetQuantiles(2, y, q)
        SignalHeight = y[1]/y[0]-1.
        self.VerbosePrint('\nApproximated Signal Amplitude: {0:0.0f}% - ({1:0.0f}%/{2:0.0f}% Quantiles approximation)\n'.format(100.*(SignalHeight), max_percent, min_percent))
        self.extremaResults[channel]['SignalHeight'] = SignalHeight
        return SignalHeight

    def ShowTotalSignalHistogram(self, channel, save = False, scale = False, showfit=False):
        '''
        The same as self.ShowSignalHisto(), but concerning only data
        with tracking information (i.e. the data that is stored in
        self.Pads[channel]).
        In addition,  this method also has the Langaus fit implemented.
        :param channel: 
        :param save: 
        :param scale: 
        :param showfit: produce a 'Langaus' fit
        :return:
        '''
        if  hasattr(self, "Pads"):
            self.Pads[channel].CreateTotalSignalHistogram(saveplot=save, scale=scale, showfit=showfit)
            # if showfit:
            #     self.GetSignalHistoFitResults() # Get the Results from self.signalHistoFitResults to self.signalHistoFitResults
        else:
            self.LoadTrackData()
            self.Pads[channel].CreateTotalSignalHistogram(saveplot=save, scale=scale, showfit=showfit)

    def GetSignalHistoFitResults(self, channel, show=False):
        '''
        Returns fit results of the signal histogram fit ('Langau'). If
        the fits of the signal histogram does not exist, it will create
        it and show the histogram containing the fit.
        :param channel: 
        :return: FitFunction, Peak, FWHM, Chi2, NDF
        '''
        if self.signalHistoFitResults[channel]["Peak"] != None:
            FitFunction = self.signalHistoFitResults[channel]["FitFunction"]
            Peak = self.signalHistoFitResults[channel]["Peak"]
            FWHM = self.signalHistoFitResults[channel]["FWHM"]
            Chi2 = self.signalHistoFitResults[channel]["Chi2"]
            NDF = self.signalHistoFitResults[channel]["NDF"]
            if show: self.ShowTotalSignalHistogram(channel=channel, save=False, showfit=True)
            return FitFunction, Peak, FWHM, Chi2, NDF
        else:
            self.ShowTotalSignalHistogram(channel=channel, save=False, showfit=True)
            if (self.signalHistoFitResults[channel]["Peak"] != None) or (hasattr(self, "ExtremeAnalysis") and self.signalHistoFitResults[channel]["Peak"] != None):
                self.GetSignalHistoFitResults()
            else:
                assert(False), "BAD SignalHistogram Fit, Stop program due to possible infinity loop"
    
    def SignalTimeEvolution(self, channel, Mode="Mean", show=True, time_spacing = 3, save = True, binnumber = None, RateTimeEvolution=False, nameExtension=None): # CUTS!! not show: save evolution data / Comment more
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
        False # MAKE FASTER USING self.GetEventAtTime()
        # assert(Mode in ["Mean", "MPV"]), "Wrong Mode, Mode has to be `Mean` or `MPV`"
        # if not hasattr(self, "Pads"):
        #     self.LoadTrackData()
        # if True: #not self.run.IsMonteCarlo:
        #
        #     # Names
        #     if binnumber != None:
        #         binCenter = self.Pads[channel].GetBinCenter(binnumber)
        #         if nameExtension is None:
        #             nameExtension = "_Bin{0:.3f}_{1:.3f}".format(*binCenter)
        #     else:
        #         binCenter = 0,0
        #         if nameExtension is None:
        #             nameExtension = "OverAll"
        #     if RateTimeEvolution:
        #         type_ = "Rate"
        #     else:
        #         type_ = "Signal"
        #
        #     results = {}
        #     self.SignalEvolution = ROOT.TGraphErrors()
        #     self.SignalEvolution.SetNameTitle(type_+"Evolution", type_+" Time Evolution "+nameExtension)
        #
        #     self.run.tree.GetEntry(1)
        #     starttime = self.run.tree.time
        #     self.run.tree.GetEntry(self.run.tree.GetEntries())
        #     # endtime = self.run.tree.time
        #     ok_deltatime = 0
        #     TimeERROR = False
        #     for i in xrange(self.run.tree.GetEntries()-1):
        #         i += 1
        #         # read the ROOT TTree
        #         self.run.tree.GetEntry(i)
        #         signal_ = self.run.tree.sig_spread[0]
        #         time_ = self.run.tree.time
        #         deltatime_min = (time_-starttime)/60000.
        #         if deltatime_min < -1: # possible time_stamp reset
        #             deltatime_min = ok_deltatime + time_/60000.
        #             TimeERROR = True
        #         else:
        #             ok_deltatime = deltatime_min
        #         time_bucket = int(deltatime_min)/int(time_spacing)*int(time_spacing)
        #         pulser = self.run.tree.pulser
        #
        #         # check if hit is inside bin (ONLY FOR BIN SIGNAL TIME EVOLUTION)
        #         if binnumber != None:
        #             x_ = self.run.tree.diam1_track_x
        #             y_ = self.run.tree.diam1_track_y
        #             binnumber_ = self.Pads[channel].GetBinNumber(x_, y_)
        #         else:
        #             binnumber_ = None
        #
        #         if not pulser and (binnumber == None or binnumber == binnumber_ ):
        #             try:
        #                 results[time_bucket].append(signal_)
        #             except KeyError:
        #                 results[time_bucket] = [signal_]
        #     if TimeERROR:
        #         print "\n\n\n\n\nWARNING: Error occured in time flow of Run "+str(self.run.run_number)+". Possible reset of timestamp during data taking!\n\n\n\n\n"
        #
        #     time_buckets = results.keys()
        #     time_buckets.sort()
        #
        #     # drop last bucket if Rate Time Scan (last bucket may be overshooting the data time window)
        #     if RateTimeEvolution:
        #         time_buckets = time_buckets[:-1]
        #
        #     count = 0
        #     _, xmin, xmax, _, ymin, ymax =  self.Pads[channel].Get2DAttributes()
        #     area = (xmax-xmin)*(ymax-ymin)
        #     for t in time_buckets:
        #         histo = ROOT.TH1D(type_+"Evolution_time_"+str(t), type_+" Time Histogram for t = {0:0.0f}-{1:0.0f}min".format(t, t+int(time_spacing)), 500, 0, 500)
        #         for i in xrange(len(results[t])):
        #             histo.Fill(results[t][i])
        #         if Mode == "Mean" or RateTimeEvolution:
        #             if RateTimeEvolution:
        #                 if binnumber != None:
        #                     c = 1./(60.*time_spacing*(self.config_object[channel].config["2DHist"]["binsize"])**2)
        #                     N = histo.GetEntries()
        #                     self.SignalEvolution.SetPoint(count, t, c*N) # Rate/cm^2 = Entries/(seconds*(binsize)**2)
        #                 else:
        #                     c = 1./(60.*time_spacing*area)
        #                     N = histo.GetEntries()
        #                     self.SignalEvolution.SetPoint(count, t, c*N) # Rate/cm^2 = Entries/(seconds*Area)
        #                 self.SignalEvolution.SetPointError(count,0,c*np.sqrt(N))
        #                 signalname = "Rate / Hz/cm^2"
        #             else:
        #                 self.SignalEvolution.SetPoint(count, t, histo.GetMean())
        #                 self.SignalEvolution.SetPointError(count, 0, histo.GetRMS()/np.sqrt(histo.GetEntries()))
        #                 signalname = "Mean of Signal Response"
        #         elif Mode == "MPV":
        #             self.SignalEvolution.SetPoint(count, t, histo.GetBinCenter(histo.GetMaximumBin()))
        #             self.SignalEvolution.SetPointError(count, 0, 0)
        #             signalname = "MPV of Signal Response"
        #         count += 1
        #         del histo
        #         ROOT.gROOT.Delete(type_+"Evolution_time_"+str(t))
        #
        #     self.SignalEvolution.GetXaxis().SetTitle("Time / min")
        #     self.SignalEvolution.GetYaxis().SetTitle(signalname)
        #     self.SignalEvolution.GetYaxis().SetTitleOffset(1.1)
        #     if RateTimeEvolution: # Start y axis from 0 when rate time evolution
        #         if binnumber != None: # just one bin
        #             self.SignalEvolution.GetYaxis().SetRangeUser(0, 4000)#self.SignalEvolution.GetYaxis().GetXmax())
        #         else: #overall
        #             self.SignalEvolution.GetYaxis().SetRangeUser(0, 1700)#self.SignalEvolution.GetYaxis().GetXmax())
        #     canvas = ROOT.TCanvas("signaltimeevolution", "signaltimeevolution")
        #
        #     self.SignalEvolution.Draw("ALP*")
        #     # Draw a second axis if it is a RateTimeEvolution:
        #     if RateTimeEvolution:
        #         canvas.Update()
        #         rightaxis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUymin()/c, ROOT.gPad.GetUymax()/c, 20210, '+L')
        #         rightaxis.SetTitle('# Hits')
        #         rightaxis.SetTitleOffset(1.2)
        #         ROOT.SetOwnership(rightaxis, False)
        #         rightaxis.Draw('SAME')
        #     if save:
        #         self.SavePlots(type_+"TimeEvolution"+Mode+nameExtension+".png", canvas=canvas)
        #     if TimeERROR:
        #         self.SignalEvolution.SaveAs(self.SaveDirectory+"ERROR_"+type_+"TimeEvolution"+Mode+nameExtension+".root")
        #     self.IfWait("Showing "+type_+" Time Evolution..")
        #     canvas.Close()
        # else:
        #     print "Run is Monte Carlo. Signal- and Rate Time Evolution cannot be created."

    def GetEventAtTime(self, dt):
        '''
        Returns the eventnunmber at time dt from beginning of the run.

        Accuracy: +- 2 Events

        The event number is evaluated using a newton's method for
        finding roots, i.e.
            f(e) := t(e) - t  -->  f(e) == 0

            ==> iteration: e = e - f(e)/f'(e)

            where t(e) is the time evaluated at event e and
            t := t_0 + dt
            break if |e_old - e_new| < 2
        :param time: time in seconds from start
        :return: event_number
        '''
        time = dt*1000 # convert to milliseconds

        if time == 0: return 0

        #get t0 and tmax
        maxevent = self.run.tree.GetEntries()
        if time < 0: return maxevent
        t_0 = self.GetTimeAtEvent(0)
        t_max = self.GetTimeAtEvent(maxevent-1)

        time = t_0 + time
        if time>t_max: return maxevent

        seedEvent = int( (1.*(time - t_0) * maxevent) / (t_max - t_0) )

        def slope_f(self, event):
            if event < 0: event = 0
            time_high = self.GetTimeAtEvent(event+10)
            time_low = self.GetTimeAtEvent(event-10)
            return 1.*(time_high-time_low)/21.

        count = 0
        goOn = True
        event = seedEvent
        while goOn and count<20:
            old_event = event
            f_event = self.GetTimeAtEvent(event) - time
            #print "f_event = {ftime} - {time} = ".format(ftime=self.run.tree.time, time=time), f_event
            #print "slope_f(self, event) = ", slope_f(self, event)
            event = int(event - 1*f_event/slope_f(self, event))
            if abs(event-old_event)<2:
                goOn = False
            count += 1
        self.run.tree.GetEntry(event)
        return event

    def _ShowPreAnalysisOverview(self, channel = None, savePlot=True):

        if channel == None:
            channels = self.run.GetChannels()
        else:
            channels = [channel]

        for ch in channels:
            self.pAOverviewCanv = ROOT.TCanvas("PAOverviewCanvas", "PAOverviewCanvas", 1500, 900)
            self.pAOverviewCanv.Divide(2,1)

            PApad = self.pAOverviewCanv.cd(1)
            rightPad = self.pAOverviewCanv.cd(2)
            rightPad.Divide(1,3)

            PApad.cd()
            self.MakePreAnalysis(channel=ch, savePlot=False, canvas=PApad)

            upperpad = rightPad.cd(1)
            upperpad.Divide(2,1)
            pulserPad = upperpad.cd(1)
            self.ShowPulserRate(canvas=pulserPad, savePlot=False)
            fftPad = upperpad.cd(2)
            self.ShowFFT("mc", cut=True, channel=ch, savePlots=False, canvas=fftPad)

            self.ResetColorPalette()
            middlepad = rightPad.cd(2)
            middlepad.Divide(2,1)
            spreadPad = middlepad.cd(1)
            self._ShowHisto(signaldef="sig_spread[{channel}]", channel=ch, canvas=spreadPad, infoid="Spread", drawruninfo=True, savePlots=False, logy=True, gridx=True, binning=150, xmin=0, xmax=150)
            medianPad = middlepad.cd(2)
            self.ShowMedianHisto(channel=ch, canvas=medianPad)

            peakPosPad = rightPad.cd(3)
            self.ShowPeakPosition(channel=ch, canvas=peakPosPad, savePlot=False)

            if savePlot:
                self.SavePlots(savename="Run{run}_PreAnalysisOverview_{dia}.png".format(run=self.run.run_number, dia=self.run.diamondname[ch]), subDir=self.run.diamondname[ch], canvas=self.pAOverviewCanv)

    def SetIndividualCuts(self, showOverview=True, savePlot=False):
        '''
        With this method the run-specific cuts can be configured.
        The cuts can be set in the terminal window, where a canvas pops-
        up in order to decide on the cut configurations.
        The individual cut configuration is stored in
            Configuration/Individual_Configs/
        and will be applied the next time an Analysis object with the
        same run number is created.
        In order to ignore the individual cuts, the file has to be
        removed from the directory.
        :param showOverview:
        :param savePlot:
        :return:
        '''
        default_dict_ = { #
            "EventRange":           None,             # [1234, 123456]
            "ExcludeFirst":         None,              # 50000 events
            "noPulser":             None,              # 1: nopulser, 0: pulser, -1: no cut
            "notSaturated":         None,
            "noBeamInter":          None,
            "FFT":                  None,
            "Tracks":               None,
            "peakPos_high":         None,
            "spread_low":           None,
            "absMedian_high":       None
        }
        self.individualCuts = {
            0: copy.deepcopy(default_dict_),
            3: copy.deepcopy(default_dict_)
        }

        try:
            for ch in [0,3]:
                if showOverview:
                    self._ShowPreAnalysisOverview(channel=ch, savePlot=savePlot)
                range_min = raw_input("{dia} - Event Range Cut. Enter LOWER Event Number: ".format(dia=self.run.diamondname[ch]))
                range_max = raw_input("{dia} - Event Range Cut. Enter UPPER Event Number: ".format(dia=self.run.diamondname[ch]))
                peakPos_high = raw_input("{dia} - Peak Position Cut. Enter maximum Peak Position Sample Point: ".format(dia=self.run.diamondname[ch]))
                spread_low = raw_input("{dia} - Spread Cut. Enter minimum Spread (max-min): ".format(dia=self.run.diamondname[ch]))
                absMedian_high = raw_input("{dia} - Median Cut. Enter maximum abs(median) value: ".format(dia=self.run.diamondname[ch]))

                if range_max != "" and range_min != "":
                    self.individualCuts[ch]["EventRange"] = [int(range_min), int(range_max)]
                elif range_max != "":
                    self.individualCuts[ch]["EventRange"] = [0, int(range_max)]
                elif range_min != "":
                    self.individualCuts[ch]["ExcludeFirst"] = int(range_min)

                if peakPos_high != "":
                    self.individualCuts[ch]["peakPos_high"] = int(peakPos_high)

                if spread_low != "":
                    self.individualCuts[ch]["spread_low"] = int(spread_low)

                if absMedian_high != "":
                    self.individualCuts[ch]["absMedian_high"] = int(absMedian_high)
        except:
            print "Aborted."
        else:
            print "Creating run-specific config file:"
            path = "Configuration/Individual_Configs/"
            filename = "{testcp}_Run{run}.json".format(testcp=self.TESTCAMPAIGN, run=self.run.run_number)
            print "\t"+path+filename

            f = open(path+filename, "w")
            json.dump(self.individualCuts, f, indent=2, sort_keys=True)
            f.close()
            print "done."

    def GetTimeAtEvent(self, event):
        '''
        Returns the time stamp at event number 'event'. For negative
        event numbers it will return the time stamp at the startevent.
        :param event: integer event number
        :return: timestamp for event
        '''
        maxevent = self.run.tree.GetEntries()
        if event < 0: event = 0
        if event >= maxevent: event = maxevent - 1
        self.run.tree.GetEntry(event)
        return self.run.tree.time

    def _checkWFChannels(self):
        nWFChannels = 0
        wf_exist = {
            0: False,
            1: False,
            2: False,
            3: False
        }
        check = self.run.tree.GetBranch("wf0")
        if check:
            nWFChannels += 1
            wf_exist[0] = True
        check = self.run.tree.GetBranch("wf1")
        if check:
            nWFChannels += 1
            wf_exist[1] = True
        check = self.run.tree.GetBranch("wf2")
        if check:
            nWFChannels += 1
            wf_exist[2] = True
        check = self.run.tree.GetBranch("wf3")
        if check:
            nWFChannels += 1
            wf_exist[3] = True
        return nWFChannels, wf_exist

    def ShowWaveForms(self, nevents=1000, cut="", startevent=None, channels=None, canvas=None, infoid="", drawoption=""):
        '''
        Draws stacked waveforms for the channels activated to be
        analyzed. Of course, this requires the all waveforms in the
        ROOT file, which due to the possible large sizes often is
        skipped.
        :param nevents:
        :param cut:
        :param startevent:
        :param channels:
        :param canvas:
        :param infoid:
        :param drawoption:
        :return:
        '''
        if startevent != None: assert(startevent>=0), "startevent as to be >= 0"
        maxevent = self.run.tree.GetEntries()
        if startevent > maxevent: return False
        # check number of wf in root file:
        nWFChannels, draw_waveforms = self._checkWFChannels()
        if nWFChannels == 0: return False

        # if userchannels:
        if channels != None:
            nWFChannels = 0
            for i in xrange(4):
                if self._GetBit(channels, i) and draw_waveforms[i]:
                    nWFChannels += 1
                else:
                    draw_waveforms[i] = False

        # check if external canvas:
        if canvas == None:
            self.waveFormCanvas = ROOT.TCanvas("waveformcanvas{run}".format(run=self.run.run_number),"WaveFormCanvas", 750, nWFChannels*300)
        else:
            self.waveFormCanvas = canvas
            self.waveFormCanvas.Clear()

        self.waveFormCanvas.Divide(1, nWFChannels)
        self.waveformplots = {}

        def drawWF(self, channel, events=1000, startevent=50000, cut=""):
            histoname = "R{run}wf{wfch}{infoID}".format(wfch=channel, run=self.run.run_number, infoID=infoid)
            # self.waveformplots[histoname] = ROOT.TH2D(histoname, self.run.GetChannelName(channel)+" {"+cut.format(channel=channel)+"}", 1024, 0, 1023, 1000, -500, 500)
            # ROOT.SetOwnership(self.waveformplots[histoname], False)
            # self.waveformplots[histoname].SetStats(0)
            print "DRAW: wf{wfch}:Iteration$>>{histoname}".format(histoname=histoname, wfch=channel)
            print "cut: ", cut.format(channel=channel), " events: ", events, " startevent: ", startevent
            n = self.run.tree.Draw("wf{wfch}:Iteration$>>{histoname}(1024, 0, 1023, 1000, -500, 500)".format(histoname=histoname, wfch=channel), cut.format(channel=channel), drawoption, events, startevent)
            self.waveformplots[histoname] = ROOT.gROOT.FindObject(histoname)
            ROOT.SetOwnership(self.waveformplots[histoname], False)
            self.waveformplots[histoname].SetStats(0)
            if cut == "":
                self.DrawRunInfo(channel=channel, comment="{nwf} Wave Forms".format(nwf=n/1024), infoid=("wf{wf}"+infoid).format(wf=channel), userWidth=0.15, userHeight=0.15)
            else:
                self.DrawRunInfo(channel=channel, comment="{nwf}/{totnwf} Wave Forms".format(nwf=n/1024, totnwf=events), infoid=("wf{wf}"+infoid).format(wf=channel), userWidth=0.18, userHeight=0.15)
            if n<=0: print "No event to draw in range. Change cut settings or increase nevents"

        start = int(self.run.tree.GetEntries()/2) if startevent==None else int(startevent)
        index = 1
        for i in xrange(4):
            self.waveFormCanvas.cd(index)
            if draw_waveforms[i]:
                drawWF(self, channel=i, events=nevents, startevent=start, cut=cut)
                self.waveFormCanvas.Update()
                index += 1
            else:
                print "Wave Form of channel ", i, " not in root file"

    def ShowWaveFormsPulser(self, nevents=1000, startevent=None, channels=None):
        nWFChannels, draw_waveforms = self._checkWFChannels()
        self.pulserWaveformCanvas = ROOT.TCanvas("pulserWaveformCanvas", "Pulser Waveform Canvas", 1500, nWFChannels*300)
        self.pulserWaveformCanvas.cd()
        self.pulserWaveformCanvas.Divide(2, 1)
        pad1 = self.pulserWaveformCanvas.cd(1)
        pad2 = self.pulserWaveformCanvas.cd(2)
        self.ShowWaveForms(nevents=nevents, cut="!pulser", startevent=startevent, channels=channels, canvas=pad1, infoid="notpulser")
        # self.pulserWaveformCanvas.cd(2)
        # pad2 = self.pulserWaveformCanvas.cd(2)
        # pad2.Clear()
        self.ShowWaveForms(nevents=nevents, cut="pulser", startevent=startevent, channels=channels, canvas=pad2, infoid="pulser")
        raw_input("asdf")

    def GetUserCutString(self, channel=None):
        '''
        Returns a short, more user-friendly cut string, which can be
        used to display the cut configuration as terminal prompt or to
        mention in canvases.
        :param channel:
        :return:
        '''
        if channel == None:
            opt0 = self.cut[0].GetUserCutString()
            opt3 = self.cut[3].GetUserCutString()
            assert(opt0==opt3), "GetUserCutString Not the same for both cut channels. Choose a specific channel."
            return opt0
        else:
            assert(channel in [0,3])
            return self.cut[channel].GetUserCutString()

    def ShowPeakPosition(self, channel=None, cut="", canvas=None, savePlot=True):
        '''
        Generates a plot, which shows the locations of the peak
        positions in the waveforms. In the 2-dimensional plot the pulse-
        height vs sample point of the peak is drawn.
        :param channel:
        :param cut:
        :param canvas:
        :return:
        '''
        if channel == None:
            channels = self.run.GetChannels()
            namesuffix = ""
        else:
            channels = [channel]
            namesuffix = "_ch{ch}".format(ch=channel)

        min_ = self.run.signalregion_low - 10
        max_ = self.run.signalregion_high +10
        binning = max_-min_

        if canvas == None:
            canvas = ROOT.TCanvas("{run}signalpositioncanvas".format(run=self.run.run_number),"{run}signalpositioncanvas".format(run=self.run.run_number), 800, len(channels)*300)
        canvas.Divide(1, len(channels))

        ROOT.gStyle.SetPalette(55) # rainbow palette
        ROOT.gStyle.SetNumberContours(200)

        for ch in channels:
            pad = canvas.cd(channels.index(ch)+1)
            if cut == "":
                thiscut = self.GetCut(ch)
            else:
                thiscut = cut.format(channel=ch)

            self.Draw(("("+self.signaldefinition+"):sig_time[{channel}]>>signalposition{run}{channel}({bins}, {low}, {high}, 600, -100, 500)").format(channel=ch, run=self.run.run_number, bins=binning, low=min_, high=max_), thiscut, "colz")
            hist = ROOT.gROOT.FindObject("signalposition{run}{channel}".format(channel=ch, run=self.run.run_number))
            pad.SetLogz()
            pad.SetGridx()
            if hist:
                hist.SetStats(0)
                hist.SetTitle("Peak Position {"+self.cut[ch].GetUserCutString()+"}")
                hist.GetXaxis().SetTitle("Sample Point of Peak")
                hist.GetXaxis().SetTitleSize(0.05)
                hist.GetXaxis().SetLabelSize(0.05)
                hist.GetXaxis().SetNdivisions(20)
                hist.GetYaxis().SetTitle("Pulse Height ({sigdef})".format(sigdef=self.signaldefinition.format(channel=ch)))
                hist.GetYaxis().SetTitleSize(0.05)
                hist.GetYaxis().SetLabelSize(0.05)
            self.DrawRunInfo(channel=ch, canvas=pad, infoid="peakpos{run}{ch}".format(run=self.run.run_number, ch=ch), userWidth=0.15, userHeight=0.15)

        canvas.Update()
        if savePlot:
            self.SavePlots("Run{run}_PeakPosition{ns}.png".format(run=self.run.run_number, ns=namesuffix), canvas=canvas, subDir="PeakPosition")
        self.IfWait("Peak Position shown")

    def ShowSignalSpread(self, channel=None, cut=""):
        if channel == None:
            channels = self.run.GetChannels()
            namesuffix = ""
        else:
            channels = [channel]
            namesuffix = "_ch{ch}".format(ch=channel)

        canvas = ROOT.TCanvas("{run}signalspreadcanvas".format(run=self.run.run_number),"{run}signalspreadcanvas".format(run=self.run.run_number), 800, len(channels)*300)
        canvas.Divide(1, len(channels))

        for ch in channels:
            canvas.cd(channels.index(ch)+1)
            if cut == "":
                thiscut = self.GetCut(ch)
            else:
                thiscut = cut.format(channel=ch)

            self.Draw(("sig_spread[{channel}]>>signalspread{run}{channel}(400, 0, 400)").format(channel=ch, run=self.run.run_number), thiscut)
            hist = ROOT.gROOT.FindObject("signalspread{run}{channel}".format(channel=ch, run=self.run.run_number))
            if hist: hist.SetStats(0)

        canvas.Update()
        self.SavePlots("Run{run}_SignalSpread{ns}.png".format(run=self.run.run_number, ns=namesuffix), canvas=canvas, subDir="Cuts")
        self.IfWait("Peak Position shown")

    def Draw(self, varexp, selection="", drawoption="", nentries=1000000000, firstentry=0):
        '''
        The ROOT TTree Draw method as a shortcut for
        self.run.tree.Draw(varexp, selection, drawoption, nentries,
        firstentry)
        :param varexp:
        :param selection:
        :param drawoption:
        :param nentries:
        :param firstentry:
        :return:
        '''
        self.run.tree.Draw(varexp, selection, drawoption, nentries, firstentry)

    def CalculateSNR(self, channel=0, signaldefinition="sig_integral2[{channel}]", pedestalname="ped_integral2", logy=False, cut="", name="SNR"):

        self.snr_canvas = ROOT.TCanvas("SNR_canvas", "SNR Canvas")

        self.ResetColorPalette()
        self._ShowHisto(pedestalname+"[{channel}]", channel=channel, canvas=self.snr_canvas, drawoption="", color=None, cut=cut, normalized=False, infoid="SNRPedestalHisto", drawruninfo=False, binning=500, xmin=-100, xmax=400,savePlots=False, logy=logy)
        self._ShowHisto(signaldefinition, channel=channel, canvas=self.snr_canvas, drawoption="sames", color=None, cut=cut, normalized=False, infoid="SNRSignalHisto", drawruninfo=False, binning=500, xmin=-100, xmax=400,savePlots=False, logy=logy)

        pedestalhisto = ROOT.gROOT.FindObject("{dia}_SNRPedestalHisto{run}".format(dia=self.run.diamondname[channel], run=self.run.run_number))
        signalhisto = ROOT.gROOT.FindObject("II6-97_SNRSignalHisto445")

        pedestalhisto.Fit("gaus", "","",-20,20)
        fitfunc = pedestalhisto.GetFunction("gaus")
        self.pedestalFitMean = fitfunc.GetParameter(1)
        pedestalSigma = fitfunc.GetParameter(2)

        signalmean = signalhisto.GetMean()

        SNR = signalmean/pedestalSigma
        print "SNR = ", SNR

        self.DrawRunInfo(channel=channel, canvas=self.snr_canvas, comment="SNR: "+str(SNR))

        self.SavePlots(savename=name+".pdf", saveDir="SNR", canvas=self.snr_canvas)

        return  SNR

    def MakeSNRAnalyis(self, channel=0):


        SNRs = {}

        windownames = ["a", "b", "c"]
        integralnames = ["1", "2", "3"]
        signaldefs = {
            "a1": "sig_a1[{channel}]",
            "a2": "sig_a2[{channel}]",
            "a3": "sig_a3[{channel}]",
            "b1": "sig_b1[{channel}]",
            "b2": "sig_b2[{channel}]",
            "b3": "sig_b3[{channel}]",
            "c1": "sig_integral1[{channel}]",
            "c2": "sig_integral2[{channel}]",
            "c3": "sig_integral3[{channel}]",
        }

        for windowname in windownames:

            for integralname in integralnames:

                SNRs[windowname+integralname] = self.CalculateSNR(signaldefinition=signaldefs[windowname+integralname], pedestalname="ped_integral"+integralname, name="SNR_"+windowname+integralname)

                if windowname=="c":
                    cut = self.GetCut(channel=channel)
                    SNRs[windowname+integralname+"b"] = self.CalculateSNR(signaldefinition=signaldefs[windowname+integralname], pedestalname="ped_integral"+integralname, cut=cut+"&&sig_peak[{channel}]<250", name="SNR_"+windowname+integralname+"b")

        print SNRs

        for key in SNRs.keys():
            print key, " - ", SNRs[key]