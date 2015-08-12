import ROOT
from ROOT import gROOT
from RunClass import Run
import numpy as np
from AbstractClasses.Elementary import Elementary
from AbstractClasses.PreAnalysisPlot import PreAnalysisPlot
from AbstractClasses.ConfigClass import *
from AbstractClasses.RunClass import Run
from AbstractClasses.Langau import Langau
import os
import copy, collections, numpy, array
import sys
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

    def __init__(self, run, config_object = Config(), verbose = False):
        '''
        Initializes the Analysis object.
        :param run: run object of type "Run"
        :param config_object: config object of type "Pad2DHistConfig"
        :return: -
        '''
        Elementary.__init__(self, verbose=verbose)
        if not isinstance(run, Run):
            assert (type(run) is t.IntType), "run has to be either an instance of Run or run number (int)"
            run = Run(run)
        else:
            assert(run.run_number > 0), "No run selected, choose run.SetRun(run_nr) before you pass the run object"
        self.run = run
        self.config_object = config_object
        self.RunInfo = copy.deepcopy(run.RunInfo)

        self.Checklist = { # True if Plot was created # -> USED?!
            "RemoveBeamInterruptions": False,
            "LoadTrackData": False,
            "MeanSignalHisto": False,
            "HitsDistribution": False,
        }
#
        #self.ExtremaResults ={
        #    "TrueNPeaks": None, # MC true number of peaks
        #    "FoundNMaxima": None, # number of Maxima found
        #    "FoundMaxima": None, # Maxima found as list [(x1, y1), (x2, y2), ...]
        #    "FoundNMinima": None, # number of Minima found
        #    "FoundMinima": None, # Minima found as list [(x1, y1), (x2, y2), ...]
        #    "Ninjas": None, # number of true peaks not found (only available when MC Run)
        #    "Ghosts": None, # nunmber of wrong Maxima found (only available when MC Run)
        #    "SignalHeight": 0.
        #}
#
        #self.SignalHistoFitResults = {
        #    "FitFunction": None,
        #    "Peak": None,
        #    "FWHM": None,
        #    "Chi2": None,
        #    "NDF": None,
#
        #}

    def LoadConfig(self):
        configfile = "Configuration/AnalysisConfig.cfg"
        parser = ConfigParser.ConfigParser()
        output = parser.read(configfile)
        print " ---- Config Parser Read ---- \n - ", output, " -\n"

        self.ShowAndWait = parser.getboolean("DISPLAY","ShowAndWait")
        self.SaveMCData = parser.getboolean("SAVE","SaveMCData")
 #       self.check_offset = parser.getboolean("DO-ANALYSIS","check_offset")
        self.signalname = parser.get("BASIC", "signalname")
        self.pedestal_correction = parser.getboolean("BASIC", "pedestal_correction")
        self.pedestalname = parser.get("BASIC", "pedestalname")
        self.cut = parser.get("BASIC", "cut")
        self.excludefirst = parser.getint("BASIC", "excludefirst")

        if not self.pedestal_correction:
            self.signaldefinition = self.signalname+"[{channel}]"
        else:
            self.signaldefinition = self.signalname+"[{channel}]-"+self.pedestalname+"[{channel}]"


    def MakePreAnalysis(self):
        '''
        Creates the Signal Time distribution (and checks for beam interruptions).
        No tracking information needed -> no spatial information provided
        :return:
        '''
        #c1 = ROOT.TCanvas("c1", "c1")
        self.preAnalysis = {}
        for ch in self.run.GetChannels():
            self.preAnalysis[ch] = PreAnalysisPlot(self, ch)
            self.preAnalysis[ch].Draw()

    def FindBeamInterruptions(self):
        '''
        Finds the beam interruptions
        :return: list of event numbers where beam interruptions occures
        '''
        print "Searching for beam interruptions.."
        nentries = self.run.tree.GetEntries()
    #    last_entry = self.run.tree.GetEntry(nentries-1)
    #    max_time = self.run.tree.time
        
        self.run.tree.Draw('time:event_number')
        
    #    graph = copy.deepcopy(ROOT.c1.FindObject('Graph'))
        histo = copy.deepcopy(ROOT.c1.FindObject('htemp'))
        
        histo.SetTitle('run %3d' %(self.run.run_number))
        histo.SetName('run %3d' %(self.run.run_number))

        # get event numbers and dt's
        dts = []
        evs = []
        i = self.excludefirst
        step = 100
        while i+step < nentries:
            self.run.tree.GetEntry(i)
            t1 = self.run.tree.time
            evs.append(self.run.tree.event_number)
            self.run.tree.GetEntry(i+step)
            t2 = self.run.tree.time
            dt = (t2 - t1)
            dts.append(dt)
            i += step
        
        self.jumps = []
        
        deq = collections.deque(dts[:100],100)
        first = True
        for i in dts[101:]:
            avg = numpy.mean(deq)
            if abs(i / avg - 1.) > 0.5: 
                if first:
                    print 'found a jump here', i, 'at event number', evs[dts.index(i)]
                    self.jumps.append(evs[dts.index(i)])
                    first = False
            else:
                if not first:
                    print 'back to normal at event', evs[dts.index(i)]
                deq.appendleft(i)
                first = True
        
        print '\n'
        print 'found %d jumps' %(len(self.jumps))
        print 'they are at event numbers', self.jumps
        
        lat = ROOT.TLatex()
        lat.SetNDC()
        lat.SetTextColor(ROOT.kRed)
        lat.DrawLatex(0.2,0.85, 'run %d' %(self.run.run_number) )
            
        if not os.path.exists("beaminterruptions"):
            os.mkdir("beaminterruptions")
        if not os.path.exists("beaminterruptions/plots"):
            os.mkdir("beaminterruptions/plots")
        if not os.path.exists("beaminterruptions/data"):
            os.mkdir("beaminterruptions/data")

        # save jump list to file
        jumpfile = open("beaminterruptions/data/Run_{run}.pickle".format(run=self.run.run_number), "wb")
        pickle.dump(self.jumps, jumpfile)
        jumpfile.close()

        if len(self.jumps):
            print 'the length of jumps is', len(self.jumps)
            jumps_array = array('d', self.jumps)
            jumps_err_array = array('d', len(self.jumps)*[histo.GetYaxis().GetXmin()])
            
            jumps_graph = ROOT.TGraph(len(self.jumps), jumps_array, jumps_err_array)
            jumps_graph.SetMarkerSize(3)
            jumps_graph.SetMarkerColor(ROOT.kRed)
            jumps_graph.SetMarkerStyle(33)
            jumps_graph.SetLineColor(ROOT.kRed)
            jumps_graph.Draw('p')
        
            outfile = open('beaminterruptions/jumps.txt','r+a')
            # check if the run is already in the file
            runInFile = False
            lines = outfile.readlines()
            for i in lines:
                if len(i.split()) > 0 and i.split()[0] == str(self.run.run_number):
                    runInFile = True
            if not runInFile:
                outfile.write(str(self.run.run_number)+'\t\t')
        
            lat.SetTextColor(ROOT.kBlack)
            for i in self.jumps:
                ind = self.jumps.index(i)
                lat.DrawLatex(0.2, 0.80-ind*0.05, '#%d at %d' %(ind, i) )
                if not runInFile:
                    outfile.write(str(i)+'\t')
            if not runInFile:
                outfile.write('\n')
            outfile.close()
                
        
        ROOT.c1.SaveAs('beaminterruptions/plots/jumpSearch_run%d.png' %(self.run.run_number))

        return self.jumps
        
    def GetBeamInterruptions(self):
        '''
        If there is beam interruption data, it will load them - otherwise it will run the beam interruption analysis
        it will create the attribute self.jumps, which is a list of event numbers, where a jump occures
        :return: list of events where beam interruptions occures
        '''
        if not hasattr(self, "jumps"):
            picklepath = "beaminterruptions/data/Run_{run}.pickle".format(run=self.run.run_number)
            if os.path.exists(picklepath):
                print "Loading beam interruption data from pickle file: \n\t"+picklepath
                jumpfile = open(picklepath, "rb")
                self.jumps = pickle.load(jumpfile)
                jumpfile.close()
            else:
                self.FindBeamInterruptions()

        return self.jumps

    def RemoveBeamInterruptions(self):
        '''
        This adds the restrictions to the cut string such that beam interruptions are excluded each time the
        cut is applied.
        :return:
        '''
        if not self.Checklist["RemoveBeamInterruptions"]:
            self.GetBeamInterruptions()
            # events cut away in [jump-2000, jump+10000] :
            highcut = 20000
            lowcut = 5000

            njumps = len(self.jumps)
            for i in xrange(njumps):
                self.cut += "&&!(event_number<={jump}+{high}&&event_number>={jump}-{low})".format(jump=self.jumps[i], high=highcut, low=lowcut)
            self.Checklist["RemoveBeamInterruptions"] = True

        return self.cut

    def GetIncludedEvents(self, maxevent=2000000):
        if maxevent>=2000000:
            maxevent = self.run.tree.GetEntries()
        self.GetBeamInterruptions()
        excluded = [i for i in np.arange(0, self.excludefirst+1)] # first events
        for jump in self.jumps:
            excluded += [i for i in np.arange(jump-5000, jump+20001)] # events around jumps
        excluded.sort()
        all_events = np.arange(0, maxevent)
        included = np.delete(all_events, excluded)
        return included

    def GetLandau(self, channel, canvas=None, drawoption="", color=ROOT.kBlue, normalized=True):
        '''

        :param channel:
        :param canvas:
        :param drawoption:
        :param color:
        :param normalized:
        :return:
        '''
        if canvas == None:
            canvas = ROOT.TCanvas("landaucanvas")
        canvas.cd()
        print "making Landau using\nSignal def:\n\t{signal}\nCut:\n\t{cut}".format(signal=self.signaldefinition, cut=self.cut)
        self.run.tree.Draw((self.signaldefinition+">>landau{run}{channel}(600, -100, 500)").format(channel=channel, run=self.run.run_number), self.cut.format(channel=channel), drawoption, 10000000, self.excludefirst)
        canvas.Update()
        histo = ROOT.gROOT.FindObject("landau{run}{channel}".format(channel=channel, run=self.run.run_number))

        histo.SetLineColor(color)
        stats = histo.FindObject("stats")
        stats.SetTextColor(color)
        canvas.Modified()

        if normalized:
            histo.Scale(1./histo.GetMaximum())
            histo.Draw(drawoption)
            canvas.Update()

    def LoadTrackData(self, minimum_bincontent = 1): # min_bincontent in config file
        '''
        Create a bin collection object as self.Pad and load data from ROOT TTree
        into the Pad object. Then get the 2-dim signal distribution from self.Pad
        :param minimum_bincontent: Bins with less hits are ignored
        :return: -
        '''
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer" # bins with less hits are ignored
        self.minimum_bincontent = minimum_bincontent

        # create a bin collection object:
        self.Pad = {}
        for ch in self.run.GetChannels():
            self.Pad[ch] = BinCollection(*self.config_object.Get2DAttributes(), parent_analysis_obj=self, channel=ch)

        # fill two 2-dim histograms to collect the hits and signal strength
        x_ = {}
        y_ = {}
        channels = self.run.GetChannels()
        included = self.GetIncludedEvents(maxevent=100000) # all event numbers without jump events and initial cut

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

                if (not pulser and not is_saturated and fft_mean>50 and fft_mean<500 and INVfft_max>1e-4):
                    self.Pad[ch].Fill(x_[ch], y_[ch], signal_)



        #self.Pad[channel].MakeFits()
        self.Signal2DDistribution = {}
        for ch in channels:
            self.Signal2DDistribution[ch] = self.Pad[ch].GetMeanSignalDistribution(self.minimum_bincontent)
            self.Signal2DDistribution[ch].SetDirectory(0)
            self.Signal2DDistribution[ch].SetStats(False)
            self.Signal2DDistribution[ch].GetXaxis().SetTitle("pos x / cm")
            self.Signal2DDistribution[ch].GetYaxis().SetTitle("pos y / cm")
            self.Signal2DDistribution[ch].GetYaxis().SetTitleOffset(1.4)

        self.Checklist["LoadTrackData"] = True

    def CreateHitMap(self, channel, saveplot = False, drawoption = "colz", RemoveLowStatBins = 0):
        assert (channel in self.run.GetChannels()), "channel not selected"
        if not self.Checklist["LoadTrackData"]:
            self.LoadTrackData()

        if RemoveLowStatBins > 0:
            if type(RemoveLowStatBins) == t.BooleanType:
                RemoveLowStatBins = 10
            bins = (self.Pad[channel].counthisto.GetNbinsX() + 2)*(self.Pad[channel].counthisto.GetNbinsY() + 2)
            for bin in xrange(bins):
                if self.Pad[channel].counthisto.GetBinContent(bin) < RemoveLowStatBins:
                    coordinates = self.Pad[channel].GetBinCenter(bin)
                    content = self.Pad[channel].counthisto.GetBinContent(bin)
                    self.Pad[channel].counthisto.Fill(coordinates[0], coordinates[1], -content)
            extension = "_min"+str(RemoveLowStatBins)
        else:
            extension = ""

        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        canvas = ROOT.TCanvas("canvas", "Hits", 500, 500) # adjust the width slightly
        canvas.cd()
        self.Pad[channel].counthisto.SetStats(False)
        self.Pad[channel].counthisto.Draw(drawoption)#"surf2")
        #self.Pad[channel].counthisto.Draw("CONT1 SAME")
        if saveplot:
            self.SavePlots("HitMap"+extension, "png")
        canvas.Update()
        self.IfWait("Hits Distribution shown")
        self.Checklist["HitsDistribution"] = True

    def CreateSignalMaps(self, saveplots = False, savename = "2DSignalDistribution",ending="png",saveDir = "Results/", show3d = False):
        '''
        Creates 2D Signal Distribution plot
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
        self.signal_canvas.cd(1)
        ROOT.SetOwnership(self.signal_canvas, False)

        ROOT.gStyle.SetPalette(53) # Dark Body Radiator palette
        ROOT.gStyle.SetNumberContours(999)

        for channel in self.run.GetChannels():
            self.signal_canvas.cd(channels.index(channel)+1)
            # Plot the Signal2D TH2D histogram

            if show3d:
                self.Signal2DDistribution[channel].Draw("SPEC dm(2,10) pa(1,1,1) ci(1,1,1) a(15,45,0) s(1,1)")
                savename += "3D"
            else:
                self.Signal2DDistribution[channel].Draw("colz")

        self.signal_canvas.Update()
        self.IfWait("2d drawn")
        if saveplots:
            self.SavePlots(savename, ending, saveDir)

    def CreateMeanSignalHistogram(self, channel, saveplots = False, savename = "MeanSignalDistribution",ending="png",saveDir = "Results/", show = False):
        '''
        Creates the Mean Signal Response Histogramm of the two-dimensional Mean Signal Map (each bin
        of the 2D dist. is one entry in the histogram)
        :param channel:
        :param saveplots:
        :param savename:
        :param ending:
        :param saveDir:
        :param show:
        :return:
        '''
        if not hasattr(self, "Pad"):
            print "Performing auto analysis with minimum binhits 1"
            self.LoadTrackData(minimum_bincontent=1)
        if saveplots:
            show = True

        #self.signal_canvas = ROOT.TCanvas()
        #ROOT.SetOwnership(self, False)
        if show:
            if hasattr(self, "signal_canvas"):
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
        self.MeanSignalHisto = {}
        self.MeanSignalHisto[channel] = ROOT.TH1D(self.MeanSignalHisto_name, "Run{run}: {diamond} Mean Signal Histogram".format(run=self.run.run_number, diamond=self.run.diamondname[channel]), 100, minimum,maximum)
        ROOT.SetOwnership(self.MeanSignalHisto[channel], False)

        nbins = (self.Signal2DDistribution[channel].GetNbinsX()+2)*(self.Signal2DDistribution[channel].GetNbinsY()+2)
        for i in xrange(nbins):
            bincontent = self.Signal2DDistribution[channel].GetBinContent(i)
            binhits = self.Pad[channel].counthisto.GetBinContent(i)
            if binhits >= self.minimum_bincontent:
                self.MeanSignalHisto[channel].Fill(bincontent)
        if show:
            self.MeanSignalHisto[channel].Draw()
            self.signal_canvas.Update()
        else:
            print "INFO: MeanSignalHisto created as attribute: self.MeanSignalHisto (ROOT.TH1D)"

        if saveplots:
            self.SavePlots(savename, ending, saveDir)

        self.Checklist["MeanSignalHisto"] = True

    def CreateExtendedSignalMaps(self,saveplots = False,savename = "SignalDistribution",ending="png",saveDir = None, PS=False, test=""):
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

            self.combined_canvas.cd(1+channels.index(channel)*2)
            #ROOT.gStyle.SetPalette(55)
            # ROOT.gStyle.SetNumberContours(999)
            ROOT.gStyle.SetPalette(53)
            ROOT.gStyle.SetNumberContours(999)
            self.Signal2DDistribution[channel].Draw("colz")#"CONT1Z")#)'colz')
            self.combined_canvas.cd(2+channels.index(channel)*2)

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
            self.SavePlots(savename, ending, saveDir)
            self.SavePlots(savename, "root", saveDir)
        if PS:
            ROOT.gStyle.SetHistFillColor(7)
            ROOT.gStyle.SetHistFillStyle(3003)
        self.IfWait("Combined 2D Signal DistributionsShown")


