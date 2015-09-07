from AbstractClasses.Elementary import Elementary
import ROOT
import copy, collections, os
from array import array
import numpy as np
import pickle
import ConfigParser


class Cut(Elementary):

    def __init__(self, parent_analysis, verbose=True):
        print "init cut"
        self.analysis = parent_analysis
        self._checklist = {
            "RemoveBeamInterruptions": {
                0: False,
                3: False
            },
            "GenerateCutString": False
        }
        self.userCutTypes = { # readable cut types
            "EventRange":           "Evts.50k-1433k",
            "ExcludeFirst":         "50k+",
            "noPulser":             "!pulser",
            "notSaturated":         "!saturated",
            "noBeamInter":          "!BeamInter.",
            "FFT":                  "FFT",
            "Tracks":               "TrackRec",
        }
        self._cutTypes = { # default values
            "EventRange":             [],             # [1234, 123456]
            "ExcludeFirst":         0,              # 50000 events
            "noPulser":             1,              # 1: nopulser, 0: pulser, -1: no cut
            "notSaturated":         True,
            "noBeamInter":          True,
            "FFT":                  False,
            "Tracks":               True,
        }
        self.excludefirst = 0
        self.cut = {
            0: "", # cut diamond 1
            3: ""  # cut diamond 2
        }
        Elementary.__init__(self, verbose=verbose)
        self.GenerateCutString()

    def LoadConfig(self):
        print "config"
        configfile = "Configuration/AnalysisConfig_"+self.TESTCAMPAIGN+".cfg"
        parser = ConfigParser.ConfigParser()
        parser.read(configfile)

        # individual additional cuts:
        self.cut[0] = parser.get("CUT", "cut1") if not parser.get("CUT", "cut1") in ["-1", "", "True", "False"] else ""
        self.cut[3] = parser.get("CUT", "cut2") if not parser.get("CUT", "cut2") in ["-1", "", "True", "False"] else ""

        # pulser cut:
        self._cutTypes["noPulser"] = parser.getint("CUT", "notPulser")

        # exclude first: (negative: time in minutes, positive: nevents)
        excludefirst = parser.getint("CUT", "excludefirst")
        self.SetExcludeFirst(excludefirst)

        # event range cuts:
        EventRange_min = parser.getint("CUT", "EventRange_min")
        EventRange_max = parser.getint("CUT", "EventRange_max")
        self.SetEventRange(min_event=EventRange_min, max_event=EventRange_max)

        # not saturated cut:
        self._cutTypes["notSaturated"] = parser.getboolean("CUT", "notSaturated")

        # not beam interruption cut:
        self._cutTypes["noBeamInter"] = parser.getboolean("CUT", "noBeamInter")
        self.excludeBeforeJump = parser.getint("CUT", "excludeBeforeJump")
        self.excludeAfterJump = parser.getint("CUT", "excludeAfterJump")

        # FFT cut:
        self._cutTypes["FFT"] = parser.getboolean("CUT", "FFT")

        # has tracks cut:
        self._cutTypes["Tracks"] = parser.getboolean("CUT", "hasTracks")
    
    def SetEventRange(self, min_event=-1, max_event=-1):
        if min_event > 0 and max_event > 0:
            self._cutTypes["EventRange"] = [min_event, max_event]
        elif min_event > 0:
            maxevent = self.analysis.GetEventAtTime(-1)
            self._cutTypes["EventRange"] = [min_event, maxevent]
        elif max_event > 0:
            self._cutTypes["EventRange"] = [self.excludefirst, max_event]
        else:
            self._cutTypes["EventRange"] = []

    def GetIncludedEvents(self, maxevent=None):
        '''
        Get List Of all event numbers, which are neither excluded by beaminerruptions nor
        events from the very beginnning
        :return: list of included event numbers
        '''
        minevent = self.GetMinEvent()
        if maxevent == None:
            maxevent = self.GetMaxEvent()

        excluded = [i for i in np.arange(0, minevent)] # first events
        if self._cutTypes["noBeamInter"]:
            self.GetBeamInterruptions()
            for i in xrange(len(self.jumpsRanges["start"])):
                excluded += [i for i in np.arange(self.jumpsRanges["start"][i], self.jumpsRanges["stop"][i]+1)] # events around jumps
        excluded.sort()
        all_events = np.arange(0, maxevent)
        included = np.delete(all_events, excluded)
        return included

    def SetExcludeFirst(self, n):
        if n >= 0:
            self._SetExcludeFirst(nevents=n)
        else:
            self._SetExcludeFirstTime(seconds=(-1)*n*60)

    def GetEventRange(self):
        return self._cutTypes["EventRange"]

    def GetMinEvent(self):
        if self._cutTypes["EventRange"] != []:
            return self._cutTypes["EventRange"][0]
        elif self._cutTypes["ExcludeFirst"] >0:
            return self._cutTypes["ExcludeFirst"]
        else:
            return 0

    def GetNEvents(self):
        totEvents = self.analysis.GetEventAtTime(-1)
        if self._cutTypes["EventRange"] != []:
            return self._cutTypes["EventRange"][1] - self._cutTypes["EventRange"][0]
        elif self._cutTypes["ExcludeFirst"] >0:
            return totEvents - self._cutTypes["ExcludeFirst"]
        else:
            return totEvents

    def GetMaxEvent(self):
        totEvents = self.analysis.GetEventAtTime(-1)
        if self._cutTypes["EventRange"] != []:
            return self._cutTypes["EventRange"][1]
        else:
            return totEvents

    def _SetExcludeFirstTime(self, seconds):
        event = self.analysis.GetEventAtTime(dt=seconds)
        self._SetExcludeFirst(nevents=event)
        if seconds>0: self.userCutTypes["ExcludeFirst"] = str(int(seconds)/60)+"min+"

    def _SetExcludeFirst(self, nevents):
        if nevents>0:
            self.userCutTypes["ExcludeFirst"] = str(int(nevents)/1000)+"k+"
        else:
            self.userCutTypes["ExcludeFirst"] = ""
        self._cutTypes["ExcludeFirst"] = nevents
        self.excludefirst = nevents

    def GenerateCutString(self, gen_PulserCut=True, gen_EventRange=True, gen_ExcludeFirst=True):
        print "generate cutstring"
        cutstring = ""
        if self._checklist["GenerateCutString"]:
            self.LoadConfig() # re-generate

        if self._cutTypes["EventRange"] != [] and gen_EventRange:
            if cutstring != "": cutstring += "&&"
            cutstring += "(event_number<={maxevent}&&event_number>={minevent})".format(minevent=self._cutTypes["EventRange"][0], maxevent=self._cutTypes["EventRange"][1])
        elif self._cutTypes["ExcludeFirst"] > 0 and gen_ExcludeFirst:
            if cutstring != "": cutstring += "&&"
            cutstring += "event_number>={minevent}".format(minevent=self._cutTypes["ExcludeFirst"])

        if self._cutTypes["noPulser"] in [0,1] and gen_PulserCut:
            if cutstring != "": cutstring += "&&"
            if self._cutTypes["noPulser"] == 1:
                cutstring += "!pulser"
            else:
                cutstring += "pulser"

        if self._cutTypes["notSaturated"] in [0,1]:
            if cutstring != "": cutstring += "&&"
            if self._cutTypes["notSaturated"] == 1:
                cutstring += "!is_saturated[{channel}]"
            else:
                cutstring += "is_saturated[{channel}]"

        if self._cutTypes["Tracks"]:
            if cutstring != "": cutstring += "&&"
            cutstring += "n_tracks"

        for channel in [0,3]:
            self.cut[channel] += cutstring
            self.cut[channel] = self.cut[channel].format(channel=channel)

        if self._cutTypes["noBeamInter"] and self._checklist["GenerateCutString"]:
            self._RemoveBeamInterruptions(justDoIt=True)
        elif self._cutTypes["noBeamInter"]:
            self._RemoveBeamInterruptions()

        self._checklist["GenerateCutString"] = True
        self._cutStringSettings = {
            "gen_PulserCut": gen_PulserCut,
            "gen_EventRange": gen_EventRange,
            "gen_ExcludeFirst": gen_ExcludeFirst
        }

    def _checkCutStringSettings(self, gen_PulserCut, gen_EventRange, gen_ExcludeFirst):
        if self._cutStringSettings["gen_PulserCut"]==gen_PulserCut and self._cutStringSettings["gen_EventRange"]==gen_EventRange and self._cutStringSettings["gen_ExcludeFirst"]==gen_ExcludeFirst:
            return True
        else:
            return False

    def _FindBeamInterruptions(self):
        '''
        Finds the beam interruptions
        :return: list of event numbers where beam interruptions occures
        '''
        print "Searching for beam interruptions.."
        nentries = self.analysis.run.tree.GetEntries()
    #    last_entry = self.analysis.run.tree.GetEntry(nentries-1)
    #    max_time = self.analysis.run.tree.time

        canvas = ROOT.TCanvas("beaminterruptioncanvas", "beaminterruptioncanvas")
        self.analysis.run.tree.Draw('time:event_number')

    #    graph = copy.deepcopy(ROOT.c1.FindObject('Graph'))
        histo = copy.deepcopy(canvas.FindObject('htemp'))

        histo.SetTitle('run %3d' %(self.analysis.run.run_number))
        histo.SetName('run %3d' %(self.analysis.run.run_number))

        # get event numbers and dt's
        dts = []
        evs = []
        i = self.excludefirst
        step = 100
        while i+step < nentries:
            self.analysis.run.tree.GetEntry(i)
            t1 = self.analysis.run.tree.time
            evs.append(self.analysis.run.tree.event_number)
            self.analysis.run.tree.GetEntry(i+step)
            t2 = self.analysis.run.tree.time
            dt = (t2 - t1)
            dts.append(dt)
            i += step

        self.jumps = []

        deq = collections.deque(dts[:100],100)
        first = True
        for i in dts[101:]:
            avg = numpy.mean(deq)
            if abs(i / avg - 1.) > 0.3:
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
        lat.DrawLatex(0.2,0.85, 'run %d' %(self.analysis.run.run_number) )

        if not os.path.exists("beaminterruptions"):
            os.mkdir("beaminterruptions")
        if not os.path.exists("beaminterruptions/plots"):
            os.mkdir("beaminterruptions/plots")
        if not os.path.exists("beaminterruptions/data"):
            os.mkdir("beaminterruptions/data")

        # save jump list to file
        jumpfile = open("beaminterruptions/data/{testcampaign}Run_{run}.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number), "wb")
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

            outfile = open('beaminterruptions/jumps_{testcampaign}.txt'.format(testcampaign=self.TESTCAMPAIGN),'r+a')
            # check if the run is already in the file
            runInFile = False
            lines = outfile.readlines()
            for i in lines:
                if len(i.split()) > 0 and i.split()[0] == str(self.analysis.run.run_number):
                    runInFile = True
            if not runInFile:
                outfile.write(str(self.analysis.run.run_number)+'\t\t')

            lat.SetTextColor(ROOT.kBlack)
            for i in self.jumps:
                ind = self.jumps.index(i)
                lat.DrawLatex(0.2, 0.80-ind*0.05, '#%d at %d' %(ind, i) )
                if not runInFile:
                    outfile.write(str(i)+'\t')
            if not runInFile:
                outfile.write('\n')
            outfile.close()


        ROOT.c1.SaveAs('beaminterruptions/plots/%djumpSearch_run%d.png' %(self.TESTCAMPAIGN, self.analysis.run.run_number))

        canvas.Close()
        return self.jumps

    def GetBeamInterruptions(self):
        '''
        If there is beam interruption data, it will load them - otherwise it will run the beam interruption analysis
        it will create the attribute self.jumps, which is a list of event numbers, where a jump occures
        :return: list of events where beam interruptions occures
        '''
        if not hasattr(self, "jumpsRanges"):
            picklepath = "beaminterruptions/data/{testcampaign}Run_{run}.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
            if os.path.exists(picklepath):
                print "Loading beam interruption data from pickle file: \n\t"+picklepath
                jumpfile = open(picklepath, "rb")
                self.jumps = pickle.load(jumpfile)
                self._ReduceJumps()
                jumpfile.close()
            else:
                print "No pickle file found at: ", picklepath, "\n .. analyzing beam interruptions.. "
                self._FindBeamInterruptions()
                self._ReduceJumps()

        return self.jumps

    def GetCutFunctionDef(self):
        defstring_ = "lambda pulser, is_saturated, n_tracks, fft_mean, INVfft_max: "
        def_ = ""

        if self._cutTypes["noPulser"] == 1:
            def_ += "not pulser"
        elif self._cutTypes["noPulser"] == 0:
            def_ += "pulser"

        if self._cutTypes["notSaturated"]:
            if def_ != "": def_ += " and "
            def_ += "not is_saturated"

        if self._cutTypes["Tracks"]:
            if def_ != "": def_ += " and "
            def_ += "n_tracks"

        if self._cutTypes["FFT"]:
            if def_ != "": def_ += " and "
            assert(False), "FFT cut not yet implemented in GetCutFunctionDef() method of Cut class. "
            # to do: FFT entry in _cutTypes should be dict and/or contain a TCutG instance

        return defstring_+def_

    def _ReduceJumps(self):
        if not hasattr(self, "jumpsRanges") and len(self.jumps)>0:
            self.jumps.sort()
            events = self.analysis.GetEventAtTime(-1)
            selection = events*[0]
            high = self.excludeAfterJump
            low = self.excludeBeforeJump
            reduced_jumps = []
            reduced_ends = []
            for jump in self.jumps:
                c = 1 if (jump-low)>0 else 0
                selection[c*(jump-low):(jump+high+1)] = len(selection[c*(jump-low):(jump+high+1)])*[1]

            for i in xrange(len(selection)-1):
                if selection[i] != selection[i+1]:
                    if selection[i] == 0:
                        print "jump start: ", i+1
                        reduced_jumps.append(i+1)
                    else:
                        print "jump end: ", i+1
                        reduced_ends.append(i+1)
            if reduced_ends[0]<reduced_jumps[0]:
                reduced_jumps = [0]+reduced_jumps
            if reduced_jumps[-1]>reduced_ends[-1]:
                reduced_ends = reduced_ends+[events]
            self.jumps = reduced_jumps
            self.jumpsRanges = {
                "start": reduced_jumps,
                "stop": reduced_ends
            }
        else:
            self.jumpsRanges = {
                "start": [],
                "stop": []
            }

    def _RemoveBeamInterruptions(self, justDoIt=False):
        '''
        This adds the restrictions to the cut string such that beam interruptions are excluded each time the
        cut is applied.
        :return:
        '''
        for channel in [0,3]:
            if not self._checklist["RemoveBeamInterruptions"][channel] or justDoIt:
                self.GetBeamInterruptions()

                njumps = len(self.jumpsRanges["start"])
                for i in xrange(njumps):
                    if self.cut[channel] != "": self.cut[channel] += "&&"
                    self.cut[channel] += "!(event_number<={upper}&&event_number>={lower})".format(upper=self.jumpsRanges["stop"][i], lower=self.jumpsRanges["start"][i])
                self._checklist["RemoveBeamInterruptions"][channel] = True

        return self.cut

    def AddCutString(self, cutstring, channel=None):
        pass

    def GetCut(self, channel, gen_PulserCut=True, gen_EventRange=True, gen_ExcludeFirst=True):
        if not self._checkCutStringSettings(gen_PulserCut, gen_EventRange, gen_ExcludeFirst):
            self.GenerateCutString(gen_PulserCut, gen_EventRange, gen_ExcludeFirst)
        return self.cut[channel]

    def GetUserCutString(self):
        string_ = ""
        return string_

    def ShowCuts(self):
        pass

    def SetCut(self):
        pass

    def SetFFTCut(self, channel):
        pass