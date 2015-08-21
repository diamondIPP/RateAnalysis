from Helper.Initializer import initializer
from Runinfos.RunInfo import RunInfo
from DiamondClass import Diamond
from Elementary import Elementary
import ROOT
import os
import ConfigParser
import json
import copy

default_info =  {
        "persons on shift": "-",
        "run info": "-",
        "type": "signal",
        "configuration": "signal",
        "mask": "-",
        "masked pixels": 0,
        "diamond 1": "CH_0",
        "diamond 2": "CH_3",
        "hv dia1": 0,
        "hv dia2": 0,
        "for1": 0,
        "for2": 0,
        "fs11": 0,
        "fsh13": 0,
        "quadrupole": "-",
        "analogue current": 0,
        "digital current": 0,
        "begin date": "-/-/-",
        "trim time": "-:-:-",
        "config time": "-:-:-",
        "start time": "-:-:-",
        "trig accept time": "-:-:-",
        "opening time": "-:-:-",
        "open time": "-:-:-",
        "stop time": "-:-:-",
        "raw rate": 0,
        "prescaled rate": 0,
        "to TLU rate": 0,
        "pulser accept rate": 0,
        "cmspixel events": 0,
        "drs4 events": 0,
        "datacollector events": 0,
        "aimed flux": 0,
        "measured flux": 0,
        "user comments": "-",
        "is good run": True
}

class Run(Elementary):
    '''

    '''

    current_run = {}
    operationmode = ''
    TrackingPadAnalysis = {}

    def __init__(self, run_number, channels=3, validate = False, verbose = False):
        '''

        :param run_number: number of the run
        :param channels: 0x1=ch0; 0x2=ch3
        :param validate:
        :param verbose:
        :return:
        '''
        Elementary.__init__(self, verbose=verbose)

        self.run_number = -1
        self.LoadConfig()

        if validate:
            self.ValidateRuns()

        if run_number != None:
            assert(run_number > 0), "incorrect run_number"
            self.SetRun(run_number)
        else:
            self.LoadRunInfo()
        self.diamondname = {
            0: str(self.RunInfo["diamond 1"]),
            3: str(self.RunInfo["diamond 2"])
        }
        self.bias = {
            0: self.RunInfo["hv dia1"],
            3: self.RunInfo["hv dia2"]
        }
        self.SetChannels(channels)
        self.IsMonteCarlo = False

    def LoadConfig(self):
        machineConfigParser = ConfigParser.ConfigParser()
        machineConfigParser.read('Configuration/Machineconfig.cfg')
        self.operationmode = machineConfigParser.get('EXEC-MACHINE','operationmode')
        self.ShowAndWait = False
        runConfigParser = ConfigParser.ConfigParser()
        runConfigParser.read("Configuration/RunConfig_"+self.TESTCAMPAIGN+".cfg")
        self.filename = runConfigParser.get('BASIC', 'filename')
        self.treename = runConfigParser.get('BASIC', 'treename')
        self.sshrunpath = runConfigParser.get('BASIC', 'runpath')
        self.runinfofile = runConfigParser.get('BASIC', 'runinfofile')

    def LoadRunInfo(self):
        self.RunInfo = {}
        try:
            f = open(self.runinfofile, "r")
            data = json.load(f)
            f.close()
            self.allRunKeys = copy.deepcopy(data.keys())
            loaderror = False
        except IOError:
            print "WARNING: unable to load json file:\n\t{file}".format(file=self.runinfofile)
            loaderror = True

        if self.run_number >= 0:
            if not loaderror:
                self.RunInfo = data.get(str(self.run_number)) # may:  = data.get("150800"+str(self.run_number).zfill(3))
                self.RenameRunInfoKeys()
            else:
                self.RunInfo = default_info
            self.current_run = self.RunInfo
            if self.RunInfo is None:
                self.RunInfo = {}
                self.RunInfo = default_info
                print "\nWARNING: No RunInfo could be loaded from file! \n"
                return 0
            else:
                return 1
        else:
            return 0

    def RenameRunInfoKeys(self):

        try:
            for key in default_info.keys():
                tmp = self.RunInfo[key]
        except KeyError:
            rename = True
        else:
            rename = False

        if rename:
            KeyConfigParser = ConfigParser.ConfigParser()
            KeyConfigParser.read("Configuration/RunInfoKeyConfig_"+self.TESTCAMPAIGN+".cfg")
            Persons =           KeyConfigParser.get("KEYNAMES", "Persons")
            Runinfo =           KeyConfigParser.get("KEYNAMES", "Runinfo")
            Typ =               KeyConfigParser.get("KEYNAMES", "Typ")
            Configuration =     KeyConfigParser.get("KEYNAMES", "Configuration")
            Mask =              KeyConfigParser.get("KEYNAMES", "Mask")
            Masked_pixels =     KeyConfigParser.get("KEYNAMES", "Masked_pixels")
            DiamondName1 =      KeyConfigParser.get("KEYNAMES", "DiamondName1")
            DiamondName2 =      KeyConfigParser.get("KEYNAMES", "DiamondName2")
            DiamondHV1 =        KeyConfigParser.get("KEYNAMES", "DiamondHV1")
            DiamondHV2 =        KeyConfigParser.get("KEYNAMES", "DiamondHV2")
            FOR1 =              KeyConfigParser.get("KEYNAMES", "FOR1")
            FOR2 =              KeyConfigParser.get("KEYNAMES", "FOR2")
            FS11 =              KeyConfigParser.get("KEYNAMES", "FS11")
            FSH13 =             KeyConfigParser.get("KEYNAMES", "FSH13")
            Quadrupole =        KeyConfigParser.get("KEYNAMES", "Quadrupole")
            AnalogCurrent =     KeyConfigParser.get("KEYNAMES", "AnalogCurrent")
            DigitalCurrent =    KeyConfigParser.get("KEYNAMES", "DigitalCurrent")
            BeginDate =         KeyConfigParser.get("KEYNAMES", "BeginDate")
            TrimTime =          KeyConfigParser.get("KEYNAMES", "TrimTime")
            ConfigTime =        KeyConfigParser.get("KEYNAMES", "ConfigTime")
            StartTime =         KeyConfigParser.get("KEYNAMES", "StartTime")
            TrigAcceptTime =    KeyConfigParser.get("KEYNAMES", "TrigAcceptTime")
            OpeningTime =       KeyConfigParser.get("KEYNAMES", "OpeningTime")
            OpenTime =          KeyConfigParser.get("KEYNAMES", "OpenTime")
            StopTime =          KeyConfigParser.get("KEYNAMES", "StopTime")
            RawRate =           KeyConfigParser.get("KEYNAMES", "RawRate")
            PrescaledRate =     KeyConfigParser.get("KEYNAMES", "PrescaledRate")
            ToTLURate =         KeyConfigParser.get("KEYNAMES", "ToTLURate")
            PulserAcceptedRate =KeyConfigParser.get("KEYNAMES", "PulserAcceptedRate")
            CMSPixelEvents =    KeyConfigParser.get("KEYNAMES", "CMSPixelEvents")
            DRS4Events =        KeyConfigParser.get("KEYNAMES", "DRS4Events")
            DataCollectorEvents = KeyConfigParser.get("KEYNAMES", "DataCollectorEvents")
            AimedFlux =         KeyConfigParser.get("KEYNAMES", "AimedFlux")
            MeasuredFlux =      KeyConfigParser.get("KEYNAMES", "MeasuredFlux")
            UserComment =       KeyConfigParser.get("KEYNAMES", "UserComment")
            IsGoodRun =         KeyConfigParser.get("KEYNAMES", "IsGoodRun")

            runinfo = copy.deepcopy(default_info)

            if Persons != "-1":             runinfo["persons on shift"] =   self.RunInfo[Persons]
            if Runinfo != "-1":             runinfo["run info"] =           self.RunInfo[Runinfo]
            if Typ != "-1":                 runinfo["type"] =               self.RunInfo[Typ]
            if Configuration != "-1":       runinfo["configuration"] =      self.RunInfo[Configuration]
            if Mask != "-1":                runinfo["mask"] =               self.RunInfo[Mask]
            if Masked_pixels != "-1":       runinfo["masked pixels"] =      self.RunInfo[Masked_pixels]
            if DiamondName1 != "-1":        runinfo["diamond 1"] =          self.RunInfo[DiamondName1]
            if DiamondName2 != "-1":        runinfo["diamond 2"] =          self.RunInfo[DiamondName2]
            if DiamondHV1 != "-1":          runinfo["hv dia1"] =            self.RunInfo[DiamondHV1]
            if DiamondHV2 != "-1":          runinfo["hv dia2"] =            self.RunInfo[DiamondHV2]
            if FOR1 != "-1":                runinfo["for1"] =               self.RunInfo[FOR1]
            if FOR2 != "-1":                runinfo["for2"] =               self.RunInfo[FOR2]
            if FS11 != "-1":                runinfo["fs11"] =               self.RunInfo[FS11]
            if FSH13 != "-1":               runinfo["fsh13"] =              self.RunInfo[FSH13]
            if Quadrupole != "-1":          runinfo["quadrupole"] =         self.RunInfo[Quadrupole]
            if AnalogCurrent != "-1":       runinfo["analogue current"] =   self.RunInfo[AnalogCurrent]
            if DigitalCurrent != "-1":      runinfo["digital current"] =    self.RunInfo[DigitalCurrent]
            if BeginDate != "-1":           runinfo["begin date"] =         self.RunInfo[BeginDate]
            if TrimTime != "-1":            runinfo["trim time"] =          self.RunInfo[TrimTime]
            if ConfigTime != "-1":          runinfo["config time"] =        self.RunInfo[ConfigTime]
            if StartTime != "-1":           runinfo["start time"] =         self.RunInfo[StartTime]
            if TrigAcceptTime != "-1":      runinfo["trig accept time"] =   self.RunInfo[TrigAcceptTime]
            if OpeningTime != "-1":         runinfo["opening time"] =       self.RunInfo[OpeningTime]
            if OpenTime != "-1":            runinfo["open time"] =          self.RunInfo[OpenTime]
            if StopTime != "-1":            runinfo["stop time"] =          self.RunInfo[StopTime]
            if RawRate != "-1":             runinfo["raw rate"] =           self.RunInfo[RawRate]
            if PrescaledRate != "-1":       runinfo["prescaled rate"] =     self.RunInfo[PrescaledRate]
            if ToTLURate != "-1":           runinfo["to TLU rate"] =        self.RunInfo[ToTLURate]
            if PulserAcceptedRate != "-1":  runinfo["pulser accept rate"] = self.RunInfo[PulserAcceptedRate]
            if CMSPixelEvents != "-1":      runinfo["cmspixel events"] =    self.RunInfo[CMSPixelEvents]
            if DRS4Events != "-1":          runinfo["drs4 events"] =        self.RunInfo[DRS4Events]
            if DataCollectorEvents != "-1": runinfo["datacollector events"]=self.RunInfo[DataCollectorEvents]
            if AimedFlux != "-1":           runinfo["aimed flux"] =         self.RunInfo[AimedFlux]
            if MeasuredFlux != "-1":        runinfo["measured flux"] =      self.RunInfo[MeasuredFlux]
            if UserComment != "-1":         runinfo["user comments"] =      self.RunInfo[UserComment]
            if IsGoodRun != "-1":           runinfo["is good run"] =        self.RunInfo[IsGoodRun]

            self.RunInfo = runinfo

    def ValidateRuns(self, list_of_runs = None):
        if list_of_runs != None:
            runs = list_of_runs
        else:
            runs = RunInfo.runs.keys() # list of all runs
        self.VerbosePrint("Validating runs: ",runs)
        for run in runs:
            self.ValidateRun(run)

    def ValidateRun(self,run_number):
        self.SetRun(run_number)
        if not os.path.exists(self.TrackingPadAnalysis['ROOTFile']):
            del RunInfo.runs[run_number]
            print "INFO: path of run number ",run_number, " not found."
            return False
        else:
            return True

    def ResetMC(self):
        pass

    def SetRun(self, run_number, validate = False, loadROOTFile = True):
        if validate:
            boolfunc = self.ValidateRun
        else:
            boolfunc = lambda run: True
        assert(run_number > 0), "incorrect run_number"
        if True or run_number in RunInfo.runs and boolfunc(run_number):
            self.run_number = run_number
            self.LoadRunInfo()

            if self.operationmode == "local-ssh":
                fullROOTFilePath = '/Volumes'+self.sshrunpath+'/'+self.filename+str(run_number).zfill(3)+'.root'
                self.TrackingPadAnalysis['ROOTFile'] = fullROOTFilePath

            if self.operationmode == "ssh":
                fullROOTFilePath = self.sshrunpath+'/'+self.filename+str(run_number).zfill(3)+'.root'
                self.TrackingPadAnalysis['ROOTFile'] = fullROOTFilePath

            if self.operationmode == "local":
                self.TrackingPadAnalysis['ROOTFile'] = 'runs/run_'+str(run_number)+'/'+self.filename+str(run_number).zfill(3)+'.root'


 #           self.RunInfo = self.current_run.copy()
            #a = self.current_run.__dict__
#            self.diamond = Diamond( self.current_run['diamond'])

            self.ResetMC()
            #self.diamond = diamond # diamond is of type Diamond
            #RunInfo.__init__(self,*args)

            if loadROOTFile:
                print "LOADING: ", fullROOTFilePath
                self.rootfile = ROOT.TFile(fullROOTFilePath)
                self.tree = self.rootfile.Get(self.treename) # Get TTree called "track_info"

            assert(bool(self.tree) and bool(self.rootfile)), "Could not load root file: \n\t"+fullROOTFilePath

            #self.diamond.SetName(self.runs[0]['diamond'])
            return True
        else:
            return False

    def SetChannels(self, channels):
        assert(channels>=1 and channels<=3), "invalid channel number: 0x1=ch0; 0x2=ch3"
        self.analyzeCh = {
            0: (channels & 1<<0) == 1<<0,
            3: (channels & 1<<1) == 1<<1
        }

    def GetChannels(self):
        return [i for i in self.analyzeCh.keys() if self.analyzeCh[i]]

    def GetDiamondName(self, channel):
        return self.diamondname[channel]

    def ShowRunInfo(self):
        print "RUN INFO:"
        print "\tRun Number: \t", self.run_number, " (",self.RunInfo["type"],")"
        print "\tDiamond1:   \t", self.diamondname[0], " (",self.bias[0],") | is selected: ", self.analyzeCh[0]
        print "\tDiamond2:   \t", self.diamondname[3], " (",self.bias[3],") | is selected: ", self.analyzeCh[3]

    def DrawRunInfo(self, channel=None, canvas=None, diamondinfo=True, showcut=False, comment=None):
        if canvas != None:
            canvas.cd()
        else:
            print "Draw run info in current pad"
            pad = ROOT.gROOT.GetSelectedPad()
            if pad:
                canvas = pad.GetCanvas()
                canvas.cd()
                pad.cd()
            else:
                print "ERROR: Can't access active Pad"

        lines = 1
        width = 0.25
        if diamondinfo:
            lines += 1
        if showcut and hasattr(self, "analysis"):
            lines += 1
            width = 0.6
        if comment != None:
            lines += 1
            width = max(0.4, width)

        if not hasattr(self, "runInfoLegends"):
            self.runInfoLegends = {}



        if channel != None:
            self.runInfoLegends[channel] = ROOT.TLegend(0.1, 0.86-(lines-1)*0.03, 0.1+width, 0.9)
            self.runInfoLegends[channel].SetMargin(0.05)
            self.runInfoLegends[channel].AddEntry(0, "Run{run} Ch{channel} ({rate})".format(run=self.run_number, channel=channel, rate=self.GetRateString()), "")
            if diamondinfo: self.runInfoLegends[channel].AddEntry(0, "{diamond} ({bias:+}V)".format(diamond=self.diamondname[channel], bias=self.bias[channel]), "")
            if showcut and hasattr(self, "analysis"): self.runInfoLegends[channel].AddEntry(0, "Cut: {cut}".format(cut=self.analysis.cut.format(channel=channel)), "")
            if comment != None: self.runInfoLegends[channel].AddEntry(0, comment, "")
            self.runInfoLegends[channel].Draw("same")
        else:
            if comment != None:
                lines = 2
            else:
                lines = 1
                width = 0.15
            self.runInfoLegends[-1] = ROOT.TLegend(0.1, 0.9-lines*0.05, 0.1+width, 0.9)
            self.runInfoLegends[-1].SetMargin(0.05)
            self.runInfoLegends[-1].AddEntry(0, "Run{run} ({rate})".format(run=self.run_number, rate=self.GetRateString()), "")
            if comment != None: self.runInfoLegends[-1].AddEntry(0, comment, "")
            self.runInfoLegends[-1].Draw("same")
            pad.Modified()

    def GetRateString(self):
        rate = self.RunInfo["measured flux"]
        if rate>1000:
            unit = "MHz"
            rate = round(rate/1000.,1)
        else:
            unit = "kHz"
            rate = int(round(rate,0))
        return "{rate:>3}{unit}".format(rate=rate, unit=unit)
