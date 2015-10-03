from Helper.Initializer import initializer
from Runinfos.RunInfo import RunInfo
from DiamondClass import Diamond
from Elementary import Elementary
from datetime import datetime as dt
import ROOT
import os
import ConfigParser
import json
import csv
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
        "begin date": "2999-03-14T15:26:53Z",
        "trim time": "-:-:-",
        "config time": "-:-:-",
        "start time": "2999-03-14T15:26:53Z",
        "trig accept time": "-:-:-",
        "opening time": "-:-:-",
        "open time": "-:-:-",
        "stop time": "2999-03-14T16:26:53Z",
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

    def __init__(self, run_number, diamonds=3, validate = False, verbose = False, maskfilename=""):
        '''

        :param run_number: number of the run
        :param diamonds: 0x1=ch0; 0x2=ch3
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
        self._LoadTiming()
        self.CalculateRate(maskfilename=maskfilename)
        self.diamondname = {
            0: str(self.RunInfo["diamond 1"]),
            3: str(self.RunInfo["diamond 2"])
        }
        self.bias = {
            0: self.RunInfo["hv dia1"],
            3: self.RunInfo["hv dia2"]
        }
        self.SetChannels(diamonds)
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
        self._runlogkeyprefix = runConfigParser.get('BASIC', 'runlog_key_prefix')
        self.runplaninfofile = runConfigParser.get('BASIC', 'runplaninfofile')
        self.maskfilepath = runConfigParser.get('BASIC', 'maskfilepath')
        self.createNewROOTFiles = runConfigParser.getboolean('BASIC', 'createNewROOTFiles')
        self.signalregion_low = runConfigParser.getint('BASIC', 'signalregion_low')
        self.signalregion_high = runConfigParser.getint('BASIC', 'signalregion_high')
        
        # Rootfile Generation Configuration:
        signal_range_low = runConfigParser.getint('ROOTFILE_GENERATION', 'signal_range_low')
        signal_range_high = runConfigParser.getint('ROOTFILE_GENERATION', 'signal_range_high')
        pedestal_range_low = runConfigParser.getint('ROOTFILE_GENERATION', 'pedestal_range_low')
        pedestal_range_high = runConfigParser.getint('ROOTFILE_GENERATION', 'pedestal_range_high')
        pulser_range_low = runConfigParser.getint('ROOTFILE_GENERATION', 'pulser_range_low')
        pulser_range_high = runConfigParser.getint('ROOTFILE_GENERATION', 'pulser_range_high')
        peakintegral1_range_low = runConfigParser.getint('ROOTFILE_GENERATION', 'peakintegral1_range_low')
        peakintegral1_range_high = runConfigParser.getint('ROOTFILE_GENERATION', 'peakintegral1_range_high')
        peakintegral2_range_low = runConfigParser.getint('ROOTFILE_GENERATION', 'peakintegral1_range_low')
        peakintegral2_range_high = runConfigParser.getint('ROOTFILE_GENERATION', 'peakintegral1_range_high')
        peakintegral3_range_low = runConfigParser.getint('ROOTFILE_GENERATION', 'peakintegral1_range_low')
        peakintegral3_range_high = runConfigParser.getint('ROOTFILE_GENERATION', 'peakintegral1_range_high')
        save_waveforms = runConfigParser.getint('ROOTFILE_GENERATION', 'save_waveforms')
        self.converterPrefix = runConfigParser.get('ROOTFILE_GENERATION', "converterPrefix")
        self.eudaqFolder = runConfigParser.get('ROOTFILE_GENERATION', "eudaqFolder")
        self.trackingFolder = runConfigParser.get('ROOTFILE_GENERATION', "trackingFolder")
        self.rawFolder = runConfigParser.get('ROOTFILE_GENERATION', "rawFolder")
        self.rawPrefix = runConfigParser.get('ROOTFILE_GENERATION', "rawPrefix")

        self.rootGenerationConfig = {
            "signal_range": [signal_range_low, signal_range_high],
            "pedestal_range": [pedestal_range_low, pedestal_range_high],
            "pulser_range": [pulser_range_low, pulser_range_high],
            "peakintegral1_range": [peakintegral1_range_low, peakintegral1_range_high],
            "peakintegral2_range": [peakintegral2_range_low, peakintegral2_range_high],
            "peakintegral3_range": [peakintegral3_range_low, peakintegral3_range_high],
            "save_waveforms": save_waveforms
        }

    def LoadRunInfo(self):
        self.RunInfo = {}
        try:
            f = open(self.runinfofile, "r")
            data = json.load(f)
            f.close()
            self.allRunKeys = copy.deepcopy(data.keys())
            loaderror = False
        except IOError:
            print "\n-------------------------------------------------"
            print "WARNING: unable to load json file:\n\t{file}".format(file=self.runinfofile)
            print "-------------------------------------------------\n"
            loaderror = True

        if self.run_number >= 0:
            if not loaderror:
                self.RunInfo = data.get(str(self.run_number)) # may:  = data.get("150800"+str(self.run_number).zfill(3))
                if self.RunInfo == None:
                    # try with run_log key prefix
                    self.RunInfo = data.get(self._runlogkeyprefix+str(self.run_number).zfill(3))
                if self.RunInfo == None:
                    print "INFO: Run not found in json run log file. Default run info will be used."
                    self.RunInfo = default_info
                else:
                    self.RenameRunInfoKeys()
            else:
                self.RunInfo = default_info
            self.current_run = self.RunInfo
        else:
            self.RunInfo = default_info
            return 0

    def CalculateRate(self, maskfilename=""):
        self.VerbosePrint("Calculate rate from mask file:\n\t"+self.RunInfo["mask"])
        if maskfilename != "": self.RunInfo["mask"] = maskfilename
        maskFilePath = self.maskfilepath+"/"+self.RunInfo["mask"] # CONFIG FILE !
        maskdata = {
            0: {
                "cornBot": [],
                "cornTop": []
            },
            3: {
                "cornBot": [],
                "cornTop": []
            }
        }
        try:
            f = open(maskFilePath, "r")
            infile = csv.reader(f, delimiter=" ")
            try:
                for i in xrange(8):
                    line = next(infile)
                    self.VerbosePrint(line)
                    if len(line)>=4:
                        #maskdata[int(line[1])][line[0]] = map(int, line[-2:])
                        maskdata[0][line[0]] = map(int, line[-2:])
            except StopIteration:
                pass
            x_low = maskdata[0]["cornBot"][0]
            y_low = maskdata[0]["cornBot"][1]
            x_high = maskdata[0]["cornTop"][0]
            y_high = maskdata[0]["cornTop"][1]

            self.RunInfo["masked pixels"] = abs((x_high-x_low+1)*(y_high-y_low+1))

            pixelarea = 0.01*0.015 # cm^2
            masked_area = self.RunInfo["masked pixels"]*pixelarea
            rate_Hz = 1.*self.RunInfo["raw rate"]/masked_area # in Hz/cm^2
            self.RunInfo["measured flux"] = rate_Hz/1000. # in kHz/cm^2

            f.close()
            return rate_Hz/1000.
        except IOError:
            print "\nERROR: Could not load mask file, thus not re-calculate rate..\n"
            return 0

    def GetRate(self):
        '''
        Returns the rate during this run. If the mask files are given,
        the rate is calculated by the raw rate and the area of the
        masked pixels in the silicon pixel plane.
        If no mask files are given, the rate returned will be the logged
        rate from the online logbook. (Typos!)
        The rate is given in kHz/cm^2
        :return:
        '''
        return self.RunInfo["measured flux"]

    def _SetConverterConfigFile(self):
        pol_dia1 = self.RunInfo["hv dia1"]
        pol_dia2 = self.RunInfo["hv dia2"]
        assert(pol_dia1!=0 and pol_dia2!=0)
        if pol_dia1 > 0:
            pol_dia1 = 1
        else:
            pol_dia1 = -1
        if pol_dia2 > 0:
            pol_dia2 = 1
        else:
            pol_dia2 = -1
        cparser = ConfigParser.ConfigParser()
        cparser.read(self.eudaqFolder+"/conf/converter.conf")
        print cparser.sections()
        cparser.set("Converter.drs4tree", "polarities", "[{pol1},0,0,{pol2}]".format(pol1=pol_dia1, pol2=pol_dia2))
        for key in self.rootGenerationConfig:
            cparser.set("Converter.drs4tree", key, str(self.rootGenerationConfig[key]))
        f = open(self.eudaqFolder+"/conf/converter.conf", "w")
        cparser.write(f)
        f.close()

        # remove white spaces:
        f = open(self.eudaqFolder+"/conf/converter.conf", "r")
        content = f.readlines()
        f.close()
        for i in xrange(len(content)):
            content[i] = content[i].replace(" ", "")
        f = open(self.eudaqFolder+"/conf/converter.conf", "w")
        f.writelines(content)
        f.close()

    def CreateROOTFile(self, do_tracking=True):
        # path and name of converter output file:
        noTracksROOTFile = os.getcwd()+"/{prefix}{run}.root".format(prefix=self.converterPrefix, run=str(self.run_number).zfill(4))

        if not os.path.exists(noTracksROOTFile):
            # the no-tracks root files doesn't exiist
            self._ConvertRAW()
        else:
            # continue with existing file (no tracks)
            print "noTracks ROOT File found here:"
            print "\t"+noTracksROOTFile

        if not do_tracking:
            # move to data folder:
            os.system("mv "+noTracksROOTFile+" "+self.TrackingPadAnalysis['ROOTFile'])
            self._LoadROOTFile(self.TrackingPadAnalysis['ROOTFile'])
            print "INFO ROOT File generated with NO Tracking information"
        else:
            self._AddTracking(noTracksROOTFile=noTracksROOTFile)

    def _ConvertRAW(self):
        # terminal command for converting raw to root
        converter_cmd = "{eudaq}/bin/Converter.exe -t drs4tree -c {eudaq}/conf/converter.conf {rawfolder}/{prefix}{run}.raw".format(eudaq=self.eudaqFolder, rawfolder=self.rawFolder, prefix=self.rawPrefix, run=str(self.run_number).zfill(4))

        self._SetConverterConfigFile() # prepare the config file
        print "\n\nSTART CONVERTING RAW FILE..."
        print converter_cmd
        os.system(converter_cmd) # convert

    def _AddTracking(self, noTracksROOTFile):

        if self.TESTCAMPAIGN == "201508":
            telescopeID = 9
        elif self.TESTCAMPAIGN == "201505":
            telescopeID = 7
        else:
            telescopeID = 0
            assert(False), "Error. unknown TESTCAMPAIGN"

        # change CWD to TrackingTelescope:
        old_cwd = os.getcwd()
        os.chdir(self.trackingFolder)
        tracking_cmd = "{trackingfolder}/TrackingTelescope {root} 0 {nr}".format(trackingfolder=self.trackingFolder, root=noTracksROOTFile, nr=telescopeID)
        print "\n\nSTART TRACKING..."
        print tracking_cmd
        os.system(tracking_cmd)
        os.chdir(old_cwd)

        tracksROOTFile = self.trackingFolder+"/{prefix}{run}_withTracks.root".format(prefix=self.converterPrefix, run=str(self.run_number).zfill(4))

        # move to data folder:
        os.system("mv "+tracksROOTFile+" "+self.TrackingPadAnalysis['ROOTFile'])

        self.rootfile = ROOT.TFile(self.TrackingPadAnalysis['ROOTFile'])
        self.tree = self.rootfile.Get(self.treename) # Get TTree called "track_info"

        assert(bool(self.tree) and bool(self.rootfile)), "Could not load root file: \n\t"+self.TrackingPadAnalysis['ROOTFile']

        # delete no tracks file:
        os.system("rm "+noTracksROOTFile)

    def _LoadTiming(self):
        try:
            self.logStartTime = dt.strptime(self.RunInfo["start time"][:10]+"-"+self.RunInfo["start time"][11:-1], "%Y-%m-%d-%H:%M:%S")
            self.logStopTime = dt.strptime(self.RunInfo["stop time"][:10]+"-"+self.RunInfo["stop time"][11:-1], "%Y-%m-%d-%H:%M:%S")
            self.logRunTime = self.logStopTime - self.logStartTime
            noerror = True
        except:
            try:
                self.logStartTime = dt.strptime(self.RunInfo["start time"][:10]+"-"+self.RunInfo["start time"][11:-1], "%H:%M:%S")
                self.logStopTime = dt.strptime(self.RunInfo["stop time"][:10]+"-"+self.RunInfo["stop time"][11:-1], "%H:%M:%S")
                self.logRunTime = self.logStopTime - self.logStartTime
                noerror = True
            except:
                noerror = False
        if noerror:
            self.VerbosePrint("Timing string translated successfully")
        else:
            print "INFO: The timing information string from run info couldn't be translated"

    def RenameRunInfoKeys(self):

        try:
            for key in default_info.keys():
                tmp = self.RunInfo[key]
                del tmp
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
            #del RunInfo.runs[run_number]
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
                self._LoadROOTFile(fullROOTFilePath)

            #self.diamond.SetName(self.runs[0]['diamond'])
            return True
        else:
            return False

    def SetChannels(self, diamonds):
        '''
        Set which diamonds (channels) should be activated for the
        analysis. The argument diamonds is an integer according to:
            1 -> Diamond 1,
            2 -> Diamond 2,
            3 -> Diamond 1 & 2
        :param diamonds:
        :return:
        '''
        assert(diamonds>=1 and diamonds<=3), "invalid diamonds number: 0x1=ch0; 0x2=ch3"
        self.analyzeCh = {
            0: self._GetBit(diamonds, 0),
            3: self._GetBit(diamonds, 1)
        }

    def GetChannels(self):
        '''
        Returns a list of channels, which are activated for analysis.
        e.g. [3] means only the second diamond is activated to be
        analyzed.
        :return:
        '''
        return [i for i in self.analyzeCh.keys() if self.analyzeCh[i]]

    def GetDiamondName(self, channel):
        '''
        Returns the diamond name.
        :param channel:
        :return:
        '''
        return self.diamondname[channel]

    def ShowRunInfo(self):
        '''
        Prints the most importnant run infos to the console. The infos
        printed are:
            Run number, Rate, Diamond names, Bias Voltages
        :return:
        '''
        print "RUN INFO:"
        print "\tRun Number: \t", self.run_number, " (",self.RunInfo["type"],")"
        print "\tRate: \t", int(self.GetRate()), " kHz"
        print "\tDiamond1:   \t", self.diamondname[0], " (",self.bias[0],") | is selected: ", self.analyzeCh[0]
        print "\tDiamond2:   \t", self.diamondname[3], " (",self.bias[3],") | is selected: ", self.analyzeCh[3]

    def DrawRunInfo(self, channel=None, canvas=None, diamondinfo=True, showcut=True, comment=None, infoid="", userWidth=None, userHeight=None):
        '''
        Draws the run infos inside the canvas. If no canvas is given, it
        will be drawn into the active Pad. If the channel number is
        passed, channel number and diamond name will be drawn.
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
        if userHeight!= None: assert(userHeight>=0 and userHeight<=0.8), "choose userHeight between 0 and 0.8 or set it to 'None'"
        if userWidth!= None: assert(userWidth>=0 and userWidth<=0.8), "choose userWidth between 0 and 0.8 or set it to 'None'"
        if canvas != None:
            pad = canvas.cd()
        else:
            print "Draw run info in current pad"
            pad = ROOT.gROOT.GetSelectedPad()
            if pad:
                # canvas = pad.GetCanvas()
                # canvas.cd()
                # pad.cd()
                pass
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
            width = max(0.3, width)
        height = (lines-1)*0.03

        if not hasattr(self, "_runInfoLegends"):
            self._runInfoLegends = {}

        if channel != None and channel in [0,3]:
            # user height and width:
            userheight = height if userHeight==None else userHeight - 0.04
            userwidth = width if userWidth==None else userWidth

            self._runInfoLegends[str(channel)+infoid] = ROOT.TLegend(0.1, 0.86-userheight, 0.1+userwidth, 0.9)
            self._runInfoLegends[str(channel)+infoid].SetMargin(0.05)
            self._runInfoLegends[str(channel)+infoid].AddEntry(0, "Run{run} Ch{channel} ({rate})".format(run=self.run_number, channel=channel, rate=self._GetRateString()), "")
            if diamondinfo: self._runInfoLegends[str(channel)+infoid].AddEntry(0, "{diamond} ({bias:+}V)".format(diamond=self.diamondname[channel], bias=self.bias[channel]), "")
            if showcut and hasattr(self, "analysis"): self._runInfoLegends[str(channel)+infoid].AddEntry(0, "Cut: {cut}".format(cut=self.analysis.GetUserCutString()), "")
            if comment != None: self._runInfoLegends[str(channel)+infoid].AddEntry(0, comment, "")
            self._runInfoLegends[str(channel)+infoid].Draw("same")
        else:
            if comment != None:
                lines = 2
            else:
                lines = 1
                width = 0.15
            height = lines*0.05
            # user height and width:
            userheight = height if userHeight==None else userHeight
            userwidth = width if userWidth==None else userWidth

            self._runInfoLegends["ch12"+infoid] = ROOT.TLegend(0.1, 0.9-userheight, 0.1+userwidth, 0.9)
            self._runInfoLegends["ch12"+infoid].SetMargin(0.05)
            self._runInfoLegends["ch12"+infoid].AddEntry(0, "Run{run} ({rate})".format(run=self.run_number, rate=self._GetRateString()), "")
            if comment != None: self._runInfoLegends["ch12"+infoid].AddEntry(0, comment, "")
            self._runInfoLegends["ch12"+infoid].Draw("same")
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

    def GetChannelName(self, channel):
        self.tree.GetEntry(1)
        return self.tree.sensor_name[channel]

    def _LoadROOTFile(self, fullROOTFilePath):
        print "LOADING: ", fullROOTFilePath
        self.rootfile = ROOT.TFile(fullROOTFilePath)
        self.tree = self.rootfile.Get(self.treename) # Get TTree called "track_info"
        if not (bool(self.tree) and bool(self.rootfile)) and self.createNewROOTFiles:
            self.CreateROOTFile()

            # print "\n\nCould not load root file!"
            # print "\t>> "+fullROOTFilePath
            # answer = raw_input("generate ROOT file instead? (y/n): ")
            # if answer == "y":
            #     tracking = raw_input("generate tracking information? (y/n): ")
            #     if tracking == "y":
            #         self.CreateROOTFile()
            #     else:
            #         self.CreateROOTFile(do_tracking=False)
        #assert(bool(self.tree) and bool(self.rootfile)), "Could not load root file: \n\t"+fullROOTFilePath
