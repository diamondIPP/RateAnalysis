# ==============================================
# IMPORTS
# ==============================================
import ROOT
import json
import re
from Elementary import Elementary
from Converter import Converter
from datetime import datetime as dt
from ROOT import TFile
from ConfigParser import ConfigParser, NoOptionError
from numpy import mean
from copy import deepcopy

default_info = {
    'persons on shift': '-',
    'run info': '-',
    'type': 'signal',
    'configuration': 'signal',
    'mask': '-',
    'masked pixels': [0] * 4,
    'diamond 1': 'CH_0',
    'diamond 2': 'CH_3',
    'hv dia1': 0,
    'hv dia2': 0,
    'for1': 0,
    'for2': 0,
    'fs11': 0,
    'fsh13': 0,
    'quadrupole': '-',
    'analogue current': 0,
    'digital current': 0,
    'begin date': '2999-03-14T15:26:53Z',
    'trim time': '-:-:-',
    'config time': '-:-:-',
    'start time': '2999-03-14T15:26:53Z',
    'trig accept time': '-:-:-',
    'opening time': '-:-:-',
    'open time': '-:-:-',
    'stop time': '2999-03-14T16:26:53Z',
    'raw rate': 0,
    'prescaled rate': 0,
    'to TLU rate': 0,
    'pulser accept rate': 0,
    'cmspixel events': 0,
    'drs4 events': 0,
    'datacollector events': 0,
    'aimed flux': 0,
    'measured flux': 0,
    'user comments': '-',
    'is good run': True
}


# ==============================================
# MAIN CLASS
# ==============================================
class Run(Elementary):
    """
    Run class containing all the information for a single run from the tree and the json file.
    """

    current_run = {}
    operationmode = ''
    TrackingPadAnalysis = {}

    def __init__(self, run_number, diamonds=3, verbose=False):
        """
        :param run_number: number of the run
        :param diamonds: 0x1=ch0; 0x2=ch3
        :param verbose:
        :return:
        """
        Elementary.__init__(self, verbose=verbose)
        self.run_number = -1

        # configuration
        self.run_config_parser = self.load_parser()
        self.ShowAndWait = False
        self.filename = self.run_config_parser.get('BASIC', 'filename')
        self.treename = self.run_config_parser.get('BASIC', 'treename')
        self.run_path = self.run_config_parser.get('BASIC', 'runpath')
        self.runinfofile = self.run_config_parser.get('BASIC', 'runinfofile')
        self._runlogkeyprefix = self.run_config_parser.get('BASIC', 'runlog_key_prefix')
        self.runplaninfofile = self.run_config_parser.get('BASIC', 'runplaninfofile')
        self.maskfilepath = self.run_config_parser.get('BASIC', 'maskfilepath')
        self.createNewROOTFiles = self.run_config_parser.getboolean('BASIC', 'createNewROOTFiles')
        self.signalregion_low = self.run_config_parser.getint('BASIC', 'signalregion_low')
        self.signalregion_high = self.run_config_parser.getint('BASIC', 'signalregion_high')

        # run info
        self.allRunKeys = None
        self.RunInfo = None

        if run_number is not None:
            self.converter = Converter(self.TESTCAMPAIGN)
            print self.run_number
            assert (run_number > 0), "incorrect run_number"
            self.set_run(run_number)

            # tree info
            self.time = self.__get_time_vec()
            self.startEvent = 0
            self.endEvent = self.tree.GetEntries() - 1
            self.startTime = self.GetTimeAtEvent(self.startEvent)
            self.endTime = self.GetTimeAtEvent(self.endEvent)
            self.totalTime = self.endTime - self.startTime
            self.totalMinutes = (self.endTime - self.startTime) / 60000
            self.n_entries = self.endEvent + 1

            # region info
            self.region_information = self.load_regions()
            self.pedestal_regions = self.get_regions('pedestal')
            self.signal_regions = self.get_regions('signal')
            self.peak_integrals = self.get_peak_integrals()

        else:
            self.load_run_info()

        # extract run info
        self.channels = self.load_channels()
        self.trigger_planes = [1, 2]
        self.flux = self.calculate_flux()
        self.__load_timing()
        self.diamondname = self.__load_diamond_name()
        self.bias = self.load_bias()
        self.SetChannels(diamonds)
        self.IsMonteCarlo = False

    def load_channels(self):
        # todo think of how to extract the dia channels!
        info = self.RunInfo
        print info
        return [0, 3]

    def load_bias(self):
        bias = {}
        for i, ch in enumerate(self.channels, 1):
            bias[ch] = self.RunInfo['hv dia{num}'.format(num=i)]
        return bias

    def load_parser(self):
        parser = ConfigParser()
        parser.read("Configuration/RunConfig_" + self.TESTCAMPAIGN + ".cfg")
        return parser

    def load_regions(self):
        root_file = TFile(self.converter.get_root_file_path(self.run_number))
        macro = root_file.Get('region_information')
        return macro.GetListOfLines()

    def show_regions(self):
        for line in self.region_information:
            print line

    def get_regions(self, string):
        ranges = []
        for line in self.region_information:
            line = str(line)
            if line.startswith(string):
                data = re.split('_|:', line)
                ranges.append(data[1])
        return ranges

    def get_peak_integrals(self):
        integrals = []
        for line in self.region_information:
            line = str(line)
            if str(line).lower().startswith('* peakintegral'):
                data = re.split('_|:', line)
                if data[0][-1].isdigit():
                    integrals.append(data[0][-1])
                else:
                    integrals.append(data[1])
        return integrals

    def set_run(self, run_number, load_root_file=True):

        assert type(run_number) is int, "incorrect run_number"

        self.run_number = run_number
        self.load_run_info()

        # check for conversion
        if load_root_file:
            location = self.converter.find_root_file(run_number)
            if not location or location == 'tracking':
                self.converter.convert_run(self.RunInfo, run_number)
            self._LoadROOTFile(run_number)

        return True

    def load_run_info(self):
        self.RunInfo = {}
        data = None
        try:
            f = open(self.runinfofile, 'r')
            data = json.load(f)
            f.close()
            self.allRunKeys = deepcopy(data.keys())
            loaderror = False
        except IOError as err:
            print '\n' + (len(str(err)) + 9) * '-'
            print 'WARNING:', err
            print 'Loading default RunInfo!'
            print (len(str(err)) + 9) * '-' + '\n'
            loaderror = True

        if self.run_number >= 0:
            if not loaderror:
                run_nr = '150500' + str(self.run_number).zfill(3) if self.TESTCAMPAIGN == '201505' else str(self.run_number)
                self.RunInfo = data.get(run_nr)
                if self.RunInfo is None:
                    # try with run_log key prefix
                    self.RunInfo = data.get(self._runlogkeyprefix + str(self.run_number).zfill(3))
                if self.RunInfo is None:
                    print "INFO: Run not found in json run log file. Default run info will be used."
                    self.RunInfo = default_info
                else:
                    self.rename_runinfo_keys()
            else:
                self.RunInfo = default_info
            self.current_run = self.RunInfo
        else:
            self.RunInfo = default_info
            return 0

    def __load_diamond_name(self):
        parser = ConfigParser()
        parser.read('Configuration/DiamondAliases.cfg')
        diamondname = {}
        for i, ch in enumerate(self.channels, 1):
            diamondname[ch] = self.RunInfo['diamond {num}'.format(num=i)]
            if diamondname[ch].lower().startswith('ch'):
                continue
            try:
                diamondname[ch] = parser.get('ALIASES', diamondname[ch])
            except NoOptionError as err:
                print err
        return diamondname

    def calculate_flux(self):
        self.VerbosePrint('Calculate rate from mask file:\n\t' + self.RunInfo['mask'])
        mask_file_path = self.maskfilepath + '/' + self.RunInfo['mask']
        maskdata = {}
        for plane in self.trigger_planes:
            maskdata[plane] = {}
        try:
            f = open(mask_file_path, 'r')
            last_i2c = None
            for line in f:
                if len(line) > 3:
                    line = line.split()
                    i2c = line[1]
                    plane = self.trigger_planes[0] if last_i2c is None or i2c == last_i2c else self.trigger_planes[1]
                    maskdata[plane][line[0]] = [int(line[2]), int(line[3])]
                    last_i2c = i2c
            f.close()
        except IOError as err:
            print '\n' + (len(str(err)) + 9) * '-'
            print 'WARNING:', err
            print 'Cannot calculate flux!'
            print (len(str(err)) + 9) * '-' + '\n'
            return None

        # check for corner method
        if not maskdata.keys()[0].keys()[0].startswith('corn'):
            return None

        # fill in the information to Run Info
        masked_pixels = {}
        for plane in self.trigger_planes:
            row = [maskdata[plane]['cornBot'][0], maskdata[plane]['cornTop'][0]]
            col = [maskdata[plane]['cornBot'][1], maskdata[plane]['cornTop'][1]]
            masked_pixels[plane] = abs((row[1] - row[0] + 1) * (col[1] - col[0] + 1))
            self.RunInfo['masked pixels'][plane] = masked_pixels[plane]

        pixel_size = 0.01 * 0.015  # cm^2
        flux = []
        for i, plane in enumerate(self.channels, 1):
            area = pixel_size * masked_pixels[plane]
            flux.append(self.RunInfo['for{num}'.format(num=i)] / area / 1000)  # in kHz/cm^2
        self.RunInfo['measured flux'] = mean(flux)
        return mean(flux)

    # todo fix
    def __load_timing(self):
        try:
            self.logStartTime = dt.strptime(self.RunInfo["start time"][:10] + "-" + self.RunInfo["start time"][11:-1], "%Y-%m-%d-%H:%M:%S")
            self.logStopTime = dt.strptime(self.RunInfo["stop time"][:10] + "-" + self.RunInfo["stop time"][11:-1], "%Y-%m-%d-%H:%M:%S")
            self.logRunTime = self.logStopTime - self.logStartTime
            noerror = True
        except ValueError:
            try:
                self.logStartTime = dt.strptime(self.RunInfo["start time"][:10] + "-" + self.RunInfo["start time"][11:-1], "%H:%M:%S")
                self.logStopTime = dt.strptime(self.RunInfo["stop time"][:10] + "-" + self.RunInfo["stop time"][11:-1], "%H:%M:%S")
                self.logRunTime = self.logStopTime - self.logStartTime
                noerror = True
            except ValueError:
                noerror = False
        if noerror:
            self.VerbosePrint("Timing string translated successfully")
        else:
            print "INFO: The timing information string from run info couldn't be translated"

    def rename_runinfo_keys(self):

        # return, if all keys from default info are in RunInfo too
        if all([key in self.RunInfo for key in default_info]):
            return

        parser = ConfigParser()
        parser.read('Configuration/KeyDict_{campaign}.cfg'.format(campaign=self.TESTCAMPAIGN))
        for new_key, old_key in parser.items('KEYNAMES'):
            if old_key in self.RunInfo:
                self.RunInfo[new_key] = self.RunInfo.pop(old_key)
        return

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
        assert (diamonds >= 1 and diamonds <= 3), "invalid diamonds number: 0x1=ch0; 0x2=ch3"
        self.analyzeCh = {
            0: self._GetBit(diamonds, 0),
            3: self._GetBit(diamonds, 1)
        }

    # GET FUNCTIONS
    # ==============================================
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

    def GetFlux(self):
        fluxes = [3, 20, 60, 600, 2000, 5000]
        range_low = [0, 10, 30, 300, 1000, 2500]
        range_high = [10, 30, 300, 1000, 2500, 5000]
        if self.TESTCAMPAIGN in ['201508', '201505']:
            for i in range(len(fluxes)):
                if range_low[i] < self.RunInfo['measured flux'] < range_high[i]:
                    return fluxes[i]

    def ShowRunInfo(self):
        '''
        Prints the most importnant run infos to the console. The infos
        printed are:
            Run number, Rate, Diamond names, Bias Voltages
        :return:
        '''
        print "RUN INFO:"
        print "\tRun Number: \t", self.run_number, " (", self.RunInfo["type"], ")"
        print "\tRate: \t", int(self.GetRate()), " kHz"
        print "\tDiamond1:   \t", self.diamondname[0], " (", self.bias[0], ") | is selected: ", self.analyzeCh[0]
        print "\tDiamond2:   \t", self.diamondname[3], " (", self.bias[3], ") | is selected: ", self.analyzeCh[3]

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
        if userHeight != None: assert (userHeight >= 0 and userHeight <= 0.8), "choose userHeight between 0 and 0.8 or set it to 'None'"
        if userWidth != None: assert (userWidth >= 0 and userWidth <= 0.8), "choose userWidth between 0 and 0.8 or set it to 'None'"
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
        height = (lines - 1) * 0.03

        if not hasattr(self, "_runInfoLegends"):
            self._runInfoLegends = {}

        if channel != None and channel in [0, 3]:
            # user height and width:
            userheight = height if userHeight == None else userHeight - 0.04
            userwidth = width if userWidth == None else userWidth

            self._runInfoLegends[str(channel) + infoid] = ROOT.TLegend(0.1, 0.86 - userheight, 0.1 + userwidth, 0.9)
            self._runInfoLegends[str(channel) + infoid].SetMargin(0.05)
            self._runInfoLegends[str(channel) + infoid].AddEntry(0, "Run{run} Ch{channel} ({rate})".format(run=self.run_number, channel=channel, rate=self.GetRateString()), "")
            if diamondinfo: self._runInfoLegends[str(channel) + infoid].AddEntry(0, "{diamond} ({bias:+}V)".format(diamond=self.diamondname[channel], bias=self.bias[channel]), "")
            if showcut and hasattr(self, "analysis"): self._runInfoLegends[str(channel) + infoid].AddEntry(0, "Cut: {cut}".format(cut=self.analysis.GetUserCutString()), "")
            if comment != None: self._runInfoLegends[str(channel) + infoid].AddEntry(0, comment, "")
            self._runInfoLegends[str(channel) + infoid].Draw("same")
        else:
            if comment != None:
                lines = 2
            else:
                lines = 1
                width = 0.15
            height = lines * 0.05
            # user height and width:
            userheight = height if userHeight == None else userHeight
            userwidth = width if userWidth == None else userWidth

            self._runInfoLegends["ch12" + infoid] = ROOT.TLegend(0.1, 0.9 - userheight, 0.1 + userwidth, 0.9)
            self._runInfoLegends["ch12" + infoid].SetMargin(0.05)
            self._runInfoLegends["ch12" + infoid].AddEntry(0, "Run{run} ({rate})".format(run=self.run_number, rate=self.GetRateString()), "")
            if comment != None: self._runInfoLegends["ch12" + infoid].AddEntry(0, comment, "")
            self._runInfoLegends["ch12" + infoid].Draw("same")
        pad.Modified()

    def GetRateString(self):
        rate = self.RunInfo["measured flux"]
        if rate > 1000:
            unit = "MHz"
            rate = round(rate / 1000., 1)
        else:
            unit = "kHz"
            rate = int(round(rate, 0))
        return "{rate:>3}{unit}".format(rate=rate, unit=unit)

    def GetChannelName(self, channel):
        self.tree.GetEntry(1)
        return self.tree.sensor_name[channel]

    def _LoadROOTFile(self, run_number):
        file_path = self.converter.get_tracking_file_path(run_number) if self.converter.do_tracking else self.converter.get_root_file_path(run_number)
        print "\nLoading infos for rootfile: ", file_path.split('/')[-1]
        self.rootfile = ROOT.TFile(file_path)
        self.tree = self.rootfile.Get(self.treename)  # Get TTree called "track_info"

    def __get_time_vec(self):
        self.tree.SetEstimate(-1)
        entries = self.tree.Draw('Entry$:time', '', 'goff')
        time = []
        for i in xrange(entries):
            time.append(self.tree.GetV2()[i])
        return time

    def GetTimeAtEvent(self, event):
        """
        Returns the time stamp at event number 'event'. For negative event numbers it will return the time stamp at the startevent.
        :param event: integer event number
        :return: timestamp for event
        """
        if event == -1:
            event = self.endEvent
        elif event < 0:
            event = 0
        elif event >= self.endEvent:
            event = self.endEvent
        return self.time[event]

    def GetEventAtTime(self, time_sec):
        """
        Returns the eventnunmber at time dt from beginning of the run. Accuracy: +- 1 Event
        :param time_sec: time in seconds from start
        :return: event_number
        """
        # return last time if input is too large
        if time_sec > self.time[-1] / 1000. or time_sec == -1:
            return self.time[-1]
        last_time = 0
        offset = self.time[0] / 1000.
        for i, time in enumerate(self.time):
            time /= 1000.
            if time >= time_sec + offset >= last_time:
                return i
            last_time = time

    def _slope_f(self, event):
        if event < 0:
            event = 0
        time_high = self.GetTimeAtEvent(event + 10)
        time_low = self.GetTimeAtEvent(event - 10)
        return 1. * (time_high - time_low) / 21.


if __name__ == "__main__":
    z = Run(None, validate=True)
