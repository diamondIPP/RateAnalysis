# ==============================================
# IMPORTS
# ==============================================
import json
import re
from Elementary import Elementary
from Converter import Converter
from datetime import datetime
from ROOT import TFile, gROOT, TLegend
from ConfigParser import ConfigParser, NoOptionError
from numpy import mean
from collections import OrderedDict
from subprocess import check_output

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
    'mean flux': None,
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
        self.run_number = None

        # configuration
        self.channels = [0, 3]
        self.trigger_planes = [1, 2]
        self.run_config_parser = self.load_parser()
        self.DUTType = self.load_dut_type()
        self.filename = self.run_config_parser.get('BASIC', 'filename')
        self.treename = self.run_config_parser.get('BASIC', 'treename')
        self.run_path = self.run_config_parser.get('BASIC', 'runpath')
        self.runinfofile = self.run_config_parser.get('BASIC', 'runinfofile')
        self._runlogkeyprefix = self.run_config_parser.get('BASIC', 'runlog_key_prefix')
        self.maskfilepath = self.run_config_parser.get('BASIC', 'maskfilepath')
        self.createNewROOTFiles = self.run_config_parser.getboolean('BASIC', 'createNewROOTFiles')
        self.signalregion_low = self.run_config_parser.getint('BASIC', 'signalregion_low')
        self.signalregion_high = self.run_config_parser.getint('BASIC', 'signalregion_high')
        
        # run info
        self.RunInfo = None

        if run_number is not None:
            self.converter = Converter(self.TESTCAMPAIGN, self.DUTType)
            assert (run_number > 0), 'incorrect run_number'
            self.set_run(run_number)

            # tree info
            self.time = self.__get_time_vec()
            self.startEvent = 0
            self.endEvent = self.tree.GetEntries() - 1
            self.startTime = self.get_time_at_event(self.startEvent)
            self.endTime = self.get_time_at_event(self.endEvent)
            self.totalTime = self.endTime - self.startTime
            self.totalMinutes = (self.endTime - self.startTime) / 60000
            self.n_entries = int(self.endEvent + 1)
            
            # region info
            if self.DUTType == 'pad':
                self.region_information = self.load_regions()
                self.pedestal_regions = self.get_regions('pedestal')
                self.signal_regions = self.get_regions('signal')
                self.peak_integrals = self.get_peak_integrals()

            self.flux = self.calculate_flux()    

        else:
            self.load_run_info()

        # extract run info
        self.analyse_ch = self.set_channels(diamonds)
        self.diamond_names = self.__load_diamond_name()
        self.bias = self.load_bias()
        self.IsMonteCarlo = False

        # times
        self.log_start = None
        self.log_stop = None
        self.duration = None
        self.__load_timing()
        # root objects
        self.run_info_legends = {}

    # ==============================================
    # region LOAD FUNCTIONS
    def load_bias(self):
        bias = {}
        for i, ch in enumerate(self.channels, 1):
            bias[ch] = self.RunInfo['hv dia{num}'.format(num=i)]
        return bias
    
    def load_dut_type(self):
        _type = self.run_config_parser.get('BASIC', 'type')
        assert _type.lower() in ["pixel", "pad"], "The DUT type {0} should be 'pixel' or 'pad'".format(_type)
        return _type

    def load_parser(self):
        parser = ConfigParser()
        parser.read("Configuration/RunConfig_" + self.TESTCAMPAIGN + ".cfg")
        return parser

    def load_regions(self):
        macro = self.rootfile.Get('region_information')
        return macro.GetListOfLines()

    def load_run_info(self):
        self.RunInfo = {}
        data = None
        try:
            f = open(self.runinfofile, 'r')
            data = json.load(f)
            f.close()
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

    def __load_timing(self):
        try:
            self.log_start = datetime.strptime(self.RunInfo['start time'], "%Y-%m-%dT%H:%M:%SZ")
            self.log_stop = datetime.strptime(self.RunInfo['stop time'], "%Y-%m-%dT%H:%M:%SZ")
        except ValueError:
            try:
                self.log_start = datetime.strptime(self.RunInfo['start time'], "%H:%M:%S")
                self.log_stop = datetime.strptime(self.RunInfo['stop time'], "%H:%M:%S")
            except ValueError as err:
                print err
                return
        self.duration = self.log_stop - self.log_start
    # endregion

    def set_run(self, run_number, load_root_file=True):

        assert type(run_number) is int, "incorrect run_number"

        self.run_number = run_number
        self.load_run_info()

        # check for conversion
        if load_root_file:
            self.converter.convert_run(self.RunInfo, run_number)
            self.__load_rootfile()

        return True

    def set_channels(self, diamonds):
        """
        Set which diamonds (channels) should be activated for the analysis.
        :param diamonds: bitwise integer (1=dia1, 2=dia2, 3=1&2)
        """
        assert 1 <= diamonds <= 3, 'invalid diamonds number: 0x1=ch0; 0x2=ch3'
        analyse_ch = {}
        for i, ch in enumerate(self.channels):
            analyse_ch[ch] = self.has_bit(diamonds, i)
        self.analyse_ch = analyse_ch
        return analyse_ch

    def calculate_flux(self):
        self.verbose_print('Calculate rate from mask file:\n\t' + self.RunInfo['mask'])
        mask_file_path = self.maskfilepath + '/' + self.RunInfo['mask']
        maskdata = {}
        for plane in self.trigger_planes:
            maskdata[plane] = {}
        try:
            f = open(mask_file_path, 'r')
            i2cs = []
            for line in f:
                if len(line) > 3:
                    line = line.split()
                    if not i2cs or i2cs[-1] != line[1]:
                        i2cs.append(line[1])
                    plane = self.trigger_planes[len(i2cs) - 1]
                    maskdata[plane][line[0]] = [int(line[2]), int(line[3])]
            f.close()
        except IOError as err:
            print '\n' + (len(str(err)) + 9) * '-'
            print 'WARNING:', err
            print 'Cannot calculate flux!'
            print (len(str(err)) + 9) * '-' + '\n'
            return None

        # check for corner method
        if not maskdata.values()[0].keys()[0].startswith('corn'):
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
        for i, plane in enumerate(self.trigger_planes, 1):
            area = pixel_size * masked_pixels[plane]
            flux.append(self.RunInfo['for{num}'.format(num=i)] / area / 1000)  # in kHz/cm^2
        self.RunInfo['mean flux'] = mean(flux)
        return mean(flux)

    def rename_runinfo_keys(self):

        # return, if all keys from default info are in RunInfo too
        if all([key in self.RunInfo for key in default_info]):
            return

        parser = ConfigParser()
        parser.read('Configuration/KeyDict_{campaign}.cfg'.format(campaign=self.TESTCAMPAIGN))
        for new_key, old_key in parser.items('KEYNAMES'):
            if old_key in self.RunInfo:
                self.RunInfo[new_key] = self.RunInfo.pop(old_key)
            else:
                self.RunInfo[new_key] = default_info[new_key]
        return

    def wf_exists(self, channel):
        wf_exist = True if self.tree.FindBranch('wf{ch}'.format(ch=channel)) else False
        if not wf_exist:
            print 'The waveform for channel {ch} is not stored in the tree'.format(ch=channel)
        return wf_exist

    # ==============================================
    # region GET FUNCTIONS
    def get_flux(self):
        return self.flux if self.flux else self.RunInfo['aimed flux']

    def get_regions(self, string):
        ranges = OrderedDict()
        for line in self.region_information:
            line = str(line)
            if line.startswith(string):
                data = re.split('_|:|-', line)
                data = [data[i].strip(' ') for i in range(len(data))]
                try:
                    ranges[data[1]] = [int(data[2]), int(data[3])]
                except IndexError:
                    ranges[data[0]] = [int(data[i]) for i in [1, 2]]
        return ranges

    def get_peak_integrals(self):
        integrals = OrderedDict()
        for line in self.region_information:
            line = str(line)
            if str(line).lower().startswith('* peakintegral'):
                data = re.split('_|:|-', line)
                if data[0][-1].isdigit():
                    if data[0][-2].isdigit():
                        integrals[data[0][-2:]] = [int(float(data[i])) for i in [1, 2]]
                    else:
                        integrals[data[0][-1]] = [int(float(data[i])) for i in [1, 2]]
                else:
                    integrals[data[1]] = [int(float(data[i])) for i in [2, 3]]
        return integrals

    def get_active_channels(self):
        """
        Returns a list of the channels, which are activated for analysis. e.g. [3] means only the channel 3 is activated for analysis.
        :return:
        """
        return [ch for ch in self.analyse_ch if self.analyse_ch[ch]]

    def get_diamond_name(self, channel):
        return self.diamond_names[channel]

    def get_channel_name(self, channel):
        self.tree.GetEntry()
        return self.tree.sensor_name[channel]

    def get_rate_string(self):
        rate = self.flux
        unit = 'MHz/cm2' if rate > 1000 else 'kHz/cm2'
        rate = round(rate / 1000., 1) if rate > 1000 else int(round(rate, 0))
        return '{rate:>3} {unit}'.format(rate=rate, unit=unit)

    def __get_time_vec(self):
        self.tree.SetEstimate(-1)
        entries = self.tree.Draw('Entry$:time', '', 'goff')
        time = []
        for i in xrange(entries):
            time.append(self.tree.GetV2()[i])
        return time

    def get_time_at_event(self, event):
        """
        For negative event numbers it will return the time stamp at the startevent.
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

    def get_event_at_time(self, time_sec):
        """
        Returns the eventnunmber at time dt from beginning of the run. Accuracy: +- 1 Event
        :param time_sec: time in seconds from start
        :return: event_number
        """
        # return time of last event if input is too large
        offset = self.time[0] / 1000.
        if time_sec > self.time[-1] / 1000. - offset or time_sec == -1:
            return self.n_entries
        last_time = 0
        for i, time in enumerate(self.time):
            time /= 1000.
            if time >= time_sec + offset >= last_time:
                return i
            last_time = time

    # endregion

    # ==============================================
    # region SHOW RUN INFO
    def show_regions(self):
        for line in self.region_information:
            print line

    def show_run_info(self):
        """
        Prints the most importnant run infos to the console. The infos printed are: Run number, Rate, Diamond names, Bias Voltages
        """
        print 'RUN INFO:'
        print '\tRun Number: \t', self.run_number, ' (', self.RunInfo['type'], ')'
        print '\tRate: \t', self.get_flux(), ' kHz'
        print '\tDiamond1:   \t', self.diamond_names[0], ' (', self.bias[0], ') | is selected: ', self.analyse_ch[0]
        print '\tDiamond2:   \t', self.diamond_names[3], ' (', self.bias[3], ') | is selected: ', self.analyse_ch[3]

    def draw_run_info(self, channel=None, canvas=None, diamondinfo=True, cut=None, comment=None, set_width=None, set_height=None, runs=None):
        """
        Draws the run infos inside the canvas. If no canvas is given, it will be drawn into the active Pad. 
        If the channel number is passed, channel number and diamond name will be drawn.
        :param channel:
        :param canvas:
        :param diamondinfo:
        :param cut:
        :param comment:
        :param runs:
        :param set_height:
        :param set_width:
        :return:
        """
        assert channel is None or channel in self.channels, 'wrong channel id "{ch}"'.format(ch=channel)
        if set_height is not None:
            assert 0 <= set_height <= 0.8, 'choose height between 0 and 0.8 or set it to "None"'
        if set_width is not None:
            assert 0 <= set_width <= 0.8, 'choose width between 0 and 0.8 or set it to "None"'
        if canvas is not None:
            pad = canvas.cd()
        else:
            print 'Draw run info in current pad'
            pad = gROOT.GetSelectedPad()
            if not pad:
                print 'ERROR: Cannot access active Pad'
                return

        lines = 2
        width = 0.4
        if diamondinfo:
            lines += 1
        if cut and hasattr(self, 'analysis'):
            lines += 1
            width = 0.6
        if comment is not None:
            lines += 1
            width = max(0.5, width)
        height = (lines - 1) * 0.03

        tc = datetime.strptime(self.TESTCAMPAIGN, '%Y%m')
        dur = '{0:02d}:{1:02.0f}'.format(int(self.totalMinutes), (self.totalMinutes - int(self.totalMinutes)) * 60)

        if not canvas.GetBottomMargin() > .105:
            canvas.SetBottomMargin(0.15)
        # user height and width:
        userheight = height if set_height is None else set_height - 0.04
        userwidth = width if set_width is None else set_width

        git_text = TLegend(.85, 0, 1, .025)
        git_text.AddEntry(0, 'git hash: {ver}'.format(ver=check_output(['git', 'describe', '--always'])), '')
        git_text.SetLineColor(0)
        legend = TLegend(.002, .00205, userwidth, userheight + 0.04)
        legend.SetName('l')
        legend.SetMargin(0.05)
        legend.AddEntry(0, 'Test Campaign: {tc}'.format(tc=tc.strftime('%b %Y')), '')
        if runs is None:
            legend.AddEntry(0, 'Run {run}: {rate}, {dur} Min ({evts} evts)'.format(run=self.run_number, rate=self.get_rate_string(), dur=dur, evts=self.n_entries), '')
        else:
            legend.AddEntry(0, 'Runs {start}-{stop} ({flux1} - {flux2})'.format(start=runs[0], stop=runs[1], flux1=runs[2].strip(' '), flux2=runs[3].strip(' ')), '')
        if channel is None:
            dias = ['{dia} @ {bias:2.0f}V'.format(dia=self.diamond_names[ch], bias=self.bias[ch]) for ch in self.channels]
            dias = str(dias).strip('[]').replace('\'', '')
            legend.AddEntry(0, 'Diamonds: {dias}'.format(dias=dias), '')
        else:
            legend.AddEntry(0, 'Diamond: {diamond} @ {bias:+}V'.format(diamond=self.diamond_names[channel], bias=self.bias[channel]), '')
        if cut and hasattr(self, 'analysis'):
            legend.AddEntry(0, 'Cut: {cut}'.format(cut=self.analysis.get_easy_cutstring()), '')
        if comment is not None:
            legend.AddEntry(0, comment, '')
        # while canvas.GetPad(n_pads + 1):
        #     n_pads += 1
        pads = [i for i in canvas.GetListOfPrimitives() if i.IsA().GetName() == 'TPad']
        if not pads:
            git_text.Draw()
            legend.Draw()
        else:
            for pad in pads:
                pad.cd()
                git_text.Draw()
                legend.Draw()
        self.run_info_legends[0] = [legend, git_text]
        pad.Modified()
        canvas.Update()

    # endregion

    def __load_rootfile(self):
        file_path = self.converter.get_final_file_path(self.run_number)
        print "\nLoading information for rootfile: ", file_path.split('/')[-1]
        self.rootfile = TFile(file_path)
        self.tree = self.rootfile.Get(self.treename)  # Get TTree called "track_info"


if __name__ == "__main__":
    z = Run(398)
