#!/usr/bin/env python
# ==============================================
# IMPORTS
# ==============================================
from json import load, loads, dumps
from ConfigParser import ConfigParser, NoOptionError
from ROOT import TFile, TTree
from argparse import ArgumentParser
from StringIO import StringIO

from Converter import Converter
from Elementary import Elementary
from Utils import *


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

    def __init__(self, run_number=None, test_campaign=None, tree=None, t_vec=None, verbose=False):
        """
        :param run_number: number of the run
        :param verbose:
        :return:
        """
        self.RunNumber = run_number
        Elementary.__init__(self, testcampaign=test_campaign, verbose=verbose)

        # configuration
        self.DUTType = self.load_dut_type()
        self.NChannels = self.load_pre_n_channels()
        self.trigger_planes = [1, 2]
        self.treename = self.run_config_parser.get('BASIC', 'treename')
        self.runinfofile = self.load_run_info_path()
        self.IrradiationFile = self.MainConfigParser.get('MISC', 'irradiation_file')
        self.maskfilepath = self.load_mask_file_dir()
        self.createNewROOTFiles = self.run_config_parser.getboolean('BASIC', 'createNewROOTFiles')
        self.Digitiser = self.run_config_parser.get('BASIC', 'digitizer') if self.DUTType == 'pad' else None

        # run info
        self.DefaultInfo = self.load_default_info()
        self.RunInfo = self.load_run_info()
        self.RootFile = None
        self.tree = None

        # general information
        self.FoundForRate = False
        self.Channels = self.load_channels()
        self.Flux = self.calculate_flux()

        # times
        self.LogStart = self.load_log_start()
        self.LogEnd = self.load_log_stop()
        self.Duration = self.LogEnd - self.LogStart
        self.LogHalfTime = time_stamp(self.LogStart + self.Duration / 2)

        # extract run info
        self.DiamondNames = self.load_diamond_names()
        self.Bias = self.load_biases()

        self.Converter = Converter(self)
        if self.set_run(run_number, tree):
            # tree info
            self.time = get_time_vec(self.tree) if t_vec is None else t_vec
            self.startEvent = 0
            self.n_entries = int(self.tree.GetEntries())
            self.endEvent = self.n_entries - 1
            self.StartTime = self.load_start_time()
            self.EndTime = self.load_end_time()
            self.TotalTime = self.load_total_time()
            self.totalMinutes = self.TotalTime / 60000.
            self.Duration = timedelta(seconds=self.TotalTime)
            self.NPlanes = self.load_n_planes()

            self.Channels = self.load_channels()

            # region info
            if self.DUTType == 'pad':
                self.TreeConfig = self.load_tree_config()
                self.DigitizerChannels = self.load_digitizer_channels()
                self.NChannels = len(self.DigitizerChannels)
                self.IntegralRegions = self.load_regions()
                self.PeakIntegrals = self.load_peak_integrals()
                self.TCal = self.load_tcal()
                self.NSamples = len(self.TCal)

    # ==============================================
    # region LOAD FUNCTIONS

    def load_digitizer_channels(self):
        return loads(self.TreeConfig.get('Sensor Names', 'Names'))

    def load_pre_n_channels(self):
        if self.DUTType == 'pad':
            return 9 if self.run_config_parser.get('BASIC', 'digitizer').lower() == 'caen' else 4
        else:
            return None

    def get_n_diamonds(self, run_number=None):
        run_info = self.load_run_info(run_number)
        return len([key for key in run_info if key.startswith('dia') and key[-1].isdigit()])

    def load_channels(self):
        if self.DUTType == 'pad':
            if hasattr(self, 'TreeConfig'):
                return [i for i, val in enumerate(self.TreeConfig.get('General', 'active regions').split()) if int(val)]
            # read it from the current config file if there is no tree loaded
            binary = self.run_config_parser.getint('ROOTFILE_GENERATION', 'active_regions')
            return [i for i in xrange(self.NChannels) if has_bit(binary, i)]
        elif self.DUTType == 'pixel':
            return [i for i in xrange(len([key for key in self.RunInfo.iterkeys() if key.startswith('dia') and key[-1].isdigit()]))]

    def load_dut_type(self):
        dut_type = self.run_config_parser.get('BASIC', 'type') if self.RunNumber is not None else None
        if dut_type not in ['pixel', 'pad', None]:
            log_critical("The DUT type {0} has to be either 'pixel' or 'pad'".format(dut_type))
        return dut_type

    def load_tree_config(self):
        macro = self.RootFile.Get('region_information')
        info = '\n'.join(str(word) for word in macro.GetListOfLines())
        config = ConfigParser()
        if not info.startswith('['):
            info = []
            for word in macro.GetListOfLines():
                word = str(word).strip('\n *:')
                if word.startswith('active'):
                    info.append('[General]')
                    data = word.replace('_', ' ').split(':')
                    word = '{0} = {1}'.format(data[0], ' '.join(bin(int(data[1]))[2:]))
                elif word.startswith('Signal') or word.startswith('Sensor'):
                    word = '[{}]'.format(word)
                elif word.startswith('tcal'):
                    info.append('[Time Calibration]')
                    word = word.replace('tcal', 'tcal =').replace(', \b\b', '')
                elif word and word[-1].isdigit() and 'pulser,' not in word.lower():
                    data = word.split(':')
                    word = '{0} = {1}'.format(data[0], str([int(num) for num in data[1].split('-')]))
                elif 'pulser' in word.lower():
                    word = 'Names = {}'.format(dumps(word.split(',')))
                info.append(word)
            info = '\n'.join(info)
        config.readfp(StringIO(info))
        return config

    def load_ranges(self, section):
        ranges = []
        for i, channel in enumerate(self.Channels):
            ranges.append(OrderedDict())
            this_section = '{} {}'.format(section, channel)
            if any('Signal' in sec for sec in self.TreeConfig.sections()):
                this_section = 'Signal Windows' if 'Region' in section else 'Signal Definitions'
            for option in self.TreeConfig.options(this_section):
                ranges[i][option.replace('peakin', 'PeakIn')] = loads(self.TreeConfig.get(this_section, option))
        return ranges

    def load_regions(self):
        return self.load_ranges('Integral Regions')

    def load_peak_integrals(self):
        return self.load_ranges('Integral Ranges')

    def load_tcal(self):
        tcal = loads(self.TreeConfig.get('Time Calibration', 'tcal').replace('nan', '0'))
        if len(tcal) < 1024:
            tcal.append(2 * tcal[-1] - tcal[-2])
        return tcal[:1024]

    def load_default_info(self):
        f = open(join(self.Dir, 'Runinfos', 'defaultInfo.json'), 'r')
        data = load(f)
        f.close()
        return data

    def load_run_info(self, run_number=None):
        data = self.load_run_info_file()

        run_number = self.RunNumber if run_number is None else run_number
        if run_number >= 0:
            run_info = data.get(str(run_number))
            # try with run_log key prefix if the number was not found
            if run_info is None:
                run_info = data.get('{tc}{run}'.format(tc=self.TESTCAMPAIGN[2:], run=str(run_number).zfill(5)))
            # abort if the run is still not found
            if run_info is None:
                log_critical('Run {} not found in json run log file!'.format(run_number))
            self.RunInfo = run_info
            self.add_default_entries()
            self.translate_diamond_names()
            return run_info
        else:
            self.RunInfo = self.DefaultInfo
            return self.DefaultInfo

    def add_default_entries(self):
        self.RunInfo['masked pixels'] = [0] * 4

    def translate_diamond_names(self):
        for key in self.RunInfo:
            if not (key.startswith('dia') and key[-1].isdigit()):
                continue
            dia = self.RunInfo[key]
            try:
                self.RunInfo[key] = self.translate_dia(dia)
            except NoOptionError as err:
                log_warning(err)

    def translate_dia(self, dia):
        if dia is None or dia.lower() in ['unknown', 'none']:
            return
        parser = ConfigParser()
        parser.read(join(self.Dir, 'Configuration', 'DiamondAliases.cfg'))
        return parser.get('ALIASES', dia.lower())

    def load_run_info_file(self):
        try:
            f = open(self.runinfofile, 'r')
            data = load(f)
            f.close()
            return data
        except IOError as err:
            # the run log file is required to get any meaningful information
            log_critical('{err}\nCould not load default RunInfo! --> Using default'.format(err=err))

    def load_diamond_names(self):
        return [self.RunInfo['dia{nr}'.format(nr=i)] for i in xrange(1, self.get_n_diamonds() + 1)]

    def load_biases(self):
        return [int(self.RunInfo['dia{nr}hv'.format(nr=i)]) for i in xrange(1, self.get_n_diamonds() + 1)]

    def load_log_start(self):
        return conv_log_time(self.RunInfo['starttime0'])

    def load_log_stop(self):
        return conv_log_time(self.RunInfo['endtime'])

    def load_start_time(self):
        return int(round(self.get_time_at_event(self.startEvent)))

    def load_end_time(self):
        return int(round(self.get_time_at_event(self.endEvent)))

    def load_total_time(self):
        return (self.time[-1] - self.time[0]) / 1000

    def load_n_planes(self):
        if self.has_branch('cluster_col'):
            self.tree.Draw('@cluster_col.size()', '', 'goff', 1)
            return int(self.tree.GetV1()[0])
        else:
            return 4
    # endregion INIT

    def get_time(self):
        return make_ufloat((self.LogHalfTime, self.Duration.seconds / 2))

    def print_run_info(self):
        print 'Run information for run', self.RunNumber
        for key, value in sorted(self.RunInfo.iteritems()):
            print '{k}: {v}'.format(k=key.ljust(13), v=value)

    def reload_run_config(self, run_number):
        self.RunNumber = run_number
        self.run_config_parser = self.load_run_config()

    def set_run(self, run_number, root_tree):

        if run_number is None:
            return False
        assert (run_number > 0), 'incorrect run_number'
        assert type(run_number) is int, 'incorrect run_number'

        self.RunNumber = run_number
        self.load_run_info()
        self.Flux = self.calculate_flux()

        # check for conversion
        if root_tree is None:
            self.Converter.convert_run()
            self.__load_rootfile()
        elif root_tree:
            # TODO: check if the tree is functional and type tree...
            self.RootFile = root_tree[0]
            self.tree = root_tree[1]
        else:
            return False
        self.validate_root_files()
        return True

    def validate_root_files(self):
        if not (self.RootFile.IsA() == TFile().IsA() and self.tree.IsA() == TTree().IsA()):
            log_critical('either TFile or TTree of wrong type')
        if self.DUTType == 'pad' and 'region_information' not in [key.GetName() for key in self.RootFile.GetListOfKeys()]:
            self.Converter.remove_final_file()
            log_critical('no region information in root tree (file not propely converted)')

    def load_mask(self):
        mask_file_path = '{path}/{mask}'.format(path=self.maskfilepath, mask=self.RunInfo['maskfile'])
        dic = {}
        try:
            f = open(mask_file_path, 'r')
            for line in f:
                if len(line) > 3:
                    line = line.split()
                    if not line[1] in dic:
                        dic[line[1]] = {}
                        dic[line[1]]['row'] = [0, 0]
                        dic[line[1]]['col'] = [0, 0]
                    if line[0] == 'cornBot':
                        dic[line[1]]['row'][0] = line[3]
                        dic[line[1]]['col'][0] = line[2]
                    elif line[0] == 'cornTop':
                        dic[line[1]]['row'][1] = line[3]
                        dic[line[1]]['col'][1] = line[2]
            f.close()
        except IOError as err:
            log_warning(err)
        print dic

    def calculate_flux(self):

        if self.RunNumber is None:
            return
        fluxes = []
        self.find_for_in_comment()
        if 'for1' in self.RunInfo:
            self.FoundForRate = True
            for i, area in enumerate(self.get_unmasked_area().itervalues(), 1):
                fluxes.append(self.RunInfo['for{num}'.format(num=i)] / area / 1000)  # in kHz/cm^2
        else:
            fluxes.append(self.RunInfo['measuredflux'])
        flux = mean(fluxes)
        self.RunInfo['mean flux'] = flux
        return make_ufloat((flux, .1 * flux))

    def get_unmasked_area(self):
        if self.RunNumber is None:
            return
        mask_file = join(self.maskfilepath, basename(self.RunInfo['maskfile'])).strip("\"")
        maskdata = {}
        if basename(self.RunInfo['maskfile']).lower() in ['no mask', 'none', 'none!']:
            pass
        else:
            try:
                f = open(mask_file, 'r')
                for line in f:
                    if line.startswith('#'):
                        continue
                    if len(line) > 3:
                        line = line.split()
                        roc = int(line[1])
                        if roc not in maskdata:
                            maskdata[roc] = {}
                        maskdata[roc][line[0]] = (int(line[2]), int(line[3]))
                f.close()
            except IOError:
                log_warning('Could not find mask file {f}! Not taking any mask!'.format(f=mask_file))
        # default pixels
        unmasked_pixels = {plane: 52 * 80 for plane in maskdata} if maskdata else {0: 52 * 80, 1: 52 * 80}
        # pass for empty file
        if not maskdata:
            pass
        # check for corner method
        elif not maskdata.values()[0].keys()[0].startswith('corn'):
            log_warning('Invalid mask file: "{}". Not taking any mask!'.format(mask_file))
        else:
            for plane, dic in maskdata.iteritems():
                row = dic['cornBot'][0], dic['cornTop'][0]
                col = dic['cornBot'][1], dic['cornTop'][1]
                unmasked_pixels[plane] = abs((row[1] - row[0] + 1) * (col[1] - col[0] + 1))
                self.RunInfo['masked pixels'][plane] = unmasked_pixels[plane]

        pixel_size = 0.01 * 0.015  # cm^2
        return {key: pixels * pixel_size for key, pixels in unmasked_pixels.iteritems()}

    def find_for_in_comment(self):
        for name in ['for1', 'for2']:
            if name not in self.RunInfo:
                for cmt in self.RunInfo['comments'].split('\r\n'):
                    cmt = cmt.replace(':', '')
                    cmt = cmt.split(' ')
                    if str(cmt[0].lower()) == name:
                        self.RunInfo[name] = int(cmt[1])

    def wf_exists(self, channel):
        wf_exist = True if self.tree.FindBranch('wf{ch}'.format(ch=channel)) else False
        if not wf_exist:
            log_warning('The waveform for channel {ch} is not stored in the tree'.format(ch=channel))
        return wf_exist

    # ==============================================
    # region GET FUNCTIONS
    def get_flux(self):
        return self.Flux if self.Flux else self.RunInfo['aimed flux']

    def get_irradiations(self):
        with open(self.IrradiationFile, 'r') as f:
            data = load(f)
        for dia in self.DiamondNames:
            if dia not in data[self.TCString].keys() + ['None']:
                log_critical('Please add "{}" to the irradiation file for {}'.format(dia, self.TCString))
        return [data[self.TCString][dia] for dia in self.DiamondNames if dia not in ['None']]

    def get_attenuators(self):
        return [str(self.RunInfo['att_dia{}'.format(i)]) for i in xrange(1, self.get_n_diamonds() + 1) if 'att_dia1' in self.RunInfo]

    def get_calibrated_times(self, trigger_cell):
        t = [self.TCal[int(trigger_cell)]]
        for i in xrange(1, self.NSamples):
            t.append(self.TCal[(int(trigger_cell) + i) % self.NSamples] + t[-1])
        return t

    def get_calibrated_time(self, trigger_cell, ibin):
        v = self.TCal[int(trigger_cell)]
        for i in xrange(ibin):
            v += self.TCal[(int(trigger_cell) + i + 1) % self.NSamples]
        return v

    def get_diamond_name(self, channel):
        return self.DiamondNames[channel]

    def get_channel_name(self, channel):
        self.tree.GetEntry()
        return self.tree.sensor_name[channel]

    def get_rate_string(self):
        rate = self.Flux
        unit = 'MHz/cm^{2}' if rate > 1000 else 'kHz/cm^{2}'
        rate = round(float(rate / 1000.), 1) if rate > 1000 else int(round(float(rate), 0))
        return '{rate:>3} {unit}'.format(rate=rate, unit=unit)

    def get_time_at_event(self, event):
        """ For negative event numbers it will return the time stamp at the startevent. """
        event = self.endEvent if event == -1 else event
        event = sorted((0, event, self.endEvent))[1]  # guarantees that the event is in [0, EndEvent]
        return self.time[event] / 1000.

    def get_event_at_time(self, seconds):
        """ Returns the event nunmber at time dt from beginning of the run. Accuracy: +- 1 Event """
        # return time of last event if input is too large
        if seconds >= self.TotalTime or seconds == -1:
            return self.n_entries
        last_time = 0
        for i, t in enumerate(self.time):
            if t / 1000. >= seconds + self.time[0] / 1000 >= last_time:
                return i
            last_time = t / 1000.
        return self.n_entries  # if loop doesn't find anything it has to be the last bin

    # endregion

    # ==============================================
    # region SHOW RUN INFO
    def show_tree_config(self):
        macro = self.RootFile.Get('region_information')
        for line in macro.GetListOfLines():
            line = str(line)
            if not (line.startswith('[Time') or line.startswith('tcal')):
                print line

    def has_branch(self, name):
        return bool(self.tree.GetBranch(name))

    # endregion

    def __load_rootfile(self):
        file_path = self.Converter.get_final_file_path()
        print '\033[1A\rLoading information for rootfile: {file}'.format(file=file_path.split('/')[-1])
        self.RootFile = TFile(file_path)
        self.tree = self.RootFile.Get(self.treename)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('run', nargs='?', default=392, type=int)
    p.add_argument('-tc', '--testcampaign', nargs='?', default=None)
    p.add_argument('-t', '--tree', action='store_true')
    args = p.parse_args()
    z = Run(args.run, tree=None if args.tree else False, test_campaign=args.testcampaign)
