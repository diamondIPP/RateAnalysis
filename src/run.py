#!/usr/bin/env python
from ConfigParser import ConfigParser, NoOptionError
from json import load, loads

from ROOT import TFile
from numpy import inf, zeros

from converter import Converter
from analysis import Analysis, join, basename
from utils import *
from glob import glob


class Run:
    """ Run class containing all the information for a single run. """

    def __init__(self, run_number=None, test_campaign=None, tree=True, t_vec=None, verbose=False):
        """
        :param run_number: if None is provided it creates a dummy run
        :param test_campaign: if None is provided ...
        :param tree: root_tree object, if None is given it will start the converter
        :param t_vec: time sequence of the run, if None is provide it will generate a corrected one
        :param verbose: turn on more output
        """
        # Basics
        self.RunNumber = run_number
        self.Verbose = verbose
        self.MainConfig = Analysis.load_main_config()

        # Directories / Test Campaign
        self.Dir = get_base_dir()
        self.DataDir = self.MainConfig.get('MAIN', 'data directory')
        self.IrradiationFile = self.MainConfig.get('MISC', 'irradiation_file')
        self.TCString = self.load_test_campaign(test_campaign)
        self.TCDir = self.generate_tc_directory()

        # Configuration & Root Files
        self.Config = self.load_config()
        self.RootFileDir = self.load_rootfile_dirname()
        self.RootFilePath = self.load_rootfile_path()

        # Run Info
        self.RunInfoFile = join(self.TCDir, 'run_log.json')
        self.RunInfo = self.load_run_info()
        self.RootFile = None
        self.Tree = None
        self.TreeName = self.Config.get('BASIC', 'treename')

        # Settings
        self.PixelSize = loads(self.MainConfig.get('PIXEL', 'size'))
        self.NPixels = loads(self.MainConfig.get('PIXEL', 'amount'))
        self.TriggerPlanes = self.load_trigger_planes()

        # General Information
        self.Flux = self.calculate_flux()
        self.DUTNames = self.load_dut_names()
        self.DUTNumbers = self.load_dut_numbers()
        self.Bias = self.load_biases()
        self.Type = self.get_type()

        # Times
        self.LogStart = self.load_log_start()
        self.LogEnd = self.load_log_stop()
        self.Duration = self.LogEnd - self.LogStart
        self.LogHalfTime = time_stamp(self.LogStart + self.Duration / 2)

        self.Converter = Converter(self)
        if self.set_run(run_number, tree):
            # tree info
            self.Time = get_time_vec(self.Tree) if t_vec is None else t_vec
            self.StartEvent = 0
            self.NEntries = int(self.Tree.GetEntries())
            self.EndEvent = self.NEntries - 1
            self.StartTime = self.load_start_time()
            self.EndTime = self.load_end_time()
            self.TotalTime = self.load_total_time()
            self.TotalMinutes = self.TotalTime / 60000.
            self.Duration = timedelta(seconds=self.TotalTime)
            self.NPlanes = self.load_n_planes()

    def set_run(self, run_number, root_tree):
        if run_number is None:
            return False
        if run_number < 0 and type(run_number) is not int:
            critical('incorrect run_number')

        self.RunNumber = run_number
        self.load_run_info()
        self.Flux = self.calculate_flux()

        # check for conversion
        if root_tree and type(root_tree) is tuple:
            self.RootFile = root_tree[0]
            self.Tree = root_tree[1]
        elif root_tree:
            self.Converter.convert_run()
            self.load_rootfile()
        else:
            return False
        if not self.rootfile_is_valid():
            self.Converter.convert_run()
            self.load_rootfile()
        return True

    def get_type(self):
        return self.Config.get('BASIC', 'type') if self.RunNumber is not None else None

    # ----------------------------------------
    # region INIT
    def load_test_campaign(self, testcampaign):
        testcampaign = self.MainConfig.get('MAIN', 'default test campaign') if testcampaign is None else testcampaign
        if testcampaign not in self.get_test_campaigns() and self.RunNumber is not None:
            critical('The Testcampaign {} does not exist!'.format(testcampaign))
        return testcampaign

    def get_test_campaigns(self):
        return [basename(path).replace('_', '').strip('psi') for path in glob(join(self.DataDir, 'psi*'))]

    def load_rootfile(self):
        file_path = self.RootFilePath
        self.info('\n\033[1A\rLoading information for rootfile: {file}'.format(file=basename(file_path)))
        self.RootFile = TFile(file_path)
        self.Tree = self.RootFile.Get(self.TreeName)

    def load_config(self):
        run_number = self.RunNumber if hasattr(self, 'RunNumber') else None
        parser = ConfigParser({'excluded_runs': '[]'})  # add non default option
        base_file_name = join(get_base_dir(), 'Configuration', self.TCString, 'RunConfig.ini')
        if not file_exists(base_file_name):
            log_critical('RunConfig.ini does not exist for {0}! Please create it in Configuration/{0}!'.format(self.TCString))
        parser.read(base_file_name)  # first read the main config file with general information for all splits
        if parser.has_section('SPLIT') and run_number is not None:
            split_runs = [0] + loads(parser.get('SPLIT', 'runs')) + [inf]
            config_nr = next(i for i in xrange(1, len(split_runs)) if split_runs[i - 1] <= run_number < split_runs[i])
            parser.read(join(get_base_dir(), 'Configuration', self.TCString, 'RunConfig{nr}.ini'.format(nr=config_nr)))  # add the content of the split config
        return parser

    def load_rootfile_path(self):
        return join(self.RootFileDir, 'TrackedRun{run:03d}.root'.format(run=self.RunNumber)) if self.RunNumber is not None else None

    def load_rootfile_dirname(self):
        fdir = 'pads' if self.get_type() == 'pad' else self.get_type()
        return ensure_dir(join(self.TCDir, 'root', fdir)) if self.RunNumber is not None else None

    def generate_tc_directory(self):
        return join(self.DataDir, 'psi_{y}_{m}'.format(y=self.TCString[:4], m=self.TCString[4:]))

    def load_trigger_planes(self):
        default = self.get_unmasked_area().keys() if self.load_mask() is not None else [1, 2]
        return loads(self.Config.get('BASIC', 'trigger planes')) if self.Config.has_option('BASIC', 'trigger planes') else default

    def get_n_diamonds(self, run_number=None):
        run_info = self.load_run_info(run_number)
        return len([key for key in run_info if key.startswith('dia') and key[-1].isdigit()])

    def load_dut_numbers(self):
        return [i + 1 for i in xrange(len([key for key in self.RunInfo.iterkeys() if key.startswith('dia') and key[-1].isdigit()]))]

    def load_dut_type(self):
        dut_type = self.Config.get('BASIC', 'type') if self.RunNumber is not None else None
        if dut_type not in ['pixel', 'pad', None]:
            log_critical("The DUT type {0} has to be either 'pixel' or 'pad'".format(dut_type))
        return dut_type

    def load_default_info(self):
        with open(join(self.Dir, 'Runinfos', 'defaultInfo.json')) as f:
            return load(f)

    def load_run_info_file(self):
        if not file_exists(self.RunInfoFile):
            log_critical('Run Log File: "{f}" does not exist!'.format(f=self.RunInfoFile))
        with open(self.RunInfoFile) as f:
            return load(f)

    def load_run_info(self, run_number=None):
        data = self.load_run_info_file()

        run_number = self.RunNumber if run_number is None else run_number
        if run_number >= 0:
            run_info = data.get(str(run_number))
            if run_info is None:  # abort if the run is still not found
                log_critical('Run {} not found in json run log file!'.format(run_number))
            self.RunInfo = run_info
            self.RunInfo['masked pixels'] = [0] * 4
            self.translate_diamond_names()
            return run_info
        else:
            self.RunInfo = self.load_default_info()
            return self.RunInfo

    def load_dut_names(self):
        return [self.RunInfo['dia{nr}'.format(nr=i)] for i in xrange(1, self.get_n_diamonds() + 1)]

    def load_biases(self):
        return [int(self.RunInfo['dia{nr}hv'.format(nr=i)]) for i in xrange(1, self.get_n_diamonds() + 1)]

    def load_log_start(self):
        return conv_log_time(self.RunInfo['starttime0'])

    def load_log_stop(self):
        return conv_log_time(self.RunInfo['endtime'])

    def load_start_time(self):
        return int(round(self.get_time_at_event(self.StartEvent)))

    def load_end_time(self):
        return int(round(self.get_time_at_event(self.EndEvent)))

    def load_total_time(self):
        return (self.Time[-1] - self.Time[0]) / 1000

    def load_n_planes(self):
        if self.has_branch('cluster_col'):
            self.Tree.Draw('@cluster_col.size()', '', 'goff', 1)
            return int(self.Tree.GetV1()[0])
        else:
            return 4

    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region MASK
    def load_mask_file_path(self):
        mask_dir = self.MainConfig.get('MAIN', 'maskfile directory') if self.MainConfig.has_option('MAIN', 'maskfile directory') else join(self.DataDir, self.TCDir, 'masks')
        if not dir_exists(mask_dir):
            log_warning('Mask file directory does not exist ({})!'.format(mask_dir))
        return join(mask_dir, basename(self.RunInfo['maskfile']))

    def load_mask(self):
        mask_file = self.load_mask_file_path()
        if basename(mask_file) in ['no mask', 'none', 'none!']:
            return
        mask_data = {}
        # format: cornBot roc x1 y1
        #         cornTop roc x2 y2
        try:
            with open(mask_file) as f:
                for line in f:
                    if len(line) > 3 and not line.startswith('#'):
                        if not line.startswith('corn'):
                            log_warning('Invalid mask file: "{}". Not taking any mask!'.format(mask_file))
                            return
                        data = line.split()
                        plane = int(data[1])
                        if plane not in mask_data:
                            mask_data[plane] = zeros(4)
                        if line.startswith('cornBot'):
                            mask_data[plane][:2] = [int(data[i]) for i in [2, 3]]
                        if line.startswith('cornTop'):
                            mask_data[plane][2:] = [int(data[i]) for i in [2, 3]]
        except Exception as err:
            log_warning(err)
            log_warning('Could not read mask file... not taking any mask!')
        return mask_data

    def get_unmasked_area(self):
        if self.RunNumber is None:
            return
        mask = self.load_mask()
        # default pixel
        pix_x, pix_y = self.PixelSize
        n_pix_x, n_pix_y = self.NPixels
        pixel_area = pix_x * pix_y
        full_area = n_pix_x * n_pix_y * pixel_area
        if mask is None:
            return {plane: full_area for plane in self.TriggerPlanes}
        # format {plane: [x1, y1, x2, y2]}
        return {plane: pixel_area * (v[2] - v[0] + 1) * (v[3] - v[1] + 1) for plane, v in mask.iteritems()}

    def find_for_in_comment(self):
        for name in ['for1', 'for2']:
            if name not in self.RunInfo:
                for cmt in self.RunInfo['comments'].split('\r\n'):
                    cmt = cmt.replace(':', '')
                    cmt = cmt.split(' ')
                    if str(cmt[0].lower()) == name:
                        self.RunInfo[name] = int(cmt[1])
        return 'for1' in self.RunInfo
    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region HELPERS
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
        return parser.get('ALIASES', dia.lower()) if dia.lower() in parser.options('ALIASES') else log_critical('Please add {} to the diamond aliases!'.format(dia.encode()))

    def reload_run_config(self, run_number):
        self.RunNumber = run_number
        self.Config = self.load_config()
        self.RunInfo = self.load_run_info()

    def rootfile_is_valid(self, file_path=None):
        tfile = self.RootFile if file_path is None else TFile(file_path)
        ttree = self.Tree if file_path is None else tfile.Get(self.TreeName)
        is_valid = not tfile.IsZombie() and tfile.ClassName() == 'TFile' and ttree.ClassName() == 'TTree'
        if not is_valid:
            warning('Invalid TFile or TTree! Deleting file {}'.format(tfile.GetName()))
            remove_file(tfile.GetName())
        return is_valid

    def calculate_flux(self):
        if self.RunNumber is None:
            return
        fluxes = []
        if self.find_for_in_comment():
            for i, area in enumerate(self.get_unmasked_area().itervalues(), 1):
                fluxes.append(self.RunInfo['for{num}'.format(num=i)] / area / 1000)  # in kHz/cm^2
        else:
            fluxes.append(self.RunInfo['measuredflux'])
        flux = mean(fluxes)
        self.RunInfo['mean flux'] = flux
        return make_ufloat((flux, .1 * flux))

    def find_n_events(self, n, cut, start):
        total_events = self.Tree.Draw('event_number', cut, 'goff', self.NEntries, start)
        evt_numbers = [self.Tree.GetV1()[i] for i in xrange(total_events)]
        return int(evt_numbers[:n][-1] + 1 - start)

    # ----------------------------------------
    # endregion

    # ----------------------------------------
    # region GET
    def get_flux(self, rel_error=0.):
        return ufloat(self.Flux.n, self.Flux.s + self.Flux.n * rel_error) if self.Flux else self.RunInfo['aimed flux']

    def get_time(self):
        return make_ufloat((time_stamp(self.LogStart + self.Duration / 2), self.Duration.seconds / 2))

    def get_irradiations(self):
        with open(join(self.Dir, self.IrradiationFile), 'r') as f:
            data = load(f)
        for dia in self.DUTNames:
            if self.TCString not in data:
                log_critical('Please add "{}" to the {}'.format(self.TCString, self.IrradiationFile))
            if dia not in data[self.TCString].keys() + ['None']:
                log_critical('Please add "{}" to the irradiation file for {}'.format(dia, self.TCString))
        return [data[self.TCString][dia] for dia in self.DUTNames if dia not in ['None']]

    def get_attenuators(self):
        return [str(self.RunInfo['att_dia{}'.format(i)]) for i in xrange(1, self.get_n_diamonds() + 1) if 'att_dia1' in self.RunInfo]

    def get_diamond_name(self, channel):
        return self.DUTNames[channel]

    def get_channel_name(self, channel):
        self.Tree.GetEntry()
        return self.Tree.sensor_name[channel]

    def get_rate_string(self):
        rate = self.Flux
        unit = 'MHz/cm^{2}' if rate > 1000 else 'kHz/cm^{2}'
        rate = round(float(rate / 1000.), 1) if rate > 1000 else int(round(float(rate), 0))
        return '{rate:>3} {unit}'.format(rate=rate, unit=unit)

    def get_time_at_event(self, event):
        """ For negative event numbers it will return the time stamp at the startevent. """
        event = self.EndEvent if event == -1 else event
        event = sorted((0, event, self.EndEvent))[1]  # guarantees that the event is in [0, EndEvent]
        return self.Time[event] / 1000.

    def get_event_at_time(self, seconds):
        """ Returns the event nunmber at time dt from beginning of the run. Accuracy: +- 1 Event """
        # return time of last event if input is too large
        if seconds >= self.TotalTime or seconds == -1:
            return self.NEntries
        last_time = 0
        for i, t in enumerate(self.Time):
            if t / 1000. >= seconds + self.Time[0] / 1000 >= last_time:
                return i
            last_time = t / 1000.
        return self.NEntries  # if loop doesn't find anything it has to be the last bin

    def get_root_vec(self, n, ind=1):
        vec = getattr(self.Tree, 'GetV{}'.format(ind))()
        vec.SetSize(n)
        return array(vec)

    # endregion
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW
    def show_run_info(self):
        print 'Run information for run', self.RunNumber
        for key, value in sorted(self.RunInfo.iteritems()):
            print '{k}: {v}'.format(k=key.ljust(13), v=value)

    def has_branch(self, name):
        return bool(self.Tree.GetBranch(name))

    # endregion
    # ----------------------------------------

    def __load_rootfile(self):
        file_path = self.Converter.get_final_file_path()
        print '\033[1A\rLoading information for rootfile: {file}'.format(file=file_path.split('/')[-1])
        self.RootFile = TFile(file_path)
        self.Tree = self.RootFile.Get(self.Config.get('BASIC', 'treename'))


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('run', nargs='?', default=392, type=int)
    p.add_argument('-tc', '--testcampaign', nargs='?', default='201510')
    p.add_argument('-t', '--tree', action='store_true')
    args = p.parse_args()
    z = Run(args.run, tree=None if not args.tree else False, test_campaign=args.testcampaign)
