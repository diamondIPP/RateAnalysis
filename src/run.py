#!/usr/bin/env python
from ROOT import TFile, TTree
from numpy import inf, log
from src.analysis import Analysis
from src.converter import *
from src.dut import DUT, Plane


class Run(Analysis):
    """ Run class containing all the information for a single run. """

    def __init__(self, number=None, testcampaign=None, load_tree=True, verbose=None):
        """
        :param number: if None is provided it creates a dummy run
        :param testcampaign: if None is provided ...
        :param load_tree: load the ROOT TTree
        :param verbose: turn on more output
        """
        # Basics
        super(Run, self).__init__(testcampaign, verbose=verbose, pickle_dir='Run')
        self.Number = number

        # Directories / Test Campaign
        self.IrradiationFile = join(self.Dir, self.MainConfig.get('MISC', 'irradiation file'))

        # Configuration & Root Files
        self.Config = self.load_run_config()
        self.RootFileDir = self.load_rootfile_dirname()
        self.RootFilePath = self.load_rootfile_path()

        # Run Info
        self.InfoFile = join(self.TCDir, 'run_log.json')
        self.Info = self.load_run_info()
        self.RootFile = None
        self.Tree = TTree()
        self.TreeName = self.Config.get('BASIC', 'treename')
        self.DUTs = [DUT(i + 1, self.Info) for i in range(self.get_n_diamonds())] if self.Number is not None else None

        # Settings
        self.Plane = Plane()
        self.TriggerPlanes = self.load_trigger_planes()

        # General Information
        self.Flux = self.get_flux()
        self.Type = self.get_type()

        # Times
        self.LogStart = self.load_log_start()
        self.LogEnd = self.load_log_stop()
        self.Duration = self.LogEnd - self.LogStart

        self.Converter = Converter(self)
        if self.set_run(number, load_tree):
            # tree info
            self.TimeOffset = None
            self.Time = self.load_time_vec()
            self.StartEvent = 0
            self.NEvents = int(self.Tree.GetEntries())
            self.EndEvent = self.NEvents - 1
            self.StartTime = self.get_time_at_event(self.StartEvent)
            self.EndTime = self.get_time_at_event(self.EndEvent)
            self.TotalTime = self.load_total_time()
            self.TotalMinutes = self.TotalTime / 60000.
            self.Duration = timedelta(seconds=self.TotalTime)
            self.LogEnd = self.LogStart + self.Duration  # overwrite if we know exact duration
            self.NPlanes = self.load_n_planes()

        self.TInit = time() - self.InitTime

    def __str__(self):
        return f'{self.__class__.__name__} {self.Number} ({self.TCString})'

    def __repr__(self):
        return self.__str__()

    def __call__(self, number, load_tree=False):
        self.set_run(number, load_tree)
        return self

    def set_run(self, number, load_tree):
        if number is None:
            return False
        if number < 0 and type(number) is not int:
            critical('incorrect run number')

        self.Number = number
        self.load_run_info()
        self.Flux = self.get_flux()

        # check for conversion
        if load_tree:
            self.Converter.convert_run()
            self.load_rootfile()
        else:
            return False
        if not self.rootfile_is_valid():
            self.Converter.convert_run()
            self.load_rootfile()
        return True

    def get_type(self):
        return self.Config.get('BASIC', 'type') if self.Number is not None else None

    # ----------------------------------------
    # region INIT
    def load_rootfile(self, prnt=True):
        self.info('Loading information for rootfile: {file}'.format(file=basename(self.RootFilePath)), endl=False, prnt=prnt)
        self.RootFile = TFile(self.RootFilePath)
        self.Tree = self.RootFile.Get(self.TreeName)
        return self.Tree

    def load_run_config(self):
        base_file_name = join(get_base_dir(), 'config', self.TCString, 'RunConfig.ini')
        if not file_exists(base_file_name):
            critical('RunConfig.ini does not exist for {0}! Please create it in config/{0}!'.format(self.TCString))
        parser = Config(base_file_name)  # first read the main config file with general information for all splits
        if parser.has_section('SPLIT') and self.Number is not None:
            split_runs = [0] + loads(parser.get('SPLIT', 'runs')) + [inf]
            config_nr = next(i for i in range(1, len(split_runs)) if split_runs[i - 1] <= self.Number < split_runs[i])
            parser.read(join(get_base_dir(), 'config', self.TCString, 'RunConfig{nr}.ini'.format(nr=config_nr)))  # add the content of the split config
        return parser

    @staticmethod
    def make_root_filename(run):
        return f'TrackedRun{run:0>3}.root'

    def make_root_subdir(self):
        return join('root', 'pads' if self.get_type() == 'pad' else self.get_type())

    def load_rootfile_path(self, run=None):
        run = choose(run, self.Number)
        return None if run is None else join(self.RootFileDir, self.make_root_filename(run))

    def load_rootfile_dirname(self):
        return ensure_dir(join(self.TCDir, self.make_root_subdir())) if self.Number is not None else None

    def load_trigger_planes(self):
        default = list(self.get_unmasked_area().keys()) if self.load_mask() else [1, 2]
        return array(loads(self.Config.get('BASIC', 'trigger planes')) if self.Config.has_option('BASIC', 'trigger planes') else default)

    def get_n_diamonds(self, run_number=None):
        run_info = self.load_run_info(run_number)
        return len([key for key in run_info if key.startswith('dia') and key[-1].isdigit()])

    def load_dut_numbers(self):
        return [i + 1 for i in range(len([key for key in self.Info.keys() if key.startswith('dia') and key[-1].isdigit()]))]

    def load_dut_type(self):
        dut_type = self.Config.get('BASIC', 'type') if self.Number is not None else None
        if dut_type not in ['pixel', 'pad', None]:
            critical("The DUT type {0} has to be either 'pixel' or 'pad'".format(dut_type))
        return dut_type

    def load_default_info(self):
        with open(join(self.Dir, 'Runinfos', 'defaultInfo.json')) as f:
            return load(f)

    def load_run_info_file(self):
        if not file_exists(self.InfoFile):
            critical('Run Log File: "{f}" does not exist!'.format(f=self.InfoFile))
        with open(self.InfoFile) as f:
            return load(f)

    def load_run_info(self, run_number=None):
        data = self.load_run_info_file()

        run_number = self.Number if run_number is None else run_number
        if run_number is not None:
            run_info = data.get(str(run_number))
            if run_info is None:  # abort if the run is still not found
                critical('Run {} not found in json run log file!'.format(run_number))
            self.Info = run_info
            self.Info['masked pixels'] = [0] * 4
            self.translate_diamond_names()
            return run_info
        else:
            self.Info = self.load_default_info()
            return self.Info

    def load_dut_names(self):
        return [self.Info['dia{nr}'.format(nr=i)] for i in range(1, self.get_n_diamonds() + 1)]

    def load_biases(self):
        return [int(self.Info['dia{nr}hv'.format(nr=i)]) for i in range(1, self.get_n_diamonds() + 1)]

    def load_log_start(self):
        return conv_log_time(self.Info['starttime0'])

    def load_log_stop(self):
        return conv_log_time(self.Info['endtime'])

    def load_total_time(self):
        return (self.Time[-1] - self.Time[0]) / 1000

    def load_n_planes(self):
        if self.has_branch('cluster_col'):
            self.Tree.Draw('@cluster_col.size()', '', 'goff', 1)
            return int(self.Tree.GetV1()[0])
        else:
            return 4

    def load_time_vec(self):
        t = get_time_vec(self.Tree)
        t0 = datetime.fromtimestamp(t[0] / 1000)
        self.TimeOffset = None if t0.year > 2000 and t0.day == self.LogStart.day else t[0] - time_stamp(self.LogStart) * 1000
        return t if self.TimeOffset is None else t - self.TimeOffset

    def load_plane_efficiency(self, plane):
        return self.load_plane_efficiencies()[plane - 1]

    def load_plane_efficiencies(self):
        return [ufloat(e, .03) for e in self.Config.get_list('BASIC', 'plane efficiencies')]
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region MASK
    def load_mask_file_path(self):
        mask_dir = self.MainConfig.get('MAIN', 'maskfile directory') if self.MainConfig.has_option('MAIN', 'maskfile directory') else join(self.DataDir, self.TCDir, 'masks')
        if not dir_exists(mask_dir):
            warning('Mask file directory does not exist ({})!'.format(mask_dir))
        return join(mask_dir, basename(self.Info['maskfile']))

    def load_mask(self):
        mask_file = self.load_mask_file_path()
        if basename(mask_file) in ['no mask', 'none', 'none!', ''] or self.Number is None:
            return
        mask_data = {}
        # format: cornBot roc x1 y1
        #         cornTop roc x2 y2
        try:
            with open(mask_file) as f:
                for line in f:
                    if len(line) > 3 and not line.startswith('#'):
                        if not line.startswith('corn'):
                            warning('Invalid mask file: "{}". Not taking any mask!'.format(mask_file))
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
            warning(err)
            warning('Could not read mask file... not taking any mask!')
        return mask_data

    def get_mask_dim(self, mm=False):
        return {plane: array([Plane.PX * (v[2] - v[0] + 1), Plane.PY * (v[3] - v[1] + 1)]) / (1 if mm else 10) for plane, v in self.load_mask().items()}

    def get_unmasked_area(self):
        if self.Number is None:
            return
        mask = self.load_mask()
        if mask is None:
            return {plane: Plane.Area / 100 for plane in self.TriggerPlanes}
        # format {plane: [x1, y1, x2, y2]}
        return {plane: Plane.PixArea / 100 * (v[2] - v[0] + 1) * (v[3] - v[1] + 1) for plane, v in mask.items()}

    def find_for_in_comment(self):
        for name in ['for1', 'for2']:
            if name not in self.Info:
                for cmt in self.Info['comments'].split('\r\n'):
                    cmt = cmt.replace(':', '')
                    cmt = cmt.split(' ')
                    if str(cmt[0].lower()) == name:
                        self.Info[name] = int(cmt[1])
        return 'for1' in self.Info
    # endregion MASK
    # ----------------------------------------

    # ----------------------------------------
    # region HELPERS
    def translate_diamond_names(self):
        for key, value in [(key, value) for key, value in self.Info.items() if key.startswith('dia') and key[-1].isdigit()]:
            self.Info[key] = self.translate_dia(value)

    def translate_dia(self, dia):
        if dia is None or dia.lower() in ['unknown', 'none']:
            return
        parser = Config(join(self.Dir, 'config', 'DiamondAliases.ini'))
        if dia in parser.get_values('ALIASES'):
            return dia
        return parser.get('ALIASES', dia.lower()) if dia.lower() in parser.options('ALIASES') else critical(f'Please add {dia} to config/DiamondAliases.ini!')

    def reload_run_config(self, run_number):
        self.Number = run_number
        self.Config = self.load_run_config()
        self.Info = self.load_run_info()
        self.RootFileDir = self.load_rootfile_dirname()
        self.RootFilePath = self.load_rootfile_path()
        return self.Config

    def rootfile_is_valid(self, file_path=None):
        tfile = self.RootFile if file_path is None else TFile(file_path)
        ttree = self.Tree if file_path is None else tfile.Get(self.TreeName)
        is_valid = not tfile.IsZombie() and tfile.ClassName() == 'TFile' and ttree.ClassName() == 'TTree'
        if not is_valid:
            warning('Invalid TFile or TTree! Deleting file {}'.format(tfile.GetName()))
            remove_file(tfile.GetName())
        return is_valid

    def calculate_plane_flux(self, plane=1, corr=True):
        """estimate the flux [kHz/cm²] through a trigger plane based on Poisson statistics."""
        rate, eff, area = self.Info[f'for{plane}'], self.load_plane_efficiency(plane), self.get_unmasked_area()[plane]
        return -log(1 - rate / Plane.Frequency) * Plane.Frequency / area / 1000 / (eff if corr else 1)  # count zero hits of Poisson

    def find_n_events(self, n, cut, start=0):
        evt_numbers = self.get_tree_vec(var='Entry$', cut=cut, nentries=self.NEvents, firstentry=start)
        return int(evt_numbers[:n][-1] + 1 - start)
    # endregion HELPERS
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_flux(self, plane=1, corr=True):
        if self.Number is None:
            return
        if not self.find_for_in_comment():
            warning('no plane rates in the data...')
            return self.Info['measuredflux'] / (ufloat(mean(self.load_plane_efficiencies()), 0.03) if corr else 1)
        return self.get_mean_flux(corr) if plane is None else self.calculate_plane_flux(plane, corr)

    def get_mean_flux(self, corr=True):
        return mean([self.get_flux(pl, corr) for pl in [1, 2]])

    def get_time(self):
        return ufloat(time_stamp(self.LogStart + self.Duration / 2), self.Duration.seconds / 2)

    def get_channel_name(self, channel):
        self.Tree.GetEntry()
        return self.Tree.sensor_name[channel]

    def get_time_at_event(self, event):
        """ For negative event numbers it will return the time stamp at the startevent. """
        return self.Time[min(event, self.EndEvent)] / 1000.

    def get_event_at_time(self, seconds, rel=False):
        """ Returns the event nunmber at time dt from beginning of the run. Accuracy: +- 1 Event """
        # return time of last event if input is too large
        if seconds - (self.StartTime if rel else 0) >= self.TotalTime or seconds == -1:
            return self.NEvents - 1
        return where(self.Time <= (seconds + (0 if rel else self.StartTime)) * 1000)[0][-1]

    def get_tree_vec(self, var, cut='', dtype=None, nentries=None, firstentry=0):
        return get_tree_vec(self.Tree, var, cut, dtype, nentries, firstentry)

    def get_tree_tuple(self):
        return (self.Tree, self.RootFile) if self.Tree is not None else False

    def get_time_vec(self):
        return self.Time if hasattr(self, 'Time') else None

    def get_bias_strings(self):
        return [str(b) for b in self.load_biases()]

    def get_high_rate_run(self):
        def f():
            from src.run_selection import RunSelector
            return int(RunSelector(testcampaign=self.TCString).get_high_rate_run(self.Number))
        return do_pickle(self.make_simple_pickle_path('HR', run=self.Number), f)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW
    def show_run_info(self):
        print('Run information for run', self.Number)
        for key, value in sorted(self.Info.items()):
            print('{k}: {v}'.format(k=key.ljust(13), v=value))

    def has_branch(self, name):
        return bool(self.Tree.GetBranch(name))

    def info(self, msg, endl=True, blank_lines=0, prnt=True):
        return info(msg, endl, prnt=self.Verbose and prnt, blank_lines=blank_lines)

    def add_to_info(self, t, txt='Done'):
        return add_to_info(t, txt, prnt=self.Verbose)
    # endregion SHOW
    # ----------------------------------------


if __name__ == '__main__':

    args = init_argparser(run=88, tc=None, tree=True)
    z = Run(args.run, testcampaign=args.testcampaign, load_tree=args.tree)
