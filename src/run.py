#!/usr/bin/env python
from json import dump
from ROOT import TFile, TTree
from numpy import inf, log

from src.analysis import Analysis
from src.converter import *
from src.dut import DUT, Plane


class Run(Analysis):
    """ Run class containing all the information for a single run. """

    NTelPlanes = 4

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
        self.DUTs = [self.dut(i + 1, self.Info) for i in range(self.get_n_diamonds())] if self.Number is not None else None

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

    def set_estimate(self, n=None):
        self.Tree.SetEstimate(choose(n, self.NEvents))

    def is_volt_scan(self):
        return any(name in self.Info['runtype'] for name in ['voltage', 'hv'])

    # ----------------------------------------
    # region INIT
    @property
    def dut(self):
        return DUT

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
        return array(self.Config.get_list('BASIC', 'trigger planes', [1, 2]))

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
        t0 = datetime.fromtimestamp(t[0] / 1000) if t[0] < 1e12 else None
        self.TimeOffset = None if t0 is None or t0.year > 2000 and t0.day == self.LogStart.day else t[0] - time_stamp(self.LogStart) * 1000
        return t if self.TimeOffset is None else t - self.TimeOffset

    def load_plane_efficiency(self, plane):
        return self.load_plane_efficiencies()[plane - 1]

    def load_plane_efficiencies(self):
        return [ufloat(e, .03) for e in self.Config.get_list('BASIC', 'plane efficiencies', default=[.95, .95])]
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region MASK
    def load_mask_file_path(self):
        mask_dir = self.MainConfig.get('MAIN', 'maskfile directory') if self.MainConfig.has_option('MAIN', 'maskfile directory') else join(self.DataDir, self.TCDir, 'masks')
        if not dir_exists(mask_dir):
            warning('Mask file directory does not exist ({})!'.format(mask_dir))
        return join(mask_dir, basename(self.Info['maskfile']))

    def load_mask(self, plane=None):
        mask_file = self.load_mask_file_path()
        if basename(mask_file).lower() in ['no mask', 'none', 'none!', ''] or self.Number is None:
            return
        try:
            data = genfromtxt(mask_file, [('id', 'U10'), ('pl', 'i'), ('x', 'i'), ('y', 'i')])
            if 'cornBot' not in data['id']:
                warning('Invalid mask file: "{}". Not taking any mask!'.format(mask_file))
            mask = [[data[where((data['pl'] == pl) & (data['id'] == n))][0][i] for n in ['cornBot', 'cornTop'] for i in [2, 3]] for pl in sorted(set(data['pl']))]
            mask = [[max(1, m[0]), max(1, m[1]), min(self.Plane.NCols - 2, m[2]), min(self.Plane.NRows - 2, m[3])] for m in mask]  # outer pixels are ignored
            return mask if plane is None else mask[plane - 1] if plane - 1 < len(mask) else None
        except Exception as err:
            warning(err)
            warning('Could not read mask file... not taking any mask!')

    def get_mask_dim(self, plane=1, mm=True):
        return Plane.get_mask_dim(self.load_mask(plane), mm)

    def get_mask_dims(self, mm=True):
        return array([self.get_mask_dim(pl, mm) for pl in [1, 2]])

    def get_unmasked_area(self, plane):
        return None if self.Number is None else Plane.get_area(self.load_mask(plane))

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

    def register_new_dut(self):
        if input('Do you want to add a new diamond? [y,n] ').lower() in ['y', 'yes']:
            dut_type = int(input('Enter the DUT type (1 for pCVD, 2 for scCVD, 3 for silicon): ')) - 1
            dut_name = input('Enter the name of the DUT (no "_"): ')
            alias = input(f'Enter the alias (no "_", for default {dut_name.lower()} press enter): ')
            self.add_alias(alias, dut_name, dut_type)
            self.add_dut_info(dut_name)
            return True
        else:
            return False

    @staticmethod
    def add_alias(alias, dut_name, dut_type):
        alias_file = join(Dir, 'config', 'DiamondAliases.ini')
        with open(alias_file, 'r+') as f:
            lines = [line.strip(' \n') for line in f.readlines()]
            i0 = lines.index(['# pCVD', '# scCVD', '# Silicon'][dut_type])
            i = next(i for i, line in enumerate(lines[i0:], i0) if line.strip() == '')
            lines.insert(i, f'{(alias if alias else dut_name).lower()} = {dut_name}')
            f.seek(0)
            f.writelines([f'{line}\n' for line in lines])
            info(f'added entry: {(alias if alias else dut_name).lower()} = {dut_name} in {alias_file}')

    def add_dut_info(self, dut_name):
        dia_info_file = join(Dir, 'Runinfos', 'dia_info.json')
        data = load_json(dia_info_file)
        if dut_name in data:
            return warning('The entered DUT name already exists!')
        tc = get_input(f'Enter the beam test [YYYYMM]', self.TCString)
        data[dut_name] = {'irradiation': {tc: get_input(f'Enter the irradiation for {tc}', '0')},
                          'boardnumber': {tc: get_input(f'Enter the board number for {tc}')},
                          'thickness': get_input('Enter the thickness'),
                          'size': get_input('Enter the lateral size ([x, y])'),
                          'manufacturer': get_input('Enter the manufacturer')}
        with open(dia_info_file, 'w') as f:
            dump(data, f, indent=2)
        info(f'added {dut_name} to {dia_info_file}')

    def translate_dia(self, dia):
        name, suf = dia.split('_')[0].lower(), '_'.join(dia.split('_')[1:])
        if name not in Config(join(self.Dir, 'config', 'DiamondAliases.ini')).options('ALIASES'):
            warning(f'{dia} was not found in config/DiamondAliases.ini!')
            if not self.register_new_dut():
                critical(f'unknown diamond {dia}')
        parser = Config(join(self.Dir, 'config', 'DiamondAliases.ini'))
        return '_'.join([parser.get('ALIASES', name)] + ([suf] if suf else []))

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
        is_valid = not tfile.IsZombie() and tfile.ClassName() == 'TFile' and ttree and ttree.ClassName() == 'TTree'
        if not is_valid:
            warning('Invalid TFile or TTree! Deleting file {}'.format(tfile.GetName()))
            remove_file(tfile.GetName())
        return is_valid

    def calculate_plane_flux(self, plane=1, corr=True):
        """estimate the flux [kHz/cm²] through a trigger plane based on Poisson statistics."""
        rate, eff, area = self.Info[f'for{plane}'], self.load_plane_efficiency(plane), self.get_unmasked_area(plane)
        return -log(1 - rate / Plane.Frequency) * Plane.Frequency / area / 1000 / (eff if corr else ufloat(1, .05))  # count zero hits of Poisson

    def find_n_events(self, n, cut, start=0):
        evt_numbers = self.get_tree_vec(var='Entry$', cut=cut, nentries=self.NEvents, firstentry=start)
        return int(evt_numbers[:n][-1] + 1 - start)

    def get_max_run(self):
        return int(max(self.load_run_info_file(), key=int))
    # endregion HELPERS
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_flux(self, plane=None, corr=True):
        if self.Number is None:
            return
        if not self.find_for_in_comment():
            # warning('no plane rates in the data...')
            return self.Info['measuredflux'] / (mean(self.load_plane_efficiencies()) if corr else 1)
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

    def get_event_at_time(self, seconds, rel=True):
        """ Returns the event nunmber at time dt from beginning of the run. Accuracy: +- 1 Event """
        if seconds - (0 if rel else self.StartTime) >= self.TotalTime or seconds == -1:  # return time of last event if input is too large
            return self.NEvents - 1
        return where(self.Time <= 1000 * (seconds + (self.StartTime if rel else 0)))[0][-1]

    def get_tree_vec(self, var, cut='', dtype=None, nentries=None, firstentry=0):
        return get_tree_vec(self.Tree, var, cut, dtype, nentries, firstentry)

    def get_tree_tuple(self):
        return (self.Tree, self.RootFile) if self.Tree is not None else False

    def get_time_vec(self):
        return self.Time if hasattr(self, 'Time') else None

    def get_bias_strings(self):
        return [str(b) for b in self.load_biases()]

    @save_pickle('HR', suf_args=0)
    def get_high_rate_run(self, high=True):
        from src.run_selection import RunSelector
        return int(RunSelector(testcampaign=self.TCString).get_high_rate_run(self.Number, high))

    def get_low_rate_run(self):
        return self.get_high_rate_run(high=False)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW
    def show_info(self):
        print('Run information for', self)
        for key, value in sorted(self.Info.items()):
            print(f'{key:<13}: {value}')

    def has_branch(self, name):
        return bool(self.Tree.GetBranch(name))

    def info(self, msg, endl=True, blank_lines=0, prnt=True):
        return info(msg, endl, prnt=self.Verbose and prnt, blank_lines=blank_lines)

    def add_to_info(self, t, txt='Done', prnt=True):
        return add_to_info(t, txt, prnt=self.Verbose and prnt)
    # endregion SHOW
    # ----------------------------------------


if __name__ == '__main__':

    args = init_argparser(run=88, tc=None, tree=True)
    z = Run(args.run, testcampaign=args.testcampaign, load_tree=args.tree)
