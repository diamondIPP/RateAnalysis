# ==============================================
# IMPORTS
# ==============================================
from ConfigParser import ConfigParser
from numpy import sign
from re import sub
from shutil import move
from os import getcwd, chdir, rename, system
from os.path import realpath
from glob import glob
from Utils import *
from PixAlignment import PixAlignment
from PadAlignment import PadAlignment
from ROOT import TProfile, TTree
from subprocess import check_call, CalledProcessError
from json import loads, load


__author__ = 'micha'


# ==============================================
# CLASS DEFINITION
# ==============================================
class Converter:
    def __init__(self, run):

        # main
        self.Run = run
        self.TestCampaign = run.TESTCAMPAIGN
        self.RunConfig = run.RunConfig
        self.SoftConfig = self.Run.MainConfigParser
        self.RunNumber = run.RunNumber
        self.RunInfo = run.RunInfo
        self.Type = self.load_type()

        # digitizer
        self.ConverterTree = '{0}tree'.format(self.RunConfig.get('BASIC', 'digitizer').lower() if self.Type == 'pad' and self.RunConfig.has_option('BASIC', 'digitizer') else 'telescope')
        self.NChannels = 9 if self.ConverterTree.startswith('caen') else 4

        # directories
        self.DataDir = run.DataDir
        self.TcDir = run.TCDir
        self.RawFileDir = self.load_raw_file()
        self.RootFileDir = self.load_root_file_dir()
        self.SoftwareDir = self.load_soft_dir()
        self.EudaqDir = '{soft}/{d}'.format(soft=self.SoftwareDir, d=self.SoftConfig.get('Converter', 'eudaqfolder'))
        self.AlignDir = '{soft}/{d}'.format(soft=self.SoftwareDir, d=self.SoftConfig.get('Converter', 'alignfolder'))

        # tracking
        self.TelescopeID = self.RunConfig.getint('BASIC', 'telescopeID') if run.RunNumber else None
        self.TrackingDir = join(self.SoftwareDir, self.SoftConfig.get('Converter', 'trackingfolder'))

        # files paths
        self.ConverterConfigFile = self.SoftConfig.get('Converter', 'converterFile')
        self.NewConfigFile = join(self.EudaqDir, 'conf', '{}.ini'.format(self.RunNumber))
        self.run_info_path = run.load_run_info_path()

        # configuration
        self.config = self.get_config()

        # alignment
        self.ErrorFile = join(self.RootFileDir, 'Errors{:03d}.txt'.format(self.RunNumber if self.RunNumber is not None else 0))
        self.DecodingErrors = self.read_errors()

    def load_type(self):
        return self.RunConfig.get('BASIC', 'type') if self.RunConfig.has_option('BASIC', 'type') else None

    def load_soft_dir(self):
        file_dir = self.SoftConfig.get('Converter', 'softwaredir')
        if not dir_exists(file_dir):
            log_warning('Could not find the software directory: {d}\nPlease set it correctly in Configuration/soft.ini'.format(d=file_dir))
        return file_dir

    def load_root_file_dir(self):
        if self.RunConfig.has_option('BASIC', 'runpath'):
            path = self.RunConfig.get('BASIC', 'runpath')
        else:
            path = join(self.DataDir, self.TcDir, 'root', '{dut}'.format(dut='pads' if self.Type == 'pad' else 'pixel'))
        ensure_dir(path)
        return path

    def load_raw_file(self):
        file_dir = join(self.DataDir, self.TcDir, self.RunConfig.get('BASIC', 'raw directory') if self.RunConfig.has_option('BASIC', 'raw directory') else 'raw')
        if not dir_exists(file_dir):
            log_warning('Could not find the raw file directory: {d}'.format(d=file_dir))
        return file_dir

    def read_errors(self):
        if not file_exists(self.ErrorFile):
            return []
        with open(self.ErrorFile) as f:
            return [int(value) for value in f.readlines()]

    def set_run(self, run_number):
        self.Run.set_run(run_number, root_tree=False)
        self.RunNumber = run_number
        self.RunInfo = self.Run.RunInfo
        self.RunConfig = self.Run.load_run_config()
        self.Type = self.load_type()
        self.RootFileDir = self.load_root_file_dir()

    @staticmethod
    def load_soft_config():
        conf = ConfigParser()
        main_dir = '/'.join(dirname(realpath(__file__)).split('/')[:-1])
        conf.read('{d}/Configuration/soft.ini'.format(d=main_dir))
        return conf

    def get_config(self):
        config = {}
        if self.RunNumber is None:
            return config
        options = self.RunConfig.options('ROOTFILE_GENERATION')
        for opt in options:
            if any(opt.endswith(ending) for ending in ['_range', '_region', '_range_drs4']) or 'polarities' in opt:
                config[opt] = loads(self.RunConfig.get('ROOTFILE_GENERATION', opt))
            elif opt not in ['excluded_runs']:
                config[opt] = self.RunConfig.getint('ROOTFILE_GENERATION', opt)
        return OrderedDict(sorted(config.iteritems()))

    def get_run_info(self, run_number):
        try:
            f = open(self.run_info_path, "r")
            data = load(f)
            f.close()
        except IOError:
            print 'could not read', self.run_info_path
            return
        run_infos = data.get(str(run_number))
        try:
            run_infos['hv dia1'] = run_infos.pop('dia1hv')
            run_infos['hv dia2'] = run_infos.pop('dia2hv')
        except KeyError:
            pass
        return run_infos
    
    def make_raw_file_path(self, run=None):
        return join(self.RawFileDir, 'run{run:06d}.raw'.format(run=self.RunNumber if run is None else run))

    def get_raw_file_path(self, crit=True):
        file_path = self.make_raw_file_path()
        if not file_exists(file_path) and crit:
            log_critical('The raw file {} does not exist ...'.format(file_path))
        return file_path

    def get_root_file_path(self):
        name = glob(join(self.RootFileDir, '*{run:05d}.root'.format(run=self.RunNumber)))
        return name[0] if name else ''

    def get_tracking_file_path(self):
        return self.get_root_file_path().replace('.root', '_withTracks.root')

    def get_final_file_path(self):
        return join(self.RootFileDir, 'TrackedRun{run:03d}.root'.format(run=self.RunNumber))

    def get_alignment_file_path(self):
        return join(self.TrackingDir, 'ALIGNMENT', 'telescope{}.dat'.format(self.TelescopeID))

    def convert_run(self):
        # check whether the root file w/ our w/o tracks already exist
        if file_exists(self.get_final_file_path()):
            self.find_plane_errors()
            return
        elif file_exists(self.get_tracking_file_path()):
            self.rename_tracking_file()
            return
        elif file_exists(self.get_root_file_path()):
            self.Run.log_info('did not find tracking file --> need conversion')
        else:
            self.Run.log_info('did not find any matching root file --> need conversion')
            self.convert_raw_to_root()
        if not self.validate_root_file():
            self.convert_raw_to_root()
        self.find_plane_errors()
        self.align_run()
        self.add_tracking()
        remove(self.get_root_file_path())

    def convert_raw_to_root(self):
        raw_file_path = self.get_raw_file_path()
        self.remove_pickle_files()
        curr_dir = getcwd()
        chdir(self.RootFileDir)  # go to root directory
        # prepare converter command
        cmd_list = [join(self.EudaqDir, 'bin', 'Converter.exe'), '-t', self.ConverterTree, '-c', join(self.EudaqDir, 'conf', self.NewConfigFile), raw_file_path]
        self.set_converter_configfile()
        print_banner('START CONVERTING RAW FILE FOR RUN {0}'.format(self.RunNumber))
        print ' '.join(cmd_list)
        max_tries = 30
        tries = 0
        while tries < max_tries:  # the command crashes randomly...
            try:
                check_call(cmd_list)
                break
            except CalledProcessError:
                tries += 1
        self.remove_new_configfile()
        chdir(curr_dir)

    def align_run(self):

        if self.Type == 'pad':
            pad_align = PadAlignment(self)
            pad_align.run()
        elif self.Type == 'pixel':
            print_banner('STARTING PIXEL EVENT ALIGNMENT OF RUN {r}'.format(r=self.RunNumber))
            pix_align = PixAlignment(self)
            if not pix_align.IsAligned and not pix_align.check_alignment():
                pix_align.write_aligned_tree()

    def align_telescope(self):
        # TODO implement!
        pass

    def find_plane_errors(self):
        if self.Run.MainConfigParser.getboolean('MISC', 'plane errors'):
            with open(self.get_alignment_file_path()) as f:
                if len(f.readlines()[3].split()) == 8:  # check if errors are already in the alignment file
                    self.Run.log_info('Plane errors already added')
                    return
                print_banner('START FINDING PLANE ERRORS FOR RUN {}'.format(self.RunNumber))
                self.tracking_tel(action='2')

    def remove_pickle_files(self):
        self.Run.log_info('Removing all pickle files for run {}'.format(self.RunNumber))
        files = glob(join(self.Run.Dir, 'Configuration', 'Individual_Configs', '*', '*{tc}*_{run}_*'.format(run=self.RunNumber, tc=self.Run.generate_tc_str())))
        for f in files:
            remove(f)

    def __rename_rootfile(self):
        rename(self.get_root_file_path(), self.get_final_file_path())

    def rename_tracking_file(self):
        rename(self.get_tracking_file_path(), self.get_final_file_path())

    def tracking_tel(self, action='0'):
        root_file = self.get_final_file_path() if file_exists(self.get_final_file_path()) else self.get_root_file_path()
        cmd_list = [join(self.TrackingDir, 'TrackingTelescope'), root_file, str(action), str(self.TelescopeID), '' if self.Type == 'pad' else '1']
        print ' '.join(cmd_list)
        curr_dir = getcwd()
        chdir(self.TrackingDir)
        check_call(cmd_list)
        chdir(curr_dir)

    def add_tracking(self):
        print_banner('START TRACKING FOR RUN {}'.format(self.RunNumber))
        self.tracking_tel()
        # move file from tracking directory to data directory
        move(join(self.TrackingDir, basename(self.get_tracking_file_path())), self.RootFileDir)
        self.rename_tracking_file()

    def load_polarities(self, pulser=False):
        option = '{}polarities'.format('pulser_' if pulser else '')
        if self.RunConfig.has_option('ROOTFILE_GENERATION', option):
            return self.RunConfig.get('ROOTFILE_GENERATION', option)
        fac = 1
        if self.RunConfig.has_option('ROOTFILE_GENERATION', 'inverted_polarities'):
            fac = -1 if int(self.RunConfig.get('ROOTFILE_GENERATION', 'inverted_polarities')) else 1
        active_regions = self.RunConfig.getint('ROOTFILE_GENERATION', 'active_regions')
        biases = deepcopy(self.Run.Bias)
        polarities = [sign(biases.pop(0)) * fac if has_bit(active_regions, i) else 0 for i in xrange(self.NChannels)]
        return str([(1 if not pol and has_bit(active_regions, i) else pol) for i, pol in enumerate(polarities)])  # pol cannot be 0, just take 1 for 0V

    def copy_raw_file(self, redo=False):
        if not file_exists(self.make_raw_file_path()) or redo:
            main_data_path = join('isg:', 'home', 'ipp', self.TcDir, 'raw', basename(self.make_raw_file_path()))
            log_message('Trying to copy {}'.format(basename(self.make_raw_file_path())))
            system('rsync -aPv {} {}'.format(main_data_path, self.RawFileDir))

    def remove_raw_file(self):
        remove_file(self.make_raw_file_path())

    def remove_final_file(self):
        remove_file(self.get_root_file_path())
        remove_file(self.get_final_file_path())

    def remove_new_configfile(self):
        remove_file(self.NewConfigFile)

    def validate_root_file(self):
        f = TFile(self.get_root_file_path())
        t = f.Get(self.Run.treename)
        if not t or t.IsA() != TTree().IsA():
            log_warning('corrupted root file: '.format(self.get_root_file_path()))
            self.remove_final_file()
            return False
        return True

    def set_converter_configfile(self):
        parser = ConfigParser()
        config_file = join(self.EudaqDir, 'conf', self.ConverterConfigFile)
        if not file_exists(config_file):
            log_critical('EUDAQ config file: "{}" does not exist!'.format(config_file))
        parser.read(config_file)
        section = 'Converter.{}'.format(self.ConverterTree)
        if self.Type == 'pad':
            parser.set(section, 'polarities', self.load_polarities())
            parser.set(section, 'pulser_polarities', self.load_polarities(pulser=True))

        # remove unset ranges and regions
        new_options = self.RunConfig.options('ROOTFILE_GENERATION')
        for opt in parser.options(section):
            if (opt.endswith('_range') or opt.endswith('_region')) and opt not in new_options:
                parser.remove_option(section, opt)
        # set the new settings
        for key, value in self.config.iteritems():
            parser.set('Converter.telescopetree' if key.startswith('decoding') else section, key, value)

        # write changes
        f = open(self.NewConfigFile, 'w')
        parser.write(f)
        f.close()

        self.format_converter_configfile()

    def format_converter_configfile(self):
        """ remove whitespaces, correct for capitalisation and sort options"""

        f = open(self.NewConfigFile, 'r+')
        content = f.readlines()
        for i, line in enumerate(content):
            line = line.replace('peaki', 'PeakI')
            line = sub('[)(\' ]', '', line)
            if len(line) > 3 and line[-2] == ',':
                line = line[:-2] + '\n'
            content[i] = line

        section_indices = [i for i, line in enumerate(content) if line.startswith('[')] + [1000]
        sorted_content = []
        for i in xrange(len(section_indices) - 1):
            sorted_content.append(content[section_indices[i]])
            sorted_content += sorted(content[section_indices[i]+1:section_indices[i+1]])[1:]
            sorted_content.append('\n')

        f.seek(0)
        f.writelines(sorted_content)
        f.truncate()
        f.close()

    def check_alignment(self, tree, binning=5000):

        n_entries = tree.GetEntries()
        nbins = n_entries / binning
        h = TProfile('h', 'Pulser Rate', nbins, 0, n_entries)
        tree.Draw('(@col.size()>1)*100:Entry$>>h', 'pulser', 'goff')
        is_aligned = self.__check_alignment_histo(h)
        if not is_aligned:
            log_warning('The events of RUN {run} are not aligned!'.format(run=self.RunNumber))
        return is_aligned

    @staticmethod
    def __check_alignment_histo(histo):
        h = histo
        for bin_ in xrange(h.FindBin(20000), h.GetNbinsX()):
            if h.GetBinContent(bin_) > 40:
                return False
        return True


if __name__ == '__main__':

    from argparse import ArgumentParser
    from Run import Run

    p = ArgumentParser()

    p.add_argument('run', nargs='?', default=398, type=int)
    p.add_argument('-tc', '--testcampaign', nargs='?', default='201510')
    args = p.parse_args()

    zrun = Run(args.run, test_campaign=args.testcampaign, tree=False, verbose=True)
    z = Converter(zrun)
