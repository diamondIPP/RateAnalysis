from glob import glob
from os import getcwd, chdir, rename, system
from os.path import expanduser, basename
from re import sub
from subprocess import check_call
from numpy import genfromtxt

from src.pad_alignment import PadAlignment
from src.pix_alignment import *


class Converter(object):
    """ The Converter checks if the root files exist and starts the conversion if neccessary.
        The conversion sequence is: raw file -> EUDAQ Converter -> Event Aligner -> Tracking Telescope.
        Every step adds information to the root file. If any of the intermediate root files exist it will start converting from this point."""

    def __init__(self, run=None):

        # Configuration
        self.RunConfig = run.Config
        self.MainConfig = run.MainConfig

        # Basics
        self.Run = run
        self.TestCampaign = run.TCString
        self.RunInfo = run.Info
        self.Type = run.Type

        # EUDAQ Converter
        self.ConverterTree = self.load_converter_tree()
        self.NChannels = 9 if self.ConverterTree.startswith('caen') else 4

        # Directories
        self.DataDir = run.DataDir
        self.TCDir = run.TCDir
        self.RawFileDir = self.load_rawfile_dirname()
        self.SoftwareDir = self.load_dirname(base='', option='software')
        self.EudaqDir = self.load_dirname(self.SoftwareDir, 'eudaq')
        self.TrackingDir = self.load_dirname(self.SoftwareDir, 'tracking')

        if self.Run.Number is not None:
            # Tracking Software
            self.TelescopeID = self.RunConfig.getint('BASIC', 'telescopeID')

            # Files
            self.RawFilePath = self.get_raw_file_path()
            self.NewConfigFile = join(self.EudaqDir, 'conf', 'tmp', f'{self.Run.Number}.ini')
            self.ErrorFile = join(self.Run.RootFileDir, 'Errors{:03d}.txt'.format(self.Run.Number if self.Run.Number is not None else 0))

            # Event Alignment
            self.DecodingErrors = self.read_errors()

    def load_dirname(self, base, option):
        path = join(base, expanduser(self.MainConfig.get('Directories', option)))
        if not dir_exists(path) and not file_exists(path):
            critical('"{}" does not exist. Please set it correctly in config/main.ini.'.format(path))
        return path

    def load_rawfile_dirname(self):
        file_dir = join(self.TCDir, self.MainConfig.get_value('Directories', 'raw', default='raw'))
        if not dir_exists(file_dir):
            warning('Could not find the raw file directory: {d}'.format(d=file_dir))
        return file_dir

    def load_filename(self, base, option):
        return self.load_dirname(base, option)

    def read_errors(self):
        if not file_exists(self.ErrorFile):
            return array([])
        with open(self.ErrorFile) as f:
            return array([value for value in f.readlines()], 'i4')

    def load_converter_tree(self):
        self.ConverterTree = '{}tree'.format(self.RunConfig.get('BASIC', 'digitizer').lower() if self.Type == 'pad' else 'telescope')
        return self.ConverterTree

    def set_run(self, run_number):
        self.Run.set_run(run_number, load_tree=False)
        self.RunInfo = self.Run.Info
        self.RunConfig = self.Run.reload_run_config(run_number)
        self.Type = self.Run.get_type()
        self.RawFilePath = self.get_raw_file_path()

    def get_raw_file_path(self, run=None):
        return join(self.RawFileDir, f'run{choose(run, self.Run.Number):0>6}.raw')

    def get_eudaqfile_path(self):
        file_names = glob(join(self.Run.RootFileDir, '{}*{:03d}.root'.format(self.MainConfig.get('Directories', 'eudaq prefix'), self.Run.Number)))
        if file_names:
            return file_names[0]
        return join(self.Run.RootFileDir, '{prefix}{run:06d}.root'.format(prefix=self.MainConfig.get('Directories', 'eudaq prefix'), run=self.Run.Number))

    def get_trackingfile_path(self):
        return self.get_eudaqfile_path().replace('.root', '_withTracks.root')

    def get_alignment_file_path(self):
        return join(self.TrackingDir, 'data', 'alignments.txt')

    def file_is_valid(self, file_path):
        return False if not file_exists(file_path) else self.Run.rootfile_is_valid(file_path)

    def convert_run(self):
        if self.file_is_valid(self.Run.RootFilePath):  # check if final root file exists
            self.add_plane_errors()  # add plane errors if set in config file
            return
        elif self.file_is_valid(self.get_trackingfile_path()):  # check if the file after tracking exists
            self.rename_tracking_file()
            return
        elif self.file_is_valid(self.get_eudaqfile_path()):  # check if eudaq file exists
            self.Run.info('Found eudaq root file --> starting conversion')
        else:
            self.Run.info('did not find any matching root file --> starting conversion')
            self.convert_raw_to_root()
        if not self.has_alignment():
            self.align_telescope()
        self.add_plane_errors()
        self.align_run()
        self.add_tracking()
        remove(self.get_eudaqfile_path())
        print_banner(f'finished {self.__class__.__name__} of run {self.Run.Number} in {get_elapsed_time(self.Run.InitTime)}', color='green')

    def convert_raw_to_root(self, tree=None, max_events=None):
        if not file_exists(self.RawFilePath):
            critical('The raw file {} does not exist ...'.format(self.RawFilePath))
        self.remove_pickle_files()
        curr_dir = getcwd()
        chdir(self.Run.RootFileDir)  # go to root directory
        # prepare converter command
        cmd_list = [join(self.EudaqDir, 'bin', 'Converter.exe'), '-t', choose(tree, self.ConverterTree), '-c', join(self.EudaqDir, 'conf', self.NewConfigFile), self.RawFilePath]
        self.set_converter_configfile(tree, max_events)
        print_banner('START CONVERTING RAW FILE FOR RUN {0}'.format(self.Run.Number))
        info('{}\n'.format(' '.join(cmd_list)))
        check_call(cmd_list)
        self.remove_new_configfile()
        self.remove_decodingfile()
        chdir(curr_dir)

    def convert_shadow_run(self):
        self.convert_raw_to_root(tree='telescopetree')
        self.add_tracking()
        remove(self.get_eudaqfile_path())

    def align_run(self):
        aligner = PadAlignment(self) if self.Type == 'pad' else PixAlignment(self)
        aligner.run()

    def align_telescope(self):
        has_eudaq_file = file_exists(self.get_eudaqfile_path())
        if not has_eudaq_file:  # convert raw file with telescope tree
            self.RunConfig.set('ROOTFILE_GENERATION', 'max_event_number', '100000')
            self.convert_raw_to_root('telescopetree')
            self.RunConfig.set('ROOTFILE_GENERATION', 'max_event_number', '0')
        if not self.has_alignment():
            self.create_new_tel_id()
        self.tracking_tel(1)
        if not has_eudaq_file:
            remove_file(self.get_eudaqfile_path())

    def has_alignment(self):
        return self.TelescopeID in genfromtxt(join(self.TrackingDir, 'config', 'telescopes.txt'), usecols=[0])

    def create_new_tel_id(self, tel_id=None):
        root_file = self.Run.RootFilePath if file_exists(self.Run.RootFilePath) else self.get_eudaqfile_path()
        cmd_list = [join(self.TrackingDir, 'python', 'add_new_telescope.py'), root_file, str(choose(tel_id, self.TelescopeID))]
        check_call(cmd_list)

    def add_plane_errors(self):
        if self.MainConfig.getboolean('MISC', 'plane errors'):
            with open(self.get_alignment_file_path()) as f:
                if len(f.readlines()[3].split()) == 8:  # check if errors are already in the alignment file
                    self.Run.info('Plane errors already added')
                    return
                print_banner('START FINDING PLANE ERRORS FOR RUN {}'.format(self.Run.Number))
                self.tracking_tel(action='2')

    def remove_pickle_files(self):
        files = glob(join(self.Run.Dir, 'metadata', '*', '*{tc}*_{run}*'.format(run=self.Run.Number, tc=self.Run.TCString)))
        self.Run.info('Removing {} pickle files for run {}'.format(len(files), self.Run.Number))
        for f in files:
            remove_file(f)

    def rename_tracking_file(self):
        rename(self.get_trackingfile_path(), self.Run.RootFilePath)

    def tracking_tel(self, action=0):
        root_file = self.Run.RootFilePath if file_exists(self.Run.RootFilePath) else self.get_eudaqfile_path()
        cmd_list = [join(self.TrackingDir, 'TrackingTelescope'), root_file, str(action), str(self.TelescopeID), '1']
        print(' '.join(cmd_list))
        check_call(cmd_list)

    def add_tracking(self):
        print_banner('START TRACKING FOR RUN {}'.format(self.Run.Number))
        curr_dir = getcwd()
        chdir(self.Run.RootFileDir)
        self.tracking_tel()
        self.rename_tracking_file()
        chdir(curr_dir)

    def load_polarities(self, pulser=False):
        option = '{}polarities'.format('pulser_' if pulser else '')
        if self.RunConfig.has_option('ROOTFILE_GENERATION', option):
            return self.RunConfig.get('ROOTFILE_GENERATION', option)
        fac = 1
        if self.RunConfig.has_option('ROOTFILE_GENERATION', 'inverted_polarities'):
            fac = -1 if int(self.RunConfig.get('ROOTFILE_GENERATION', 'inverted_polarities')) else 1
        active_regions = self.RunConfig.getint('ROOTFILE_GENERATION', 'active_regions')
        biases = self.Run.load_biases()
        polarities = [sign(biases.pop(0)) * fac if has_bit(active_regions, i) else 0 for i in range(self.NChannels)]
        return str([(1 if not pol and has_bit(active_regions, i) else pol) for i, pol in enumerate(polarities)])  # pol cannot be 0, just take 1 for 0V

    def copy_raw_file(self, redo=False, raw_dir='raw'):
        if not file_exists(self.RawFilePath) or redo:
            main_data_path = join('isg:', 'home', 'ipp', basename(self.TCDir), raw_dir, basename(self.RawFilePath))
            self.Run.info('Trying to copy {}'.format(basename(self.RawFilePath)))
            system('rsync -aPv {} {}'.format(main_data_path, self.RawFileDir))

    def remove_raw_file(self):
        remove_file(self.RawFilePath)

    def remove_final_file(self):
        remove_file(self.get_eudaqfile_path())
        remove_file(self.Run.RootFilePath)

    def remove_new_configfile(self):
        remove_file(self.NewConfigFile)

    def remove_decodingfile(self):
        remove_file(f'decoding_{self.Run.Number}.root')

    def set_converter_configfile(self, tree=None, max_events=None):
        filename = self.NewConfigFile if file_exists(self.NewConfigFile) else self.load_filename(join(self.EudaqDir, 'conf'), 'converter config')
        if not file_exists(filename):
            critical(f'EUDAQ config file: "{filename}" does not exist!')
        parser = Config(filename)
        section = 'Converter.{}'.format(choose(tree, self.ConverterTree))
        if self.Type == 'pad':
            parser.set(section, 'polarities', self.load_polarities())
            parser.set(section, 'pulser_polarities', self.load_polarities(pulser=True))
            parser.set(section, 'save_waveforms', self.RunConfig.get('ROOTFILE_GENERATION', 'active_regions') if self.MainConfig.get_value('SAVE', 'waveforms', bool) else '0')

        # remove unset ranges and regions
        new_options = self.RunConfig.options('ROOTFILE_GENERATION')
        for opt in parser.options(section):
            if (opt.endswith('_range') or opt.endswith('_region')) and opt not in new_options:
                parser.remove_option(section, opt)
        # set the new settings
        for key, value in self.RunConfig.items('ROOTFILE_GENERATION'):
            parser.set('Converter.telescopetree' if 'decoding' in key else section, key, value)
        if max_events is not None:
            parser.set(section, 'max_event_number', str(max_events))

        # write changes
        with open(self.NewConfigFile, 'w') as f:
            parser.write(f)
        self.format_converter_configfile()

    def format_converter_configfile(self):
        """ remove whitespaces, correct for capitalisation and sort options"""
        with open(self.NewConfigFile, 'r+') as f:
            content = f.readlines()
            for i, line in enumerate(content):
                line = line.replace('peaki', 'PeakI')
                line = sub('[)(\' ]', '', line)
                if len(line) > 3 and line[-2] == ',':
                    line = line[:-2] + '\n'
                content[i] = line

            section_indices = [i for i, line in enumerate(content) if line.startswith('[')] + [1000]
            sorted_content = []
            for i in range(len(section_indices) - 1):
                sorted_content.append(content[section_indices[i]])
                sorted_content += sorted(content[section_indices[i]+1:section_indices[i+1]])[1:]
                sorted_content.append('\n')

            f.seek(0)
            f.writelines(sorted_content)
            f.truncate()
