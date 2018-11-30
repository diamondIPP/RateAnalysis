# ==============================================
# IMPORTS
# ==============================================
from ConfigParser import ConfigParser
from math import copysign
from re import sub
from shutil import move
from os import getcwd, chdir, rename
from os.path import dirname, realpath
from glob import glob
from Utils import *
from PixAlignment import PixAlignment
from PadAlignment import PadAlignment
from ROOT import TProfile
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
        self.RunParser = run.run_config_parser
        self.SoftConfig = self.load_soft_config()
        self.RunNumber = run.RunNumber
        self.RunInfo = run.RunInfo
        self.Type = self.RunParser.get('BASIC', 'type') if self.RunParser.has_option('BASIC', 'type') else None

        # digitizer
        self.ConverterTree = '{0}tree'.format(self.RunParser.get('BASIC', 'digitizer').lower() if self.Type == 'pad' and self.RunParser.has_option('BASIC', 'digitizer') else 'telescope')
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
        self.TelescopeID = self.RunParser.getint('BASIC', 'telescopeID') if run.RunNumber else None
        self.TrackingDir = join(self.SoftwareDir, self.SoftConfig.get('Converter', 'trackingfolder'))

        # files paths
        self.ConverterConfigFile = self.SoftConfig.get('Converter', 'converterFile')
        self.run_info_path = run.load_run_info_path()

        # configuration
        self.config = self.get_config()

        # alignment
        self.ErrorFile = join(self.RootFileDir, 'Errors{:03d}.txt'.format(self.RunNumber if self.RunNumber is not None else 0))
        self.DecodingErrors = self.read_errors()

    def load_soft_dir(self):
        file_dir = self.SoftConfig.get('Converter', 'softwaredir')
        if not dir_exists(file_dir):
            log_warning('Could not find the software directory: {d}\nPlease set it correctly in Configuration/soft.ini'.format(d=file_dir))
        return file_dir

    def load_root_file_dir(self):
        if self.RunParser.has_option('BASIC', 'runpath'):
            path = self.RunParser.get('BASIC', 'runpath')
        else:
            path = join(self.DataDir, self.TcDir, 'root', '{dut}'.format(dut='pads' if self.Type == 'pad' else 'pixel'))
        ensure_dir(path)
        return path

    def load_raw_file(self):
        if self.RunParser.has_option('ConverterFolders', 'rawfolder'):
            file_dir = self.RunParser.get('ConverterFolders', 'rawfolder')
        else:
            file_dir = join(self.DataDir, self.TcDir, 'raw')
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
        options = self.RunParser.options('ROOTFILE_GENERATION')
        for opt in options:
            if any(opt.endswith(ending) for ending in ['_range', '_region', '_range_drs4']) or opt == 'polarities':
                config[opt] = loads(self.RunParser.get('ROOTFILE_GENERATION', opt))
            elif opt not in ['excluded_runs']:
                config[opt] = self.RunParser.getint('ROOTFILE_GENERATION', opt)
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

    def get_raw_file_path(self):
        file_path = join(self.RawFileDir, 'run{run:06d}.raw'.format(run=self.RunNumber))
        if not file_exists(file_path):
            log_critical('The raw file {} does not exist ...')
        return file_path

    def get_root_file_path(self):
        return join(self.RootFileDir, 'test{run:06d}.root'.format(run=self.RunNumber))

    def get_tracking_file_path(self):
        return self.get_root_file_path().replace('.root', '_withTracks.root')

    def get_final_file_path(self):
        return join(self.RootFileDir, 'TrackedRun{run:03d}.root'.format(run=self.RunNumber))

    def convert_run(self):
        # check whether the root file w/ our w/o tracks already exist
        if file_exists(self.get_final_file_path()):
            return
        elif file_exists(self.get_tracking_file_path()):
            self.rename_tracking_file()
            return
        elif file_exists(self.get_root_file_path()):
            self.Run.log_info('did not find tracking file --> need conversion')
        else:
            self.Run.log_info('did not find any matching root file --> need conversion')
            self.convert_raw_to_root()
        self.align_run()
        self.add_tracking()
        remove(self.get_root_file_path())

    def convert_raw_to_root(self):
        raw_file_path = self.get_raw_file_path()
        self.remove_pickle_files()
        curr_dir = getcwd()
        chdir(self.RootFileDir)  # go to root directory
        # prepare converter command
        cmd_list = [join(self.EudaqDir, 'bin', 'Converter.exe'), '-t', self.ConverterTree, '-c', join(self.EudaqDir, 'conf', self.ConverterConfigFile), raw_file_path]
        self.set_converter_configfile()
        print_banner('START CONVERTING RAW FILE FOR RUN {0}'.format(self.RunNumber))
        print ' '.join(cmd_list)
        while True:  # the command crashes randomly...
            try:
                check_call(cmd_list)
                break
            except CalledProcessError:
                continue
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

    def remove_pickle_files(self):
        self.Run.log_info('Removing all pickle files for run {}'.format(self.RunNumber))
        files = glob(join(self.Run.Dir, 'Configuration', 'Individual_Configs', '*', '*{tc}*_{run}_*'.format(run=self.RunNumber, tc=self.Run.generate_tc_str())))
        for f in files:
            remove(f)

    def __rename_rootfile(self):
        rename(self.get_root_file_path(), self.get_final_file_path())

    def rename_tracking_file(self):
        rename(self.get_tracking_file_path(), self.get_final_file_path())

    def add_tracking(self):
        print_banner('START TRACKING FOR RUN {}'.format(self.RunNumber))
        curr_dir = getcwd()
        chdir(self.TrackingDir)
        cmd_list = [join(self.TrackingDir, 'TrackingTelescope'), self.get_root_file_path(), '0', str(self.TelescopeID), '' if self.Type == 'pad' else '1']
        print ' '.join(cmd_list)
        check_call(cmd_list)
        chdir(curr_dir)
        # move file from tracking directory to data directory
        move(join(self.TrackingDir, basename(self.get_tracking_file_path())), self.RootFileDir)
        self.rename_tracking_file()

    def load_polarities(self, info):
        if self.RunParser.has_option('ROOTFILE_GENERATION', 'polarities'):
            return self.RunParser.get('ROOTFILE_GENERATION', 'polarities')
        active_regions = self.RunParser.getint('ROOTFILE_GENERATION', 'active_regions')
        pols = []
        i = 1
        for j in xrange(self.NChannels):
            if has_bit(active_regions, j):
                pols.append(int(copysign(1, info['dia{0}hv'.format(i)])))
                i += 1
            else:
                pols.append(0)
        return str(pols)

    def set_converter_configfile(self):

        parser = ConfigParser()
        config_file = join(self.EudaqDir, 'conf', self.ConverterConfigFile)
        parser.read(config_file)
        section = 'Converter.{}'.format(self.ConverterTree)
        if self.Type == 'pad':
            parser.set(section, 'polarities', self.load_polarities(self.RunInfo))
            parser.set(section, 'pulser_polarities', self.load_polarities(self.RunInfo))

        # remove unset ranges and regions
        new_options = self.RunParser.options('ROOTFILE_GENERATION')
        for opt in parser.options(section):
            if (opt.endswith('_range') or opt.endswith('_region')) and opt not in new_options:
                parser.remove_option(section, opt)
        # set the new settings
        for key, value in self.config.iteritems():
            parser.set('Converter.telescopetree' if key.startswith('decoding') else section, key, value)

        # write changes
        f = open(config_file, 'w')
        parser.write(f)
        f.close()

        self.format_converter_configfile()

    def format_converter_configfile(self):
        """ remove whitespaces, correct for capitalisation and sort options"""

        f = open(join(self.EudaqDir, 'conf', self.ConverterConfigFile), 'r+')
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
    from Run import Run
    zrun = Run(343, test_campaign=None, tree=False, verbose=True)
    z = Converter(zrun)
