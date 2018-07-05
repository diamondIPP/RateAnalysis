# ==============================================
# IMPORTS
# ==============================================
import json

from ConfigParser import ConfigParser
from math import copysign
from re import sub
from shutil import move
from os import remove, getcwd, chdir, system, rename
from os.path import dirname, realpath, getsize
from glob import glob
from Utils import *
from PixAlignment import PixAlignment
from PadAlignment import PadAlignment
from ROOT import TProfile
do_gui = False
if do_gui:
    from Tkinter import *


__author__ = 'micha'


def print_banner(message):
    print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(message)) * '=', msg=message)


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
        self.Type = self.RunParser.get('BASIC', 'type')

        # digitizer
        self.ConverterTree = '{0}tree'.format(self.RunParser.get('BASIC', 'digitizer').lower() if self.Type == 'pad' else 'telescope')
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
        self.TrackingDir = '{soft}/{d}'.format(soft=self.SoftwareDir, d=self.SoftConfig.get('Converter', 'trackingfolder'))

        # files paths
        self.ConverterConfigFile = self.SoftConfig.get('Converter', 'converterFile')
        self.run_info_path = run.load_run_info_path()
        # prefixes
        self.raw_prefix = self.load_prefix()
        self.root_prefix = self.raw_prefix.replace('run', 'test')

        # configuration
        self.config = self.get_config()
        self.ErrorEvents = self.load_errors()

        # gui
        if do_gui:
            self.root = Tk()
            self.root.withdraw()
            self.frame = Frame(self.root, bd=5, relief=GROOVE)
            self.__stop_conversion = False

            self.spinboxes = self.create_spinboxes()
            self.labels = self.create_labels()
            self.buttons = self.create_buttons()
            self.intvars = self.create_intvars()

            self.set_start_values()
            self.make_gui()

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

    def load_prefix(self):
        pref = self.RunParser.get('ConverterFolders', 'rawprefix')
        if pref == 'long':
            return 'run{y}{m}0'.format(y=self.TestCampaign[2:4], m=self.TestCampaign[-2:])
        elif pref == 'short':
            return 'run00'
        else:
            return pref

    def load_raw_file(self):
        if self.RunParser.has_option('ConverterFolders', 'rawfolder'):
            file_dir = self.RunParser.get('ConverterFolders', 'rawfolder')
        else:
            file_dir = join(self.DataDir, self.TcDir, 'raw')
        if not dir_exists(file_dir):
            log_warning('Could not find the raw file directory: {d}'.format(d=file_dir))
        return file_dir

    def load_errors(self):
        path = join(self.RootFileDir, 'Errors{r}.txt'.format(r=str(self.RunNumber).zfill(3)))
        if file_exists(path):
            f = open(path, 'r')
            return [int(line.strip('\n')) for line in f.readlines()]
        else:
            return []

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
            if any(opt.endswith(ending) for ending in ['_range', '_region', '_range_drs4']):
                config[opt] = json.loads(self.RunParser.get('ROOTFILE_GENERATION', opt))
            elif opt not in ['excluded_runs']:
                config[opt] = self.RunParser.getint('ROOTFILE_GENERATION', opt)
        return OrderedDict(sorted(config.iteritems()))

    def get_run_info(self, run_number):
        try:
            f = open(self.run_info_path, "r")
            data = json.load(f)
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

    def find_raw_file(self, run_number, prnt=True):
        file_path = join(self.RawFileDir, '{pref}{run}.raw'.format(pref=self.raw_prefix, run=str(run_number).zfill(4)))
        if file_exists(file_path):
            return file_path
        else:
            if prnt:
                print file_path, 'does not exist!'
            return False

    def get_root_file_path(self):
        file_name = '{prefix}{run}.root'.format(prefix=self.root_prefix, run=str(self.RunNumber).zfill(4))
        return self.RootFileDir + '/' + file_name

    def get_tracking_file_path(self, run_number):
        file_name = '{prefix}{run}_withTracks.root'.format(prefix=self.root_prefix, run=str(run_number).zfill(4))
        return self.RootFileDir + '/' + file_name

    def get_final_file_path(self, run_number):
        return join(self.RootFileDir, 'TrackedRun{run:03d}.root'.format(run=run_number))

    def find_root_file(self, run_number):
        old_track_file = self.get_tracking_file_path(run_number)
        track_file = self.get_final_file_path(run_number)
        final_file = self.get_root_file_path()
        # print 'looking for:\n ', track_file, '\n ', final_file
        if file_exists(track_file):
            return 'found_file'
        elif file_exists(old_track_file):
            return 'found_old'
        elif file_exists(final_file):
            log_message('did not find tracking file --> need conversion')
            return 'found_untracked'
        else:
            log_message('did not find any matching root file --> need conversion')
            return False

    def convert_run(self, run_infos, run_number):
        if do_gui and self.Type == 'pad':
            self.__stop_conversion = False
            self.root.deiconify()
            self.root.mainloop()
            if self.__stop_conversion:
                return
        # check whether the root file w/ our w/o tracks already exist
        found_root_file = self.find_root_file(run_number)
        if found_root_file == 'found_file':
            return
        if found_root_file == 'found_old':
            self.__rename_tracking_file(run_number)
            return
        if not found_root_file:
            # remove all old pickle files for a new conversion
            self.remove_pickle_files(run_number)
            curr_dir = getcwd()
            # check if raw file exists
            raw_file_path = self.find_raw_file(run_number)
            assert raw_file_path
            # go to root directory
            chdir(self.RootFileDir)
            # prepare converter command
            conf_string = '-c {eudaq}/conf/{file}'.format(eudaq=self.EudaqDir, file=self.ConverterConfigFile)
            tree_string = '-t {tree}'.format(tree=self.ConverterTree)
            converter_cmd = '{eudaq}/bin/Converter.exe {tree} {conf} {raw}'.format(eudaq=self.EudaqDir, raw=raw_file_path, tree=tree_string, conf=conf_string)
            self.set_converter_configfile(run_infos)
            print_banner('START CONVERTING RAW FILE FOR RUN {0}'.format(run_number))
            print converter_cmd
            system(converter_cmd)
            sleep(1)
            while getsize(self.get_root_file_path()) < 500:
                remove(self.get_root_file_path())
                system(converter_cmd)
            chdir(curr_dir)
        self.align_run()
        self.__add_tracking(run_number)
        self.__rename_tracking_file(run_number)
        remove(self.get_root_file_path())

    def align_run(self):

        if self.Type == 'pad':
            print_banner('STARTING PAD EVENT ALIGNMENT OF RUN {r}'.format(r=self.RunNumber))
            pad_align = PadAlignment(self)
            if not pad_align.IsAligned:
                pad_align.write_aligned_tree()
        elif self.Type == 'pixel':
            print_banner('STARTING PIXEL EVENT ALIGNMENT OF RUN {r}'.format(r=self.RunNumber))
            pix_align = PixAlignment(self)
            if not pix_align.IsAligned and not pix_align.check_alignment():
                pix_align.write_aligned_tree()

    def remove_pickle_files(self, run_number):
        self.Run.log_info('Removing all pickle files for run {}'.format(run_number))
        files = glob(join(self.Run.Dir, 'Configuration', 'Individual_Configs', '*', '*{tc}*_{run}_*'.format(run=run_number, tc=self.Run.generate_tc_str())))
        for f in files:
            remove(f)

    def __rename_rootfile(self, run_number):
        rename(self.get_root_file_path(), self.get_final_file_path(run_number))

    def __rename_tracking_file(self, run_number):
        rename(self.get_tracking_file_path(run_number), self.get_final_file_path(run_number))

    def __add_tracking(self, run_number):
        root_file_path = self.get_root_file_path()
        curr_dir = getcwd()
        chdir(self.TrackingDir)
        tracking_cmd = '{dir}/TrackingTelescope {root} 0 {nr}{tr}'.format(dir=self.TrackingDir, root=root_file_path, nr=self.TelescopeID, tr='' if self.Type == 'pad' else ' 1')
        print '\nSTART TRACKING FOR RUN', run_number, '\n'
        print tracking_cmd
        system(tracking_cmd)
        chdir(curr_dir)
        # move file to data folder
        file_name = '/{prefix}{run}_withTracks.root'.format(prefix=self.root_prefix, run=str(run_number).zfill(4))
        path = self.TrackingDir + file_name
        move(path, self.RootFileDir)

    def load_polarities(self, info):
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

    def set_converter_configfile(self, run_infos):

        parser = ConfigParser()
        config_file = join(self.EudaqDir, 'conf', self.ConverterConfigFile)
        parser.read(config_file)
        section = 'Converter.{}'.format(self.ConverterTree)
        if self.Type == 'pad':
            parser.set(section, 'polarities', self.load_polarities(run_infos))
            parser.set(section, 'pulser_polarities', self.load_polarities(run_infos))

        # remove unset ranges and regions
        new_options = self.RunParser.options('ROOTFILE_GENERATION')
        for opt in parser.options(section):
            if (opt.endswith('_range') or opt.endswith('_region')) and opt not in new_options:
                parser.remove_option(section, opt)
        # set the new settings
        for key, value in self.config.iteritems():
            parser.set(section, key, value)

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

    # ============================================
    # WIDGETS
    def create_spinboxes(self):
        dic = OrderedDict()
        dic['min'] = OrderedDict()
        dic['max'] = OrderedDict()
        dic['single'] = OrderedDict()
        for key, value in self.config.iteritems():
            if type(value[0]) is str or type(value) is str:
                dic['single'][key] = Spinbox(self.frame, width=8, justify=CENTER, from_=0, to=16, increment=1)
            else:
                dic['min'][key] = Spinbox(self.frame, width=8, justify=CENTER, from_=0, to=1024, increment=10)
                dic['max'][key] = Spinbox(self.frame, width=8, justify=CENTER, from_=0, to=1024, increment=10)
        return dic

    def create_labels(self):
        dic = OrderedDict()
        dic['main'] = Label(self.frame, text='Converter', font='font/Font 16 bold')
        dic['config'] = OrderedDict()
        for key in self.config.iterkeys():
            dic['config'][key] = Label(self.frame, text=key, width=20, anchor=W)
        return dic

    def create_buttons(self):
        dic = OrderedDict()
        dic['start'] = Button(self.frame, text='Start', width=12, command=self.root.destroy)
        dic['save'] = Button(self.frame, text='Save Config', width=12, command=self.save_config_values)
        dic['stop'] = Button(self.frame, text='Stop', width=12, command=self.__stop_conversion)
        return dic

    @staticmethod
    def create_intvars():
        dic = OrderedDict()
        dic['tracking'] = IntVar()
        dic['tracking'].set(1)
        return dic

    def set_start_values(self):
        for key, value in self.config.iteritems():
            if type(value[0]) is str or type(value) is str:
                value = value if type(value) is str else value[0]
                self.spinboxes['single'][key].delete(0, 'end')
                self.spinboxes['single'][key].insert(0, value)
            else:
                value = value if type(value) is list else value[0]
                self.spinboxes['min'][key].delete(0, 'end')
                self.spinboxes['min'][key].insert(0, value[0])
                self.spinboxes['max'][key].delete(0, 'end')
                self.spinboxes['max'][key].insert(0, value[1])

    def save_config_values(self):
        for key, value in self.config.iteritems():
            if type(value[0]) is str or type(value) is str:
                self.config[key] = self.spinboxes['single'][key].get()
            else:
                lst = [int(self.spinboxes['min'][key].get()), int(self.spinboxes['max'][key].get())]
                self.config[key] = lst

        parser = ConfigParser()
        conf_file = 'Configuration/RunConfig_' + self.TestCampaign + '.cfg'
        parser.read(conf_file)
        for key, value in self.config.iteritems():
            value = value[0] if type(value) is tuple else value
            parser.set('ROOTFILE_GENERATION', key, value)
        # write changes
        f = open(conf_file, 'w')
        parser.write(f)
        f.close()

    def __stop_conversion(self):
        self.__stop_conversion = True
        self.root.destroy()
        return

    def make_gui(self):
        self.frame.pack()
        self.labels['main'].grid(columnspan=3)
        j = 0
        k = 0
        for i, label in enumerate(self.labels['config'].values()):
            label.grid(row=i + 1)
        for i, box in enumerate(self.spinboxes['min'].values()):
            box.grid(row=i + 1, column=1)
        for i, box in enumerate(self.spinboxes['max'].values()):
            box.grid(row=i + 1, column=2)
            j = i
        for i, box in enumerate(self.spinboxes['single'].values(), j + 2):
            box.grid(row=i, column=1, columnspan=2)
            k = i
        self.buttons['save'].grid(row=k + 2, pady=3, column=1, columnspan=2)
        self.buttons['stop'].grid(row=k + 3)
        self.buttons['start'].grid(row=k + 3, columnspan=2, column=1)

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
    zrun = Run(143, test_campaign='201707-2', tree=False, verbose=True)
    z = Converter(zrun)
    # run_info = z.get_run_info(run_number=393)
    # z.convert_run(run_info)
    # z.root.deiconify()
    # z.root.mainloop()
