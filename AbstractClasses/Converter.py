# ==============================================
# IMPORTS
# ==============================================
import os
import json

from ConfigParser import ConfigParser
from math import copysign
from collections import OrderedDict
from re import sub
from shutil import move
from os import remove
from glob import glob
do_gui = False
if do_gui:
    from tkinter import *


__author__ = 'micha'


def print_banner(message):
    print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(message)) * '=', msg=message)


# ==============================================
# CLASS DEFINITION
# ==============================================
class Converter:
    def __init__(self, test_campaign):

        # main
        self.test_campaign = test_campaign
        self.parser = self.setup_configparser()
        self.Type = self.parser.get('BASIC', 'type')

        # tracking
        self.telescope_id = self.parser.getint('BASIC', 'telescopeID')
        self.tracking_dir = self.parser.get('ConverterFolders', 'trackingfolder')

        # directories
        self.raw_file_dir = self.parser.get('ConverterFolders', 'rawfolder')
        self.root_file_dir = self.parser.get('BASIC', 'runpath')
        self.eudaq_dir = self.parser.get('ConverterFolders', 'eudaqfolder')
        # files paths
        self.converter_config_path = self.parser.get('ConverterFolders', 'converterFile')
        self.run_info_path = self.parser.get('BASIC', 'runinfofile')
        # prefixes
        self.root_prefix = self.parser.get('ConverterFolders', "converterPrefix")
        self.raw_prefix = self.parser.get('ConverterFolders', "rawprefix")

        # configuration for pad
        self.config = self.get_config() if self.Type == 'pad' else None

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

    def get_config(self):
        config = OrderedDict()
        options = self.parser.options('ROOTFILE_GENERATION')
        for opt in options:
            if opt.endswith('_range') or opt.endswith('_region'):
                config[opt] = json.loads(self.parser.get('ROOTFILE_GENERATION', opt))
        config['pulser_range_drs4'] = json.loads(self.parser.get('ROOTFILE_GENERATION', 'pulser_range_drs4')),
        config['save_waveforms'] = self.parser.get('ROOTFILE_GENERATION', 'save_waveforms'),
        config['pulser_drs4_threshold'] = self.parser.get('ROOTFILE_GENERATION', 'pulser_drs4_threshold'),
        config['pulser_channel'] = self.parser.get('ROOTFILE_GENERATION', 'pulser_channel'),
        config['trigger_channel'] = self.parser.get('ROOTFILE_GENERATION', 'trigger_channel')
        return config

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

    def find_raw_file(self, run_number):
        file_path = self.raw_file_dir + '/{pref}{run}.raw'.format(pref=self.raw_prefix, run=str(run_number).zfill(4))
        if os.path.exists(file_path):
            return file_path
        else:
            print file_path, 'does not exist!'
            return False

    def get_root_file_path(self, run_number):
        file_name = '{prefix}{run}.root'.format(prefix=self.root_prefix, run=str(run_number).zfill(4))
        return self.root_file_dir + '/' + file_name

    def get_tracking_file_path(self, run_number):
        file_name = '{prefix}{run}_withTracks.root'.format(prefix=self.root_prefix, run=str(run_number).zfill(4))
        return self.root_file_dir + '/' + file_name

    def get_final_file_path(self, run_number):
        return '{dir}/TrackedRun{run:03d}.root'.format(run=run_number, dir=self.root_file_dir)

    def find_root_file(self, run_number):
        old_track_file = self.get_tracking_file_path(run_number)
        track_file = self.get_final_file_path(run_number)
        final_file = self.get_root_file_path(run_number)
        # print 'looking for:\n ', track_file, '\n ', final_file
        if os.path.exists(track_file):
            return 'found_file'
        elif os.path.exists(old_track_file):
            return 'found_old'
        elif os.path.exists(final_file):
            print 'did not find tracking file --> need conversion'
            return 'found_untracked'
        else:
            print 'did not find any matching root file --> need conversion'
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
            curr_dir = os.getcwd()
            # check if raw file exists
            raw_file_path = self.find_raw_file(run_number)
            assert raw_file_path
            # go to root directory
            os.chdir(self.root_file_dir)
            # prepare converter command
            conf_string = '-c {eudaq}/conf/{file}'.format(eudaq=self.eudaq_dir, file=self.converter_config_path)
            tree_string = '-t {tree}'.format(tree='drs4tree' if self.Type == 'pad' else 'telescopetree')
            converter_cmd = '{eudaq}/bin/Converter.exe {tree} {conf} {raw}'.format(eudaq=self.eudaq_dir, raw=raw_file_path, tree=tree_string, conf=conf_string)
            if self.Type == 'pad':
                self.__set_converter_configfile(run_infos)
            print_banner('START CONVERTING RAW FILE FOR RUN {0}'.format(run_number))
            print converter_cmd
            os.system(converter_cmd)
            os.chdir(curr_dir)
        self.remove_pickle_files(run_number)
        self.__add_tracking(run_number)
        self.__rename_tracking_file(run_number)
        os.remove(self.get_root_file_path(run_number))

    @staticmethod
    def remove_pickle_files(run_number):
        program_dir = ''
        for i in __file__.split('/')[:-2]:
            program_dir += i + '/'
        files = glob('{prog}Configuration/Individual_Configs/*/*{run}*'.format(prog=program_dir, run=run_number))
        for _file in files:
            remove(_file)

    def __rename_tracking_file(self, run_number):
        os.rename(self.get_tracking_file_path(run_number), self.get_final_file_path(run_number))

    def __add_tracking(self, run_number):
        root_file_path = self.get_root_file_path(run_number)
        curr_dir = os.getcwd()
        os.chdir(self.tracking_dir)
        tracking_cmd = "{dir}/TrackingTelescope {root} 0 {nr}".format(dir=self.tracking_dir, root=root_file_path, nr=self.telescope_id)
        print '\nSTART TRACKING FOR RUN', run_number, '\n'
        print tracking_cmd
        os.system(tracking_cmd)
        os.chdir(curr_dir)
        # move file to data folder
        file_name = '/{prefix}{run}_withTracks.root'.format(prefix=self.root_prefix, run=str(run_number).zfill(4))
        path = self.tracking_dir + file_name
        shutil.move(path, self.root_file_dir)

    def __set_converter_configfile(self, run_infos):
        pol_dia1 = int(copysign(1, run_infos['hv dia1']))
        pol_dia2 = int(copysign(1, run_infos['hv dia2']))

        parser = ConfigParser()
        conf_file = '{eudaq}/conf/{file}'.format(eudaq=self.eudaq_dir, file=self.converter_config_path)
        parser.read(conf_file)
        parser.set('Converter.drs4tree', 'polarities', '[{pol1},0,0,{pol2}]'.format(pol1=pol_dia1, pol2=pol_dia2))
        for key, value in self.config.iteritems():
            parser.set('Converter.drs4tree', key, value)

        # write changes
        f = open(conf_file, 'w')
        parser.write(f)
        f.close()

        # remove whitespaces and correct for lower mistakes
        f = open(conf_file, 'r+')
        content = f.readlines()
        for i, line in enumerate(content):
            line = line.replace('peaki', 'PeakI')
            line = sub('[)(\' ]', '', line)
            if len(line) > 3 and line[-2] == ',':
                line = line[:-2] + '\n'
            content[i] = line
        f.seek(0)
        f.writelines(content)
        f.truncate()
        f.close()

    def setup_configparser(self):
        conf = ConfigParser()
        conf.read('Configuration/RunConfig_' + self.test_campaign + '.cfg')
        return conf

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
        conf_file = 'Configuration/RunConfig_' + self.test_campaign + '.cfg'
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


if __name__ == "__main__":
    z = Converter('201510')
    run_info = z.get_run_info(run_number=393)
    # z.convert_run(run_info)
    z.root.deiconify()
    z.root.mainloop()
