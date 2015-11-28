__author__ = 'micha'

# ==============================================
# IMPORTS
# ==============================================
import os
import json

from ConfigParser import ConfigParser
from math import copysign


# ==============================================
# CLASS DEFINITION
# ==============================================
class Converter:
    def __init__(self, test_campaign, run_number):

        # main
        self.test_campaign = test_campaign
        self.run_number = run_number
        self.parser = self.setup_configparser()

        # directories
        self.raw_file_dir = self.parser.get('ROOTFILE_GENERATION', 'rawfolder')
        self.root_file_dir = self.parser.get('BASIC', 'runpath')
        self.eudaq_dir = self.parser.get('ROOTFILE_GENERATION', 'eudaqfolder')
        # files pathes
        self.converter_config_path = self.parser.get('ROOTFILE_GENERATION', 'converterFile')
        self.run_info_path = self.parser.get('BASIC', 'runinfofile')
        # prefixes
        self.root_prefix = self.parser.get('ROOTFILE_GENERATION', "converterPrefix")
        self.raw_prefix = self.parser.get('ROOTFILE_GENERATION', "rawprefix")

        # configuration
        self.config = self.get_config()

    # todo make small gui to enter parameters
    def get_config(self):
        config = {
            'signal_range': json.loads(self.parser.get('ROOTFILE_GENERATION', 'signal_range')),
            'signal_a_region': json.loads(self.parser.get('ROOTFILE_GENERATION', 'signal_a_region')),
            'signal_b_region': json.loads(self.parser.get('ROOTFILE_GENERATION', 'signal_b_region')),
            "signal_c_region": json.loads(self.parser.get('ROOTFILE_GENERATION', 'signal_c_region')),
            "signal_d_region": json.loads(self.parser.get('ROOTFILE_GENERATION', 'signal_d_region')),
            "pedestal_range": json.loads(self.parser.get('ROOTFILE_GENERATION', 'pedestal_range')),
            "pedestal_a_region": json.loads(self.parser.get('ROOTFILE_GENERATION', 'pedestal_a_region')),
            "pedestal_aa_region": json.loads(self.parser.get('ROOTFILE_GENERATION', 'pedestal_aa_region')),
            "pedestal_b_region": json.loads(self.parser.get('ROOTFILE_GENERATION', 'pedestal_b_region')),
            "pedestal_c_region": json.loads(self.parser.get('ROOTFILE_GENERATION', 'pedestal_c_region')),
            "pulser_range": json.loads(self.parser.get('ROOTFILE_GENERATION', 'pulser_range')),
            "peakintegral1_range": json.loads(self.parser.get('ROOTFILE_GENERATION', 'peakintegral1_range')),
            "peakintegral2_range": json.loads(self.parser.get('ROOTFILE_GENERATION', 'peakintegral2_range')),
            "peakintegral3_range": json.loads(self.parser.get('ROOTFILE_GENERATION', 'peakintegral3_range')),
            "pulser_range_drs4": json.loads(self.parser.get('ROOTFILE_GENERATION', 'pulser_range_drs4')),
            "save_waveforms": self.parser.get('ROOTFILE_GENERATION', 'save_waveforms'),
            "pulser_drs4_threshold": self.parser.get('ROOTFILE_GENERATION', 'pulser_drs4_threshold'),
            "pulser_channel": self.parser.get('ROOTFILE_GENERATION', 'pulser_channel'),
            "trigger_channel": self.parser.get('ROOTFILE_GENERATION', 'trigger_channel'),
        }
        return config

    def get_run_info(self):
        try:
            f = open(self.run_info_path, "r")
            data = json.load(f)
            f.close()
        except IOError:
            print 'could not read', self.run_info_path
            return
        run_infos = data.get(str(self.run_number))
        try:
            run_infos['hv dia1'] = run_infos.pop('dia1hv')
            run_infos['hv dia2'] = run_infos.pop('dia2hv')
        except KeyError:
            pass
        return run_infos

    def find_raw_file(self):
        file_path = self.raw_file_dir + '/{pref}{run}.raw'.format(pref=self.raw_prefix, run=str(self.run_number).zfill(4))
        if os.path.exists(file_path):
            return file_path
        else:
            print file_path, 'does not exist!'
            return False

    def find_root_file(self):
        file_name = '{prefix}{run}.root'.format(prefix=self.root_prefix, run=str(self.run_number).zfill(4))
        local_file = os.getcwd() + '/' + file_name
        final_file = self.root_file_dir + '/' + file_name
        print 'looking for:\n ', local_file, '\n ', final_file
        if os.path.exists(local_file):
            print 'found file in local folder'
            return 'local'
        elif os.path.exists(final_file):
            print 'found file in final folder'
            return 'final'
        else:
            print 'did not find any matching root file --> need conversion'
            return False

    def convert_run(self, run_infos):
        file_path = self.find_raw_file()
        if not file_path:
            return
        curr_dir = os.getcwd()
        # goto root directory
        os.chdir(self.root_file_dir)
        converter_cmd = '{eudaq}/bin/Converter.exe -t drs4tree -c {eudaq}/conf/{file} {raw}'.format(eudaq=self.eudaq_dir, file=self.converter_config_path, raw=file_path)
        self.__set_converter_configfile(run_infos)
        print '\n========================================'
        print 'START CONVERTING RAW FILE FOR RUN', self.run_number, '\n'
        print converter_cmd
        os.system(converter_cmd)
        os.chdir(curr_dir)

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
        f = open(conf_file, "w")
        parser.write(f)
        f.close()

    def setup_configparser(self):
        conf = ConfigParser()
        conf.read('Configuration/RunConfig_' + self.test_campaign + '.cfg')
        return conf


if __name__ == "__main__":
    z = Converter('201510', 393)
    run_info = z.get_run_info()
    z.convert_run(run_info)
