#!/usr/bin/env python
# --------------------------------------------------------
#       pad run class
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ConfigParser import ConfigParser
from StringIO import StringIO
from collections import OrderedDict
from json import dumps, loads

from numpy import array, concatenate, cumsum

from run import Run, join
from utils import has_bit, critical, warning, ensure_dir, init_argparser


class PadRun(Run):
    """ Run class containing all the information for a single run. """

    def __init__(self, number=None, test_campaign=None, tree=True, t_vec=None, verbose=False):
        """
        :param number: if None is provided it creates a dummy run
        :param test_campaign: if None is provided ...
        :param tree: root_tree object, if None is given it will start the converter
        :param t_vec: time sequence of the run, if None is provide it will generate a corrected one
        :param verbose: turn on more output
        """

        Run.__init__(self, number, test_campaign, tree, t_vec, verbose)

        # Settings
        self.NChannels = self.load_pre_n_channels()
        self.Digitiser = self.Config.get('BASIC', 'digitizer').lower()
        self.Channels = self.load_channels()

        if self.Tree is not None:
            self.validate_region_information()
            self.TreeConfig = self.load_tree_config()
            self.DigitizerChannels = self.load_digitizer_channels()
            self.NChannels = len(self.DigitizerChannels)
            self.IntegralRegions = self.load_regions()
            self.PeakIntegrals = self.load_peak_integrals()
            self.TCal = self.load_tcal()
            self.TCalSum = cumsum(concatenate([[0], self.TCal, self.TCal])).astype('f4')
            self.NSamples = len(self.TCal)
            self.Channels = self.load_channels()

    def load_rootfile_dirname(self):
        return ensure_dir(join(self.TCDir, 'root', 'pads'))

    def load_pre_n_channels(self):
        return 9 if self.Config.get('BASIC', 'digitizer').lower() == 'caen' else 4

    def load_channels(self):
        if hasattr(self, 'TreeConfig'):
            return [i for i, val in enumerate(self.TreeConfig.get('General', 'active regions').split()) if int(val)]
        # read it from the current config file if there is no tree loaded
        binary = self.Config.getint('ROOTFILE_GENERATION', 'active_regions')
        return [i for i in xrange(self.NChannels) if has_bit(binary, i)]

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

    def load_digitizer_channels(self):
        return loads(self.TreeConfig.get('Sensor Names', 'Names'))

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
        return array(tcal[:1024], dtype='f4')

    def get_calibrated_time(self, trigger_cell, ibin):
        v = self.TCal[int(trigger_cell)]
        for i in xrange(ibin):
            v += self.TCal[(int(trigger_cell) + i + 1) % self.NSamples]
        return v

    def validate_region_information(self):
        if 'region_information' not in [key.GetName() for key in self.RootFile.GetListOfKeys()]:
            self.Converter.remove_final_file()
            critical('no region information in root tree (file not propely converted)')

    def wf_exists(self, channel):
        wf_exist = True if self.Tree.FindBranch('wf{ch}'.format(ch=channel)) else False
        if not wf_exist:
            warning('The waveform for channel {ch} is not stored in the tree'.format(ch=channel))
        return wf_exist

    def show_tree_config(self):
        macro = self.RootFile.Get('region_information')
        for line in macro.GetListOfLines():
            line = str(line)
            print 'tcal = {}'.format(array(loads(line.strip('tcal =')))) if line.startswith('tcal') else line


if __name__ == '__main__':
    args = init_argparser(run=23, tc='201908', tree=True)
    z = PadRun(args.run, tree=args.tree, test_campaign=args.testcampaign)
