#!/usr/bin/env python
# --------------------------------------------------------
#       pad run class
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from collections import OrderedDict
from json import loads
from numpy import array, concatenate, cumsum
from helpers.utils import has_bit, critical, warning, ensure_dir, init_argparser, Config
from src.run import Run, join


class PadRun(Run):
    """ Run class containing all the information for a single run. """

    def __init__(self, number=None, testcampaign=None, load_tree=True, verbose=False):
        """
        :param number: if None is provided it creates a dummy run
        :param testcampaign: if None is provided ...
        :param load_tree: root_tree object, if None is given it will start the converter
        :param verbose: turn on more output
        """

        Run.__init__(self, number, testcampaign, load_tree, verbose)

        # Settings
        self.NChannels = self.load_pre_n_channels()
        self.Digitiser = self.Config.get('BASIC', 'digitizer').lower()
        self.Channels = self.load_channels()

        if self.Tree.Hash():
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
        return [i for i in range(self.NChannels) if has_bit(binary, i)]

    def load_tree_config(self):
        info = [str(line).strip('\n') for line in self.RootFile.Get('region_information').GetListOfLines()]
        if not info[0].startswith('['):
            critical('very old data! rerun conversion!')
        return Config(info)

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
        for i in range(ibin):
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
            print('tcal = {}'.format(array(loads(line.strip('tcal =')))) if line.startswith('tcal') else line)


if __name__ == '__main__':
    args = init_argparser(run=88, tc='201908', tree=True, has_verbose=True)
    z = PadRun(args.run, testcampaign=args.testcampaign, load_tree=args.tree, verbose=args.verbose)
