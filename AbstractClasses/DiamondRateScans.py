# --------------------------------------------------------
#       DIAMOND RATE SCANS
# created on June 24th 2016 by M. Reichmann
# --------------------------------------------------------

from AnalysisCollection import AnalysisCollection
from Elementary import Elementary
from ConfigParser import ConfigParser, NoOptionError
from Utils import *
from time import time
from argparse import ArgumentParser
from RunSelection import RunSelection
from json import load, dump
from collections import OrderedDict


class DiaScans(Elementary):
    def __init__(self, diamond, testcampaigns=None, verbose=False):
        Elementary.__init__(self, verbose=verbose)

        self.Parser = self.load_diamond_parser()

        # information
        self.DiamondName = self.load_diamond(diamond)
        self.TestCampaigns = self.load_tcs(testcampaigns)
        self.RunInfos = self.load_runinfos()
        self.AllRunPlans = self.load_all_runplans()
        self.RunPlans = self.find_diamond_runplans()

        # run plan selection
        self.Selections = self.load_selections()
        self.Selection = self.load_selection()


    # ==========================================================================
    # region INIT

    def load_selections(self):
        file_path = self.get_program_dir() + self.MainConfigParser.get('MISC', 'runplan_selection_file')
        f = open(file_path)
        selections = load(f)
        f.close()
        return selections

    def load_selection(self):
        return self.Selections[self.DiamondName] if self.DiamondName in self.Selections else {}

    def load_diamond_parser(self):
        parser = ConfigParser()
        parser.read('{0}/Configuration/DiamondAliases.cfg'.format(self.get_program_dir()))
        return parser

    def load_diamond(self, dia):
        try:
            return self.Parser.get('ALIASES', dia)
        except NoOptionError:
            log_warning('{0} is not a known diamond name! Please choose one from \n{1}'.format(dia, self.Parser.options('ALIASES')))
            exit()

    def load_tcs(self, tcs):
        if tcs is None:
            return ['201508', '201510']
        valid_tcs = self.find_test_campaigns()
        tcs = [tcs] if type(tcs) is not list else tcs
        if not all(tc in valid_tcs for tc in tcs):
            log_warning('You entered and invalid test campaign! Aborting!')
            exit()
        else:
            return tcs

    def load_all_runplans(self):
        runplan_path = self.get_program_dir() + self.MainConfigParser.get('MISC', 'runplan_file')
        f = open(runplan_path, 'r')
        runplans = load(f)
        f.close()
        return runplans

    def load_runinfos(self):
        run_infos = {}
        for tc in self.TestCampaigns:
            self.TESTCAMPAIGN = tc
            parser = self.load_run_configs(None)
            file_path = parser.get('BASIC', 'runinfofile')
            f = open(file_path)
            run_infos[tc] = load(f)
            f.close()
        return run_infos

    def find_diamond_runplans(self):
        runplans = {}
        for tc in self.TestCampaigns:
            runplans[tc] = {}
            for rp, runs in self.AllRunPlans[tc]['rate_scan'].iteritems():
                for ch in [1, 2]:
                    if all(self.DiamondName == self.load_diamond(self.RunInfos[tc][str(run)]['dia{0}'.format(ch)]) for run in runs):
                        bias = self.RunInfos[tc][str(runs[0])]['dia{0}hv'.format(ch)]
                        if all(self.RunInfos[tc][str(run)]['dia{0}hv'.format(ch)] == bias for run in runs):
                            if bias not in runplans[tc]:
                                runplans[tc][bias] = {}
                            runplans[tc][bias][rp] = ch
        return runplans

    # endregion

    # ==========================================================================
    # region SHOW

    def show_runplans(self):
        for tc, vals in self.RunPlans.iteritems():
            print_small_banner(tc.rjust(15))
            for bias, val1s in vals.iteritems():
                print '{0} V:'.format(int(bias))
                for rp, ch in val1s.iteritems():
                    print ' ', rp.ljust(5), ch
                print

    def show_selection(self):
        if self.Selection:
            for tc, rps in sorted(self.Selection.iteritems()):
                print_small_banner((tc + ':').ljust(15))
                for rp, ch in sorted(rps.iteritems()):
                    runs = self.get_runs(rp, tc)
                    print rp.ljust(5), '{0}-{1}'.format(str(runs[0]).zfill(3), str(runs[-1]).zfill(3))
        else:
            log_warning('Selection is empty!')

    def get_runs(self, rp, tc):
        return self.AllRunPlans[tc]['rate_scan'][rp]

    def show_all_runplans(self):
        for tc in self.TestCampaigns:
            print_small_banner(tc)
            for rp, runs in sorted(self.AllRunPlans[tc]['rate_scan'].iteritems()):
                dias = [self.load_diamond(self.RunInfos[tc][str(runs[0])]['dia{0}'.format(ch)]) for ch in [1, 2]]
                print rp.ljust(5), '{0}-{1}'.format(str(runs[0]).zfill(3), str(runs[-1]).zfill(3)), dias[0].ljust(11), dias[1].ljust(11)

    # endregion

    # ==========================================================================
    # region SELECTION
    def select_runplan(self, runplan, ch=1, testcampaign=None):
        rp = self.make_runplan_string(runplan)
        tc = str(testcampaign) if testcampaign is not None else self.TestCampaigns[-1]
        if rp in self.AllRunPlans[tc]['rate_scan'].keys():
            if tc not in self.Selection:
                self.Selection[tc] = {}
            self.Selection[tc][rp] = ch
        else:
            log_warning('The runplan {0} does not exist!'.format(rp))

    def select_runplans_by_bias(self, value):
        self.clear_selection()
        for tc, vals in self.RunPlans.iteritems():
            for bias, rps in vals.iteritems():
                if abs(bias) == value:
                    for rp, ch in rps.iteritems():
                        self.select_runplan(rp, ch, tc)

    def clear_selection(self):
        self.Selection = {}

    def save_selection(self, name):
        file_path = self.get_program_dir() + self.MainConfigParser.get('MISC', 'runplan_selection_file')
        f = open(file_path, 'r+')
        selections = load(f)
        if self.Selection:
            selections[name] = self.Selection
        else:
            log_warning('Selection is empty!')
        f.seek(0)
        dump(selections, f, indent=2, sort_keys=True)
        f.truncate()
        f.close()

    # endregion

    @staticmethod
    def make_runplan_string(nr):
        nr = str(nr)
        return nr.zfill(2) if len(nr) <= 2 else nr.zfill(4)

    def draw_rate_scans(self):
        if not self.Selection:
            log_warning('Selection is empty!')
            return
        run_selections = self.load_run_selections()
        print run_selections

    def load_run_selections(self):
        run_selections = OrderedDict()
        for tc, rps in sorted(self.Selection.iteritems()):
            for rp, ch in sorted(rps.iteritems()):
                self.log_info('Loading runplan {0} of testcampaign {1}'.format(rp.rjust(4), datetime.strptime(tc, '%Y%m').strftime('%b %Y')))
                sel = RunSelection(tc)
                sel.select_runs_from_runplan(rp)
                run_selections[sel] = ch
        return run_selections

if __name__ == "__main__":
    main_parser = ArgumentParser()
    main_parser.add_argument('dia', nargs='?', default='S129')
    main_parser.add_argument('-tcs', nargs='?', default=None)
    args = main_parser.parse_args()
    print_banner('STARTING DIAMOND RATE SCAN COLLECTION OF DIAMOND {0}'.format(args.dia))

    z = DiaScans(args.dia, args.tcs, verbose=True)
