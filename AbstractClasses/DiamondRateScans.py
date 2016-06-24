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


class DiaScans(Elementary):

    def __init__(self, diamond, testcampaigns=None, verbose=False):
        Elementary.__init__(self, verbose=verbose)

        self.DiamondName = self.load_diamond(diamond)
        self.TestCampaigns = self.load_tcs(testcampaigns)

    def load_diamond(self, dia):
        parser = ConfigParser()
        parser.read('{0}/Configuration/DiamondAliases.cfg'.format(self.get_program_dir()))
        try:
            return parser.get('ALIASES', dia)
        except NoOptionError:
            log_warning('{0} is not a known diamond name! Please choose one from \n{1}'.format(dia, parser.options('ALIASES')))
            dia = raw_input('Enter diamond name: ')
            try:
                return parser.get('ALIASES', dia)
            except NoOptionError:
                log_warning('Invalid entry! Aborting!')
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

if __name__ == "__main__":
    st = time()
    main_parser = ArgumentParser()
    main_parser.add_argument('dia', nargs='?', default='S129')
    main_parser.add_argument('-tcs', nargs='?', default=None)
    args = main_parser.parse_args()
    print_banner('STARTING DIAMOND RATE SCAN COLLECTION OF DIAMOND {0}'.format(args.dia))

    z = DiaScans(args.dia, args.tcs, verbose=True)
    z.print_elapsed_time(st, 'Instantiation')
