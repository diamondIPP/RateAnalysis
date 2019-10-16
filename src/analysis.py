from ConfigParser import ConfigParser
from argparse import ArgumentParser
from glob import glob
from sys import stdout

from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar

from draw import *
from utils import *

# global test campaign
g_test_campaign = None


class Analysis(Draw):
    """ This class provides default behaviour objects in the analysis framework and is the parent of all analyses.
        It contains, among other things, the root drawing methods, the main config and information about the directory structure. """

    def __init__(self, testcampaign=None, verbose=False):

        # Configuration
        self.Verbose = verbose
        self.ConfigDir = join(get_base_dir(), 'Configuration')
        self.MainConfig = self.load_main_config()

        # Directories
        self.Dir = get_base_dir()
        self.PickleDir = join(self.Dir, self.MainConfig.get('SAVE', 'pickle directory'))
        self.DataDir = self.MainConfig.get('MAIN', 'data directory')

        # Test Campaign
        self.TCString = self.load_test_campaign(testcampaign)
        self.TestCampaign = datetime.strptime(self.TCString.split('-')[0], '%Y%M')

        # Analysis Config
        self.Config = self.load_config()

        # Drawing Class
        Draw.__init__(self, self.TCString, verbose=verbose, config=self.MainConfig)

        # Progress Bar
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
        self.ProgressBar = None

    @staticmethod
    def load_main_config():
        parser = ConfigParser()
        file_name = join(get_base_dir(), 'Configuration', 'main.ini')
        if not file_exists(file_name):
            log_critical('{} does not exist. Please copy it from the main.default and adapt it to your purpose!'.format(file_name))
        parser.read(file_name)
        return parser

    def load_config(self):
        parser = ConfigParser()
        file_name = join(self.ConfigDir, self.TCString, 'AnalysisConfig.ini')
        if not file_exists(file_name):
            log_critical('AnalysisConfig.ini does not exist for {0}! Please create it in Configuration/{0}!'.format(self.TestCampaign))
        parser.read(file_name)
        return parser

    def load_test_campaign(self, testcampaign):
        global g_test_campaign
        if g_test_campaign is None and testcampaign is None:
            g_test_campaign = self.MainConfig.get('MAIN', 'default test campaign')
        elif testcampaign is not None:
            g_test_campaign = testcampaign
        if g_test_campaign not in self.get_test_campaigns():
            critical('The Testcampaign {} does not exist!'.format(g_test_campaign))
        return g_test_campaign

    def set_test_campaign(self, testcampaign):
        self.TestCampaign = testcampaign

    def print_testcampaign(self, pr=True):
        info = self.TCString.split('-')
        out = self.TestCampaign.strftime('%b %Y')
        subset = ' Part {}'.format(info[-1]) if len(info) > 1 else ''
        if pr:
            print '\nTESTCAMPAIGN: {}{}'.format(out, subset)
        return out

    def get_test_campaigns(self):
        return [basename(path).replace('_', '').strip('psi') for path in glob(join(self.DataDir, 'psi*'))]

    # endregion

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=self.Widgets, maxval=n)
        self.ProgressBar.start()

    def verbose_print(self, *args):
        """ Print command if verbose is activated. Print each argument separated with a comma """
        if self.Verbose:
            print ', '.join(args)

    def log_info(self, msg, next_line=True, prnt=True):
        t1 = time()
        if prnt and self.Verbose:
            t = datetime.now().strftime('%H:%M:%S')
            print '{head} {t} --> {msg}'.format(head=colored('INFO:', 'cyan', attrs=['dark']), t=t, msg=msg),
            stdout.flush()
            if next_line:
                print
        return t1

    def add_info(self, t, msg='Done'):
        if self.Verbose:
            print '{m} ({t:2.2f} s)'.format(m=msg, t=time() - t)

    def make_pickle_path(self, sub_dir, name=None, run=None, ch=None, suf=None, camp=None):
        ensure_dir(join(self.PickleDir, sub_dir))
        campaign = self.TestCampaign if camp is None else camp
        run = '_{r}'.format(r=run) if run is not None else ''
        ch = '_{c}'.format(c=ch) if ch is not None else ''
        suf = '_{s}'.format(s=suf) if suf is not None else ''
        name = '{n}_'.format(n=name) if name is not None else ''
        return join(self.PickleDir, sub_dir, '{name}{tc}{run}{ch}{suf}.pickle'.format(name=name, tc=campaign, run=run, ch=ch, suf=suf))


if __name__ == '__main__':
    arg_parser = ArgumentParser(description='Basic analysis class')
    arg_parser.add_argument('-tc', '--testcampaign', nargs='?', default=None)
    arg_parser.add_argument('-v', '--verbose', action='store_true')
    pargs = arg_parser.parse_args()
    z = Analysis(pargs.testcampaign, verbose=pargs.verbose)
