from draw import *
from ConfigParser import ConfigParser
from glob import glob


# global test campaign
g_test_campaign = None


class Analysis(Draw):
    """ This class provides default behaviour objects in the analysis framework and is the parent of all analyses.
        It contains, among other things, the root drawing methods, the main config and information about the directory structure. """

    def __init__(self, testcampaign=None, verbose=False):

        self.InitTime = time()

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
        self.PBar = PBar()

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

    def print_testcampaign(self, pr=True):
        data = self.TCString.split('-')
        out = self.TestCampaign.strftime('%b %Y')
        subset = ' Part {}'.format(data[-1]) if len(data) > 1 else ''
        if pr:
            print '\nTESTCAMPAIGN: {}{}'.format(out, subset)
        return out

    def get_test_campaigns(self):
        return [basename(path).replace('_', '').strip('psi') for path in glob(join(self.DataDir, 'psi*'))]

    def verbose_print(self, *args):
        """ Print command if verbose is activated. Print each argument separated with a comma """
        if self.Verbose:
            print ', '.join(args)

    def info(self, msg, next_line=True, prnt=True):
        return info(msg, next_line, prnt=self.Verbose and prnt)

    def add_to_info(self, t, txt='Done'):
        return add_to_info(t, txt, prnt=self.Verbose)

    def print_start(self, run=None, prnt=True, tc=True, dut=None):
        if prnt:
            ana_name = self.__class__.__name__.replace('Analysis', '')
            run = ' FOR RUN{} {}'.format('PLAN' if 'Coll' in ana_name else '', run) if run is not None else ''
            tc = ' OF {}'.format(self.TCString) if tc else ''
            dia = '{} '.format(dut) if dut is not None else ''
            print_banner('STARTING {}{} ANALYSIS{}{}'.format(dia, ana_name.upper(), run, tc), symbol='~', color='green')

    def print_finished(self, prnt=True):
        if prnt:
            print_banner('Finished Instantiation in {}'.format(get_elapsed_time(self.InitTime)), color='green')

    def make_pickle_path(self, sub_dir, name=None, run=None, ch=None, suf=None, camp=None):
        ensure_dir(join(self.PickleDir, sub_dir))
        campaign = self.TCString if camp is None else camp
        run = '_{r}'.format(r=run) if run is not None else ''
        ch = '_{c}'.format(c=ch) if ch is not None else ''
        suf = '_{s}'.format(s=suf) if suf is not None else ''
        name = '{n}_'.format(n=name) if name is not None else ''
        return join(self.PickleDir, sub_dir, '{name}{tc}{run}{ch}{suf}.pickle'.format(name=name, tc=campaign, run=run, ch=ch, suf=suf))


if __name__ == '__main__':

    pargs = init_argparser(has_verbose=True)
    z = Analysis(pargs.testcampaign, verbose=pargs.verbose)
