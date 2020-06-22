from draw import *
from ConfigParser import ConfigParser
from glob import glob
from numpy import deg2rad, rad2deg, arange, round_


class Analysis(Draw):
    """ This class provides default behaviour objects in the analysis framework and is the parent of all analyses.
        It contains, among other things, the root drawing methods, the main config and information about the directory structure. """

    def __init__(self, testcampaign=None, verbose=False):

        self.InitTime = time()

        # Configuration
        self.Verbose = verbose
        self.ConfigDir = join(get_base_dir(), 'config')
        self.MainConfig = self.load_main_config()
        self.Momentum = self.MainConfig.getfloat('BEAM', 'momentum')
        self.PathLength = self.MainConfig.getfloat('BEAM', 'path length')
        self.BeamFrequency = self.MainConfig.getfloat('BEAM', 'frequency') * 1e6  # [HZ]
        self.BunchSpacing = 1 / self.BeamFrequency * 1e9  # [ns]

        # Directories
        self.Dir = get_base_dir()
        self.PickleDir = join(self.Dir, self.MainConfig.get('SAVE', 'pickle directory'))
        self.DataDir = self.MainConfig.get('MAIN', 'data directory')
        self.PickleSubDir = ''

        # Test Campaign
        self.TCString = self.load_test_campaign(testcampaign)
        self.TestCampaign = datetime.strptime(self.TCString.split('-')[0], '%Y%m')

        # Analysis Config
        self.Config = self.load_config()

        # Drawing Class
        Draw.__init__(self, self.TCString, verbose=verbose, config=self.MainConfig)

        # Progress Bar
        self.PBar = PBar()

    @staticmethod
    def load_main_config():
        parser = ConfigParser()
        file_name = join(get_base_dir(), 'config', 'main.ini')
        if not file_exists(file_name):
            log_critical('{} does not exist. Please copy it from the main.default and adapt it to your purpose!'.format(file_name))
        parser.read(file_name)
        return parser

    def load_config(self):
        parser = ConfigParser()
        file_name = join(self.ConfigDir, self.TCString, 'AnalysisConfig.ini')
        if not file_exists(file_name):
            log_critical('AnalysisConfig.ini does not exist for {0}! Please create it in config/{0}!'.format(self.TestCampaign))
        parser.read(file_name)
        return parser

    def load_test_campaign(self, testcampaign):
        tc = testcampaign if testcampaign is not None else self.MainConfig.get('MAIN', 'default test campaign')
        if tc not in self.get_test_campaigns():
            critical('The Testcampaign {} does not exist!'.format(tc))
        return tc

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

    def set_pickle_sub_dir(self, name):
        self.PickleSubDir = name

    def make_pickle_path(self, sub_dir=None, name=None, run=None, ch=None, suf=None, camp=None):
        ensure_dir(join(self.PickleDir, sub_dir))
        campaign = self.TCString if camp is None else camp
        run = '_{r}'.format(r=run) if run is not None else ''
        ch = '_{c}'.format(c=ch) if ch is not None else ''
        suf = '_{s}'.format(s=suf) if suf is not None else ''
        name = '{n}_'.format(n=name) if name is not None else ''
        return join(self.PickleDir, sub_dir, '{name}{tc}{run}{ch}{suf}.pickle'.format(name=name, tc=campaign, run=run, ch=ch, suf=suf))

    def make_simple_pickle_path(self, name='', suf='', sub_dir=None, run=None, dut=None, camp=None):
        directory = join(self.PickleDir, self.PickleSubDir if sub_dir is None else sub_dir)
        ensure_dir(directory)
        campaign = self.TCString if camp is None else camp
        run_str = str(run) if run is not None else self.RunPlan if hasattr(self, 'RunPlan') else ''
        # noinspection PyUnresolvedReferences
        run_str = run_str if run is not None or run_str else str(self.Run.Number) if hasattr(self.Run, 'Number') else ''
        # noinspection PyUnresolvedReferences
        dut = str(dut if dut is not None else self.DUT.Number if hasattr(self, 'DUT') and hasattr(self.DUT, 'Number') else '')
        return join(directory, '{}.pickle'.format('_'.join([v for v in [name, campaign, run_str, dut, str(suf)] if v])))

    def make_hdf5_path(self, sub_dir, name=None, run=None, ch=None, suf=None, camp=None):
        return self.make_pickle_path(sub_dir, name, run, ch, suf, camp).replace('pickle', 'hdf5')

    def make_simple_hdf5_path(self, *args, **kwargs):
        return self.make_simple_pickle_path(*args, **kwargs).replace('pickle', 'hdf5')

    # TODO: move to higher analysis
    def calc_time_difference(self, m1, m2, p=None):
        return t_diff(self.PathLength, self.Momentum if p is None else p, m1, m2) % self.BunchSpacing

    def calc_td(self, p, pars):
        return t_diff(self.PathLength, p[0], pars[0], pars[1]) % self.BunchSpacing

    def get_time_differences(self, p=None):
        t0 = array([self.calc_time_difference(M_PI, m, p) for m in [M_MU, M_E]])
        return round_(sorted(abs(concatenate([t0, t0 - self.BunchSpacing]))), 1)

    def draw_time_differences(self):
        leg = self.make_legend(w=.2)
        for m, n in zip([M_E, M_MU, M_P], ['positron', 'muon', 'proton']):
            f = TF1('f{}'.format(n), self.calc_td, 200, 300, 2)
            f.SetParameters(M_PI, m)
            format_histo(f, title='Time Differences', x_tit='Momentum [MeV/c]', y_tit='Time Difference [ns]', y_off=1.3, color=self.get_color(), lw=2, y_range=[3, 18])
            self.draw_histo(f, grid=True, lm=.12, draw_opt='' if m == M_E else 'same', canvas=None if m == M_E else get_last_canvas())
            t0 = self.calc_time_difference(M_PI, m)
            self.draw_tlatex(self.Momentum * 1.02, t0, text='{:2.1f}'.format(t0), size=.04)
            leg.AddEntry(f, n, 'l')
        get_object('fproton').SetNpx(1000)
        self.draw_vertical_line(260, 0, 20, w=2, color=2)
        leg.Draw()
        self.reset_colors()

    def get_decay_ratio(self, p=None, d=None):
        r = decay_ratio(self.Momentum if p is None else p, M_PI, self.PathLength if d is None else d, TAU_PI)
        print('{:1.1f}% of the particles are left...'.format(r * 100))
        return r

    def draw_decay_angle(self, p=None, show=True):
        def f(a, pars):
            return rad2deg(decay_angle(deg2rad(a[0]), m=M_PI, m1=M_MU, p=pars[0]))
        f = TF1('fda', f, 0, 180, 1)
        f.SetParameter(0, self.Momentum if p is None else p)
        format_histo(f, x_tit='Decay Angle [deg]', y_tit='Boosted Angle [deg]', y_off=1.2, color=self.get_color(), lw=2)
        self.draw_histo(f, show=show)
        return f

    def draw_decay_angles(self, momenta=None):
        momenta = arange(200, 301, 20) if momenta is None else array(momenta, dtype='d')
        leg = self.make_legend(nentries=momenta.size)
        graphs = [(p, self.draw_decay_angle(p, show=False)) for p in momenta]
        c = self.draw_histo(graphs[0][1], grid=True, lm=.121)
        leg.AddEntry(graphs[0][1], 'p = {} MeV/c'.format(momenta[0]), 'l')
        for p, f in graphs[1:]:
            self.draw_histo(f, draw_opt='same', canvas=c)
            leg.AddEntry(f, str(p), 'l')
        leg.Draw()
        self.reset_colors()

    def draw_decay_ratio(self, p=None, show=True):
        def f(d, pars):
            return decay_ratio(pars[0], M_PI, d[0], TAU_PI) * 100
        f = TF1('fdr', f, 6, 7, 1)
        f.SetParameter(0, self.Momentum if p is None else p)
        format_histo(f, x_tit='Travel Distance [m]', y_tit='Decay Ratio [%]', y_off=1.2, color=self.get_color(), lw=2)
        self.draw_histo(f, show=show)
        return f

    def draw_decay_ratios(self):
        leg = self.make_legend(nentries=6)
        graphs = [(p, self.draw_decay_ratio(p, show=False)) for p in arange(200, 301, 20)]
        c = self.draw_histo(graphs[0][1], grid=True, lm=.12)
        leg.AddEntry(graphs[0][1], 'p = 200 MeV/c', 'l')
        for p, f in graphs[1:]:
            self.draw_histo(f, draw_opt='same', canvas=c)
            leg.AddEntry(f, str(p), 'l')
        leg.Draw()
        self.reset_colors()


if __name__ == '__main__':
    pargs = init_argparser(has_verbose=True)
    z = Analysis(pargs.testcampaign, verbose=pargs.verbose)
