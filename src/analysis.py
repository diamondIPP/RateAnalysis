from numpy import deg2rad, rad2deg, sort
from helpers.utils import *
from plotting.save import SaveDraw, join, Draw, format_histo, get_last_canvas
from helpers.info_legend import AnaInfo
from os.path import getsize
from os import getcwd, chdir
from subprocess import check_call


class Analysis(object):
    """ This class provides default behaviour objects in the analysis framework and is the parent of all analyses.
        It contains, among other things, the root drawing methods, the main config and information about the directory structure. """

    MainConfig = load_main_config()
    Verbose = False

    # Beam
    Momentum = MainConfig.getfloat('BEAM', 'momentum')
    PathLength = ufloat(MainConfig.getfloat('BEAM', 'path length'), 1)
    BeamFrequency = MainConfig.getfloat('BEAM', 'frequency') * 1e6  # [HZ]
    BunchSpacing = 1 / BeamFrequency * 1e9  # [ns]

    # Directories
    PickleDir = Dir.joinpath(MainConfig.get('SAVE', 'pickle directory'))
    DataDir = Path(MainConfig.get('Directories', 'data'))

    def __init__(self, testcampaign=None, results_dir=None, sub_dir='', pickle_dir='', verbose=None):

        self.InitTime = time()

        Analysis.Verbose = choose(verbose, Analysis.Verbose)
        self.PickleSubDir = self.PickleDir.joinpath(str(pickle_dir))

        # Test Campaign
        self.TCString = self.load_test_campaign(testcampaign)
        self.TestCampaign = datetime.strptime(self.TCString.split('-')[0], '%Y%m')
        self.TCDir = self.load_tc_directory()

        # Modules
        self.Config = self.load_config()
        self.InfoLegend = AnaInfo
        self.Draw = SaveDraw(self, choose(results_dir, self.TCString), sub_dir)
        self.PBar = PBar()

    def __repr__(self):
        return f'{self.__class__.__name__.strip("Analysis")} Analysis'

    @property
    def server_save_dir(self):
        return

    def load_config(self):
        file_name = Dir.joinpath('config', self.TCString, 'AnalysisConfig.ini')
        if not file_exists(file_name):
            critical('AnalysisConfig.ini does not exist for {0}! Please create it in config/{0}!'.format(self.TestCampaign))
        return Config(file_name)

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
            print('\nTESTCAMPAIGN: {}{}'.format(out, subset))
        return out

    def set_verbose(self, status: bool):
        self.Verbose = status
        for field in self.__dict__.values():
            if hasattr(field, 'Verbose'):
                field.Verbose = status
            if hasattr(field, '__dict__'):
                for subfield in field.__dict__.values():
                    if hasattr(subfield, 'Verbose'):
                        subfield.Verbose = status

    @staticmethod
    def get_test_campaigns():
        return [path.name.replace('_', '').strip('psi') for path in Analysis.DataDir.glob('psi*')]

    @staticmethod
    def find_testcampaign(default=None):
        return next((tc for tc in sorted(Analysis.get_test_campaigns(), reverse=True) if f'psi_{tc[:4]}_{tc[4:]}' in getcwd()), default) if default is None else default

    def load_tc_directory(self, data_dir=None):
        return join(choose(data_dir, default=self.DataDir), 'psi_{y}_{m}'.format(y=self.TCString[:4], m=self.TCString[4:]))

    def info(self, msg, endl=True, prnt=True):
        return info(msg, endl, prnt=self.Verbose and prnt)

    def add_to_info(self, t, txt='Done', prnt=True):
        return add_to_info(t, txt, prnt=self.Verbose and prnt)

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
        directory = self.PickleSubDir if sub_dir is None else self.PickleDir.joinpath(sub_dir)
        ensure_dir(directory)
        campaign = self.TCString if camp is None else camp
        dut = str(dut if dut is not None else self.DUT.Number if hasattr(self, 'DUT') and hasattr(self.DUT, 'Number') else '')
        return join(directory, '{}.pickle'.format('_'.join([v for v in [name, campaign, self.make_run_str(run), dut, str(suf)] if v])))

    def make_run_str(self, run=None):
        if run is not None:
            return str(getattr(self, str(run)) if hasattr(self, str(run)) else run)
        return self.RunPlan if hasattr(self, 'RunPlan') else str(self.Run.Number) if hasattr(self, 'Run') and hasattr(self.Run, 'Number') else str(self.Number) if hasattr(self, 'Number') else ''

    def make_hdf5_path(self, sub_dir, name=None, run=None, ch=None, suf=None, camp=None):
        return self.make_pickle_path(sub_dir, name, run, ch, suf, camp).replace('pickle', 'hdf5')

    def make_simple_hdf5_path(self, *args, **kwargs):
        return self.make_simple_pickle_path(*args, **kwargs).replace('pickle', 'hdf5')

    def get_meta_files(self, all_=False):
        runs = self.Runs if hasattr(self, 'Runs') else [self.Run.Number] if hasattr(self, 'Run') else []
        dut_nr = f'_{self.DUT.Number}' if hasattr(self, 'DUT') else ''
        rp_files = list(self.PickleDir.rglob(f'*{self.TCString}_{self.RunPlan}{dut_nr}*')) if hasattr(self, 'RunPlan') else []
        sel_files = [n for i in self.Info for n in self.PickleDir.joinpath('selections').rglob(f'*{i.TCString}_{i.RunPlan}_{i.DUTNr}*')] if hasattr(self, 'Info') and hasattr(self, 'Selection') else []
        return sel_files + rp_files + [f for run in runs for f in (self.PickleDir if all_ else self.PickleSubDir).rglob(f'*{self.TCString}_{run}{dut_nr}*')]

    def remove_metadata(self, all_subdirs=False):
        for f in self.get_meta_files(all_subdirs):
            remove_file(f)

    def remove_tc_metadata(self):
        files = glob(join(self.PickleDir, '*', f'*{self.TCString}*'))
        info(f'removing {len(files)} meta files with a total size of {make_byte_string(sum(getsize(f) for f in files))}')
        for f in glob(join(self.PickleDir, '*', f'*{self.TCString}*')):
            remove_file(f, prnt=False)

    def get_metadata_size(self, all_subdirs=True):
        info('total size of metadata: {}'.format(make_byte_string(sum(getsize(f) for f in self.get_meta_files(all_subdirs)))))

    # TODO: move to higher analysis
    def calc_time_difference(self, p=None, m1=M_MU, m2=M_PI):
        return t_diff(self.PathLength.n, choose(p, self.Momentum), m1, m2) % self.BunchSpacing

    def get_time_differences(self, s=None, p=None, spacing=None):
        t = array([t_diff(choose(s, self.PathLength.n), choose(p, self.Momentum), m1, m2) for m1, m2 in [(M_PI, M_E), (M_PI, M_MU), (M_E, M_PI), (M_MU, M_PI)]])
        return sort(t % choose(spacing, self.BunchSpacing))

    def draw_time_differences(self):
        masses = [(M_PI, M_E), (M_PI, M_MU), (M_E, M_PI), (M_MU, M_PI)]
        fs = [self.Draw.make_tf1(None, self.calc_time_difference, 210, 310, color=self.Draw.get_color(2), w=2, style=[1, 2][m2 == M_PI], npx=500, m1=m1, m2=m2) for m1, m2 in masses]
        for i, f in enumerate(fs):
            tit, xtit, ytit = 'TOFDiff', 'Momentum [MeV/c]', 'Difference in Time-of-Flight [ns]'
            self.Draw(f, tit, x_tit=xtit, y_tit=ytit, line_color=None, y_range=[0, round(self.BunchSpacing)], grid=True, c=None if not i else get_last_canvas(), draw_opt='same' if i else '')
            Draw.tlatex(self.Momentum * 1.02, f(self.Momentum), text=f'{f(self.Momentum):2.1f}', size=.04)
        self.Draw.legend(fs, ['e^{+}-#pi^{+}', '#mu^{+}-#pi^{+}', '#pi^{+}-e^{+}', '#pi^{+}-#mu^{+}'], 'l')
        Draw.vertical_line(self.Momentum, w=2, color=2)
        self.Draw.save_plots('TDiff')

    def get_decay_ratio(self, p=None, d=None):
        r = decay_ratio(self.Momentum if p is None else p, M_PI, self.PathLength if d is None else d, TAU_PI)
        print('{:1.1f}% of the particles are left...'.format(r * 100))
        return r

    def draw_decay_angle(self, p=None, show=True):
        def f(a, pars):
            return rad2deg(decay_angle(deg2rad(a[0]), m=M_PI, m1=M_MU, p=pars[0]))
        f = TF1('fda', f, 0, 180, 1)
        f.SetParameter(0, self.Momentum if p is None else p)
        format_histo(f, x_tit='Decay Angle [deg]', y_tit='Boosted Angle [deg]', y_off=1.2, color=self.Draw.get_color(6), lw=2)
        Draw.histo(f, show=show)
        return f

    def draw_decay_angles(self, momenta=None):
        momenta = arange(200, 301, 20) if momenta is None else array(momenta, dtype='d')
        leg = Draw.make_legend(nentries=momenta.size)
        plots = [self.draw_decay_angle(p, show=False) for p in momenta]
        c = Draw.histo(plots[0], grid=True, lm=.121)
        leg.AddEntry(plots[0], 'p = {} MeV/c'.format(momenta[0]), 'l')
        for p, f in zip(momenta[1:], plots[1:]):
            Draw.histo(f, draw_opt='same', canvas=c)
            leg.AddEntry(f, str(p), 'l')
        leg.Draw()

    def draw_decay_ratio(self, p=None, show=True):
        def f(d, pars):
            return decay_ratio(pars[0], M_PI, d[0], TAU_PI) * 100
        f = TF1('fdr', f, 6, 7, 1)
        f.SetParameter(0, self.Momentum if p is None else p)
        format_histo(f, x_tit='Travel Distance [m]', y_tit='Decay Ratio [%]', y_off=1.2, color=self.Draw.get_color(6), lw=2)
        Draw.histo(f, show=show)
        return f

    def draw_decay_ratios(self):
        leg = Draw.make_legend(nentries=6)
        graphs = [(p, self.draw_decay_ratio(p, show=False)) for p in arange(200, 301, 20)]
        c = Draw.histo(graphs[0][1], grid=True, lm=.12)
        leg.AddEntry(graphs[0][1], 'p = 200 MeV/c', 'l')
        for p, f in graphs[1:]:
            Draw.histo(f, draw_opt='same', canvas=c)
            leg.AddEntry(f, str(p), 'l')
        leg.Draw()

    @staticmethod
    def go2data(tc, pad=True):
        chdir(join(Analysis.DataDir, f'psi_{tc[:4]}_{tc[4:]}', 'root', 'pads' if pad else 'pixel'))
        check_call('/bin/bash')


if __name__ == '__main__':
    aparser = ArgumentParser()
    aparser.add_argument('tc')
    aparser.add_argument('-pix', action='store_false')  # noqa default pixel
    # aparser.add_argument('-pix', action='store_true')  # noqa default pad
    pargs = aparser.parse_args()

    Analysis.go2data(pargs.tc, not pargs.pix)
