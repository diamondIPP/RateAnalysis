from ConfigParser import ConfigParser
from ROOT import TCanvas
from glob import glob
from json import loads
from numpy import inf
from shutil import copyfile
from sys import stdout

from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar

from Draw import Draw
from Utils import *

# global test campaign and resolution
g_test_campaign = None
g_resolution = None


class Elementary(Draw):
    """
    The Elementary class provides default behaviour objects in the analysis framework and is the Mother of all myPadAnalysis objects.
    It provides, among other things, a verbose printing method or a save plot method containing a global save directory handling.
    """

    def __init__(self, testcampaign=None, verbose=False, resolution=None):

        self.verbose = verbose

        # Config
        self.Dir = self.get_program_dir()
        self.MainConfigParser = self.load_main_config()
        self.PickleDir = join(self.Dir, self.MainConfigParser.get('SAVE', 'pickle directory'))
        self.DataDir = self.MainConfigParser.get('MAIN', 'data_dir')

        # screen resolution
        self.Res = self.load_resolution(resolution)

        # test campaign
        self.TESTCAMPAIGN = self.load_test_campaign(testcampaign)
        self.SubSet = None
        self.TCString = self.generate_tc_str()
        self.ResultsDir = self.generate_results_directory()

        # Drawing Class
        Draw.__init__(self, self)

        # read configuration files
        self.RunConfig = self.load_run_config()
        self.AnaConfig = self.load_ana_config()
        self.TCDir = self.generate_tc_directory()

        # progress bar
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
        self.ProgressBar = None

        # container for the ROOT objects
        self.ROOTObjects = []

    # ============================================

    def load_main_config(self):
        parser = ConfigParser()
        file_name = join(self.Dir, 'Configuration', 'main.ini')
        if not file_exists(file_name):
            log_critical('{} does not exist. Please copy it from the main.default and adapt it to your purpose!'.format(file_name))
        parser.read(file_name)
        return parser

    def load_ana_config(self):
        parser = ConfigParser()
        file_name = join(self.Dir, 'Configuration', self.TCString, 'AnalysisConfig.ini')
        if not file_exists(file_name):
            log_critical('AnalysisConfig.ini does not exist for {0}! Please create it in Configuration/{0}!'.format(self.TESTCAMPAIGN))
        parser.read(file_name)
        # also read dut config
        if hasattr(self, 'RunNumber') and self.RunNumber is not None:
            dut_type = self.RunConfig.get('BASIC', 'type')
            parser.read(join(self.Dir, 'Configuration', self.TCString, '{}Config.ini').format(dut_type.title()))
        return parser

    def load_run_config(self):
        run_number = self.RunNumber if hasattr(self, 'RunNumber') else None
        parser = ConfigParser({'excluded_runs': '[]'})  # add non default option
        base_file_name = join(self.Dir, 'Configuration', self.TCString, 'RunConfig.ini')
        if not file_exists(base_file_name):
            log_critical('RunConfig.ini does not exist for {0}! Please create it in Configuration/{0}!'.format(self.TESTCAMPAIGN))
        parser.read(base_file_name)  # first read the main config file with general information for all splits
        if parser.has_section('SPLIT') and run_number is not None:
            split_runs = [0] + loads(parser.get('SPLIT', 'runs')) + [inf]
            config_nr = next(i for i in xrange(1, len(split_runs)) if split_runs[i - 1] <= run_number < split_runs[i])
            parser.read(join(self.Dir, 'Configuration', self.TCString, 'RunConfig{nr}.ini'.format(nr=config_nr)))  # add the content of the split config
        return parser

    @staticmethod
    def load_resolution(resolution):
        global g_resolution
        if g_resolution is None:
            g_resolution = get_resolution() if resolution is None else round_down_to(resolution, 500)
        return g_resolution

    def load_mask_file_dir(self):
        if self.RunConfig.has_option('BASIC', 'maskfilepath'):
            file_path = self.RunConfig.get('BASIC', 'maskfilepath')
        else:
            file_path = join(self.DataDir, self.TCDir, 'masks')
        if not dir_exists(file_path):
            log_warning('Mask file directory does not exist ({})!'.format(file_path))
        return file_path

    def load_run_info_path(self):
        file_path = join(self.DataDir, self.TCDir, 'run_log.json')
        if not file_exists(file_path):
            log_critical('Run Log File: "{f}" does not exist!'.format(f=file_path))
        return file_path

    def generate_sub_set_str(self):
        return '-{0}'.format(self.SubSet) if self.SubSet is not None else ''

    def generate_tc_str(self):
        return '{tc}{s}'.format(tc=self.TESTCAMPAIGN, s=self.generate_sub_set_str())

    def generate_tc_directory(self):
        return 'psi_{y}_{m}'.format(y=self.TESTCAMPAIGN[:4], m=self.TESTCAMPAIGN[4:])

    def load_test_campaign(self, testcampaign):
        global g_test_campaign
        if g_test_campaign is None and testcampaign is None:
            g_test_campaign = self.MainConfigParser.get('MAIN', 'default test campaign')
        elif testcampaign is not None:
            g_test_campaign = testcampaign
        if g_test_campaign not in self.get_test_campaigns():
            critical('The Testcampaign {} does not exist!'.format(g_test_campaign))
        return g_test_campaign

    def set_test_campaign(self, testcampaign):
        self.TESTCAMPAIGN = testcampaign

    def print_testcampaign(self, pr=True):
        info = self.TESTCAMPAIGN.split('-')
        out = datetime.strptime(info[0], '%Y%m').strftime('%b %Y')
        subset = ' Part {}'.format(info[-1]) if len(info) > 1 else ''
        if pr:
            print '\nTESTCAMPAIGN: {}{}'.format(out, subset)
        return out

    # def get_test_campaigns(self):
    #     return loads(self.MainConfigParser.get('MAIN', 'test_campaigns'))

    def get_test_campaigns(self):
        return [basename(path).replace('_', '').strip('psi') for path in glob(join(self.DataDir, 'psi*'))]

    # endregion

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=self.Widgets, maxval=n)
        self.ProgressBar.start()

    def verbose_print(self, *args):
        """
        Print command if verbose is activated.
        :param args:
        """
        if self.verbose:
            # Print each argument separately so caller doesn't need to put everything to be printed into a single string.
            for arg in args:
                print arg,
            print

    def log_info(self, msg, next_line=True, prnt=True):
        t1 = time()
        if prnt and self.verbose:
            t = datetime.now().strftime('%H:%M:%S')
            print '{head} {t} --> {msg}'.format(head=colored('INFO:', 'cyan', attrs=['dark']), t=t, msg=msg),
            stdout.flush()
            if next_line:
                print
        return t1

    def add_info(self, t, msg='Done'):
        if self.verbose:
            print '{m} ({t:2.2f} s)'.format(m=msg, t=time() - t)

    def make_pickle_path(self, sub_dir, name=None, run=None, ch=None, suf=None, camp=None):
        ensure_dir(join(self.PickleDir, sub_dir))
        campaign = '{tc}{s}'.format(tc=self.TESTCAMPAIGN, s=self.generate_sub_set_str()) if camp is None else camp
        run = '_{r}'.format(r=run) if run is not None else ''
        ch = '_{c}'.format(c=ch) if ch is not None else ''
        suf = '_{s}'.format(s=suf) if suf is not None else ''
        name = '{n}_'.format(n=name) if name is not None else ''
        return '{dir}/{sdir}/{name}{tc}{run}{ch}{suf}.pickle'.format(dir=self.PickleDir, sdir=sub_dir, name=name, tc=campaign, run=run, ch=ch, suf=suf)

    def create_new_testcampaign(self):
        year = raw_input('Enter the year of the test campgaign (YYYY): ')
        month = raw_input('Enter the month of the testcampaign: ').zfill(2)
        if year + month in self.get_test_campaigns():
            print 'This test campaign already exists! --> returning'
            return
        new_tc = year + month
        new_tc_cfg = year + month + '.cfg'
        conf_dir = self.get_program_dir() + 'Configuration/'
        names = []
        old_tc_cfg = ''
        for f in glob(conf_dir + '*'):
            name = f.split('/')[-1].split('_')
            if len(name) > 1 and name[1].startswith('20') and name[0] not in names:
                if name[1] > old_tc_cfg:
                    old_tc_cfg = name[1]
                names.append(name[0])
        old_tc = old_tc_cfg.split('.')[0]
        for name in names:
            file_name = conf_dir + name + '_'
            copyfile(file_name + old_tc_cfg, file_name + new_tc_cfg)
            f = open(file_name + new_tc_cfg, 'r+')
            lines = []
            for line in f.readlines():
                print line
                print old_tc[2:], new_tc[2:], old_tc[:4] + '_' + old_tc[4:], year + '_' + month
                line = line.replace(old_tc[2:], new_tc[2:])
                old = old_tc[:4] + '_' + old_tc[4:]
                lines.append(line.replace(old, year + '_' + month))
            f.seek(0)
            f.writelines(lines)
            f.close()

    @staticmethod
    def calc_fwhm(histo):
        h = histo
        max_ = h.GetMaximum()
        bin1 = h.FindFirstBinAbove(max_ / 2)
        bin2 = h.FindLastBinAbove(max_ / 2)
        fwhm = h.GetBinCenter(bin2) - h.GetBinCenter(bin1)
        return fwhm

    @staticmethod
    def get_program_dir():
        return dirname(dirname(__file__))

    @staticmethod
    def adj_length(value):
        string = str(value)
        num = len(string) / 4 * 4 + 4
        return string.ljust(num)

    @staticmethod
    def fit_fwhm(histo, fitfunc='gaus', do_fwhm=True, draw=False):
        h = histo
        if do_fwhm:
            peak_pos = h.GetBinCenter(h.GetMaximumBin())
            bin1 = h.FindFirstBinAbove(h.GetMaximum() / 2)
            bin2 = h.FindLastBinAbove(h.GetMaximum() / 2)
            fwhm = h.GetBinLowEdge(bin2 + 2) - h.GetBinLowEdge(bin1 - 1)
            option = 'qs' if draw else 'qs0'
            fit = h.Fit(fitfunc, option, '', peak_pos - fwhm / 2, peak_pos + fwhm / 2)
        else:
            fit = h.Fit(fitfunc, 'qs')
        return FitRes(fit)

    @staticmethod
    def del_rootobj(obj):
        if obj is None:
            return
        try:
            if obj.IsA().GetName() != 'TCanvas':
                obj.Delete()
        except AttributeError:
            pass

    @staticmethod
    def make_class_from_instance(instance):
        copy = deepcopy(instance.__dict__)
        instance_factory = type('instance_factory', (instance.__class__, ), {})
        instance_factory.__init__ = lambda self, *args, **kwargs: self.__dict__.update(copy)
        return instance_factory

    def save_combined_pulse_heights(self, mg, mg1, mg_y, show=True, name=None, pulser_leg=None,
                                    x_range=None, y_range=None, rel_y_range=None, draw_objects=None):
        set_root_output(show)
        c = TCanvas('c', 'c', int(self.Res * 10 / 11.), self.Res)
        make_transparent(c)
        bm = .11
        scale = 1.5
        pm = bm + (1 - bm - .1) / 5.

        # set unified x-range:
        mg1.GetXaxis().SetLimits(1, 3e4) if x_range is None else do_nothing()
        mg.GetXaxis().SetLimits(1, 3e4) if x_range is None else do_nothing()

        # bottom pad with 20%
        p0 = self.draw_tpad('p0', 'p0', pos=[0, 0, 1, pm], margins=[.14, .03, bm / pm, 0], transparent=True, logx=True, gridy=True)
        scale_multigraph(mg1)
        rel_y_range = [.7, 1.3] if rel_y_range is None else rel_y_range
        self.format_histo(mg1, title='', y_range=rel_y_range, y_tit='Rel. ph [au]' if not scale > 1 else ' ', y_off=66, tit_size=.1 * scale, x_off=99, lab_size=.1 * scale)
        mg1.GetYaxis().SetNdivisions(3)
        hide_axis(mg1.GetXaxis())
        mg1.Draw('alp')
        x_range = [mg1.GetXaxis().GetXmin(), mg1.GetXaxis().GetXmax()] if x_range is None else x_range
        self.draw_x_axis(1.3, x_range[0], x_range[1], mg1.GetXaxis().GetTitle() + ' ', opt='SG+-=', tit_size=.1, lab_size=.1 * scale, off=99, tick_size=.1, l_off=0)
        c.cd()

        # top pad with zero suppression
        self.draw_tpad('p1', 'p1', pos=[0, pm, 1, 1], margins=[.14, .03, 0, .1], transparent=True, logx=True)
        mg.Draw('alp')
        hide_axis(mg.GetXaxis())
        if pulser_leg:
            pulser_leg()
        if y_range:
            mg.SetMinimum(y_range[0])
            mg.SetMaximum(y_range[1])
        self.format_histo(mg, tit_size=.04 * scale, y_off=1.75 / scale, lab_size=.04 * scale)
        self.draw_x_axis(mg_y, x_range[0], x_range[1], mg1.GetXaxis().GetTitle() + ' ', opt='SG=', tit_size=.035 * scale, lab_size=0, off=1, l_off=99)
        l = mg.GetListOfFunctions()[0]
        move_legend(l, .17, .03)
        l.Draw()
        if draw_objects is not None:
            for obj, opt in draw_objects:
                obj.Draw(opt)

        if hasattr(self, 'InfoLegend'):
            run_info = self.InfoLegend.draw(p0, all_pads=False)
            scale_legend(run_info[0], txt_size=.09, height=0.098 / pm)
            run_info[1].SetTextSize(.05)

        for obj in p0.GetListOfPrimitives():
            if obj.GetName() == 'title':
                obj.SetTextColor(0)
        self.save_canvas(c, name='CombinedPulseHeights' if name is None else name, show=show)

        self.ROOTObjects.append([c, draw_objects])
        set_root_output(True)


if __name__ == '__main__':
    z = Elementary()