import os
import pickle
import re
from copy import deepcopy
from glob import glob
from shutil import copyfile
from time import time
from ConfigParser import ConfigParser
from json import loads
from Utils import *
from screeninfo import get_monitors

from ROOT import gROOT, TGraphErrors, TGaxis, TLatex, TGraphAsymmErrors, TSpectrum, TF1, TMath, TCanvas, gStyle, TLegend, TLine
# global test campaign and resolution
tc = None
res = None
default_tc = '201510'


class Elementary(object):
    """
    The Elementary class provides default behaviour objects in the analysis framework and is the Mother of all myPadAnalysis objects.
    It provides, among other things, a verbose printing method or a save plot method containing a global save directory handling.
    """

    def __init__(self, testcampaign=None, verbose=False, resolution=None):
        self.verbose = verbose

        self.TESTCAMPAIGN = None
        self.set_global_testcampaign(testcampaign)
        self.results_directory = '{dir}/Results{tc}/'.format(dir=self.get_program_dir(), tc=self.TESTCAMPAIGN)

        # read configuration files
        self.MainConfigParser = self.load_main_config()
        self.run_config_parser = self.load_run_config()
        self.ana_config_parser = self.load_ana_config()

        self.Felix = self.MainConfigParser.get('SAVE', 'felix')

        self.Stuff = []

        # screen resolution
        self.Res = self.load_resolution(resolution)

        # container for the ROOT objects
        self.ROOTObjects = []

        # colors
        self.count = 0
        self.colors = self.create_colorlist()
        self.FillColor = 821
        gStyle.SetLegendFont(42)

    # ============================================
    # region CONFIG
    def load_main_config(self):
        parser = ConfigParser()
        parser.read('{dir}/Configuration/main.cfg'.format(dir=self.get_program_dir()))
        return parser

    def load_run_config(self):
        run_parser = ConfigParser({'excluded_runs': '[]'})
        run_parser.read('Configuration/RunConfig_{tc}.cfg'.format(tc=self.TESTCAMPAIGN))
        return run_parser

    def load_run_configs(self, run_number):
        run_parser = ConfigParser({'excluded_runs': '[]'})
        # set run_number to zero if none is given to prevent crash
        run_number = 0 if run_number is None else run_number
        if self.MainConfigParser.has_section(self.TESTCAMPAIGN):
            n_splits = self.MainConfigParser.getint(self.TESTCAMPAIGN, 'n_splits')
            split_runs = [0] + loads(self.MainConfigParser.get(self.TESTCAMPAIGN, 'split_runs')) + [int(1e10)]
            for i in xrange(1, n_splits + 1):
                if split_runs[i - 1] <= run_number < split_runs[i]:
                    run_parser.read('{dir}/Configuration/RunConfig_{tc}_pad{i}.cfg'.format(dir=self.get_program_dir(), tc=self.TESTCAMPAIGN, i=i))
                    break
        else:
            run_parser.read('Configuration/RunConfig_{tc}.cfg'.format(tc=self.TESTCAMPAIGN))
        return run_parser

    def load_ana_config(self):
        ana_parser = ConfigParser()
        ana_parser.read('Configuration/AnalysisConfig_{tc}.cfg'.format(tc=self.TESTCAMPAIGN))
        return ana_parser

    @staticmethod
    def load_resolution(resolution):
        if resolution is not None:
            global res
            res = resolution
        if res is not None:
            return round_down_to(res, 500)
        else:
            try:
                m = get_monitors()
                return round_down_to(m[0].height, 500)
            except Exception as err:
                log_warning(err)
                return 1000

    def set_save_directory(self, name):
        self.results_directory = '{dir}/{nam}/'.format(dir=self.get_program_dir(), nam=name)

    def set_global_testcampaign(self, testcampaign):
        if testcampaign is not None:
            global tc
            tc = testcampaign
        if tc is not None:
            self.set_test_campaign(tc)
        else:
            self.TESTCAMPAIGN = default_tc

    def set_test_campaign(self, campaign='201508'):
        campaigns = self.find_test_campaigns()
        if not str(campaign) in campaigns:
            print 'This Testcampaign does not exist yet! Use create_new_testcampaign!\nExisting campaigns: {camp}'.format(camp=campaigns)
            return
        self.TESTCAMPAIGN = str(campaign)

    def print_testcampaign(self, pr=True):
        out = datetime.strptime(self.TESTCAMPAIGN, '%Y%m')
        out = 'TESTCAMPAIGN: {0}'.format(out.strftime('%b %Y'))
        if pr:
            print out
        return out

    @classmethod
    def find_test_campaigns(cls):
        conf_dir = cls.get_program_dir() + 'Configuration/'
        f_location = conf_dir + 'RunConfig_*'
        names = glob(f_location)
        campaigns = [re.split('_|\.', name)[1] for name in names]
        campaigns = [camp for i, camp in enumerate(campaigns) if camp not in campaigns[i + 1:]]
        return sorted(campaigns)

    # endregion

    @staticmethod
    def create_colorlist():
        col_names = [ROOT.kGreen, ROOT.kOrange, ROOT.kViolet, ROOT.kYellow, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta, ROOT.kAzure, ROOT.kCyan, ROOT.kTeal]
        colors = []
        for color in col_names:
            colors.append(color + 1)
        for color in col_names:
            colors.append(color + 3)
        return colors

    @staticmethod
    def ensure_dir(f):
        d = os.path.dirname(f)
        if not os.path.exists(d):
            os.makedirs(d)

    def get_color(self):
        self.count %= 20
        color = self.colors[self.count]
        self.count += 1
        return color

    def reset_colors(self):
        self.count = 0

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

    def log_info(self, msg):
        if self.verbose:
            t = datetime.now().strftime('%H:%M:%S')
            print 'INFO: {t} --> {msg}'.format(t=t, msg=msg)

    @staticmethod
    def log_warning(msg):
        t = datetime.now().strftime('%H:%M:%S')
        print '{head} {t} --> {msg}'.format(t=t, msg=msg, head=colored('WARNING:', 'red'))

    @staticmethod
    def has_bit(num, bit):
        assert (num >= 0 and type(num) is int), 'num has to be non negative int'
        return bool(num & 1 << bit)

    def make_bias_string(self):
        if hasattr(self, 'bias'):
            bias = self.bias
            pol = 'm' if bias < 0 else 'p'
            return '_{pol}{bias:04d}'.format(pol=pol, bias=int(abs(bias)))
        else:
            return ''

    def make_info_string(self):
        info = ''
        if not self.MainConfigParser.getboolean('SAVE', 'short_name'):
            info = '_{dia}'.format(dia=self.diamond_name) if hasattr(self, 'diamond_name') else ''
            info += self.make_bias_string()
            info += '_{tc}'.format(tc=self.TESTCAMPAIGN)
            info = info.replace('-', '')
        return info

    def save_canvas(self, canvas, sub_dir=None, name=None, print_names=True):
        sub_dir = self.save_dir if hasattr(self, 'save_dir') and sub_dir is None else '{subdir}/'.format(subdir=sub_dir)
        canvas.Update()
        file_name = canvas.GetName() if name is None else name
        file_path = '{save_dir}{res}/{{typ}}/{file}'.format(res=sub_dir, file=file_name, save_dir=self.results_directory)
        ftypes = ['root', 'png', 'pdf', 'eps']
        out = 'Saving plots: {nam}'.format(nam=name)
        run_number = self.run_number if hasattr(self, 'run_number') else None
        run_number = 'rp{nr}'.format(nr=self.run_plan) if hasattr(self, 'run_plan') else run_number
        file_path += self.make_info_string()
        gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
        for f in ftypes:
            ext = '.{typ}'.format(typ=f)
            if not f == 'png' and run_number is not None:
                ext = '_{run}.{typ}'.format(run=run_number, typ=f)
            self.ensure_dir(file_path.format(typ=f))
            out_file = '{fname}{ext}'.format(fname=file_path, ext=ext)
            out_file = out_file.format(typ=f)
            canvas.SaveAs(out_file)
        if print_names:
            log_message(out)
        gROOT.ProcessLine("gErrorIgnoreLevel = kError;")

    def save_plots(self, savename, sub_dir=None, canvas=None, ind=0, ch='dia', x=1, y=1, prnt=True):
        """
        Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
        such as in SSl terminals, it is recommended to pass the canvas to the method.
        :param savename:
        # :param file_type:
        # :param save_dir:
        :param sub_dir:
        :param canvas:
        :param ind: index of the collection
        :param ch: if None print both dias (dirty fix)
        """
        # save_dir = self.save_directory if save_dir is None else save_dir
        # file_type = '.png' if file_type is None else '.{end}'.format(end=file_type)

        if canvas is None:
            try:
                c = gROOT.GetListOfCanvases()
                # pad = gROOT.GetSelectedPad()
                canvas = c[-1]
            except Exception as inst:
                print '\n\n{delim}\nERROR in get canvas!\n{msg}\n{delim}\n\n'.format(delim=len(str(inst)) * '-', msg=inst)
                return
        channel = self.channel if hasattr(self, 'channel') else None
        if hasattr(self, 'run'):
            self.run.draw_run_info(channel=ch if ch is None else channel, canvas=canvas, x=x, y=y)
        if hasattr(self, 'Run'):
            self.Run.draw_run_info(channel=ch if ch is None else channel, canvas=canvas, x=x, y=y)
        elif hasattr(self, 'analysis'):
            try:
                self.analysis.run.draw_run_info(channel=ch if ch is None else channel, canvas=canvas, x=x, y=y)
            except AttributeError as err:
                self.log_warning(err)
        elif hasattr(self, 'collection'):
            runs = [self.collection.keys()[0], self.collection.keys()[-1], self.collection.values()[0].run.get_rate_string(), self.collection.values()[-1].run.get_rate_string()]
            if not ind:
                self.collection.values()[ind].run.draw_run_info(channel=ch if ch is None else self.collection.values()[ind].channel, canvas=canvas, runs=runs, x=x, y=y)
            else:
                self.collection.values()[ind].run.draw_run_info(channel=ch if ch is None else self.collection.values()[ind].channel, canvas=canvas, x=x, y=y)
        canvas.Update()

        try:
            self.save_canvas(canvas, sub_dir=sub_dir, name=savename, print_names=prnt)
        except Exception as inst:
            print self.print_banner('ERROR in save plots!\n{0}'.format(inst), '-')

    def create_new_testcampaign(self):
        year = raw_input('Enter the year of the test campgaign (YYYY): ')
        month = raw_input('Enter the month of the testcampaign: ').zfill(2)
        if year + month in self.find_test_campaigns():
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

    def print_elapsed_time(self, start, what='This', show=True):
        string = '{1} took {0:2.2f} seconds'.format(time() - start, what)
        self.print_banner(string) if show else self.do_nothing()
        return string

    @staticmethod
    def do_pickle(path, function, value=None):
        if value is not None:
            f = open(path, 'w')
            pickle.dump(value, f)
            f.close()
            return value
        try:
            f = open(path, 'r')
            ret_val = pickle.load(f)
            f.close()
        except IOError:
            ret_val = function()
            f = open(path, 'w')
            pickle.dump(ret_val, f)
            f.close()
        return ret_val

    @staticmethod
    def set_root_output(status=True):
        if status:
            gROOT.SetBatch(0)
            gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        else:
            gROOT.SetBatch(1)
            gROOT.ProcessLine("gErrorIgnoreLevel = kError;")

    @staticmethod
    def make_tgrapherrors(name, title, color=1, marker=20, marker_size=1, width=1, asym_err=False):
        gr = TGraphErrors() if not asym_err else TGraphAsymmErrors()
        gr.SetTitle(title)
        gr.SetName(name)
        gr.SetMarkerStyle(marker)
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)
        gr.SetMarkerSize(marker_size)
        gr.SetLineWidth(width)
        return gr

    @staticmethod
    def make_tgaxis(x1, x2, y1, y2, title, color=1, width=1, offset=.15, tit_size=.035, lab_size=0.035, line=True, opt='+SU'):
        ran = [y1, y2] if x1 == x2 else [x1, x2]
        a = TGaxis(x1, y1, x2, y2, ran[0], ran[1], 510, opt)
        a.SetName('ax')
        a.SetLineColor(color)
        a.SetLineWidth(width)
        a.SetLabelSize(lab_size)
        if line:
            a.SetTickSize(0)
            a.SetLabelSize(0)
        a.SetTitleSize(tit_size)
        a.SetTitleOffset(offset)
        a.SetTitle(title)
        a.SetTitleColor(color)
        a.SetLabelColor(color)
        return a

    def draw_axis(self, xmin, xmax, ymin, ymax, tit, col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035):
        a1 = self.make_tgaxis(xmin, xmax, ymin, ymax, tit, color=col, offset=off, line=False, opt=opt, width=w, tit_size=tit_size, lab_size=lab_size)
        a1.SetLabelOffset(0.01)
        a1.SetLabelFont(42)
        a1.SetTitleFont(42)
        a1.Draw()
        self.Stuff.append(a1)

    def draw_y_axis(self, x, ymin, ymax, tit, col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035):
        self.draw_axis(x, x, ymin, ymax, tit, col=col, off=off, opt=opt, w=w, tit_size=tit_size, lab_size=lab_size)

    def draw_x_axis(self, y, xmin, xmax, tit, col=1, off=1, w=1, opt='+L', tit_size=.035, lab_size=0.035):
        self.draw_axis(xmin, xmax, y, y, tit, col=col, off=off, opt=opt, w=w, tit_size=tit_size, lab_size=lab_size)

    def make_legend(self, x1=.58, y2=.88, nentries=2, w=.3, scale=1, name='l', felix=True):
        x2 = x1 + w
        y1 = y2 - nentries * .05 * scale
        l = TLegend(x1, y1, x2, y2)
        l.SetName(name)
        l.SetTextFont(42)
        l.SetTextSize(0.03 * scale)
        if self.Felix and felix:
            l.SetLineWidth(2)
            l.SetBorderSize(0)
            l.SetFillColor(0)
            l.SetFillStyle(0)
            l.SetTextAlign(12)
        return l

    def format_histo(self, histo, name='', title='', x_tit='', y_tit='', z_tit='', marker=20, color=1, markersize=1, x_off=1, y_off=1, z_off=1, lw=1, fill_color=0, stats=True,
                     tit_size=.04, draw_first=False, x_range=None, y_range=None):
        h = histo
        if draw_first:
            gROOT.SetBatch(1)
            h.Draw('a')
            gROOT.SetBatch(0)
        h.SetTitle(title) if title else h.SetTitle(h.GetTitle())
        h.SetName(name) if name else h.SetName(h.GetName())
        try:
            h.SetStats(stats)
        except AttributeError or ReferenceError:
            pass
        # markers
        try:
            h.SetMarkerStyle(marker)
            h.SetMarkerColor(color) if color is not None else h.SetMarkerColor(h.GetMarkerColor())
            h.SetMarkerSize(markersize)
        except AttributeError or ReferenceError:
            pass
        # lines/fill
        try:
            h.SetLineColor(color) if color is not None else h.SetLineColor(h.GetLineColor())
            h.SetFillColor(fill_color)
            h.SetLineWidth(lw)
        except AttributeError or ReferenceError:
            pass
        # axis titles
        try:
            x_tit = untitle(x_tit) if self.Felix else x_tit
            y_tit = untitle(y_tit) if self.Felix else y_tit
            z_tit = untitle(z_tit) if self.Felix else z_tit
            # x-axis
            x_axis = h.GetXaxis()
            x_axis.SetTitle(x_tit) if x_tit else h.GetXaxis().GetTitle()
            x_axis.SetTitleOffset(x_off)
            x_axis.SetTitleSize(tit_size)
            x_axis.SetRangeUser(x_range[0], x_range[1]) if x_range is not None else do_nothing()
            # y-axis
            y_axis = h.GetYaxis()
            y_axis.SetTitle(y_tit) if y_tit else y_axis.GetTitle()
            y_axis.SetTitleOffset(y_off)
            y_axis.SetTitleSize(tit_size)
            y_axis.SetRangeUser(y_range[0], y_range[1]) if y_range is not None else do_nothing()
            # z-axis
            h.GetZaxis().SetTitle(z_tit) if z_tit else h.GetZaxis().GetTitle()
            h.GetZaxis().SetTitleOffset(z_off)
            h.GetZaxis().SetTitleSize(tit_size)
        except AttributeError or ReferenceError:
            pass

    def save_histo(self, histo, save_name='test', show=True, sub_dir=None, lm=.1, rm=0.1, bm=.15, tm=.1, draw_opt='', x_fac=None, y_fac=None,
                   l=None, logy=False, logx=False, logz=False, canvas=None, gridx=False, gridy=False, save=True, ch='dia', prnt=True):
        x_fac = self.Res if x_fac is None else int(x_fac * self.Res)
        y_fac = self.Res if y_fac is None else int(y_fac * self.Res)
        h = histo
        if not show:
            gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
            gROOT.SetBatch(1)
        c = TCanvas('c_{0}'.format(h.GetName()), h.GetTitle().split(';')[0], x_fac, y_fac) if canvas is None else canvas
        c.SetMargin(lm, rm, bm, tm)
        c.SetLogx() if logx else self.do_nothing()
        c.SetLogy() if logy else self.do_nothing()
        c.SetLogz() if logz else self.do_nothing()
        c.SetGridx() if gridx else self.do_nothing()
        c.SetGridy() if gridy else self.do_nothing()
        h.Draw(draw_opt)
        l.Draw() if l is not None else self.do_nothing()
        if save:
            self.save_plots(save_name, sub_dir=sub_dir, x=x_fac, y=y_fac, ch=ch, prnt=prnt)
        gROOT.SetBatch(0)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        lst = [c, h, l] if l is not None else [c, h]
        self.ROOTObjects.append(lst)
        return lst

    def draw_histo(self, histo, save_name='', show=True, sub_dir=None, lm=.1, rm=0.1, bm=.15, tm=.1, draw_opt='', x=None, y=None,
                   l=None, logy=False, logx=False, logz=False, canvas=None, gridy=False, gridx=False, ch='dia', prnt=True):
        return self.save_histo(histo, save_name, show, sub_dir, lm, rm, bm, tm, draw_opt, x, y, l, logy, logx, logz, canvas, gridx, gridy, False, ch, prnt)

    @staticmethod
    def make_tlatex(x, y, text, align=20, color=1, size=.05):
        l = TLatex(x, y, text)
        l.SetTextAlign(align)
        l.SetTextColor(color)
        l.SetTextSize(size)
        return l

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
        ret_val = ''
        for i in __file__.split('/')[:-2]:
            ret_val += i + '/'
        return ret_val

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
        return fit

    @staticmethod
    def del_rootobj(obj):
        if obj is None:
            return
        if obj.IsA().GetName() != 'TCanvas':
            obj.Delete()

    @staticmethod
    def normalise_histo(histo, to100=False):
        h = histo
        fac = 100 if to100 else 1
        h.Scale(fac / h.Integral(1, h.GetNbinsX()))
        return h

    @staticmethod
    def do_nothing():
        pass

    @staticmethod
    def triple_gauss_fit(histo, show=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        h = histo
        fit = TF1('fit', 'gaus(0) + gaus(3) + gaus(6)')
        s = TSpectrum(2)
        s.Search(h)
        fit.SetParLimits(0, .8 * s.GetPositionY()[1], 1.2 * s.GetPositionY()[1])
        fit.SetParLimits(1, s.GetPositionX()[1] - 10, s.GetPositionX()[1] + 10)
        fit.SetParLimits(2, 5, 50)
        fit.SetParLimits(3, .8 * s.GetPositionY()[0], 1.2 * s.GetPositionY()[0])
        fit.SetParLimits(4, s.GetPositionX()[0] - 5, s.GetPositionX()[0] + 5)
        fit.SetParLimits(5, 1, 10)
        fit.SetParLimits(6, 10, s.GetPositionY()[1])
        fit.SetParLimits(7, s.GetPositionX()[0], s.GetPositionX()[1])
        fit.SetParLimits(8, 1, 10)
        for i in xrange(5):
            h.Fit(fit, 'qs{0}'.format('' if show else '0'), '', -50, s.GetPositionX()[1])
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        return fit

    @staticmethod
    def make_class_from_instance(instance):
        copy = deepcopy(instance.__dict__)
        instance_factory = type('instance_factory', (instance.__class__, ), {})
        instance_factory.__init__ = lambda self, *args, **kwargs: self.__dict__.update(copy)
        return instance_factory

    @staticmethod
    def isfloat(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

    @staticmethod
    def find_graph_margins(graphs):
        extrema = [max([TMath.MaxElement(gr.GetN(), gr.GetY()) for gr in graphs]), min([TMath.MinElement(gr.GetN(), gr.GetY()) for gr in graphs])]
        return [extrema[1] - (extrema[0] - extrema[1]) * .1, extrema[0] + (extrema[0] - extrema[1]) * .1]

    @staticmethod
    def print_banner(msg, symbol='='):
        print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(msg)) * symbol, msg=msg)

if __name__ == "__main__":
    z = Elementary()
