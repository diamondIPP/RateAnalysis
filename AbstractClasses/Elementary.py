import os
import ROOT
from time import time
from ROOT import gROOT, TGraphErrors, TGaxis
import pickle
import sys
from glob import glob
import re


class Elementary(object):
    """
    The Elementary class provides default behaviour objects in the analysis framework and is the Mother of all myPadAnalysis objects.
    It provides, among other things, a verbose printing method or a save plot method containing a global save directory handling.
    """

    default_testcampaign = '201510'

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.save_directory = self.get_program_dir() + 'Results/'

        self.TESTCAMPAIGN = None
        self.set_test_campaign(self.default_testcampaign)

        self.load_config()
        self.aimedFluxes = [3, 20, 60, 600, 2000, 5000]
        # colors
        self.count = 0
        self.colors = self.create_colorlist()

    def load_config(self):
        pass

    @staticmethod
    def create_colorlist():
        col_names = [ROOT.kGreen, ROOT.kOrange, ROOT.kViolet, ROOT.kYellow, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta, ROOT.kAzure, ROOT.kCyan, ROOT.kTeal]
        colors = []
        for color in col_names:
            colors.append(color + 1)
        for color in col_names:
            colors.append(color + 3)
        return colors

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

    @staticmethod
    def has_bit(num, bit):
        assert (num >= 0 and type(num) is int), 'num has to be non negative int'
        return bool(num & 1 << bit)

    def set_save_directory(self, directory="Results/"):
        if not directory[-1] == "/":
            directory += "/"
        self.save_directory = directory

    def save_plots(self, savename, file_type=None, save_dir=None, sub_dir=None, canvas=None):
        """
        Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
        such as in SSl terminals, it is recommended to pass the canvas to the method.
        :param savename:
        :param file_type:
        :param save_dir:
        :param sub_dir:
        :param canvas:
        """
        save_dir = self.save_directory if save_dir is None else save_dir
        file_type = '.png' if file_type is None else '.{end}'.format(end=file_type)
        sub_dir = '' if sub_dir is None else '{subdir}/'.format(subdir=sub_dir)
        resultsdir = save_dir + sub_dir
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)
        if canvas is None:
            try:
                pad = gROOT.GetSelectedPad()
                canvas = pad.GetCanvas()
            except Exception as inst:
                print '\n\n{delim}\nERROR in save plots!\n{msg}\n{delim}\n\n'.format(delim=len(str(inst)) * '-', msg=inst)
                return
        try:
            canvas.SaveAs(resultsdir + savename + file_type)
        except Exception as inst:
            print '\n\n{delim}\nERROR in save plots!\n{msg}\n{delim}\n\n'.format(delim=len(str(inst)) * '-', msg=inst)

    def create_new_testcampaign(self):
        # todo query date and create all new necessary configfiles etc
        pass

    def set_test_campaign(self, campaign='201508'):
        campaigns = self.find_test_campaigns()
        if not str(campaign) in campaigns:
            print 'This Testcampaign does not exist yet! Use create_new_testcampaign!\nExisting campaigns: {camp}'.format(camp=campaigns)
            return
        self.TESTCAMPAIGN = str(campaign)
        print 'Testcampaign set to: {tc} '.format(tc=campaign)

    def find_test_campaigns(self):
        conf_dir = self.get_program_dir() + 'Configuration/'
        names = glob(conf_dir + 'RunConfig_*')
        campaigns = [re.split('_|\.', name)[1] for name in names]
        campaigns = [camp for i, camp in enumerate(campaigns) if camp not in campaigns[i + 1:]]
        return sorted(campaigns)

    @staticmethod
    def elapsed_time(start):
        string = str('{0:2.2f}'.format(time() - start)) + ' seconds'
        return string

    @staticmethod
    def do_pickle(path, function):
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
    def make_tgrapherrors(name, title, color=1, marker=20):
        gr = TGraphErrors()
        gr.SetTitle(title)
        gr.SetName(name)
        gr.SetMarkerStyle(marker)
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)
        return gr

    @staticmethod
    def make_tgaxis(x, y1, y2, title, color=1):
        a = TGaxis(x, y1, x, y2, y1, y2, 510, '+SU')
        a.SetLineColor(color)
        a.SetTickSize(0)
        a.SetLabelSize(0)
        a.SetTitleSize(0.04)
        a.SetTitleOffset(0.15)
        a.SetTitle(title + '  ')
        a.SetTitleColor(color)
        return a

    @staticmethod
    def get_program_dir():
        arg = 2 if len(sys.argv) > 2 else 0
        path = os.path.dirname(os.path.realpath(sys.argv[arg])).split('/')
        ret_val = ''
        for i in range(len(path) - 1):
            ret_val += path[i] + '/'
        return ret_val

if __name__ == "__main__":
    z = Elementary()
