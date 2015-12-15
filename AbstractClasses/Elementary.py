import os
import ROOT
import types as t
from time import time
from ROOT import gROOT
import pickle


class Elementary(object):
    """
    The Elementary class provides default behaviour objects in the analysis framework and is the Mother of all myPadAnalysis objects.
    It provides, among other things, a verbose printing method or a save plot method containing a global save directory handling.
    """
    GLOBAL_COUNT = 0
    colors = [2, 4, 3, 6, 7, 8, 9, 1, 30, 40, 28, 38, 39, 20, 41, 42, 43, 44, 45, 46, 47, 48, 49, 12, 13, 14, 15]
    GLOBAL_COLOR_INDEX = 0
    SaveDirectory = "Results/"
    TESTCAMPAIGN = ""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.showAndWait = False
        if self.TESTCAMPAIGN == "":
            Elementary.SetTestCampaign("201510")
            # print "No Testcampaign was set. Testcampaign is now set to: 201508"
            # print "To change Testcampaign: cls.SetTestCampaign(namestring)"
        self.LoadConfig()
        self.aimedFluxes = [3, 20, 60, 600, 2000, 5000]
        self.count = 0
        self.colors = self.create_colorlist()
        # self.colors = [2, 3, 4, ROOT.kOrange - 3, 28, 30, 41, 46, 44, ROOT.kGreen - 1, ROOT.kViolet + 4]

    def LoadConfig(self):
        pass

    def create_colorlist(self):
        col_names = [ROOT.kGreen, ROOT.kYellow, ROOT.kOrange, ROOT.kRed, ROOT.kMagenta, ROOT.kViolet, ROOT.kBlue, ROOT.kAzure, ROOT.kCyan, ROOT.kTeal]
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

    def VerbosePrint(self, *args):
        """
        Print command if verbose is activated.
        :param args:
        """
        if self.verbose:
            # Print each argument separately so caller doesn't need to put everything to be printed into a single string.
            for arg in args:
                print arg,
            print

    def IfWait(self, message):
        if self.showAndWait:
            raw_input(message)

    @staticmethod
    def has_bit(num, bit):
        assert (num >= 0 and type(num) is int), 'num has to be non negative int'
        return bool(num & 1 << bit)

    @classmethod
    def SetSaveDirectory(cls, directory="Results/"):
        """
        Sets the SaveDirectory globally.
        :param directory:
        """
        if not directory[-1] == "/":
            directory += "/"
        Elementary.SaveDirectory = directory

    @classmethod
    def SavePlots(cls, savename, ending=None, saveDir=None, subDir=None, canvas=None):
        """
        Saves the canvas at the desired location. If no canvas is passed as argument, the active canvas will be saved. However for applications without graphical interface,
        such as in SSl terminals, it is recommended to pass the canvas to the method.
        :param savename:
        :param ending:
        :param saveDir:
        :param subDir:
        :param canvas:
        """
        if saveDir is None:
            saveDir = Elementary.SaveDirectory
        if ending is None:
            ending = ""
        else:
            assert (type(ending) == t.StringType), "ending has to be string type"
            ending = "." + ending
        # Results directories:
        # resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        if subDir is None:
            subDir = ""
        else:
            if subDir[-1] != "/":
                subDir += "/"
        resultsdir = saveDir + subDir  # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)
        # print "PRINT: ", [resultsdir+savename+ending]
        if canvas is None:
            pad = ROOT.gROOT.GetSelectedPad()
            canvas = pad.GetCanvas()
        try:
            canvas.SaveAs(resultsdir + savename + ending)
        except Exception as inst:
            print "\n\n\n-----------------------------------"
            print inst
            print "ERROR in SAVE PLOTs !"
            print "tried to save:\n\t", [resultsdir + savename + ending]
            print "-----------------------------------\n\n\n"

    @classmethod
    def GC(cls):
        gc = cls.GLOBAL_COUNT
        cls.GLOBAL_COUNT += 1
        return gc

    @classmethod
    def get_new_color(cls):
        """
        Returns a new color number from the global color palette, which has a length of 27 colors.
        :return: color int
        """
        cls.GLOBAL_COLOR_INDEX %= 27
        color = cls.colors[cls.GLOBAL_COLOR_INDEX]
        cls.GLOBAL_COLOR_INDEX += 1
        return color

    @classmethod
    def ResetColorPalette(cls):
        """
        Resets the color palette, such that the next color which will be
        returned by the GetNewColor() method will be the first color in
        the color palette.
        :return:
        """
        cls.GLOBAL_COLOR_INDEX = 0

    @classmethod
    def SetTestCampaign(cls, namestring="201508"):
        """
        Sets the test campaign name globally.
        :param namestring:
        :return:
        """
        Elementary.TESTCAMPAIGN = str(namestring)
        print "Testcampaign set to: ", namestring

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

