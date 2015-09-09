import os
import ROOT
import types as t


class Elementary(object):
    '''
    The Elementary class provides default behaviour objects in
    the analysis framework and is the Mother of all myPadAnalysis objects.
    It provides, among other things, a verbose printing method or a
    save plot method containing a global save directory handling.
    '''
    GLOBAL_COUNT = 0
    colors = [2,4,3,6,7,8,9,1,30,40,28,38,39,20,41,42,43,44,45,46,47,48,49,12,13,14,15]
    GLOBAL_COLOR_INDEX = 0
    SaveDirectory = "Results/"
    TESTCAMPAIGN = ""

    def __init__(self, verbose = False):
        self.verbose = verbose
        self.showAndWait = False
        if self.TESTCAMPAIGN == "":
            Elementary.SetTestCampaign("201508")
            print "No Testcampaign was set. Testcampaign is now set to: 201508"
            print "To change Testcampaign: cls.SetTestCampaign(namestring)"
        self.LoadConfig()


    def LoadConfig(self):
        pass

    def VerbosePrint(self, *args):
        if self.verbose:
            # Print each argument separately so caller doesn't need to
            # stuff everything to be printed into a single string
            for arg in args:
               print arg,
            print
        else:
            pass

    def IfWait(self, message):
        if self.showAndWait:
            raw_input(message)

    def _GetBit(self, num, bit):
        assert(num>=0 and type(num) is t.IntType), "num as to be non negative int"
        return (num & 1<<bit) == 1<<bit


    @classmethod
    def SetSaveDirectory(cls, directory = "Results/"):
        if not directory[-1] == "/": directory+="/"
        Elementary.SaveDirectory = directory

    @classmethod
    def SavePlots(cls, savename, ending=None, saveDir=None, subDir=None, canvas=None):
        if saveDir == None:
            saveDir = Elementary.SaveDirectory
        if ending == None:
            ending = ""
        else:
            assert(type(ending) == t.StringType), "ending has to be string type"
            ending = "."+ending
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        if subDir == None:
            subDir=""
        else:
            if subDir[-1] != "/": subDir = subDir+"/"
        resultsdir = saveDir+subDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)
        print "PRINT: ", [resultsdir+savename+ending]
        if canvas == None:
            pad = ROOT.gROOT.GetSelectedPad()
            canvas = pad.GetCanvas()
        try:
            canvas.SaveAs(resultsdir+savename+ending)
        except:
            print "\n\n\n-----------------------------------"
            print "ERROR in SAVE PLOTs !"
            print "tried to save:\n\t", [resultsdir+savename+ending]
            print "-----------------------------------\n\n\n"


    @classmethod
    def GC(cls):
        gc = cls.GLOBAL_COUNT
        cls.GLOBAL_COUNT += 1
        return gc

    @classmethod
    def GetNewColor(cls):
        cls.GLOBAL_COLOR_INDEX %= 27
        color = cls.colors[cls.GLOBAL_COLOR_INDEX]
        cls.GLOBAL_COLOR_INDEX += 1
        return color

    @classmethod
    def ResetColorPalette(cls):
        cls.GLOBAL_COLOR_INDEX = 0

    @classmethod
    def SetTestCampaign(cls, namestring="201508"):
        Elementary.TESTCAMPAIGN = namestring
        print "Testcampaign set to: ", namestring
