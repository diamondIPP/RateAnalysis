import os
import ROOT
import types as t


class Elementary(object):
    '''
    Mother of all myPadAnalysis Objects
    '''
    GLOBAL_COUNT = 0
    SaveDirectory = "Results/"

    def __init__(self, verbose = False):
        self.verbose = verbose
        self.ShowAndWait = False
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
        if self.ShowAndWait:
            raw_input(message)

    @classmethod
    def SetSaveDirectory(cls, directory = "Results/"):
        Elementary.SaveDirectory = directory

    @classmethod
    def SavePlots(cls, savename, ending=None, saveDir=None):
        if saveDir == None:
            saveDir = Elementary.SaveDirectory
        if ending == None:
            ending = ""
        else:
            assert(type(ending) == t.StringType), "ending has to be string type"
            ending = "."+ending
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        resultsdir = saveDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)
        ROOT.gPad.Print(resultsdir+savename+ending)

    @classmethod
    def GC(cls):
        gc = cls.GLOBAL_COUNT
        cls.GLOBAL_COUNT += 1
        return gc