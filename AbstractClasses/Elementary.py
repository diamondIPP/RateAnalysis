import os
import ROOT


class Elementary(object):
    '''
    Mother of all myPadAnalysis Objects
    '''

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

    def SavePlots(self, savename, ending, saveDir):
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        resultsdir = saveDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)
        ROOT.gPad.Print(resultsdir+savename+'.'+ending)