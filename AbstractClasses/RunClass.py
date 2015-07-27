from Helper.Initializer import initializer
from Runinfos.RunInfo import RunInfo
from DiamondClass import Diamond
from Elementary import Elementary
import ROOT
import os
import ConfigParser

class Run(Elementary):
    '''

    '''

    current_run = {}
    operationmode = ''
    TrackingPadAnalysis = {}

    def __init__(self, run_number=None, validate = False, verbose = False):
        Elementary.__init__(self, verbose=verbose)
        self.run_number = -1
        self.LoadConfig()
#        RunInfo.load(self.runinfofile)
#        self.runinfo = RunInfo.runs # used at all ?!?
        self.RunInfo = {}
        if validate:
            self.ValidateRuns()

        if run_number != None:
            assert(run_number > 0), "incorrect run_number"
            self.run_number = run_number

            self.SetRun(run_number)
        self.IsMonteCarlo = False

    def LoadConfig(self):
        machineConfigParser = ConfigParser.ConfigParser()
        machineConfigParser.read('Configuration/Machineconfig.cfg')
        self.operationmode = machineConfigParser.get('EXEC-MACHINE','operationmode')
        self.ShowAndWait = False
        runConfigParser = ConfigParser.ConfigParser()
        runConfigParser.read('Configuration/RunConfig.cfg')
        self.filename = runConfigParser.get('BASIC', 'filename')
        self.treename = runConfigParser.get('BASIC', 'treename')
        self.signalname = runConfigParser.get('BASIC', 'signalname')
        self.sshrunpath = runConfigParser.get('BASIC', 'sshrunpath')
        self.runinfofile = runConfigParser.get('BASIC', 'runinfofile')


    def ValidateRuns(self, list_of_runs = None):
        if list_of_runs != None:
            runs = list_of_runs
        else:
            runs = RunInfo.runs.keys() # list of all runs
        self.VerbosePrint("Validating runs: ",runs)
        for run in runs:
            self.ValidateRun(run)

    def ValidateRun(self,run_number):
        self.SetRun(run_number)
        if not os.path.exists(self.TrackingPadAnalysis['ROOTFile']):
            del RunInfo.runs[run_number]
            print "INFO: path of run number ",run_number, " not found."
            return False
        else:
            return True

    def ResetMC(self):
        pass

    def SetRun(self, run_number, validate = False):
        if validate:
            boolfunc = self.ValidateRun
        else:
            boolfunc = lambda run: True
        assert(run_number > 0), "incorrect run_number"
        if True or run_number in RunInfo.runs and boolfunc(run_number):
            self.run_number = run_number

            if self.operationmode == "local-ssh":
                fullROOTFilePath = '/Volumes/'+self.sshrunpath+'/'+self.filename+str(run_number).zfill(3)+'.root'
                self.TrackingPadAnalysis['ROOTFile'] = fullROOTFilePath

            if self.operationmode == "ssh":
                fullROOTFilePath = '/'+self.sshrunpath+'/'+self.filename+str(run_number).zfill(3)+'.root'
                self.TrackingPadAnalysis['ROOTFile'] = fullROOTFilePath

            if self.operationmode == "local":
                self.TrackingPadAnalysis['ROOTFile'] = 'runs/run_'+str(run_number)+'/'+self.filename+str(run_number).zfill(3)+'.root'

            self.run_number = run_number
 #           self.current_run = RunInfo.runs[run_number].__dict__ # store dict containing all run infos
 #           self.RunInfo = self.current_run.copy()
            #a = self.current_run.__dict__
#            self.diamond = Diamond( self.current_run['diamond'])

            self.ResetMC()
            #self.diamond = diamond # diamond is of type Diamond
            #RunInfo.__init__(self,*args)

            self.rootfile = ROOT.TFile(fullROOTFilePath)
            print "LOADING: ", fullROOTFilePath
            self.tree = self.rootfile.Get(self.treename) # Get TTree called "track_info"

            #self.diamond.SetName(self.runs[0]['diamond'])
            return True
        else:
            return False
