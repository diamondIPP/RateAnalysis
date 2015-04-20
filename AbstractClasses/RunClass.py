from Helper.Initializer import initializer
from Runinfos.RunInfo import RunInfo
from DiamondClass import Diamond
from Elementary import Elementary
import os
import ConfigParser

class Run(Elementary):

    current_run = {}
    operationmode = ''
    TrackingPadAnalysis = {}

    def __init__(self,validate = True, run_number=None, verbose = False):
        Elementary.__init__(self, verbose=verbose)
        self.run_number = -1
        RunInfo.load('Runinfos/runs.json')
        self.runinfo = RunInfo.runs
        if validate:
            self.ValidateRuns()

        if run_number is not None:
            assert(run_number > 0), "incorrect run_number"
            self.run_number = run_number

            self.SetRun(run_number)
        self.IsMonteCarlo = False

    def LoadConfig(self):
        parser = ConfigParser.ConfigParser()
        parser.read('Configuration/Machineconfig.cfg')
        self.operationmode = parser.get('EXEC-MACHINE','operationmode')
        self.ShowAndWait = False

    def ValidateRuns(self, list_of_runs = None):
        if list_of_runs is not None:
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
        if run_number in RunInfo.runs and boolfunc(run_number):
            self.run_number = run_number

            if self.operationmode == "local-ssh":
                self.TrackingPadAnalysis['ROOTFile'] = '/Volumes/scratch/PAD-testbeams/PSI_sept_14/software/TrackingPadAnalysis/results/runs/run_'+str(run_number)+'/track_info.root'

            if self.operationmode == "ssh":
                self.TrackingPadAnalysis['ROOTFile'] = '/scratch/PAD-testbeams/PSI_sept_14/software/TrackingPadAnalysis/results/runs/run_'+str(run_number)+'/track_info.root'

            if self.operationmode == "local":
                self.TrackingPadAnalysis['ROOTFile'] = 'runs/run_'+str(run_number)+'/track_info.root'

            self.run_number = run_number
            self.current_run = RunInfo.runs[run_number].__dict__ # store dict containing all run infos
            #a = self.current_run.__dict__
            self.diamond = Diamond( self.current_run['diamond'])

            self.ResetMC()
            #self.diamond = diamond # diamond is of type Diamond
            #RunInfo.__init__(self,*args)

            #self.diamond.SetName(self.runs[0]['diamond'])
            return True
        else:
            return False
