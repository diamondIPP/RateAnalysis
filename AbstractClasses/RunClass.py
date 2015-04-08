from Helper.Initializer import initializer
from Runinfos.RunInfo import RunInfo
from DiamondClass import Diamond
import os
import ConfigParser

class Run(object):

    current_run = {}
    operationmode = ''
    TrackingPadAnalysis = {}

    def __init__(self,validate = True, run_number=None):
        self.run_number = -1
        self.runinfo = RunInfo.load('Runinfos/runs.json')
        if validate:
            self.ValidateRuns()

        if run_number is not None:
            assert(run_number > 0), "incorrect run_number"
            self.run_number = run_number

            self.SetRun(run_number)

    def ValidateRuns(self, list_of_runs = None):
        runs = RunInfo.runs.keys() # list of all runs
        if list_of_runs is not None:
            runs = list_of_runs
        print "Validating runs: ",runs
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

    def SetRun(self,run_number):

        assert(run_number > 0), "incorrect run_number"
        if run_number in RunInfo.runs:
            self.run_number = run_number

            parser = ConfigParser.ConfigParser()
            parser.read('Configuration/Machineconfig.cfg')
            operationmode = parser.get('EXEC-MACHINE','operationmode')

            if operationmode == "local-ssh":
                self.TrackingPadAnalysis['ROOTFile'] = '/Volumes/scratch/PAD-testbeams/PSI_sept_14/software/TrackingPadAnalysis/results/runs/run_'+str(run_number)+'/track_info.root'

            if operationmode == "ssh":
                self.TrackingPadAnalysis['ROOTFile'] = '/scratch/PAD-testbeams/PSI_sept_14/software/TrackingPadAnalysis/results/runs/run_'+str(run_number)+'/track_info.root'

            if operationmode == "local":
                self.TrackingPadAnalysis['ROOTFile'] = 'runs/run_'+str(run_number)+'/track_info.root'

            self.run_number = run_number
            self.current_run = RunInfo.runs[run_number].__dict__ # store dict containing all run infos
            #a = self.current_run.__dict__
            self.diamond = Diamond( self.current_run['diamond'])



            #self.diamond = diamond # diamond is of type Diamond
            #RunInfo.__init__(self,*args)

            #self.diamond.SetName(self.runs[0]['diamond'])
            return True
        else:
            return False
