from Helper.Initializer import initializer
from Runinfos.RunInfo import RunInfo
from DiamondClass import Diamond

class Run(object):

    run_number = -1
    current_run = {}
    operationmode = ''
    #TrackingPadAnalysis.ROOTFile = '' # filepath
    TrackingPadAnalysis = {}

    def __init__(self,run_number=None, operationmode="local-ssh"):
        self.operationmode = operationmode
        self.runinfo = RunInfo.load('Runinfos/runs.json')
        if run_number is not None:
            assert(run_number > 0), "incorrect run_number"
            self.run_number = run_number

            assert(operationmode == "local" or operationmode == "ssh"), "operation mode has to be 'local', 'local-ssh' or 'ssh' "
            if operationmode == "local-ssh":
                self.TrackingPadAnalysis['ROOTFile'] = '/Volumes/scratch/PAD-testbeams/PSI_sept_14/software/TrackingPadAnalysis/results/runs/run_'+str(run_number)+'/track_info.root'

            if operationmode == "ssh":
                self.TrackingPadAnalysis['ROOTFile'] = '/scratch/PAD-testbeams/PSI_sept_14/software/TrackingPadAnalysis/results/runs/run_'+str(run_number)+'/track_info.root'

            if operationmode == "local":
                self.TrackingPadAnalysis['ROOTFile'] = 'runs/run_'+str(run_number)+'/track_info.root'

            self.SetRun(run_number)


    def SetRun(self,run_number):

        assert(run_number > 0), "incorrect run_number"
        self.run_number = run_number
        operationmode = self.operationmode
        assert(operationmode == "local" or operationmode == "ssh" or operationmode == "local-ssh"), "operation mode has to be 'local', 'local-ssh' or 'ssh'.. op.mode.is: "+operationmode
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
