from Helper.Initializer import initializer
#from RunInfo import RunInfo
from DiamondClass import Diamond
#from Helper.Initializer import initializer

class Run(object):

    #TrackingPadAnalysis.ROOTFile = '' # filepath
    TrackingPadAnalysis = {
        'ROOTFile': ''

    }
    #@initializer
    def __init__(self,run_number):
        self.run_number = run_number
        self.TrackingPadAnalysis['ROOTFile'] = '/scratch/PAD-testbeams/PSI_sept_14/software/TrackingPadAnalysis/results/runs/run_'+str(run_number)+'/track_info.root'
        #self.diamond = diamond # diamond is of type Diamond
        #RunInfo.__init__(self,*args)

        #self.diamond.SetName(self.runs[0]['diamond'])
