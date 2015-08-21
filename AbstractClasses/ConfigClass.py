#from Helper.Initializer import initializer
from math import ceil
from AbstractClasses.Elementary import Elementary
import ConfigParser

class Pad2DHistConfig(object):
    #@initializer
    def __init__(self,bins_x,bins_y=None):

        self.bins_x = bins_x
        self.binsxy = bins_x
        if bins_y == None:
            self.bins_y = bins_x
        else:
            self.bins_y = bins_y

class MeanSignalHistConfig(object):
    #@initializer
    def __init__(self,bins,min_x,max_x):
        pass

class BinCollectionConfig(Elementary):
    '''
    A config object for analysis
    e.g. binning size
    '''
    def __init__(self, binningsize=None,  **kwargs):
        '''

        :param binningsize: size of bins in microns
        :param kwargs:
        :return:
        '''

        self.config = {
            '2DHist': {
                'binsize': 100/10000., # in cm
                'binsx': 0,
                'xmin': 0.,
                'xmax': 0.,
                'binsy': 0,
                'ymin': 0.,
                'ymax': 0.
            },
            'Hist1': ''
        }

        self.LoadConfigFile()

        if binningsize != None:
            self.SetBinning(binningsize)



    def LoadConfigFile(self):
        configfile = "Configuration/AnalysisConfig_"+self.TESTCAMPAIGN+".cfg"
        print "configclass configfile: ", configfile
        parser = ConfigParser.ConfigParser()
        parser.read(configfile)
        windowXmin = parser.getfloat("TRACKING", "windowXmin")
        windowXmax = parser.getfloat("TRACKING", "windowXmax")
        windowYmin = parser.getfloat("TRACKING", "windowYmin")
        windowYmax = parser.getfloat("TRACKING", "windowYmax")
        self.binning = parser.getint("TRACKING", "padBinning")
        self.SetBinning(self.binning)
        self.SetWindow(windowXmin, windowXmax, windowYmin, windowYmax) # default window


    def Get2DAttributes(self):
        '''

        :return:
        '''
        binsx = self.config['2DHist']['binsx']
        xmin = self.config['2DHist']['xmin']
        xmax = self.config['2DHist']['xmax']
        binsy = self.config['2DHist']['binsy']
        ymin = self.config['2DHist']['ymin']
        ymax = self.config['2DHist']['ymax']

        return binsx, xmin, xmax, binsy, ymin, ymax

    def SetWindow(self, xlow, xhigh, ylow, yhigh):
        '''
        calculates the xmin, xmax, ymin, ymax for the binsize used
        :param xlow:
        :param xhigh:
        :param ylow:
        :param yhigh:
        :return:
        '''
        binsize = self.config['2DHist']['binsize']
        self.config['2DHist']['xmin'] = xlow
        self.config['2DHist']['ymin'] = ylow
        minwidth = xhigh - xlow
        minheight = yhigh - ylow
        binsx = ceil(1.*minwidth/binsize)
        binsy = ceil(1.*minheight/binsize)
        self.config['2DHist']['binsx'] = int(binsx)
        self.config['2DHist']['binsy'] = int(binsy)
        self.config['2DHist']['xmax'] = xlow + binsx * binsize
        self.config['2DHist']['ymax'] = ylow + binsy * binsize

    def SetBinning(self, binning):
        self.config['2DHist']['binsize'] = binning/10000.
        self.binning = binning

    def SetWindowFromDiamond(self, diamond):
        '''
        Uses the diamond configureation
        :param diamond:
        :return:
        '''
        xmin = diamond.Position['xmin']
        xmax = diamond.Position['xmax']
        ymin = diamond.Position['ymin']
        ymax = diamond.Position['ymax']
        self.SetWindow(xmin,xmax,ymin,ymax)
