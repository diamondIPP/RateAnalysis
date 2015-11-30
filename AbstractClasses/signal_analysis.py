__author__ = 'micha'

# ==============================================
# IMPORTS
# ==============================================
import ROOT
from AbstractClasses.Elementary import Elementary
from ROOT import gROOT, TGraphErrors, TCanvas, TH2D, gStyle, TF1
from newAnalysis import Analysis


# ==============================================
# MAIN CLASS
# ==============================================
class SignalAnalysis(Elementary):

    def __init__(self, analysis, channel, canvas=None, binning=5000):
        Elementary.__init__(self)

        # main
        self.analysis = analysis
        self.analysis = Analysis(5)
        self.channel = channel

        # stuff
        self.cut = self.analysis.GetNEventsCut(channel)