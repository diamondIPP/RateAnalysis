import ROOT
from RunClass import Run
import os
from ConfigClass import Pad2DHistConfig

class Analysis(object):

    Signal2DDistribution = {
        'Histogram':''

    }

    def __init__(self, run_object, config_object = Pad2DHistConfig(50,-1,1)):
        print('initializing....')
        self.run_object = run_object
        self.config_object = config_object
        self.TrackingPadAnalysisROOTFile = run_object.TrackingPadAnalysis['ROOTFile']
        self.Signal2DDistribution['Histogram'] = ROOT.TH2D("Signal2D",
                                                        "2D Signal distribution",
                                                        config_object.bins_x,
                                                        config_object.min_x,
                                                        config_object.max_x,
                                                        config_object.bins_y,
                                                        config_object.min_y,
                                                        config_object.max_y
                                                        )

        self.signal_canvas = ROOT.TCanvas("signal_canvas", "signal distribution")


        self.signal_sum = ROOT.TH2D("signal_sum",
                                    "Sum of signal distribution",
                                    config_object.bins_x,
                                    config_object.min_x,
                                    config_object.max_x,
                                    config_object.bins_y,
                                    config_object.min_y,
                                    config_object.max_y
                                    )
        self.signal_counts = self.signal_sum.Clone("Signal Counts")

        # loading data file
        assert (os.path.exists(self.TrackingPadAnalysisROOTFile)), 'cannot find '+filename
        self.rootfile = ROOT.TFile(self.TrackingPadAnalysisROOTFile)


    def DoAnalysis(self,minimum_bincontent = 1):
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer"
        minimum_statistics = minimum_bincontent # bins with less hits are ignored

        self.track_info = self.rootfile.Get('track_info') # Get TTree called "track_info"

        bins = self.config_object.bins_x
        xmin = self.config_object.min_x
        xmax = self.config_object.max_x
        ymin = self.config_object.min_y
        ymax = self.config_object.max_y

        # fill two 2-dim histograms to collect the hits and signal strength
        for i in xrange(self.track_info.GetEntries()):

            self.track_info.GetEntry(i)
            self.signal_sum.Fill(self.track_info.track_x, self.track_info.track_y, self.track_info.integral50)
            self.signal_counts.Fill(self.track_info.track_x, self.track_info.track_y, 1)

        # go through every bin, calculate the average signal strength and fill the main 2D hist
        binwidth_x = 1.*(xmax-xmin)/bins
        binwidth_y = 1.*(ymax-ymin)/bins
        current_pos_x = xmin + 1.*binwidth_x/2.
        for bin_x in xrange(1,bins+1):

            current_pos_y = ymin + 1.*binwidth_y/2.

            for bin_y in xrange(1,bins+1):

                binsignalsum = abs(self.signal_sum.GetBinContent(bin_x, bin_y))
                binsignalcount = self.signal_counts.GetBinContent(bin_x, bin_y)

                if binsignalcount >= minimum_statistics :
                    self.Signal2DDistribution['Histogram'].Fill(current_pos_x, current_pos_y, abs(binsignalsum/binsignalcount))

                current_pos_y += binwidth_y

            current_pos_x += binwidth_x

    def CreatePlots(self,saveplots = False,savename = '2DSignalDistribution',ending='pdf',saveDir = 'Results/'):

        # Plot the Signal2D TH2D histogram
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)
        self.Signal2DDistribution['Histogram'].SetStats(False)
        #signal_canvas.cd()
        self.Signal2DDistribution['Histogram'].Draw('colz')

        if saveplots:
            # Results directories:
            resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
            if not os.path.exists(resultsdir):
                os.makedirs(resultsdir)

            ROOT.gPad.Print(resultsdir+savename+'.'+ending)