import ROOT
import numpy as np
import os
from Bin import Bin

class BinCollection(object):

    def __init__(self, binsx, xmin, xmax, binsy, ymin, ymax):

        self.ListOfBins = [Bin(i,self) for i in xrange((binsx+2)*(binsy+2))]
        self.binnumbers = [i for i in xrange((binsx+2)*(binsy+2))]
        self.Attributes = {
            'binsx': binsx, # bins in x without frame of 1 bin
            'binsy': binsy, # bins in y without frame of 1 bin
            'XMIN': xmin,
            'XMAX': xmax,
            'YMIN': ymin,
            'YMAX': ymax,
            'binwidth_x': 1.*(xmax-xmin)/binsx,
            'binwidth_y': 1.*(ymax-ymin)/binsy
        }

        self.counthisto = ROOT.TH2D('counthisto',
                                    '2D hit distribution',
                                    self.Attributes['binsx'],
                                    self.Attributes['XMIN'],
                                    self.Attributes['XMAX'],
                                    self.Attributes['binsy'],
                                    self.Attributes['YMIN'],
                                    self.Attributes['YMAX']
                                    )
        self.totalsignal = ROOT.TH2D('totalsignal',
                                    '2D total signal distribution',
                                    self.Attributes['binsx'],
                                    self.Attributes['XMIN'],
                                    self.Attributes['XMAX'],
                                    self.Attributes['binsy'],
                                    self.Attributes['YMIN'],
                                    self.Attributes['YMAX']
                                    )

    def Fill(self,x,y,signal):
        self.counthisto.Fill(x,y)
        self.totalsignal.Fill(x,y,signal)
        self.ListOfBins[self.GetBinNumber(x,y)].AddData(signal)

    def ShowBinXYSignalHisto(self,x,y,saveplot = False):
        self.ListOfBins[self.GetBinNumber(x,y)].CreateBinSignalHisto(saveplot)

    def CalculateMeanSignalDistribution(self,minimum_bincontent = 1):
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer"
        self.meansignaldistribution = ROOT.TH2D('meansignaldistribution',
                                                "Mean Signal Distribution",
                                                self.Attributes['binsx'],
                                                self.Attributes['XMIN'],
                                                self.Attributes['XMAX'],
                                                self.Attributes['binsy'],
                                                self.Attributes['YMIN'],
                                                self.Attributes['YMAX']
                                                )
        # go through every bin, calculate the average signal strength and fill the main 2D hist
        binwidth_x = self.Attributes['binwidth_x']
        binwidth_y = self.Attributes['binwidth_y']
        x_ = self.Attributes['XMIN'] + 1.*binwidth_x/2.
        for bin_x in xrange(1,self.Attributes['binsx']+1):

            y_ = self.Attributes['YMIN'] + 1.*binwidth_y/2.

            for bin_y in xrange(1,self.Attributes['binsy']+1):

                binsignalsum = abs(self.totalsignal.GetBinContent(bin_x, bin_y))
                binsignalcount = self.counthisto.GetBinContent(bin_x, bin_y)

                if binsignalcount >= minimum_bincontent :
                    self.meansignaldistribution.Fill(x_, y_, abs(binsignalsum/binsignalcount))

                y_ += binwidth_y

            x_ += binwidth_x

    # select Bins in a rectangular region and return list of bins
    def SelectRectangularBins(self, xlow, xhigh, ylow, yhigh, activate = True):
        list_of_bins = []
        lower_left_bin = self.GetBinNumber(xlow,ylow)
        lower_right_bin = self.GetBinNumber(xhigh,ylow)
        upper_left_bin = self.GetBinNumber(xlow,yhigh)
        upper_right_bin = self.GetBinNumber(xhigh,yhigh)
        totalbinsx = self.Attributes['binsx']+2 # binsx plus two frame bins

        while lower_left_bin <= upper_left_bin:
            list_of_bins += [i for i in xrange(lower_left_bin,lower_right_bin+1)]
            lower_left_bin += totalbinsx
            lower_right_bin += totalbinsx
        assert(upper_right_bin == lower_right_bin-totalbinsx), "Bin Mismatch in SelectRectangularBins.\n\tupper right bin: "+str(upper_right_bin)+"\n\tlower right bin: "+str(lower_right_bin)

        # selection
        if activate:
            for binnr in list_of_bins:
                self.ListOfBins[binnr].selected = True

        return list_of_bins

    def UnselectAllBins(self):
        for bin in self.ListOfBins:
            bin.selected = False

    # select bins within a mean signalstrength around the signal of a reference bin
    def SelectSignalStrengthRegion(self,
                                   refBin,
                                   sensitivity = 0.1,
                                   activate = True,
                                   xlow = None,
                                   xhigh = None,
                                   ylow = None,
                                   yhigh = None):
        selected_bins = []
        if yhigh == None:
            list_of_bins = self.binnumbers
        else:
            list_of_bins = self.SelectRectangularBins(xlow,xhigh,ylow,yhigh,False)

        binnumber = refBin.Attributes['binnumber']
        assert(binnumber in list_of_bins), "Bin given is not in selected region."

        bin_avg_signal = self.meansignaldistribution.GetBinContent(binnumber)
        signal_lowerbound = bin_avg_signal*(1-sensitivity)
        signal_upperbound = bin_avg_signal*(1+sensitivity)

        for binnumber in list_of_bins:
            signal = self.meansignaldistribution.GetBinContent(binnumber)
            if signal_lowerbound <= signal <= signal_upperbound:
                selected_bins.append(binnumber)
                if activate:
                    self.ListOfBins[binnumber].selected = True

        return selected_bins

    # select a single bin with bin number binnumber
    def SelectBin(self,binnumber):
        self.ListOfBins[binnumber].selected = True

    # draw a 2d distribution which shows the selected bins
    def ShowSelectedBins(self,draw = True):
        if draw:
            ROOT.gStyle.SetPalette(51)
            ROOT.gStyle.SetNumberContours(2)
            selection_canvas = ROOT.TCanvas('selection_canvas', 'Selected Bins', 500, 500)
        binsx = self.Attributes['binsx']
        binsy = self.Attributes['binsy']
        xmin = self.Attributes['XMIN']
        xmax = self.Attributes['XMAX']
        ymin = self.Attributes['YMIN']
        ymax = self.Attributes['YMAX']

        selection_pad = ROOT.TH2D('selection_pad', "Selected Bins", binsx, xmin, xmax, binsy, ymin, ymax)
        i = 0
        for bin in self.ListOfBins:
            if bin.selected:
                x_, y_ = bin.GetBinCenter()
                selection_pad.Fill(x_, y_)
                i += 1
        selection_pad.SetTitle(str(i)+" Bins selected")
        selection_pad.SetStats(False)
        selection_pad.GetXaxis().SetTitle('pos x / cm')
        selection_pad.GetYaxis().SetTitle('pos y / cm')
        selection_pad.GetYaxis().SetTitleOffset(1.4)
        if draw:
            selection_canvas.cd()
            selection_pad.Draw("col")
            raw_input("Selected bins shown")
            ROOT.gStyle.SetPalette(53)
            ROOT.gStyle.SetNumberContours(999)
        return selection_pad


    def GetListOfSelectedBins(self):
        selected_bins = []
        for bin in self.ListOfBins:
            if bin.selected:
                selected_bins.append(bin.Attributes['binnumber'])
        return selected_bins

    # show distribution of K_i from SIGMA = K_i * sigma_i / sqrt(n) for selected bins
    def ShowKDistribution(self,draw = True):
        if draw:
            Kcanvas = ROOT.TCanvas('Kcanvas','K Canvas')

        selection = self.GetListOfSelectedBins()
        binmeans = []
        n = []
        sigma = []
        for bin_nr in selection:
            binmeans.append(self.meansignaldistribution.GetBinContent(bin_nr))
            sigma.append(self.ListOfBins[bin_nr].GetSigma())
            n.append(self.ListOfBins[bin_nr].GetEntries())
        #N = len(selection)
        SIGMA = np.std(binmeans)
        sigma = np.array(sigma)
        K = SIGMA * np.sqrt(n) / sigma
        Khisto = ROOT.TH1D('Khisto', 'K Distribution', 50, 0, int(K.max())+1)
        Khisto.GetXaxis().SetTitle('K value')
        for i in xrange(len(K)):
            Khisto.Fill(K[i])
        if draw:
            Kcanvas.cd()
            Khisto.Draw()
            raw_input("K dostribution shown..")
        return Khisto

    def ShowCombinedKDistribution(self, saveplots = False, savename = 'CombinedKDistribution', ending='png', saveDir = 'Results/'):
        ROOT.gStyle.SetPalette(51)
        ROOT.gStyle.SetNumberContours(2)

        canvas = ROOT.TCanvas('canvas', 'combined', 1000,500)
        selection = self.ShowSelectedBins(False)
        Khisto = self.ShowKDistribution(False)
        canvas.Divide(2,1)
        canvas.cd(1)
        selection.Draw('col')
        canvas.cd(2)
        Khisto.Draw()
        if saveplots:
            self.SavePlots(savename, ending, saveDir)
        raw_input('Combined K Distribution Drawn')
        ROOT.gStyle.SetPalette(53)
        ROOT.gStyle.SetNumberContours(999)

    def GenerateTotalMeanDistribution(self):
        pass

    def GetTotalCountDistribution(self):
        return self.totalsignal

    def GetMeanSignalDistribution(self, minimum_bincontent = 1):
        self.CalculateMeanSignalDistribution(minimum_bincontent)
        return self.meansignaldistribution

    # def AddBin(self):
    #     pass

    def GetBinNumber(self,x,y):
        binnumber = self.counthisto.FindBin(x,y)
        return binnumber

    def GetBinByNumber(self, bin_nunmber):
        return self.ListOfBins[bin_nunmber]

    def GetBinByCoordinates(self, x, y):
        nr = self.GetBinNumber(x,y)
        return self.ListOfBins[nr]

    def SavePlots(self, savename, ending, saveDir):
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        resultsdir = saveDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)

        ROOT.gPad.Print(resultsdir+savename+'.'+ending)

