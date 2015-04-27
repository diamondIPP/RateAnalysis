import ROOT
from ROOT import gPad
import SignalInRCClass
import numpy as np
import os
import types as t
from AbstractClasses.FindExtrema import FindMaxima, FindMinima
from Elementary import Elementary
from Bin import Bin
from Langaus import *
from array import array
import copy


class BinCollection(Elementary):
    '''
    A BinCollection Object Contains Bin-objects. It is used to store the data, as well
    as to read the data, make selections or to make Plots related to collection of many bins.

    '''

    def __init__(self, parent_analysis_obj, binsx, xmin, xmax, binsy, ymin, ymax, verbose=False):
        '''
        Constructor of a Bincollection. Since the data collection is based on ROOT.TH2D,
        the bins are ordered in a rectangular pattern inside a frame which is 1 bin thick leading
        to a total number of bins of (binsx+2)*(binsy+2)
        :param binsx: Number of bins in x direction of data collection window
        :param xmin: data collection window lower x bound
        :param xmax: data collection window upper x bound
        :param binsy: Number of bins in y direction of data collection window
        :param ymin: data collection window lower y bound
        :param ymax: data collection window upper y bound
        :return: -
        '''
        Elementary.__init__(self, verbose = verbose)
        if type(binsx) is not t.IntType or type(binsy) is not t.IntType:
            "INFO: binsx or binsy not of int type. Changing it to int..."
            binsx = int(binsx)
            binsy = int(binsy)
        self.ListOfBins = [Bin(i,self) for i in xrange((binsx+2)*(binsy+2))] # A list, containing all Bin objects
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
        self.parent_analysis_obj = parent_analysis_obj
        self.run_object = parent_analysis_obj.run_object
        self.counthisto = ROOT.TH2D('counthisto', '2D hit distribution', *self.Get2DAttributes())
        self.totalsignal = ROOT.TH2D('totalsignal', '2D total signal distribution', *self.Get2DAttributes())
        self.SignalHisto = ROOT.TH1D('SignalHisto,', 'Signal response Histogram', 500, 0, 500)

    def __del__(self):
        for bin in self.ListOfBins:
            bin.__del__()
            del bin
        ROOT.gROOT.Delete('counthisto')
        ROOT.gROOT.Delete('totalsignal')
        ROOT.gROOT.Delete('SignalHisto')
        if hasattr(self, "meansignaldistribution"):
            ROOT.gROOT.Delete('meansignaldistribution')
        if hasattr(self, "MinimaSearch"):
            self.MinimaSearch.__del__()
            del self.MinimaSearch

    def LoadConfig(self):
        self.ShowAndWait = False

    def Get2DAttributes(self):
        '''
        Returns attributes to initialize a ROOT TH2D object
        :return: binsx, xmin, xmax, binsy, ymin, ymax
        '''
        binsx = self.Attributes['binsx']
        xmin = self.Attributes['XMIN']
        xmax = self.Attributes['XMAX']
        binsy = self.Attributes['binsy']
        ymin = self.Attributes['YMIN']
        ymax = self.Attributes['YMAX']
        return binsx, xmin, xmax, binsy, ymin, ymax

    def MakeFits(self,minimum_entries = 5):
        assert(minimum_entries >= 1), "number of entries has to be greater or equal to 1"
        for i in xrange(len(self.ListOfBins)):
            if self.ListOfBins[i].GetEntries() >= minimum_entries:
                self.ListOfBins[i].FitLandau()

    # !! cannot be inherented to non rectangular
    def Fill(self,x,y,signal):
        '''
        Adds the datapoint into the corresponding bin inside the bincollection as
        well as into the two histograms counthisto and totalsignal inside this bin collection
        :param x:
        :param y:
        :param signal:
        :return:
        '''
        self.counthisto.Fill(x,y)
        self.totalsignal.Fill(x,y,signal)
        self.ListOfBins[self.GetBinNumber(x,y)].AddData(signal)
        self.SignalHisto.Fill(signal)

    def ShowBinXYSignalHisto(self,x,y,saveplot = False, show_fit=False):
        '''
        Shows a Histogram of the Signal response distribution inside the bin which
        contains the coordinates (x,y)
        :param x: coordinate x in cm which is contained in the bin of interest
        :param y: coordinate y in cm which is contained in the bin of interest
        :param saveplot: if True save plot as as Results/Bin_X0.123Y-0.123_SignalHisto.png
        :return: -
        '''
        self.ListOfBins[self.GetBinNumber(x,y)].CreateBinSignalHisto(saveplot,show_fit)

    def CalculateMeanSignalDistribution(self, minimum_bincontent = 1):
        '''

        :param minimum_bincontent:
        :return:
        '''
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer"
        self.meansignaldistribution = ROOT.TH2D('meansignaldistribution'+str(self.GLOBAL_COUNT), "Mean Signal Distribution", *self.Get2DAttributes())
        self.GLOBAL_COUNT += 1
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
        print len(list_of_bins), " bins selected (Rectangualr region)"
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
                                   yhigh = None,
                                   minimum_bincontent = 5):
        '''
        Creates and returns a list of all binnumbers in a region with bins that
        have a similar mean signal response as the mean signal response of a
        reference bin inside this region. If activate = True, the bins get selected
        (bin.selection = True). If no region is passed, all bins are considered.
        :param refBin: sets the default value of mean response
        :param sensitivity: bins are picked inside
            refSignal*(1-sensitivity) <= signal <= refSignal*(1+sensitivity)
        :param activate: if True the bins get set to bin.selected = True
        :param xlow: Window to restrict the considered bins
        :param xhigh:
        :param ylow:
        :param yhigh:
        :return:
        '''
        selected_bins = []
        if yhigh == None:
            list_of_bins = self.binnumbers
        else:
            list_of_bins = self.SelectRectangularBins(xlow,xhigh,ylow,yhigh,False)

        binnumber = refBin.GetBinNumber()
        assert(binnumber in list_of_bins), "Bin given is not in selected region."

        bin_avg_signal = self.meansignaldistribution.GetBinContent(binnumber)
        signal_lowerbound = bin_avg_signal*(1-sensitivity)
        signal_upperbound = bin_avg_signal*(1+sensitivity)

        for binnumber in list_of_bins:
            signal = self.meansignaldistribution.GetBinContent(binnumber)
            if self.ListOfBins[binnumber].GetEntries() > minimum_bincontent:
                if signal_lowerbound <= signal <= signal_upperbound:
                    selected_bins.append(binnumber)
                    if activate:
                        self.ListOfBins[binnumber].selected = True
        print len(selected_bins), " bins selected"
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
            self.IfWait("Selected bins shown")
            ROOT.gStyle.SetPalette(53)
            ROOT.gStyle.SetNumberContours(999)
        return selection_pad

    def GetSortedListOfBins(self, attribute='average', ascending = True):
        '''
        Returns list of bins (binnunmbers) in an order with respect to "attribute"
        :return: ordered_list
        '''
        # self.UpdateBinAttributes()
        # SortedListOfBins = sorted(self.ListOfBins, key = lambda bin: bin.Attributes[attribute], reverse = not ascending)
        # ordered_list = [SortedListOfBins[i].Attributes['binnumber'] for i in xrange(len(SortedListOfBins))]
        sortdata = np.ones((3,len(self.ListOfBins))) # sortdata[0,:] numbers, [1,:] means, [3,:] hits
        count = 0
        for i in xrange(len(self.ListOfBins)):
            self.ListOfBins[i].UpdateAttributes()
            if self.ListOfBins[i].Attributes['entries'] >= 5:
                sortdata[0,i] = self.ListOfBins[i].Attributes['binnumber']
                sortdata[1,i] = self.ListOfBins[i].Attributes['average']
                sortdata[2,i] = self.ListOfBins[i].Attributes['entries']
                count += 1
                #print "nr: ",sortdata[0,i]," av.: ", sortdata[1,i]," ent.: ", sortdata[2,i]
        #print "*************************************************************"
        data = list(-sortdata[1][:]) #select data to sort ([1]->by average)
        arg_sorted = np.argsort(data)
        sorted_data = sortdata[:,arg_sorted] # ?! WTF? why does sortdata[:][arg_sorted] not work??!?!
        # for i in xrange(len(sorted_data[0,:count])):
        #     print "nr: ",sorted_data[0,i]," av.: ", sorted_data[1,i]," ent.: ", sorted_data[2,i]
        means = list(sorted_data[1,:])
        ordered_list = list(sorted_data[0,:count])
        #print "ordered list:", ordered_list
        # print "means: ", means
        # print "entries : ", sorted_data[2,:]
        # print "len of sorted_data: ", len(sorted_data)
        return map(int,ordered_list)

    def GetMaximumSignalResponseBinNumber(self):
        return self.GetSortedListOfBins(ascending=False)[0]

    def GetListOfSelectedBins(self):
        '''
        Returns a List of bin numbers, which correspond to selected bins (bin.selected = True)
        :return: selected_bins (bin numbers)
        '''
        selected_bins = []
        for bin in self.ListOfBins:
            if bin.selected:
                selected_bins.append(bin.GetBinNumber())
        return selected_bins

    # show distribution of K_i from SIGMA = K_i * sigma_i / sqrt(n) for selected bins
    def ShowKDistribution(self,draw = True):
        selection = self.GetListOfSelectedBins()
        assert(len(selection) is not 0), "no Bins selected!"
        if draw:
            Kcanvas = ROOT.TCanvas('Kcanvas','K Canvas')

        binmeans = []
        n = []
        sigma = []
        for bin_nr in selection:
            binmeans.append(self.ListOfBins[bin_nr].GetMean())
            sigma.append(self.ListOfBins[bin_nr].GetSigma())
            n.append(self.ListOfBins[bin_nr].GetEntries())
        N = len(selection)
        SIGMA = np.std(binmeans)
        # print "sigmas : ",sorted(sigma)
        # print "means : ",sorted(binmeans)
        # print "n : ",sorted(n)


        sigma = np.array(sigma)
        K = SIGMA * np.sqrt(N) / sigma
        Khisto = ROOT.TH1D('Khisto', 'K Distribution', 50, 0, int(K.max())+1)
        Khisto.GetXaxis().SetTitle('K value')
        for i in xrange(len(K)):
            Khisto.Fill(K[i])
        if draw:
            Kcanvas.cd()
            Khisto.Draw()
            self.IfWait("K distribution shown..")
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
        self.IfWait('Combined K Distribution Drawn')
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

    def GetBinByNumber(self, bin_number):
        '''
        Returns the bin object with number "bin_number"
        :param bin_nunmber: the bin number of the bin to return
        :return: bin with bin number bin_number
        '''
        if not(type(bin_number) is t.IntType):
            binnumber = int(bin_number)
        else:
            binnumber = bin_number
        return self.ListOfBins[binnumber]

    def GetBinByCoordinates(self, x, y):
        nr = self.GetBinNumber(x,y)
        return self.ListOfBins[nr]

    def GetBinCenter(self, bin_numbers):
        '''
        returns the coordinates of the center of a bin with bin number 'binnumbers' or the coordinates
        of a list of bin numbers
        :param bin_numbers:
        :return:
        '''
        if type(bin_numbers) is t.ListType:
            nr = len(bin_numbers)
            coordinates = []
            for i in xrange(nr):
                coordinates.append(self.ListOfBins[bin_numbers[i]].GetBinCenter())
            return coordinates
        elif type(bin_numbers) is t.IntType or type(bin_numbers) is t.FloatType:
            bin_numbers = int(bin_numbers)
            coordinates = self.ListOfBins[bin_numbers].GetBinCenter()
            return coordinates
        else:
            assert(False), "type of bin_number not accepted. bin_number as to be a list or int or of type float"

    def GetBinsInColumn(self, position):
        '''
        Returns a list, containing all bin objects of a column at
        x position 'position'
        :param heigth:
        :return:
        '''
        start_y = self.Attributes['YMIN']+ self.Attributes['binwidth_y']/2.
        start_bin = self.GetBinByCoordinates(position, start_y)
        start_binnumber = start_bin.GetBinNumber()

        end_y = self.Attributes['YMAX']- self.Attributes['binwidth_y']/2.
        end_bin = self.GetBinByCoordinates(position, end_y)
        end_binnumber = end_bin.GetBinNumber()

        list_of_bins = []
        for i in xrange(self.Attributes['binsy']):
            list_of_bins.append( self.ListOfBins[start_binnumber + i*(self.Attributes['binsx']+2)] )
        assert (end_binnumber == start_binnumber + (self.Attributes['binsy']-1)*(self.Attributes['binsx']+2)), "Bin Mismatch in GetBinsInColumn"
        return list_of_bins

    def GetBinsInRow(self, height):
        '''
        Returns a List of bin numbers in the Row at height 'height'.
        The height doesen't have to be the height of the bin centers.
        :param heigth:
        :return:
        '''
        start_x = self.Attributes['XMIN']+ self.Attributes['binwidth_x']/2.
        start_bin = self.GetBinByCoordinates(start_x, height)
        start_binnumber = start_bin.GetBinNumber()
        list_of_bins = [self.ListOfBins[i] for i in xrange(start_binnumber, start_binnumber+self.Attributes['binsx'])]
        return list_of_bins

    def GetSignalInColumn(self, position , show = False, show_hits = True):
        '''

        :param position:
        :return:
        '''
        list_of_bins = self.GetBinsInColumn(position)
        columnposition, _ = list_of_bins[0].GetBinCenter()
        signals = []
        attributes = self.Attributes['binsy'], self.Attributes['YMIN'], self.Attributes['YMAX']
        graph = SignalInRCClass.SignalInColumnGraph(position, *attributes)
        count = 0
        for i in xrange(len(list_of_bins)):
            _ , y_ = list_of_bins[i].GetBinCenter()
            signal_ = list_of_bins[i].GetMean()
            sigma_ = list_of_bins[i].GetSigma()
            counts_ = list_of_bins[i].GetEntries()
            graph.SetHistPoint(y_, counts_)
            if counts_ > 0:
                graph.SetGraphPoint(count, y_, signal_)
                graph.SetGraphPointError(count, 0, sigma_/np.sqrt(counts_))
                signals.append(signal_)
                count += 1
            else:
                signals.append(0)
        if show:
            if show_hits:
                graph.DrawBoth()
            else:
                graph.DrawGraph()
            self.SavePlots('signals_in_column{:.3f}'.format(columnposition),'png')
            self.IfWait('show signal in column {:.3f}..'.format(columnposition))
        return signals

    def GetSignalInRow(self, height, show = False, show_hits = True):
        '''

        :param height:
        :return:
        '''
        list_of_bins = self.GetBinsInRow(height)
        _, rowheight = list_of_bins[0].GetBinCenter()
        signals = []
        attributes = self.Attributes['binsx'], self.Attributes['XMIN'], self.Attributes['XMAX']
        graph = SignalInRCClass.SignalInRowGraph(height, *attributes)
        count = 0
        for i in xrange(len(list_of_bins)):
            x_ , _ = list_of_bins[i].GetBinCenter()
            signal_ = list_of_bins[i].GetMean()
            sigma_ = list_of_bins[i].GetSigma()
            counts_ = list_of_bins[i].GetEntries()
            graph.SetHistPoint(x_, counts_)
            if counts_ > 0:
                graph.SetGraphPoint(count, x_, signal_)
                graph.SetGraphPointError(count, 0, sigma_/np.sqrt(counts_))
                signals.append(signal_)
                count += 1
            else:
                signals.append(0)
        if show:
            if show_hits:
                graph.DrawBoth()
            else:
                graph.DrawGraph()
            self.SavePlots('signals_in_row{:.3f}'.format(rowheight),'png')
            self.IfWait('show signal in row {:.3f}..'.format(rowheight))
        return signals

    def GetMPVInColumn(self, position, show=False):
        '''

        :param height:
        :return:
        '''
        list_of_bins = self.GetBinsInColumn(position)
        columnposition, _ = list_of_bins[0].GetBinCenter()
        signals = []
        graph = ROOT.TGraphErrors(self.Attributes['binsy'])
        #histo = ROOT.TH1D('histo', 'bins in row at height {:.3f}'.format(rowheight), *attributes)

        count = 0
        for i in xrange(len(list_of_bins)):
            _ , y_ = list_of_bins[i].GetBinCenter()
            MPV_ = list_of_bins[i].Fit['MPV']
            MPVerr_ = list_of_bins[i].Fit['MPVErr']
            if MPV_ is not None:
                graph.SetPoint(count, y_ , MPV_)
                graph.SetPointError(count, 0, MPVerr_)
                signals.append(MPV_)
                count += 1
            else:
                signals.append(0)
        if show:
            canvas = ROOT.TCanvas('signal_in_column', 'signals in column at position {:.3f}'.format(columnposition))
            canvas.cd()
            graph.SetTitle('Most Probable Signal response in bins in column at position {:.3f}'.format(columnposition))
            graph.GetXaxis().SetTitle("y position / cm")
            graph.GetYaxis().SetTitle("Signal")
            graph.Draw('AP')
            self.SavePlots('signals_in_column{:.3f}'.format(columnposition),'pdf')
            self.IfWait('show signal in column {:.3f}..'.format(columnposition))
        return signals

    def GetMPVInRow(self, height, show=False):
        '''

        :param height:
        :return:
        '''
        list_of_bins = self.GetBinsInRow(height)
        _, rowheight = list_of_bins[0].GetBinCenter()
        signals = []
        graph = ROOT.TGraphErrors(self.Attributes['binsx'])
        count = 0
        for i in xrange(len(list_of_bins)):
            x_ , _ = list_of_bins[i].GetBinCenter()
            MPV_ = list_of_bins[i].Fit['MPV']
            MPVerr_ = list_of_bins[i].Fit['MPVErr']
            if MPV_ is not None:
                graph.SetPoint(count, x_ , MPV_)
                graph.SetPointError(count, 0, MPVerr_)
                signals.append(MPV_)
                count += 1
            else:
                signals.append(0)
        if show:
            canvas = ROOT.TCanvas('signal_in_row', 'signals in row at height {:.3f}'.format(rowheight))
            canvas.cd()
            graph.SetTitle('Most Probable Signal response in bins in row at height {:.3f}'.format(rowheight))
            graph.GetXaxis().SetTitle("x position / cm")
            graph.GetYaxis().SetTitle("Signal")
            graph.Draw('AP')
            self.SavePlots('signals_in_row{:.3f}'.format(rowheight),'pdf')
            self.IfWait('show signal in row {:.3f}..'.format(rowheight))
        return signals

    def CreateSignalHistogram(self,saveplot = False, scale = False, showfit = True):
        canvas = ROOT.TCanvas('canvas', 'canvas')
        canvas.cd()
        self.SignalHisto.GetXaxis().SetTitle("Signal Response [ADC Units]")
        self.SignalHisto.GetYaxis().SetTitle("Counts")
        if scale or showfit:
            maximum = self.SignalHisto.GetMaximum()
            scale_val = 1./maximum
            self.SignalHisto.Scale(scale_val)
            scale_str = " (scaled to 1)"
            tmpTitle = self.SignalHisto.GetTitle()
            self.SignalHisto.SetTitle(tmpTitle+scale_str)
        self.SignalHisto.Draw()
        if saveplot:
            self.SavePlots('TotalSignalDistribution', 'png')
        if showfit:
            # Fitting SNR histo
            print "Fitting...\n"

            # Setting fit range and start values
            fr = array("d", [0,0])
            sv = array("d", [0,0,0,0])
            pllo = array("d", [0,0,0,0])
            plhi = array("d", [0,0,0,0])
            fp = array("d", [0,0,0,0])
            fpe = array("d", [0,0,0,0])

            # fitrange for fit:
            fr[0]=0.3*self.SignalHisto.GetMean()
            fr[1]=3.0*self.SignalHisto.GetMean()

            # lower parameter limits
            pllo[0]=2.0 # par[0]=Width (scale) parameter of Landau density
            pllo[1]=20.0 # par[1]=Most Probable (MP, location) parameter of Landau density
            pllo[2]=1.0 # par[2]=Total area (integral -inf to inf, normalization constant)
            pllo[3]=1.0 # par[3]=Width (sigma) of convoluted Gaussian function

            # upper parameter limits
            plhi[0]=50.0         # par[0]=Width (scale) parameter of Landau density
            plhi[1]=250.0        # par[1]=Most Probable (MP, location) parameter of Landau density
            plhi[2]=10000.0   # par[2]=Total area (integral -inf to inf, normalization constant)
            plhi[3]=50.0         # par[3]=Width (sigma) of convoluted Gaussian function

            # Startvalues for fit:
            sv[0]=11.0       # par[0]=Width (scale) parameter of Landau density
            sv[1]=self.SignalHisto.GetBinCenter(self.SignalHisto.GetMaximumBin())      # par[1]=Most Probable (MP, location) parameter of Landau density
            sv[2]=40.0   # par[2]=Total area (integral -inf to inf, normalization constant)
            sv[3]=11.0       # par[3]=Width (sigma) of convoluted Gaussian function

            chisqr = array("d", [0]) # returns the chi square
            ndf = array("d", [0]) # returns ndf

            self.LangauFitFunction = langaufit(self.SignalHisto,fr,sv,pllo,plhi,fp,fpe,chisqr,ndf)

            SNRPeak = array("d", [0])
            SNRFWHM = array("d", [0])
            langaupro(fp,SNRPeak,SNRFWHM)

            print "Fitting done\n"

            self.LangauFitFunction.Draw("lsame")
            raw_input("WAIT!")
        self.IfWait('Signal Histogram drawn')

    def FindMaxima(self, threshold = None, minimum_bincount = 5, show = False):
        '''

        :param threshold:
        :param show:
        :return: list of coordinates
        '''
        self.MaximaSearch = FindMaxima(self)
        self.MaximaSearch.Find(threshold=threshold, minimum_bincount=minimum_bincount, show=show)

    def FindMinima(self, threshold = None, minimum_bincount = 5, show = False):
        '''

        :param threshold:
        :param show:
        :return: list of coordinates
        '''
        self.MinimaSearch = FindMinima(self)
        self.MinimaSearch.Find(threshold=threshold, minimum_bincount=minimum_bincount, show=show)

    def UpdateBinAttributes(self):
        for i in xrange(len(self.ListOfBins)):
            self.ListOfBins[i].UpdateAttributes()

