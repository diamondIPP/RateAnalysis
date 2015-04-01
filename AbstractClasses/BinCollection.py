import ROOT
import numpy as np
import os
import types as t
from Bin import Bin


class BinCollection(object):
    '''
    A BinCollection Object Contains Bin-objects. It is used to store the data, as well
    as to read the data, make selections or to make Plots related to collection of many bins.

    '''

    def __init__(self, binsx, xmin, xmax, binsy, ymin, ymax):
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

        self.counthisto = ROOT.TH2D('counthisto', '2D hit distribution', *self.Get2DAttributes())
        self.totalsignal = ROOT.TH2D('totalsignal', '2D total signal distribution', *self.Get2DAttributes())

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

    def ShowBinXYSignalHisto(self,x,y,saveplot = False):
        '''
        Shows a Histogram of the Signal response distribution inside the bin which
        contains the coordinates (x,y)
        :param x: coordinate x in cm which is contained in the bin of interest
        :param y: coordinate y in cm which is contained in the bin of interest
        :param saveplot: if True save plot as as Results/Bin_X0.123Y-0.123_SignalHisto.png
        :return: -
        '''
        self.ListOfBins[self.GetBinNumber(x,y)].CreateBinSignalHisto(saveplot)

    def CalculateMeanSignalDistribution(self, minimum_bincontent = 1):
        '''

        :param minimum_bincontent:
        :return:
        '''
        assert (minimum_bincontent > 0), "minimum_bincontent has to be a positive integer"
        self.meansignaldistribution = ROOT.TH2D('meansignaldistribution', "Mean Signal Distribution", *self.Get2DAttributes())
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
            raw_input("Selected bins shown")
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

        :param heigth:
        :return:
        '''
        start_x = self.Attributes['XMIN']+ self.Attributes['binwidth_x']/2.
        start_bin = self.GetBinByCoordinates(start_x, height)
        start_binnumber = start_bin.GetBinNumber()
        list_of_bins = [self.ListOfBins[i] for i in xrange(start_binnumber, start_binnumber+self.Attributes['binsx'])]
        return list_of_bins

    def SavePlots(self, savename, ending, saveDir):
        # Results directories:
        #resultsdir = saveDir+'run_'+str(self.run_object.run_number)+'/' # eg. 'Results/run_364/'
        resultsdir = saveDir # eg. 'Results/run_364/'
        if not os.path.exists(resultsdir):
            os.makedirs(resultsdir)

        ROOT.gPad.Print(resultsdir+savename+'.'+ending)

    def GetSignalInColumn(self, position , show = False):
        '''

        :param position:
        :return:
        '''
        list_of_bins = self.GetBinsInColumn(position)
        signals = []
        attributes = self.Attributes['binsy'], self.Attributes['YMIN'], self.Attributes['YMAX']
        histo = ROOT.TH1D('histo', 'bins in column', *attributes)
        for i in xrange(len(list_of_bins)):
            _ , y_ = list_of_bins[i].GetBinCenter()
            signal_ = list_of_bins[i].GetMean()
            histo.Fill(y_, signal_)
            signals.append(signal_)
        if show:
            canvas = ROOT.TCanvas('signal_in_column', 'signals in column')
            canvas.cd()
            histo.Draw()
            raw_input('show signal in column..')
        return signals

    def GetSignalInRow(self, height, show = False):
        '''

        :param height:
        :return:
        '''
        list_of_bins = self.GetBinsInRow(height)
        signals = []
        attributes = self.Attributes['binsx'], self.Attributes['XMIN'], self.Attributes['XMAX']
        histo = ROOT.TH1D('histo', 'bins in row', *attributes)
        for i in xrange(len(list_of_bins)):
            x_ , _ = list_of_bins[i].GetBinCenter()
            signal_ = list_of_bins[i].GetMean()
            histo.Fill(x_, signal_)
            signals.append(signal_)
        if show:
            canvas = ROOT.TCanvas('signal_in_row', 'signals in row')
            canvas.cd()
            histo.Draw()
            raw_input('show signal in row..')
        return signals

    def FindMaxima(self, threshold = None, minimum_bincount = 5,show = False):
        '''

        :param threshold:
        :param show:
        :return: list of coordinates
        '''
        list_of_coordinates = []
        if threshold is None:
            threshold = self.meansignaldistribution.GetMean()
        self.voting_histo = ROOT.TH2D('voting_histo', 'Voting', *self.Get2DAttributes())
        bin_SW_coordinates = self.ListOfBins[self.Attributes['binsx']+3].GetBinCenter()
        bin_SE_coordinates = self.ListOfBins[2*self.Attributes['binsx']+2].GetBinCenter()
        binsx_, xmin_, xmax_, binsy_, ymin_, ymax_ = self.Get2DAttributes()

        # horizontal scan:
        def horizontal_scan(self):
            height = bin_SW_coordinates[1]
            # scan rows:
            while height < ymax_:

                signals = self.GetSignalInRow(height)
                bins = self.GetBinsInRow(height)

                for i in xrange(len(signals)-2):
                    if signals[i] < signals[i+1] and signals[i+2] < signals[i+1] and signals[i+1] > threshold and bins[i+1].GetEntries() >= minimum_bincount:
                        binnumber = bins[i+1].GetBinNumber()
                        FillHistoByBinnumber(self, binnumber, 1)

                height += self.Attributes['binwidth_y']

        # vertical scan:
        def vertical_scan(self):
            position = bin_SW_coordinates[0]

            # scan rows:
            while position < xmax_:

                signals = self.GetSignalInColumn(position)
                bins = self.GetBinsInColumn(position)

                for i in xrange(len(signals)-2):
                    if signals[i] < signals[i+1] and signals[i+2] < signals[i+1] and signals[i+1] > threshold and bins[i+1].GetEntries() >= minimum_bincount:
                        binnumber = bins[i+1].GetBinNumber()
                        FillHistoByBinnumber(self, binnumber, 1)

                position += self.Attributes['binwidth_x']

        # southwest to northeast scan:
        def SWNE_scan(self):
            # creating an array containing all bin numbers
            bin_numbers = np.ones((binsx_+2, binsy_+2))
            for i in xrange(binsy_+2):
                bin_numbers[:,i] = np.arange(i*(binsx_+2), (i+1)*(binsx_+2))
                window_bin_numbers = bin_numbers[1:binsx_+1,1:binsy_+1]

            def CheckBinInsideWindow(binnumber):
                '''
                checks if the bin number 'binnumber' lies inside the window of binsx * binsy
                :param binnumber: the binnumber under test
                :return: True if binnumber is contained inside window
                '''
                itemindex = np.where(window_bin_numbers == binnumber)
                if np.size(itemindex) > 0:
                    bin_in_window = True
                else:
                    bin_in_window = False
                return bin_in_window

            # creating the vertical start array of bins
            vertical_bins = self.GetBinsInColumn(position=bin_SW_coordinates[0])
            vertical_binnumbers = []
            for i in xrange(len(vertical_bins)):
                vertical_binnumbers.append(vertical_bins[i].GetBinNumber())
            vertical_starts = vertical_binnumbers[::-1][2:] # reversing the list and cutting to not to start from edge

            def scan_swne(start_list):
                for nr in start_list:
                    left_nr = nr
                    middle_nr = left_nr + binsx_ + 3
                    right_nr = middle_nr + binsx_ + 3
                    while CheckBinInsideWindow(right_nr):
                        left_signal = self.ListOfBins[left_nr].GetMean()
                        middle_signal = self.ListOfBins[middle_nr].GetMean()
                        right_signal = self.ListOfBins[right_nr].GetMean()
                        if left_signal < middle_signal and right_signal < middle_signal and middle_signal > threshold and self.ListOfBins[middle_nr].GetEntries() >= minimum_bincount:
                            FillHistoByBinnumber(self, middle_nr, 1)
                        left_nr = middle_nr
                        middle_nr = right_nr
                        right_nr = right_nr + binsx_ + 3

            scan_swne(vertical_starts)

            horizontal_bins = self.GetBinsInRow(height=bin_SW_coordinates[1])
            horizontal_binnumbers = []
            for i in xrange(len(horizontal_bins)):
                horizontal_binnumbers.append(horizontal_bins[i].GetBinNumber())
            horizontal_starts = horizontal_binnumbers[1:-2] # cutting to not to start from edge

            scan_swne(horizontal_starts)

        def SENW_scan(self):
            # creating an array containing all bin numbers
            bin_numbers = np.ones((binsx_+2, binsy_+2))
            for i in xrange(binsy_+2):
                bin_numbers[:,i] = np.arange(i*(binsx_+2), (i+1)*(binsx_+2))
                window_bin_numbers = bin_numbers[1:binsx_+1,1:binsy_+1]

            def CheckBinInsideWindow(binnumber):
                '''
                checks if the bin number 'binnumber' lies inside the window of binsx * binsy
                :param binnumber: the binnumber under test
                :return: True if binnumber is contained inside window
                '''
                itemindex = np.where(window_bin_numbers == binnumber)
                if np.size(itemindex) > 0:
                    bin_in_window = True
                else:
                    bin_in_window = False
                return bin_in_window

            # creating the vertical start array of bins
            vertical_bins = self.GetBinsInColumn(position=bin_SE_coordinates[0])
            vertical_binnumbers = []
            for i in xrange(len(vertical_bins)):
                vertical_binnumbers.append(vertical_bins[i].GetBinNumber())
            vertical_starts = vertical_binnumbers[::-1][2:] # reversing the list and cutting to not to start from edge

            def scan_senw(start_list):
                for nr in start_list:
                    right_nr = nr
                    middle_nr = right_nr + binsx_ + 1
                    left_nr = middle_nr + binsx_ + 1
                    while CheckBinInsideWindow(left_nr):
                        left_signal = self.ListOfBins[left_nr].GetMean()
                        middle_signal = self.ListOfBins[middle_nr].GetMean()
                        right_signal = self.ListOfBins[right_nr].GetMean()
                        if left_signal < middle_signal and right_signal < middle_signal and middle_signal > threshold and self.ListOfBins[middle_nr].GetEntries() >= minimum_bincount:
                            FillHistoByBinnumber(self, middle_nr, 1)
                        right_nr = middle_nr
                        middle_nr = left_nr
                        left_nr = left_nr + binsx_ + 1

            scan_senw(vertical_starts)

            horizontal_bins = self.GetBinsInRow(height=bin_SE_coordinates[1])
            horizontal_binnumbers = []
            for i in xrange(len(horizontal_bins)):
                horizontal_binnumbers.append(horizontal_bins[i].GetBinNumber())
            horizontal_starts = horizontal_binnumbers[2:-1] # cutting to not to start from edge

            scan_senw(horizontal_starts)


        def FillHistoByBinnumber(self, binnumber, weight = 1):
            bin = self.GetBinByNumber(binnumber)
            x_, y_ = bin.GetBinCenter()
            self.voting_histo.Fill(x_, y_, weight)

        def GetBinsInVoteRange(votes, maxvotes=None):
            '''
            Returns all bin numbers which contain 'votes' number of votes or returns
            all bins which contain at least 'votes' number of votes and at most
            'maxvotes' number of votes
            :param value:
            :param maxvalue:
            :return: found_maxima
            '''
            if maxvotes == None:
                maxvotes = votes
            found_maxima = []
            for i in xrange((binsx_+2)*(binsy_+2)):
                content = self.voting_histo.GetBinContent(i)
                if content >= votes and content <= maxvotes:
                    found_maxima.append(i)
            return found_maxima

        horizontal_scan(self)
        vertical_scan(self)
        SWNE_scan(self)
        SENW_scan(self)


        maxima = GetBinsInVoteRange(4)
        print len(maxima)," maxima found: "
        print self.GetBinCenter(maxima)

        if show:
            vote_canvas = ROOT.TCanvas('vote_canvas', "find Max vote-histo")
            vote_canvas.cd()
            ROOT.gStyle.SetPalette(51)
            ROOT.gStyle.SetNumberContours(5)
            self.voting_histo.SetStats(False)
            self.voting_histo.Draw('colz')
            print self.voting_histo.GetMaximumBin()
            raw_input("vote histo drawn..")
            ROOT.gStyle.SetPalette(53)
            ROOT.gStyle.SetNumberContours(999)


    def UpdateBinAttributes(self):
        for i in xrange(len(self.ListOfBins)):
            self.ListOfBins[i].UpdateAttributes()

