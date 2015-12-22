import ROOT
import numpy as np
from array import array
from Elementary import Elementary



class FindExtrema(Elementary):
    
    def __init__(self, bincollection, verbose = False):
        Elementary.__init__(self, verbose=verbose)
        self.BinCollectionObj = bincollection
        self.channel = self.BinCollectionObj.channel
        self.ExtremaType = ''
        self.SetType()

        self.voting_histo_name = 'voting_histo_'+self.ExtremaType+str(self.GLOBAL_COUNT)
        self.voting_histo = ROOT.TH2D(self.voting_histo_name, 'Voting for '+self.ExtremaType, *self.BinCollectionObj.Get2DAttributes())
        self.bin_SW_coordinates = self.BinCollectionObj.listOfBins[self.BinCollectionObj.attributes['binsx']+3].GetBinCenter()
        self.bin_SE_coordinates = self.BinCollectionObj.listOfBins[2*self.BinCollectionObj.attributes['binsx']+2].GetBinCenter()
        self.binsx_, self.xmin_, self.xmax_, self.binsy_, self.ymin_, self.ymax_ = self.BinCollectionObj.Get2DAttributes()
        self.GLOBAL_COUNT += 1

    def __del__(self):
        ROOT.gROOT.Delete(self.voting_histo_name)
        del self.voting_histo

    def load_config(self):
        self.ShowAndWait = False

    def Local1DExtrema(self, signal_1, signal_2, signal_3, entries):
        pass

    def FifthVoting(self, nbhd_mean, mean):
        pass

    def SetType(self):
        pass

    def SetThreshold(self):
        pass

    # horizontal scan:
    def horizontal_scan(self):
        height = self.bin_SW_coordinates[1]
        # scan rows:
        while height < self.ymax_:

            signals = self.BinCollectionObj.GetSignalInRow(height)
            bins = self.BinCollectionObj.GetBinsInRow(height)

            for i in xrange(len(signals)-2):
                if self.Local1DExtrema(signals[i], signals[i+1], signals[i+2], bins[i+1].GetEntries()):
                    binnumber = bins[i+1].GetBinNumber()
                    self.FillHistoByBinnumber(binnumber, 1)

            height += self.BinCollectionObj.attributes['binwidth_y']

    # vertical scan:
    def vertical_scan(self):
        position = self.bin_SW_coordinates[0]

        # scan rows:
        while position < self.xmax_:

            signals = self.BinCollectionObj.GetSignalInColumn(position)
            bins = self.BinCollectionObj.GetBinsInColumn(position)

            for i in xrange(len(signals)-2):
                if self.Local1DExtrema(signals[i], signals[i+1], signals[i+2], bins[i+1].GetEntries()):
                    binnumber = bins[i+1].GetBinNumber()
                    self.FillHistoByBinnumber(binnumber, 1)

            position += self.BinCollectionObj.attributes['binwidth_x']

    # southwest to northeast scan:
    def SWNE_scan(self):
        # creating an array containing all bin numbers
        bin_numbers = np.ones((self.binsx_+2, self.binsy_+2))
        for i in xrange(self.binsy_+2):
            bin_numbers[:,i] = np.arange(i*(self.binsx_+2), (i+1)*(self.binsx_+2))
            window_bin_numbers = bin_numbers[1:self.binsx_+1,1:self.binsy_+1]

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
        vertical_bins = self.BinCollectionObj.GetBinsInColumn(position=self.bin_SW_coordinates[0])
        vertical_binnumbers = []
        for i in xrange(len(vertical_bins)):
            vertical_binnumbers.append(vertical_bins[i].GetBinNumber())
        vertical_starts = vertical_binnumbers[::-1][2:] # reversing the list and cutting to not to start from edge

        def scan_swne(start_list):
            for nr in start_list:
                left_nr = nr
                middle_nr = left_nr + self.binsx_ + 3
                right_nr = middle_nr + self.binsx_ + 3
                while CheckBinInsideWindow(right_nr):
                    left_signal = self.BinCollectionObj.listOfBins[left_nr].GetMean()
                    middle_signal = self.BinCollectionObj.listOfBins[middle_nr].GetMean()
                    right_signal = self.BinCollectionObj.listOfBins[right_nr].GetMean()
                    if self.Local1DExtrema(left_signal, middle_signal, right_signal, self.BinCollectionObj.listOfBins[middle_nr].GetEntries()):
                        self.FillHistoByBinnumber(middle_nr, 1)
                    left_nr = middle_nr
                    middle_nr = right_nr
                    right_nr = right_nr + self.binsx_ + 3

        scan_swne(vertical_starts)

        horizontal_bins = self.BinCollectionObj.GetBinsInRow(height=self.bin_SW_coordinates[1])
        horizontal_binnumbers = []
        for i in xrange(len(horizontal_bins)):
            horizontal_binnumbers.append(horizontal_bins[i].GetBinNumber())
        horizontal_starts = horizontal_binnumbers[1:-2] # cutting to not to start from edge

        scan_swne(horizontal_starts)

    def SENW_scan(self):
        # creating an array containing all bin numbers
        bin_numbers = np.ones((self.binsx_+2, self.binsy_+2))
        for i in xrange(self.binsy_+2):
            bin_numbers[:,i] = np.arange(i*(self.binsx_+2), (i+1)*(self.binsx_+2))
            window_bin_numbers = bin_numbers[1:self.binsx_+1,1:self.binsy_+1]

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
        vertical_bins = self.BinCollectionObj.GetBinsInColumn(position=self.bin_SE_coordinates[0])
        vertical_binnumbers = []
        for i in xrange(len(vertical_bins)):
            vertical_binnumbers.append(vertical_bins[i].GetBinNumber())
        vertical_starts = vertical_binnumbers[::-1][2:] # reversing the list and cutting to not to start from edge

        def scan_senw(start_list):
            for nr in start_list:
                right_nr = nr
                middle_nr = right_nr + self.binsx_ + 1
                left_nr = middle_nr + self.binsx_ + 1
                while CheckBinInsideWindow(left_nr):
                    left_signal = self.BinCollectionObj.listOfBins[left_nr].GetMean()
                    middle_signal = self.BinCollectionObj.listOfBins[middle_nr].GetMean()
                    right_signal = self.BinCollectionObj.listOfBins[right_nr].GetMean()
                    if self.Local1DExtrema(left_signal, middle_signal, right_signal, self.BinCollectionObj.listOfBins[middle_nr].GetEntries()):
                        self.FillHistoByBinnumber(middle_nr, 1)
                    right_nr = middle_nr
                    middle_nr = left_nr
                    left_nr = left_nr + self.binsx_ + 1

        scan_senw(vertical_starts)

        horizontal_bins = self.BinCollectionObj.GetBinsInRow(height=self.bin_SE_coordinates[1])
        horizontal_binnumbers = []
        for i in xrange(len(horizontal_bins)):
            horizontal_binnumbers.append(horizontal_bins[i].GetBinNumber())
        horizontal_starts = horizontal_binnumbers[2:-1] # cutting to not to start from edge

        scan_senw(horizontal_starts)

    def FillHistoByBinnumber(self, binnumber, weight = 1):
        bin = self.BinCollectionObj.GetBinByNumber(binnumber)
        x_, y_ = bin.GetBinCenter()
        self.voting_histo.Fill(x_, y_, weight)

    def GetBinsInVoteRange(self, votes, maxvotes=None):
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
        for i in xrange((self.binsx_+2)*(self.binsy_+2)):
            content = self.voting_histo.GetBinContent(i)
            if content >= votes and content <= maxvotes:
                found_maxima.append(i)
        return found_maxima

    def GetBinsInNbhd(self, binnumber, include_center = False, extended = False):
        '''
        Returns the bin numbers of the surrounding bins of bin 'binnumber'
        :param binnumber:
        :param include_center:
        :return:
        '''
        binsx = self.BinCollectionObj.attributes['binsx']
        range_x = binsx + 2
        nbhd = []
        nbhd.append(binnumber+range_x-1)
        nbhd.append(binnumber+range_x)
        nbhd.append(binnumber+range_x+1)
        nbhd.append(binnumber-1)
        if include_center:
            nbhd.append(binnumber)
        nbhd.append(binnumber+1)
        nbhd.append(binnumber-range_x-1)
        nbhd.append(binnumber-range_x)
        nbhd.append(binnumber-range_x+1)
        if extended:
            nbhd.append(binnumber+2*range_x-1)
            nbhd.append(binnumber+2*range_x)
            nbhd.append(binnumber+2*range_x+1)
            nbhd.append(binnumber+range_x-2)
            nbhd.append(binnumber+range_x+2)
            nbhd.append(binnumber-2)
            nbhd.append(binnumber+2)
            nbhd.append(binnumber-range_x-2)
            nbhd.append(binnumber-range_x+2)
            nbhd.append(binnumber-2*range_x-1)
            nbhd.append(binnumber-2*range_x)
            nbhd.append(binnumber-2*range_x+1)
        return nbhd

    def GetBinsInRadius(self, center_bin, radius, include_center = True): # maybe move this to bincollection class
        binnnumbers = self.BinCollectionObj.binnumbers
        c_x, c_y = self.BinCollectionObj.GetBinCenter(center_bin)
        list_of_bins = []
        for binnumber in binnnumbers:
            x, y = self.BinCollectionObj.GetBinCenter(binnumber)
            if 0 < (c_x-x)**2 + (c_y-y)**2 <= radius**2 or ((c_x-x)**2 + (c_y-y)**2 == 0 and include_center):
                list_of_bins.append(binnumber)
        return list_of_bins

    def Find(self, threshold = None, minimum_bincount = 5, show = False):

        if threshold == None:
            self.SetThreshold()
            # self.threshold = self.BinCollectionObj.SignalHisto.GetMean()
            self.verbose_print('threshold set automatically to ', self.threshold)
        else:
            self.threshold = threshold
        self.minimum_bincount = minimum_bincount
        self.show = show

        self.horizontal_scan()
        self.vertical_scan()
        self.SWNE_scan()
        self.SENW_scan()
    
        Vote4Extrema = self.GetBinsInVoteRange(4)
        self.verbose_print(len(Vote4Extrema), " " + self.ExtremaType + " found containing 4 votings")
    
        # mean = self.BinCollectionObj.SignalHisto.GetMean()
        if hasattr(self.BinCollectionObj.parent_analysis_obj, "MeanSignalHisto"):
            mean = self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetMean()
        else:
            self.BinCollectionObj.parent_analysis_obj.CreateMeanSignalHistogram()
            mean = self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetMean()

        for i in xrange(len(Vote4Extrema)):
            center_bin = Vote4Extrema[i]
            nbhd = self.GetBinsInNbhd(center_bin)
            MeanNbhdSignals = []
            for binnr in nbhd:
                bin_mean = self.BinCollectionObj.listOfBins[binnr].GetMean()
                bin_entries = self.BinCollectionObj.listOfBins[binnr].GetEntries()
                if bin_entries >= self.minimum_bincount:
                    MeanNbhdSignals.append(bin_mean)
            if len(MeanNbhdSignals)>0:
                nbhd_mean = np.mean(MeanNbhdSignals)
                if self.FifthVoting(nbhd_mean, mean):
                    self.FillHistoByBinnumber(center_bin, 1)
    
        Vote5Extrema = self.GetBinsInVoteRange(5)
        self.verbose_print(len(Vote5Extrema), " " + self.ExtremaType + " found containing 5 votings: ")
        self.verbose_print(self.BinCollectionObj.GetBinCenter(Vote5Extrema))
    
    
        # If Monte Carlo, match with real peak positions:
        if self.BinCollectionObj.run_object.IsMonteCarlo and self.ExtremaType == 'Maxima':
            npeaks = int(self.BinCollectionObj.run_object.SignalParameters[0])
            #peak_height = self.BinCollectionObj.run_object.SignalParameters[2]
            peaks_x = []
            peaks_y = []
            Ghosts = []
            Ninjas = []
            Peak_nbhd = [] # all bins around eff peaks
            Maxima_nbhd = [] # all bins around found maxima
            self.real_peaks = ROOT.TGraphErrors()
            for i in xrange(npeaks):
                x = self.BinCollectionObj.run_object.SignalParameters[3+4*i]
                y = self.BinCollectionObj.run_object.SignalParameters[4+4*i]
                # sigma_x = self.BinCollectionObj.run_object.SignalParameters[5+4*i]
                # sigma_y = self.BinCollectionObj.run_object.SignalParameters[6+4*i]
                # sigma = 1.*(sigma_x+sigma_y)/2.
                peaks_x.append(x)
                peaks_y.append(y)
                # Peak_nbhd_temp = self.GetBinsInNbhd(self.BinCollectionObj.GetBinNumber(x, y), include_center=True, extended=True)
                Peak_nbhd_temp = self.GetBinsInRadius(self.BinCollectionObj.GetBinNumber(x, y), radius=0.02, include_center=True)
                Peak_nbhd += Peak_nbhd_temp
                self.real_peaks.SetPoint(i, x, y)
                self.real_peaks.SetPointError(i, self.BinCollectionObj.run_object.SignalParameters[5+4*i], self.BinCollectionObj.run_object.SignalParameters[6+4*i])
            for bin_nr in Vote5Extrema:
                # Maxima_nbhd_temp = self.GetBinsInNbhd(bin_nr, include_center=True, extended=True)
                Maxima_nbhd_temp = self.GetBinsInRadius(bin_nr, radius=0.02, include_center=True)
                Maxima_nbhd += Maxima_nbhd_temp
            # Looking for Ghosts:
            for bin_nr in Vote5Extrema:
                if not bin_nr in Peak_nbhd:
                    self.verbose_print('Ghost peak found at position ({0:.3f}/{1:.3f})'.format(*self.BinCollectionObj.GetBinCenter(bin_nr)))
                    Ghosts.append(bin_nr)
            # Looking for Ninjas:
            for i in xrange(npeaks):
                bin_nr = self.BinCollectionObj.GetBinNumber(peaks_x[i],peaks_y[i])
                if not bin_nr in Maxima_nbhd:
                    self.verbose_print('Ninja peak found at position ({0:.3f}/{1:.3f})'.format(peaks_x[i], peaks_y[i]))
                    Ninjas.append(bin_nr)
            if npeaks > 0:
                print "\n{0:.1f}% of generated peaks found.".format(100.*(npeaks-len(Ninjas))/npeaks)
            print len(Ninjas)," Ninjas."
            print len(Ghosts)," Ghosts.\n"
            self.BinCollectionObj.parent_analysis_obj.extremaResults[self.channel]['TrueNPeaks'] = npeaks
            self.BinCollectionObj.parent_analysis_obj.extremaResults[self.channel]['Ninjas'] = len(Ninjas)
            self.BinCollectionObj.parent_analysis_obj.extremaResults[self.channel]['Ghosts'] = len(Ghosts)
        self.BinCollectionObj.parent_analysis_obj.extremaResults[self.channel]['FoundN'+self.ExtremaType] = len(Vote5Extrema)
        self.BinCollectionObj.parent_analysis_obj.extremaResults[self.channel]['Found'+self.ExtremaType] = self.BinCollectionObj.GetBinCenter(Vote5Extrema) # store REAL parent analysis, not parent.ExtremeAnalysis

        self.found_extrema = ROOT.TGraph()
        for i in xrange(len(Vote5Extrema)):
            self.found_extrema.SetPoint(i, *(self.BinCollectionObj.GetBinCenter(Vote5Extrema)[i]))
    
        if self.show:
            vote_canvas = ROOT.TCanvas('vote_canvas_'+self.ExtremaType, "find "+self.ExtremaType+" vote-histo", 500, 500)
            vote_canvas.cd()
            ROOT.gStyle.SetPalette(51)
            ROOT.gStyle.SetNumberContours(6)
            self.voting_histo.SetStats(False)
            self.voting_histo.Draw('colz')
            self.found_extrema.SetMarkerStyle(29)
            self.found_extrema.SetMarkerSize(2)
            self.found_extrema.SetMarkerColor(ROOT.kMagenta)
            self.found_extrema.Draw('SAME P0')
            if self.BinCollectionObj.run_object.IsMonteCarlo and self.ExtremaType == 'Maxima':
                self.real_peaks.SetMarkerStyle(20)
                self.real_peaks.SetMarkerSize(2)
                self.real_peaks.SetMarkerColor(ROOT.kRed)
                self.real_peaks.SetLineColor(ROOT.kRed)
                self.real_peaks.Draw('SAME P0')
            vote_canvas.Update()
            self.if_wait("Vote Histogram drawn")
            self.BinCollectionObj.save_plots('vote_histo_' + self.ExtremaType, 'png')
            # ROOT.gStyle.SetPalette(53)
            # ROOT.gStyle.SetNumberContours(999)


class FindMaxima(FindExtrema):
    
    def Local1DExtrema(self, signal_1, signal_2, signal_3, entries):
        if signal_1 < signal_2 and signal_3 < signal_2 and signal_2 > self.threshold and entries >= self.minimum_bincount:
            return True
        else:
            return False

    def FifthVoting(self, nbhd_mean, mean):
        if nbhd_mean > 1.04*mean:
            return True
        else:
            return False
        
    def SetType(self):
        self.ExtremaType = 'Maxima'

    def SetThreshold(self):
        percent = 55
        q = array('d', [1.*percent/100.])
        y = array('d', [0])

        if hasattr(self.BinCollectionObj.parent_analysis_obj, "MeanSignalHisto"):
            if self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto.has_key(self.channel):
                self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetQuantiles(1, y, q)
            else:
                self.BinCollectionObj.parent_analysis_obj.CreateMeanSignalHistogram(channel=self.channel)
                self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetQuantiles(1, y, q)
        else:
            self.BinCollectionObj.parent_analysis_obj.CreateMeanSignalHistogram(channel=self.channel)
            self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetQuantiles(1, y, q)

        # self.BinCollectionObj.SignalHisto.GetQuantiles(1, y, q)
        self.threshold = y[0]
        
class FindMinima(FindExtrema):
    
    def Local1DExtrema(self, signal_1, signal_2, signal_3, entries):
        if signal_1 > signal_2 and signal_3 > signal_2 and signal_2 < self.threshold and entries >= self.minimum_bincount:
            return True
        else:
            return False

    def FifthVoting(self, nbhd_mean, mean):
        if nbhd_mean < 0.97*mean:
            return True
        else:
            return False
        
    def SetType(self):
        self.ExtremaType = 'Minima'

    def SetThreshold(self):
        percent = 45
        q = array('d', [1.*percent/100.])
        y = array('d', [0])

        if hasattr(self.BinCollectionObj.parent_analysis_obj, "MeanSignalHisto"):
            if self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto.has_key(self.channel):
                self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetQuantiles(1, y, q)
            else:
                self.BinCollectionObj.parent_analysis_obj.CreateMeanSignalHistogram(channel=self.channel)
                self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetQuantiles(1, y, q)
        else:
            self.BinCollectionObj.parent_analysis_obj.CreateMeanSignalHistogram(channel=self.channel)
            self.BinCollectionObj.parent_analysis_obj.MeanSignalHisto[self.channel].GetQuantiles(1, y, q)

        # self.BinCollectionObj.SignalHisto.GetQuantiles(1, y, q)
        self.threshold = y[0]