import ROOT
from Bin import Bin

class BinCollection(object):

    def __init__(self, binsx, xmin, xmax, binsy, ymin, ymax):

        self.ListOfBins = [Bin(i,self) for i in xrange((binsx+2)*(binsy+2))]
        self.binnumbers = [i for i in xrange((binsx+2)*(binsy+2))]
        self.Attributes = {
            'binsx': binsx,
            'binsy': binsy,
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

    def ShowBinXYSignalHisto(self,x,y):
        self.ListOfBins[self.GetBinNumber(x,y)].CreateBinSignalHisto()

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


