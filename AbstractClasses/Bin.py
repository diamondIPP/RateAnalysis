import ROOT
from AbstractClasses.Elementary import Elementary


class Bin(Elementary):
    '''
    A bin contains localized signal data.
    '''
    def __init__(self, binnumber, bincollectionobject, verbose = False):
        '''
        Constructor of a Bin object
        :param binnumber: indistinguishable number
        :param bincollectionobject: the bincollection object which this bin belongs to
        :return: -
        '''
        Elementary.__init__(self, verbose=verbose)
        self.signals = []
        self.bincollectionobject = bincollectionobject
        self.selected = False
        self.Attributes = {
            'binnumber': binnumber,
            'average': 0.,
            'entries': 0.,
            'sigma': 0.,
            'IsActive': True,# not yet implemented
            'minx': 0.,
            'maxx': 0.,
            'miny': 0.,
            'maxy': 0.,
            'bincenter_x': 0.,
            'bincenter_y': 0.,
            'coordinate_x': 0, # int
            'coordinate_y': 0 # int
        }
        self.Fit = {
            'Function': '',
            'FunctionName': 'Landau',
            'Constant': None,
            'ConstantErr': None,
            'MPV': None,
            'MPVErr': None,
            'Sigma': None,
            'SigmaErr': None,
            'Chi2': None
        }
        self.BinSignalHisto = ROOT.TH1D('BinSignalHisto'+str(binnumber), 'Signal Distribution',500,0,500)

    def __del__(self):
        ROOT.gROOT.Delete('BinSignalHisto'+str(self.Attributes['binnumber']))

    def AddData(self, signal, update_bin_attributes = False):
        '''
        Fill data into this bin
        :param signal: data signal
        :param update_bin_attributes: if True updates the bin attributes after filling
        :return: -
        '''
        self.signals.append(signal)
        self.BinSignalHisto.Fill(signal)
        if update_bin_attributes: self.UpdateAttributes()

    def FitLandau(self):
        '''
        asdf
        :return:
        '''
        self.BinSignalHisto.Fit('landau', 'L0Q') # 'L': Log Likelihood, '0': do not draw
        self.Fit['Function'] = self.BinSignalHisto.GetFunction('landau')
        self.Fit['Constant'] = self.Fit['Function'].GetParameter(0)
        self.Fit['ConstantErr'] = self.Fit['Function'].GetParError(0)
        self.Fit['MPV'] = self.Fit['Function'].GetParameter(1)
        self.Fit['MPVErr'] = self.Fit['Function'].GetParError(1)
        self.Fit['Sigma'] = self.Fit['Function'].GetParameter(2)
        self.Fit['SigmaErr'] = self.Fit['Function'].GetParError(2)
        self.Fit['Chi2'] = self.Fit['Function'].GetChisquare()


    def CreateBinSignalHisto(self, saveplot = False, savedir = None, show_fit = False):
        '''
        Creates and draws a histogram of the signal response contained in this bin
        :param saveplot: if True, saves the plot as Results/Bin_X0.123Y-0.123_SignalHisto.png
        :return: -
        '''
        x_, y_ = self.GetBinCenter()
        canvasname = 'Bin {0:.3f} / {1:.3f} Signal Histo'.format(x_, y_)
        canvas = ROOT.TCanvas('canvas',canvasname)
        canvas.cd()
        self.BinSignalHisto.SetTitle(canvasname)
        self.BinSignalHisto.GetXaxis().SetTitle('Signal response')
        self.BinSignalHisto.Draw()
        if self.BinSignalHisto.GetEntries() >= 5 and show_fit:
            if self.Fit["MPV"] is None:
                self.FitLandau()
            self.VerbosePrint("Most Probable Signal Response: {0:.2f} +- {1:.2f}".format(self.Fit['MPV'],self.Fit['MPVErr']))
            self.VerbosePrint( "Sigma of Landau dist: {0:.3f} +- {1:.3f}".format(self.Fit['Sigma'],self.Fit['SigmaErr']))
            kNotDraw = 1<<9 # bit 9
            self.BinSignalHisto.GetFunction("landau").ResetBit(kNotDraw)
            self.BinSignalHisto.Draw()
        if saveplot:
            self.SavePlots('Bin_X{:.3f}Y{:.3f}_SignalHisto'.format(x_,y_),'png', saveDir=savedir)
        self.IfWait('Bin signal histo drawn')

    def GetBinCenter(self):
        '''
        Returns the x and y coordinates in cm of this bin
        :return: x_, y_ in cm
        '''
        self.UpdateAttributes()
        x_ = self.Attributes['bincenter_x']
        y_ = self.Attributes['bincenter_y']
        return x_, y_

    def GetBinCoordinates(self):
        '''
        Returns the x and y coordinates in number of bins from the left and from the bottom
        :return: x_, y_
        '''
        self.UpdateAttributes()
        x_ = self.Attributes['coordinate_x']
        y_ = self.Attributes['coordinate_y']
        return x_, y_

    def GetEntries(self):
        '''
        Returns the number of signals in this bin
        :return: self.Attributes['entries']
        '''
        self.UpdateAttributes()
        if int(self.Attributes['entries']) != len(self.signals):
            raw_input("WARNING: Error in number of entries.. "+str(self.Attributes['entries'])+" != "+str(len(self.signals)))
        return self.Attributes['entries']

    def GetMean(self):
        '''
        Returns the mean value of the signals contained in this bin
        :return:
        '''
        self.UpdateAttributes()
        return self.Attributes['average']

    def GetBinNumber(self):
        '''
        Returns the bin number of this bin
        :return: self.Attributes['binnumber']
        '''
        return self.Attributes['binnumber']

    def GetSigma(self):
        '''
        Returns the standard deviation of the signals contained in this bin.
        :return:
        '''
        self.UpdateAttributes()
        return self.Attributes['sigma']

    def UpdateAttributes(self):
        '''
        Updates the Attributes dict
        :return:
        '''
        self.Attributes['average'] = self.BinSignalHisto.GetMean()
        self.Attributes['sigma'] = self.BinSignalHisto.GetRMS()
        self.Attributes['entries'] = self.BinSignalHisto.GetEntries()

        collection_attributes = self.bincollectionobject.Attributes
        binsx_ = collection_attributes['binsx']+2
        binsy_ = collection_attributes['binsy']+2
        xmin_ = collection_attributes['XMIN']
        xmax_ = collection_attributes['XMAX']
        ymin_ = collection_attributes['YMIN']
        ymax_ = collection_attributes['YMAX']
        binwidthx_ = 1.*(xmax_-xmin_)/(binsx_-2)
        binwidthy_ = 1.*(ymax_-ymin_)/(binsy_-2)

        self.Attributes['coordinate_x'] = int(self.Attributes['binnumber'])%int(binsx_)
        self.Attributes['coordinate_y'] = int(self.Attributes['binnumber'])/int(binsx_)
        self.Attributes['bincenter_x'] = xmin_ - binwidthx_/2. + self.Attributes['coordinate_x']*binwidthx_
        self.Attributes['bincenter_y'] = ymin_ - binwidthy_/2. + self.Attributes['coordinate_y']*binwidthy_
        self.Attributes['minx'] = self.Attributes['bincenter_x'] - binwidthx_/2.
        self.Attributes['maxx'] = self.Attributes['bincenter_x'] + binwidthx_/2.
        self.Attributes['miny'] = self.Attributes['bincenter_y'] - binwidthy_/2.
        self.Attributes['maxy'] = self.Attributes['bincenter_y'] + binwidthy_/2.
