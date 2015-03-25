import ROOT


class Bin(object):

    def __init__(self, binnumber, bincollectionobject):

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
            'bincenter_x': 0,
            'bincenter_y': 0.,
            'coordinate_x': 0, # int
            'coordinate_y': 0 # int
        }
        self.BinSignalHisto = ROOT.TH1D('BinSignalHisto'+str(binnumber), ' Signal Distribution',500,0,500)

    def UpdateAttributes(self):
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

    def GetBinCenter(self):
        self.UpdateAttributes()
        x_ = self.Attributes['bincenter_x']
        y_ = self.Attributes['bincenter_y']
        return x_, y_

    def GetBinCoordinates(self):
        self.UpdateAttributes()
        x_ = self.Attributes['coordinate_x']
        y_ = self.Attributes['coordinate_y']
        return x_, y_

    def AddData(self, signal, update_bin_attributes = False):
        self.signals.append(signal)
        self.BinSignalHisto.Fill(signal)
        if update_bin_attributes: self.UpdateAttributes()

    def GetEntries(self):
        self.UpdateAttributes()
        if int(self.Attributes['entries']) != len(self.signals):
            raw_input("WARNING: Error in number of entries.. "+str(self.Attributes['entries'])+" != "+str(len(self.signals)))
        return self.Attributes['entries']

    def GetMean(self):
        self.UpdateAttributes()
        return self.Attributes['average']

    def GetSigma(self):
        self.UpdateAttributes()
        return self.Attributes['sigma']

    def CreateBinSignalHisto(self, saveplot = False):
        x_, y_ = self.GetBinCenter()
        canvasname = 'Bin '+str(x_)+' / '+str(y_)
        canvas = ROOT.TCanvas('canvas',canvasname+' Signal Histo')
        canvas.cd()
        title = canvasname+self.BinSignalHisto.GetTitle()
        self.BinSignalHisto.SetTitle(title)
        self.BinSignalHisto.GetXaxis().SetTitle('Signal response')
        self.BinSignalHisto.Draw()
        if saveplot:
            self.bincollectionobject.SavePlots('Bin_X{:.3f}Y{:.3f}_SignalHisto'.format(x_,y_),'png','Results/')
        raw_input('Bin signal histo drawn')