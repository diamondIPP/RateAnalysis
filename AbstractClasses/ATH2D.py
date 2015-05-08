from ROOT import TH2D
import types as t

class ATH2D(TH2D):

    def __init__(self, name, title, binsx, xmin, xmax, binsy, ymin, ymax):
        TH2D.__init__(self, name, title, binsx, xmin, xmax, binsy, ymin, ymax)

    def GetBinCenter(self, binnumber):
        assert(type(binnumber) == t.IntType), "binnumber has to be int type"

        xmin = self.GetXaxis().GetXmin()
        xmax = self.GetXaxis().GetXmax()
        nbinsx = self.GetXaxis().GetNbins()
        ymin = self.GetYaxis().GetXmin()
        ymax = self.GetYaxis().GetXmax()
        nbinsy = self.GetYaxis().GetNbins()

        binwidthx_ = 1.*(xmax-xmin)/nbinsx
        binwidthy_ = 1.*(ymax-ymin)/nbinsy

        coordinate_x = binnumber%int(nbinsx+2)
        coordinate_y = binnumber/int(nbinsx+2)
        bincenter_x = xmin - binwidthx_/2. + coordinate_x*binwidthx_
        bincenter_y = ymin - binwidthy_/2. + coordinate_y*binwidthy_

        return bincenter_x, bincenter_y

    def GetMaximumPosition(self):
        return self.GetBinCenter(self.GetMaximumBin())

    def GetBinsInNbhd(self, binnumber, include_center = False, extended = False):
        '''
        Returns the bin numbers of the surrounding bins of bin 'binnumber'
        :param binnumber:
        :param include_center:
        :return:
        '''
        binsx = self.GetXaxis().GetNbins()
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

