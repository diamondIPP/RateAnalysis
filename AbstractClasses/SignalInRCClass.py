from ROOT import TGraphErrors, TH1D, TCanvas, gPad, TGaxis
import ROOT

# Color for hit distribution (2nd graph):
kColor = 15

class SignalInRCGraph(object):

    def __init__(self, position, NBinsX, XMin, XMax):
        self.graph = TGraphErrors()
        self.hit_histo = TH1D('hit_histo', 'hit_histo', NBinsX, XMin, XMax)
        self.RowOrColumn = ''
        self.HeightOrPosition = ''
        self.direction = ''
        self.positionvalue = position
        self.SetOrientation()

    def SetOrientation(self):
        pass

    def SetGraphPoint(self, nr, x, y):
        self.graph.SetPoint(nr, x, y)

    def SetGraphPointError(self, nr, x, y):
        self.graph.SetPointError(nr, x, y)

    def SetHistPoint(self, x, counts):
        self.hit_histo.Fill(x, counts)

    def DrawGraph(self, drawoption='AP'):
        self.canvas = TCanvas('canvas', 'Signal Response Canvas')
        self.canvas.cd()
        self.graph.SetTitle('Signal response in bin '+self.RowOrColumn+' at '+self.HeightOrPosition+' {:.3f}'.format(self.positionvalue))
        self.graph.GetXaxis().SetTitle(self.direction+' position / cm')
        self.graph.GetYaxis().SetTitle('Mean Bin Signal Response')
        self.graph.Draw(drawoption)

    def DrawBoth(self, drawoption='AP'):
        self.DrawGraph(drawoption)
        self.canvas.Update()

        # scale hit histo to pad:
        self.graph.GetYaxis().SetRangeUser(min(gPad.GetUymin(), self.hit_histo.GetMinimum()),gPad.GetUymax())
        TitleFont = self.graph.GetYaxis().GetTitleFont()
        TitleSize = self.graph.GetYaxis().GetTitleSize()
        LabelFont = self.graph.GetYaxis().GetLabelFont()
        LabelSize = self.graph.GetYaxis().GetLabelSize()
        self.canvas.Update()
        rightmax = 1.1*self.hit_histo.GetMaximum()
        scale = gPad.GetUymax()/rightmax

        self.hit_histo.SetLineColor(kColor)
        self.hit_histo.Scale(scale)
        self.hit_histo.Draw('HIST SAME')
        self.canvas.Update()

        #draw axis on the right side:
        rightaxis = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(), 0, rightmax, 20210, '+L')
        rightaxis.SetLineColor(kColor)
        rightaxis.SetTextColor(kColor)
        rightaxis.SetTitle('# Hits in Bin')
        rightaxis.SetTitleFont(TitleFont)
        rightaxis.SetTitleSize(TitleSize)
        rightaxis.SetLabelFont(LabelFont)
        rightaxis.SetLabelSize(LabelSize)
        rightaxis.SetLabelColor(kColor)
        ROOT.SetOwnership(rightaxis, False)
        rightaxis.Draw('SAME')

class SignalInRowGraph(SignalInRCGraph):

    def __init__(self, height, NBinsX, XMin, XMax):
        SignalInRCGraph.__init__(self, position=height, NBinsX=NBinsX, XMin=XMin, XMax=XMax)

    def SetOrientation(self):
        self.RowOrColumn = 'Row'
        self.HeightOrPosition = 'Height'
        self.direction = 'x'

class SignalInColumnGraph(SignalInRCGraph):

    def __init__(self, height, NBinsY, YMin, YMax):
        SignalInRCGraph.__init__(self, position=height, NBinsX=NBinsY, XMin=YMin, XMax=YMax)

    def SetOrientation(self):
        self.RowOrColumn = 'Column'
        self.HeightOrPosition = 'Position'
        self.direction = 'y'