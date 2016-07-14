# ==============================================
# IMPORTS
# ==============================================
from ROOT import TGraphErrors, TCanvas, TH1D, TH2D, gStyle, TH1F, gROOT, TLegend, TCut, TGraph, TProfile2D, TH2F, TProfile, TCutG, kRed, kBlack, kBlue, kMagenta, kGreen, TF1, TPie
# from TelescopeAnalysis import Analysis
# from CurrentInfo import Currents
# from numpy import array
from math import sqrt, ceil, log
# from argparse import ArgumentParser
# from Extrema import Extrema2D
# from ChannelCut import ChannelCut
# from time import time, sleep
# from collections import OrderedDict
# from sys import stdout
from copy import deepcopy

__author__ = 'DA'


# ==============================================
# MAIN CLASS
# ==============================================

class Plots:
    def __init__(self, num_entries, run=312, do_hit_map=True, do_pulse_height=True, num_devices=7, binning=-1, runinfo=None):
        gStyle.SetPalette(53)
        gStyle.SetNumberContours(999)
        self.do_hit_map = do_hit_map
        self.do_pulse_height = do_pulse_height
        self.run = run
        self.runinfo = None
        self.binning = binning
        self.num_devices = num_devices
        self.num_entries = num_entries
        self.plot_settings = {
            'ph1Dbins': 100,
            'ph1Dmin': 0,
            'ph1Dmax': 50000,
            'maxphplots': int(ceil(8*self.num_entries/100)),
            'nBinsX': 80,
            'xmin': -6,
            'xmax': 6,
            'nBinsY': 120,
            'ymin': -6,
            'ymax': 6,
            'nBinCol': 52,
            'minCol': 0,
            'maxCol': 52,
            'nBinRow': 80,
            'minRow': 0,
            'maxRow': 80,
            'binWidthTGraph': 200000,
            'num_diff_cluster_sizes': 4
        }
        self.plot_settings['deltaX'] = float(self.plot_settings['xmax']-self.plot_settings['xmin'])/self.plot_settings['nBinsX']
        self.plot_settings['deltaY'] = float(self.plot_settings['ymax']-self.plot_settings['ymin'])/self.plot_settings['nBinsY']

    def create_TGraphErrors(self, title='tgraph', xTitle='X', yTitle='Y', linecolor=kBlack, markercolor=kBlack):
        graph = TGraphErrors()
        graph.SetNameTitle(title,title)
        graph.SetLineColor(linecolor)
        graph.SetMarkerColor(markercolor)
        graph.GetXaxis().SetTitle(xTitle)
        graph.GetYaxis().SetTitle(yTitle)
        return deepcopy(graph)

    def create_1D_histogram(self, name='histo', title='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0):
        ph1Dbins = self.plot_settings['ph1Dbins']
        ph1Dmin = self.plot_settings['ph1Dmin']
        ph1Dmax = self.plot_settings['ph1Dmax']
        histo1D = TH1D(name, title, int(ph1Dbins + 1), ph1Dmin - float(ph1Dmax - ph1Dmin)/(2*ph1Dbins), ph1Dmax + float(ph1Dmax - ph1Dmin)/(2*ph1Dbins))
        self.set_1D_options(histo1D, xTitle, yTitle, color, min_val)
        return deepcopy(histo1D)

    def set_1D_options(self, histo, xTitle='X', yTitle='Y', color=kBlack, min_val=0):
        histo.GetXaxis().SetTitle(xTitle)
        histo.GetYaxis().SetTitle(yTitle)
        histo.SetMaximum(self.plot_settings['maxphplots'])
        histo.SetMinimum(min_val)
        histo.SetLineColor(color)

    def create_2D_histogram(self, type='spatial', name='histo', title='histo', xTitle='X', yTitle='Y', min_val=0, max_val=-1):
        xbins = self.plot_settings['nBinsX'] if type is 'spatial' else self.plot_settings['nBinCol']
        xmin = self.plot_settings['xmin'] if type is 'spatial' else self.plot_settings['minCol']
        xmax = self.plot_settings['xmax'] if type is 'spatial' else self.plot_settings['maxCol']
        ybins = self.plot_settings['nBinsY'] if type is 'spatial' else self.plot_settings['nBinRow']
        ymin = self.plot_settings['ymin'] if type is 'spatial' else self.plot_settings['minRow']
        ymax = self.plot_settings['ymax'] if type is 'spatial' else self.plot_settings['maxRow']
        histo2D = TH2D(name, title, int(xbins + 1), xmin - float(xmax-xmin)/(2*xbins), xmax + float(xmax-xmin)/(2*xbins),
                       int(ybins + 1), ymin - float(ymax-ymin)/(2*ybins), ymax + float(ymax-ymin)/(2*ybins))
        self.set_2D_options(histo2D, xTitle, yTitle, min_val, max_val)
        return deepcopy(histo2D)

    def set_2D_options(self, histo, xTitle='X', yTitle='Y', min_val=0, max_val=-1):
        histo.GetXaxis().SetTitle(xTitle)
        histo.GetYaxis().SetTitle(yTitle)
        histo.SetMinimum(min_val)
        if max_val is not -1: histo.SetMaximum(max_val)

    def create_histograms(self):
        if self.do_hit_map:
            pass
        if self.do_pulse_height:
            # 1D
            self.phROC_all = {i: self.create_1D_histogram('phROC{n}_all'.format(n=i),
                                                          'Pulse Height ROC {n} all cluster sizes'.format(n=i), 'Charge (e)',
                                                          'Num Clusters', kBlack) for i in xrange(self.num_devices)}
            self.phROC_1cl = {i: self.create_1D_histogram('phROC{n}_1cl'.format(n=i),
                                                          'Pulse Height ROC {n} 1pix cluster'.format(n=i), 'Charge (e)',
                                                          'Num Clusters', kBlue) for i in xrange(self.num_devices)}
            self.phROC_2cl = {i: self.create_1D_histogram('phROC{n}_2cl'.format(n=i),
                                                          'Pulse Height ROC {n} 2pix cluster'.format(n=i), 'Charge (e)',
                                                          'Num Clusters', kGreen) for i in xrange(self.num_devices)}
            self.phROC_3cl = {i: self.create_1D_histogram('phROC{n}_3cl'.format(n=i),
                                                          'Pulse Height ROC {n} 3pix cluster'.format(n=i), 'Charge (e)',
                                                          'Num Clusters', kRed) for i in xrange(self.num_devices)}
            self.phROC_M4cl = {i: self.create_1D_histogram('phROC{n}_M4cl'.format(n=i),
                                                           'Pulse Height ROC {n} 4 or more pix cluster'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kMagenta) for i in xrange(self.num_devices)}
            # 2D
            self.hitMap = {i: self.create_2D_histogram('pixel', 'hitMapROC{n}'.format(n=i), 'Hit Map ROC {n}'.format(n=i),
                                                       'Column', 'Row', 0, -1) for i in xrange(self.num_devices)}
            self.avPhROC_local_1cl = {i: self.create_2D_histogram('spatial', 'avPh_ROC{n}_local_1cl'.format(n=i),
                                                                  'Average Pulse Height ROC {n} Local Coord. 1pix cluster'.format(n=i),
                                                                  'x (mm)', 'y (mm)', 0, -1) for i in xrange(self.num_devices)}
            self.phROC_hitMap_local_1cl = {i: self.create_2D_histogram('spatial', 'ph_ROC{n}_hitMap_local_1cl'.format(n=i),
                                                                       'Pulse Height ROC {n} Hit Map Local Coord. 1pix cluster'.format(n=i),
                                                                       'x (mm)', 'y (mm)', 0, -1) for i in xrange(self.num_devices)}
            self.avPhROC_telescope_1cl = {i: self.create_2D_histogram('spatial', 'avPh_ROC{n}_telescope_1cl',
                                                                      'Average Pulse Height ROC {n} telescope Coord. 1 pix cluster'.format(n=i),
                                                                      'x (mm)', 'y (mm)', 0, -1) for i in xrange(self.num_devices)}
            self.phROC_hitMap_telescope_1cl = {i: self.create_2D_histogram('spatial', 'ph_ROC{n}_hitMap_telescope_1cl'.format(n=i),
                                                                       'Pulse Height ROC {n} Hit Map telescope Coord. 1pix cluster'.format(n=i),
                                                                       'x (mm)', 'y (mm)', 0, -1) for i in xrange(self.num_devices)}
            self.avPhROC_pixelated_1cl = {i: self.create_2D_histogram('pixel', 'avPh_ROC{n}_pixelated_1cl',
                                                                      'Average Pulse Height ROC {n} pixelated Coord. 1 pix cluster'.format(n=i),
                                                                      'Column', 'Row', 0, -1) for i in xrange(self.num_devices)}
            self.phROC_hitMap_pixelated_1cl = {i: self.create_2D_histogram('pixel', 'ph_ROC{n}_hitMap_pixelated_1cl'.format(n=i),
                                                                       'Pulse Height ROC {n} Hit Map pixelated Coord. 1pix cluster'.format(n=i),
                                                                       'Column', 'Row', 0, -1) for i in xrange(self.num_devices)}
            # TGraphErrors
            self.meanPhROC_all = {i: self.create_TGraphErrors('meanPHROC{n}_all'.format(n=i), 'Event', 'Charge (e)',
                                                              kBlack, kBlack) for i in xrange(self.num_devices)}
            self.meanPhROC_1cl = {i: self.create_TGraphErrors('meanPHROC{n}_1cl'.format(n=i), 'Event', 'Charge (e)',
                                                              kBlack, kBlack) for i in xrange(self.num_devices)}
            self.meanPhROC_2cl = {i: self.create_TGraphErrors('meanPHROC{n}_2cl'.format(n=i), 'Event', 'Charge (e)',
                                                              kBlack, kBlack) for i in xrange(self.num_devices)}
            self.meanPhROC_3cl = {i: self.create_TGraphErrors('meanPHROC{n}_3cl'.format(n=i), 'Event', 'Charge (e)',
                                                              kBlack, kBlack) for i in xrange(self.num_devices)}
            self.meanPhROC_M4cl = {i: self.create_TGraphErrors('meanPHROC{n}_M4cl'.format(n=i), 'Event', 'Charge (e)',
                                                              kBlack, kBlack) for i in xrange(self.num_devices)}
            self.nPointsTGraph = 0

    def AverageBinHistogram(self, histToAv, histHits, binx, biny):
        temp = float(histToAv.GetBinContent(binx, biny))/float(histHits.GetBinContent(binx,biny))
        histToAv.SetBinContent(binx, biny, temp)

    def DoAverageHistogramDUT(self, histToAvDUT1, histHitsDUT1, histToAvDUT2, histHitsDUT2, histToAvDUT3, histHitsDUT3,
                              xbins, ybins):
        for i in xrange(1, xbins + 1):
            for j in xrange(1, ybins + 1):
                if histHitsDUT1.GetBinContent(i,j) >= 1: self.AverageBinHistogram(histToAvDUT1, histHitsDUT1, i, j)
                if histHitsDUT2.GetBinContent(i,j) >= 1: self.AverageBinHistogram(histToAvDUT2, histHitsDUT2, i, j)
                if histHitsDUT3.GetBinContent(i,j) >= 1: self.AverageBinHistogram(histToAvDUT3, histHitsDUT3, i, j)

    def DoAverageHistogramTPlanes(self, histToAvPlane0, histHitsPlane0, histToAvPlane1, histHitsPlane1, histToAvPlane2,
                                  histHitsPlane2, histToAvPlane3, histHitsPlane3, xbins, ybins):
        for i in xrange(1,xbins + 1):
            for j in xrange(1,ybins + 1):
                if histHitsPlane0.GetBinContent(i,j) >= 1: self.AverageBinHistogram(histToAvPlane0, histHitsPlane0, i, j)
                if histHitsPlane1.GetBinContent(i,j) >= 1: self.AverageBinHistogram(histToAvPlane1, histHitsPlane1, i, j)
                if histHitsPlane2.GetBinContent(i,j) >= 1: self.AverageBinHistogram(histToAvPlane2, histHitsPlane2, i, j)
                if histHitsPlane3.GetBinContent(i,j) >= 1: self.AverageBinHistogram(histToAvPlane3, histHitsPlane3, i, j)
