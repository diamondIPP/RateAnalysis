# ==============================================
# IMPORTS
# ==============================================
from ROOT import TGraphErrors, TCanvas, TH1D, TH2D, gStyle, TH1F, gROOT, gErrorIgnoreLevel, kError,TLegend, TCut, TGraph, TProfile2D, TH2F, TProfile, TCutG, kRed, kBlack, kPink, kBlue, kViolet ,kMagenta, kTeal, kGreen, kOrange, TF1, TPie, gPad, TLatex, THStack
# from TelescopeAnalysis import Analysis
# from CurrentInfo import Currents
# from numpy import array
from math import sqrt, ceil, log
from Elementary import Elementary
import os
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

class Plots(Elementary):
    def __init__(self, num_entries, run=None, num_devices=7, binning=-1, roc_tel=[0,1,2,3], roc_d1=4, roc_d2=5, roc_si=6):
        Elementary.__init__(self)
        # gStyle.SetPalette(53)  # kDarkBodyRadiator
        gStyle.SetPalette(55)  # kRainBow
        gStyle.SetNumberContours(999)
        self.run = run
        self.runinfo = self.run.RunInfo
        self.binning = binning
        self.num_devices = num_devices
        self.num_entries = num_entries
        self.plot_settings = {
            'ph1DbinsD4': 80,
            'ph1DminD4': 0,
            'ph1DmaxD4': 30000,
            'ph1DbinsD5': 80,
            'ph1DminD5': 0,
            'ph1DmaxD5': 60000,
            'ph1DbinsSi': 160,
            'ph1DminSi': 0,
            'ph1DmaxSi': 90000,
            'nEventsAv': 20000,
            'event_bins': max(int(ceil(float(self.num_entries)/10)), 200),
            'event_min': 0,
            'event_max': self.num_entries,
            'maxphplots': int(ceil(8*self.num_entries/100)),  ## for landau histograms histograms
            'nBinsX': 240,  #80 # 277 #240
            'xmin': -6,
            'xmax': 6,
            'nBinsY': 360,  #120 # 415 #360
            'ymin': -6,
            'ymax': 6,
            'nBinCol': 51,
            'minCol': 0,
            'maxCol': 51,
            'nBinRow': 79,
            'minRow': 0,
            'maxRow': 79,
            'num_diff_cluster_sizes': 4,
            'chi2_1Dbins': 60,
            'chi2_1Dmin': 0,
            'chi2_1Dmax': 30,
            'angle_1Dbins': 60,
            'angle_1Dmin': -3,
            'angle_1Dmax': 3,
            'rhit_1Dbins': 100,
            'rhit_1Dmin': 0,
            'rhit_1Dmax': 10
        }
        self.plot_settings['event_bins'] = int(ceil(float(self.num_entries)/5000)) if self.num_entries <= 100000 else \
            int(ceil(float(self.num_entries)/100)) if self.num_entries <= 500000 else int(ceil(float(self.num_entries)/self.plot_settings['nEventsAv']))
        self.plot_settings['deltaX'] = float(self.plot_settings['xmax']-self.plot_settings['xmin'])/self.plot_settings['nBinsX']
        self.plot_settings['deltaY'] = float(self.plot_settings['ymax']-self.plot_settings['ymin'])/self.plot_settings['nBinsY']
        self.roc_tel, self.roc_d1, self.roc_d2, self.roc_si = roc_tel, roc_d1, roc_d2, roc_si
        self.save_dir = './'

    def create_TGraphErrors(self, title='tgraph', xTitle='X', yTitle='Y', linecolor=kBlack, markercolor=kBlack):
        graph = TGraphErrors()
        graph.SetNameTitle(title,title)
        graph.SetLineColor(linecolor)
        graph.SetMarkerColor(markercolor)
        graph.GetXaxis().SetTitle(xTitle)
        graph.GetYaxis().SetTitle(yTitle)
        return (graph)

    def create_1D_histogram(self, type='landau', name='histo', title='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0, roc=4):
        if type is 'landau' or type is 'landaus':
            ph1Dbins = self.plot_settings['ph1DbinsD4'] if roc is self.roc_d1 else self.plot_settings['ph1DbinsD5'] if roc is self.roc_d2 else self.plot_settings['ph1DbinsSi']
            ph1Dmin = self.plot_settings['ph1DminD4'] if roc is self.roc_d1 else self.plot_settings['ph1DminD5'] if roc is self.roc_d2 else self.plot_settings['ph1DminSi']
            ph1Dmax = self.plot_settings['ph1DmaxD4'] if roc is self.roc_d1 else self.plot_settings['ph1DmaxD5'] if roc is self.roc_d2 else self.plot_settings['ph1DmaxSi']
        elif type is 'chi2':
            ph1Dbins = self.plot_settings['chi2_1Dbins']
            ph1Dmin = self.plot_settings['chi2_1Dmin']
            ph1Dmax = self.plot_settings['chi2_1Dmax']
        elif type is 'angle':
            ph1Dbins = self.plot_settings['angle_1Dbins']
            ph1Dmin = self.plot_settings['angle_1Dmin']
            ph1Dmax = self.plot_settings['angle_1Dmax']
        elif type is 'rhit':
            ph1Dbins = self.plot_settings['rhit_1Dbins']
            ph1Dmin = self.plot_settings['rhit_1Dmin']
            ph1Dmax = self.plot_settings['rhit_1Dmax']
        histo1D = TH1D(name, title, int(ph1Dbins + 1), ph1Dmin - float(ph1Dmax - ph1Dmin)/(2*ph1Dbins),
                       ph1Dmax + float(ph1Dmax - ph1Dmin)/(2*ph1Dbins))
        self.set_1D_options(type, histo1D, xTitle, yTitle, color, min_val)
        return (histo1D)

    def create_1D_profile(self, type='event', name='histo', title='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0, roc=4):
        nbins = int(ceil(float(self.plot_settings['event_max'] - self.plot_settings['event_min'])/self.plot_settings['nEventsAv']))
        xmin = self.plot_settings['event_min']
        xmax = self.plot_settings['event_max']
        histo1D = TProfile(name, title, int(nbins + 1), xmin - float(xmax - xmin)/(2*nbins), xmax + float(xmax - xmin)/(2*nbins))
        self.set_1D_options(type, histo1D, xTitle, yTitle, color, min_val, roc)
        return (histo1D)

    def set_1D_options(self, type='event', histo='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0, roc=4):
        histo.GetXaxis().SetTitle(xTitle)
        histo.GetYaxis().SetTitle(yTitle)
        histo.GetYaxis().SetTitleOffset(1.3)
        if type is 'event':
            if roc is self.roc_d1:
                histo.SetMaximum(self.plot_settings['ph1DmaxD4'])
            elif roc is self.roc_d2:
                histo.SetMaximum(self.plot_settings['ph1DmaxD5'])
            else:
                histo.SetMaximum(self.plot_settings['ph1DmaxSi'])
        histo.SetMinimum(min_val)
        histo.SetLineColor(color)
        histo.SetLineWidth(3*gStyle.GetLineWidth())
        if type is 'landaus':
            histo.SetFillColor(color)

    def create_2D_profile(self, type='spatial', name='histo', title='histo', xTitle='X', yTitle='Y', zTitle='Z', min_val=0, max_val=-1):
        xbins = self.plot_settings['nBinsX'] if type is 'spatial' else self.plot_settings['nBinCol']
        xmin = self.plot_settings['xmin'] if type is 'spatial' else self.plot_settings['minCol']
        xmax = self.plot_settings['xmax'] if type is 'spatial' else self.plot_settings['maxCol']
        ybins = self.plot_settings['nBinsY'] if type is 'spatial' else self.plot_settings['nBinRow']
        ymin = self.plot_settings['ymin'] if type is 'spatial' else self.plot_settings['minRow']
        ymax = self.plot_settings['ymax'] if type is 'spatial' else self.plot_settings['maxRow']
        histo2D = TProfile2D(name, title, int(xbins + 1), xmin - float(xmax-xmin)/(2*xbins),
                             xmax + float(xmax-xmin)/(2*xbins), int(ybins + 1), ymin - float(ymax-ymin)/(2*ybins),
                             ymax + float(ymax-ymin)/(2*ybins))
        self.set_2D_options(histo2D, xTitle, yTitle, zTitle, min_val, max_val)
        return histo2D

    def create_2D_histogram(self, type='spatial', name='histo', title='histo', xTitle='X', yTitle='Y', zTitle='Z', min_val=0, max_val=-1, roc=4):
        if type is 'spatial':
            xbins = self.plot_settings['nBinsX']
            xmin = self.plot_settings['xmin']
            xmax = self.plot_settings['xmax']
            ybins = self.plot_settings['nBinsY']
            ymin = self.plot_settings['ymin']
            ymax = self.plot_settings['ymax']
        elif type is 'pixel':
            xbins = self.plot_settings['nBinCol']
            xmin = self.plot_settings['minCol']
            xmax = self.plot_settings['maxCol']
            ybins = self.plot_settings['nBinRow']
            ymin = self.plot_settings['minRow']
            ymax = self.plot_settings['maxRow']
        elif type is 'correlpixcol':
            xbins = self.plot_settings['nBinCol']
            xmin = self.plot_settings['minCol']
            xmax = self.plot_settings['maxCol']
            ybins = self.plot_settings['nBinCol']
            ymin = self.plot_settings['minCol']
            ymax = self.plot_settings['maxCol']
        elif type is 'correlpixx':
            xbins = self.plot_settings['nBinsX']
            xmin = self.plot_settings['xmin']
            xmax = self.plot_settings['xmax']
            ybins = self.plot_settings['nBinsX']
            ymin = self.plot_settings['xmin']
            ymax = self.plot_settings['xmax']
        elif type is 'correlpixrow':
            xbins = self.plot_settings['nBinRow']
            xmin = self.plot_settings['minRow']
            xmax = self.plot_settings['maxRow']
            ybins = self.plot_settings['nBinRow']
            ymin = self.plot_settings['minRow']
            ymax = self.plot_settings['maxRow']
        elif type is 'correlpixy':
            xbins = self.plot_settings['nBinsY']
            xmin = self.plot_settings['ymin']
            xmax = self.plot_settings['ymax']
            ybins = self.plot_settings['nBinsY']
            ymin = self.plot_settings['ymin']
            ymax = self.plot_settings['ymax']
        else:
            xbins = self.plot_settings['event_bins']
            xmin = self.plot_settings['event_min']
            xmax = self.plot_settings['event_max']
            ybins = self.plot_settings['ph1DbinsD4'] if roc is self.roc_d1 else self.plot_settings['ph1DbinsD5'] if roc is self.roc_d2 else self.plot_settings['ph1DbinsSi']
            ymin = self.plot_settings['ph1DminD4'] if roc is self.roc_d1 else self.plot_settings['ph1DminD5'] if roc is self.roc_d2 else self.plot_settings['ph1DminSi']
            ymax = self.plot_settings['ph1DmaxD4'] if roc is self.roc_d1 else self.plot_settings['ph1DmaxD5'] if roc is self.roc_d2 else self.plot_settings['ph1DmaxSi']
        histo2D = TH2D(name, title, int(xbins + 1), xmin - float(xmax-xmin)/(2*xbins), xmax + float(xmax-xmin)/(2*xbins),
                       int(ybins + 1), ymin - float(ymax-ymin)/(2*ybins), ymax + float(ymax-ymin)/(2*ybins))
        self.set_2D_options(histo2D, xTitle, yTitle, zTitle, min_val, max_val)
        return histo2D

    def set_2D_options(self, histo, xTitle='X', yTitle='Y', zTitle='Z', min_val=0, max_val=-1):
        histo.GetXaxis().SetTitle(xTitle)
        histo.GetYaxis().SetTitle(yTitle)
        histo.GetZaxis().SetTitle(zTitle)
        histo.GetYaxis().SetTitleOffset(1.3)
        histo.GetZaxis().SetTitleOffset(1.4)
        histo.GetZaxis().CenterTitle(True)
        histo.SetMinimum(min_val)
        if max_val is not -1: histo.SetMaximum(max_val)

    def create_histograms(self, doTlscp=False):
        # 1D Histograms
        devini = 0 if doTlscp else self.roc_d1
        self.print_banner('Creating 1D histograms...')

        self.chi2_x = self.create_1D_histogram('chi2', 'chi2_x', 'Chi2 in X', 'chi2_x', 'Num Entries', kBlue)
        self.chi2_x_cut = self.create_1D_histogram('chi2', 'chi2_x_cut', 'Chi2 in X after cut', 'chi2_x', 'Num Entries', kRed)
        self.chi2_y = self.create_1D_histogram('chi2', 'chi2_y', 'Chi2 in Y', 'chi2_y', 'Num Entries', kBlue)
        self.chi2_y_cut = self.create_1D_histogram('chi2', 'chi2_y_cut', 'Chi2 in Y after cut', 'chi2_y', 'Num Entries', kRed)
        self.angle_x = self.create_1D_histogram('angle', 'angle_x', 'angle in X', 'angle_x', 'Num Entries', kBlue)
        self.angle_x_cut = self.create_1D_histogram('angle', 'angle_x_cut', 'angle in X after cut', 'angle_x', 'Num Entries', kRed)
        self.angle_y = self.create_1D_histogram('angle', 'angle_y', 'angle in Y', 'angle_y', 'Num Entries', kBlue)
        self.angle_y_cut = self.create_1D_histogram('angle', 'angle_y_cut', 'angle in Y after cut', 'angle_y', 'Num Entries', kRed)
        self.rhit = {i: self.create_1D_histogram('rhit', 'rhit_ROC{n}'.format(n=i), 'Cluster - track positions distance ROC {n}'.format(n=i),
                                                 'R Hit (mm)', 'Num Events', kBlue, roc=i) for i in xrange(devini, self.num_devices)}

        self.colors = [kBlack, kBlue, kRed, kOrange, kGreen, kMagenta, kViolet, kTeal]
        self.cuts = ['no', 'mask', 'fid', 'beam', 'tracks', 'angle', 'chi2', 'rhit']
        self.landaus0 = {}
        for j in xrange(len(self.cuts)):
            self.landaus0[self.cuts[j]] = {i: self.create_1D_histogram('landaus', 'phROC{n}_all_{cut}'.format(n=i, cut=self.cuts[j]),
                                                                       'Pulse Height ROC {n} all cluster sizes after {c} cut'.format(n=i, c=self.cuts[j]),
                                                                       'Charge (e)', 'Num Clusters', self.colors[j], 0.1, i) for i in xrange(devini, self.num_devices)}
        self.landaus1 = {}
        for j in xrange(len(self.cuts)):
            self.landaus1[self.cuts[j]] = {i: self.create_1D_histogram('landaus', 'phROC{n}_1cl_{cut}'.format(n=i, cut=self.cuts[j]),
                                                                       'Pulse Height ROC {n} 1 pix cluster size after {c} cut'.format(n=i, c=self.cuts[j]),
                                                                       'Charge (e)', 'Num Clusters', self.colors[j], 0.1, i) for i in xrange(devini, self.num_devices)}

        self.landaus2 = {}
        for j in xrange(len(self.cuts)):
            self.landaus2[self.cuts[j]] = {i: self.create_1D_histogram('landaus', 'phROC{n}_2cl_{cut}'.format(n=i, cut=self.cuts[j]),
                                                                       'Pulse Height ROC {n} 2 pix cluster size after {c} cut'.format(n=i, c=self.cuts[j]),
                                                                       'Charge (e)', 'Num Clusters', self.colors[j], 0.1, i) for i in xrange(devini, self.num_devices)}

        self.landaus3 = {}
        for j in xrange(len(self.cuts)):
            self.landaus3[self.cuts[j]] = {i: self.create_1D_histogram('landaus', 'phROC{n}_3cl_{cut}'.format(n=i, cut=self.cuts[j]),
                                                                       'Pulse Height ROC {n} 3 pix cluster size after {c} cut'.format(n=i, c=self.cuts[j]),
                                                                       'Charge (e)', 'Num Clusters', self.colors[j], 0.1, i) for i in xrange(devini, self.num_devices)}

        self.landausM4 = {}
        for j in xrange(len(self.cuts)):
            self.landausM4[self.cuts[j]] = {i: self.create_1D_histogram('landaus', 'phROC{n}_M4cl_{cut}'.format(n=i, cut=self.cuts[j]),
                                                                       'Pulse Height ROC {n} 4 or more pix cluster sizes after {c} cut'.format(n=i, c=self.cuts[j]),
                                                                       'Charge (e)', 'Num Clusters', self.colors[j], 0.1, i) for i in xrange(devini, self.num_devices)}

        self.rhit_cut = {i: self.create_1D_histogram('rhit', 'rhit_ROC{n}_cut'.format(n=i), 'Cluster - track positions distance ROC {n} after cut'.format(n=i),
                                                     'R Hit (mm)', 'Num Events', kRed, roc=i) for i in xrange(devini, self.num_devices)}

        self.phROC_all = {i: self.create_1D_histogram('landau', 'phROC{n}_all'.format(n=i),
                                                      'Pulse Height ROC {n} all cluster sizes'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kBlack, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_all_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_all_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} all cluster sizes after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kBlack, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_1cl = {i: self.create_1D_histogram('landau', 'phROC{n}_1cl'.format(n=i),
                                                      'Pulse Height ROC {n} 1pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kBlue, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_1cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_1cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 1pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kBlue, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_2cl = {i: self.create_1D_histogram('landau', 'phROC{n}_2cl'.format(n=i),
                                                      'Pulse Height ROC {n} 2pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kGreen, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_2cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_2cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 2pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kGreen, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_3cl = {i: self.create_1D_histogram('landau', 'phROC{n}_3cl'.format(n=i),
                                                      'Pulse Height ROC {n} 3pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kRed, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_3cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_3cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 3pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kRed, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_M4cl = {i: self.create_1D_histogram('landau', 'phROC{n}_M4cl'.format(n=i),
                                                       'Pulse Height ROC {n} 4 or more pix cluster'.format(n=i), 'Charge (e)',
                                                       'Num Clusters', kMagenta, roc=i) for i in xrange(devini, self.num_devices)}
        self.phROC_M4cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_M4cl_cuts'.format(n=i),
                                                            'Pulse Height ROC {n} 4 or more pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                            'Num Clusters', kMagenta, roc=i) for i in xrange(devini, self.num_devices)}
        
        self.smallChROC_all = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_all'.format(n=i), 'Smallest Charge ROC {n} all cluster sizes'.format(n=i),
                                                           'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_all_cuts = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_all_cuts'.format(n=i), 'Smallest Charge ROC {n} all cluster sizes after cuts'.format(n=i),
                                                           'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_1cl = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_1cl'.format(n=i), 'Smallest Charge ROC {n} 1pix cluster size'.format(n=i),
                                                           'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_1cl_cuts = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_1cl_cuts'.format(n=i), 'Smallest Charge ROC {n} 1pix cluster size after cuts'.format(n=i),
                                                           'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_2cl = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_2cl'.format(n=i), 'Smallest Charge ROC {n} 2pix cluster size'.format(n=i),
                                                           'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_2cl_cuts = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_2cl_cuts'.format(n=i), 'Smallest Charge ROC {n} 2pix cluster size after cuts'.format(n=i),
                                                           'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_3cl = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_3cl'.format(n=i), 'Smallest Charge ROC {n} 3pix cluster size'.format(n=i),
                                                           'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_3cl_cuts = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_3cl_cuts'.format(n=i), 'Smallest Charge ROC {n} 3pix cluster size after cuts'.format(n=i),
                                                                'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_M4cl = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_M4cl'.format(n=i), 'Smallest Charge ROC {n} 4 or mor pix clusters'.format(n=i),
                                                            'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.smallChROC_M4cl_cuts = {i: self.create_1D_histogram('landau', 'small_charge_ROC{n}_M4cl_cuts'.format(n=i), 'Smallest Charge ROC {n} 4 or mor pix clusters after cuts'.format(n=i),
                                                                 'Charge (e)', 'Num Clusters', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        
        self.print_banner('1D histograms creation -> Done')
        # 2D Histograms
        self.print_banner('Creating 2D histograms...')
        self.hitMap = {i: self.create_2D_histogram('pixel', 'hitMapROC{n}'.format(n=i), 'Hit Map ROC {n}'.format(n=i),
                                                   'Column', 'Row', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.hitMap_cuts = {i: self.create_2D_histogram('pixel', 'hitMapROC{n}_cuts'.format(n=i), 'Hit Map ROC {n} after cuts'.format(n=i),
                                                        'Column', 'Row', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.ph1cl_vs_event = {i: self.create_2D_histogram('event', 'phCl1VsEventROC{n}'.format(n=i),
                                                           'PH 1 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.ph1cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl1VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 1 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.ph2cl_vs_event = {i: self.create_2D_histogram('event', 'phCl2VsEventROC{n}'.format(n=i),
                                                           'PH 2 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.ph2cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl2VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 2 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.ph3cl_vs_event = {i: self.create_2D_histogram('event', 'phCl3VsEventROC{n}'.format(n=i),
                                                           'PH 3 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.ph3cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl3VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 3 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.phM4cl_vs_event = {i: self.create_2D_histogram('event', 'phClM4VsEventROC{n}'.format(n=i),
                                                            'PH 4 or more pix cluster Vs Event ROC {n}'.format(n=i),
                                                            'Event', 'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.phM4cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phClM4VsEventROC{n}_cuts'.format(n=i),
                                                                 'PH 4 or more pix cluster Vs Event ROC {n} after cuts'.format(n=i),
                                                                 'Event', 'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.phAll_vs_event = {i: self.create_2D_histogram('event', 'phAllVsEventROC{n}'.format(n=i),
                                                           'PH all cluster sizes Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}
        self.phAll_vs_event_cuts = {i: self.create_2D_histogram('event', 'phAllVsEventROC{n}_cuts'.format(n=i),
                                                                'PH all cluster sizes Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1, roc=i) for i in xrange(devini, self.num_devices)}

        self.correl_col, self.correl_row = {}, {}
        self.correl_col[self.roc_tel[1]], self.correl_col[self.roc_d1], self.correl_col[self.roc_si], self.correl_col[self.roc_tel[2]] = {}, {}, {}, {}
        self.correl_row[self.roc_tel[1]], self.correl_row[self.roc_d1], self.correl_row[self.roc_si], self.correl_row[self.roc_tel[2]] = {}, {}, {}, {}
        self.correl_col[self.roc_tel[1]][self.roc_d1] = self.create_2D_histogram('correlpixcol', 'corr_{rx}_{ry}_col'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'Col correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'col ROC {rx}'.format(rx=self.roc_tel[1]), 'col ROC {ry}'.format(ry=self.roc_d1), 'Entries', 0, -1)
        self.correl_col[self.roc_tel[2]][self.roc_si] = self.create_2D_histogram('correlpixcol', 'corr_{rx}_{ry}_col'.format(rx=self.roc_tel[2], ry=self.roc_si), 'Col correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[2], ry=self.roc_si), 'col ROC {rx}'.format(rx=self.roc_tel[2]), 'col ROC {ry}'.format(ry=self.roc_si), 'Entries', 0, -1)
        self.correl_row[self.roc_tel[1]][self.roc_d1] = self.create_2D_histogram('correlpixrow', 'corr_{rx}_{ry}_row'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'Row correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'row ROC {rx}'.format(rx=self.roc_tel[1]), 'row ROC {ry}'.format(ry=self.roc_d1), 'Entries', 0, -1)
        self.correl_row[self.roc_tel[2]][self.roc_si] = self.create_2D_histogram('correlpixrow', 'corr_{rx}_{ry}_row'.format(rx=self.roc_tel[2], ry=self.roc_si), 'Row correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[2], ry=self.roc_si), 'row ROC {rx}'.format(rx=self.roc_tel[2]), 'row ROC {ry}'.format(ry=self.roc_si), 'Entries', 0, -1)
        
        self.correl_x, self.correl_y = {}, {}
        self.correl_x[self.roc_tel[1]], self.correl_x[self.roc_d1], self.correl_x[self.roc_si], self.correl_x[self.roc_tel[2]] = {}, {}, {}, {}
        self.correl_y[self.roc_tel[1]], self.correl_y[self.roc_d1], self.correl_y[self.roc_si], self.correl_y[self.roc_tel[2]] = {}, {}, {}, {}
        self.correl_x[self.roc_tel[1]][self.roc_d1] = self.create_2D_histogram('correlpixx', 'corr_{rx}_{ry}_x'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'X correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'x ROC {rx}'.format(rx=self.roc_tel[1]), 'x ROC {ry}'.format(ry=self.roc_d1), 'Entries', 0, -1)
        self.correl_x[self.roc_tel[2]][self.roc_si] = self.create_2D_histogram('correlpixx', 'corr_{rx}_{ry}_x'.format(rx=self.roc_tel[2], ry=self.roc_si), 'X correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[2], ry=self.roc_si), 'x ROC {rx}'.format(rx=self.roc_tel[2]), 'x ROC {ry}'.format(ry=self.roc_si), 'Entries', 0, -1)
        self.correl_y[self.roc_tel[1]][self.roc_d1] = self.create_2D_histogram('correlpixy', 'corr_{rx}_{ry}_y'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'Y correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[1], ry=self.roc_d1), 'y ROC {rx}'.format(rx=self.roc_tel[1]), 'y ROC {ry}'.format(ry=self.roc_d1), 'Entries', 0, -1)
        self.correl_y[self.roc_tel[2]][self.roc_si] = self.create_2D_histogram('correlpixy', 'corr_{rx}_{ry}_y'.format(rx=self.roc_tel[2], ry=self.roc_si), 'Y correlation between ROCs {rx} and {ry}'.format(rx=self.roc_tel[2], ry=self.roc_si), 'y ROC {rx}'.format(rx=self.roc_tel[2]), 'y ROC {ry}'.format(ry=self.roc_si), 'Entries', 0, -1)

        self.print_banner('2D histograms creation -> Done')
        # alternative 2D
        self.print_banner('Creating 2D profiles...')
        self.avPhROC_local_all = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_local_all'.format(n=i),
                                                            'Average Pulse Height ROC {n} Local Coord. all cluster sizes'.format(n=i),
                                                            'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(devini, self.num_devices)}
        self.avPhROC_local_all_cuts = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_local_all_cuts'.format(n=i),
                                                                 'Average Pulse Height ROC {n} Local Coord. all cluster sizes after cuts'.format(n=i),
                                                                 'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(devini, self.num_devices)}
        self.avPhROC_telescope_all = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_telescope_all'.format(n=i),
                                                                'Average Pulse Height ROC {n} telescope Coord. all cluster sizes'.format(n=i),
                                                                'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(devini, self.num_devices)}
        self.avPhROC_telescope_all_cuts = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_telescope_all_cuts'.format(n=i),
                                                                     'Average Pulse Height ROC {n} telescope Coord. all cluster sizes after cuts'.format(n=i),
                                                                     'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(devini, self.num_devices)}
        self.avPhROC_pixelated_all = {i: self.create_2D_profile('pixel', 'avPh_ROC{n}_pixelated_all'.format(n=i),
                                                                'Average Pulse Height ROC {n} pixelated Coord. all cluster sizes'.format(n=i),
                                                                'Column', 'Row', 'Charge (e)', 0, -1) for i in xrange(devini, self.num_devices)}
        self.avPhROC_pixelated_all_cuts = {i: self.create_2D_profile('pixel', 'avPh_ROC{n}_pixelated_all_cuts'.format(n=i),
                                                                     'Average Pulse Height ROC {n} pixelated Coord. all cluster sizes after cuts'.format(n=i),
                                                                     'Column', 'Row', 'Charge (e)', 0, -1) for i in xrange(devini, self.num_devices)}
        self.print_banner('2D profiles creation -> Done')

        # Alternative to TGraphErrors
        self.print_banner('Creating 1D profiles...')
        self.meanPhROC_all = {i: self.create_1D_profile('event', 'meanPHROC{n}_all'.format(n=i),
                                                        'Mean PH ROC {n} all cluster sizes'.format(n=i), 'Event',
                                                        'Charge(e)', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_all_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_all_cuts'.format(n=i),
                                                             'Mean PH ROC {n} all cluster sizes after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kBlack, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_1cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_1cl'.format(n=i),
                                                        'Mean PH ROC {n} 1 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kBlue, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_1cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_1cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 1 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kBlue, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_2cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_2cl'.format(n=i),
                                                        'Mean PH ROC {n} 2 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kGreen, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_2cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_2cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 2 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kGreen, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_3cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_3cl'.format(n=i),
                                                        'Mean PH ROC {n} 3 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kRed, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_3cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_3cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 3 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kRed, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_M4cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_M4cl'.format(n=i),
                                                         'Mean PH ROC {n} 4 or more pixs cluster'.format(n=i), 'Event',
                                                         'Charge(e)', kMagenta, 0, i) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_M4cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_M4cl_cuts'.format(n=i),
                                                              'Mean PH ROC {n} 4 or more pixs cluster after cuts'.format(n=i), 'Event',
                                                              'Charge(e)', kMagenta, 0, i) for i in xrange(devini, self.num_devices)}
                
        self.print_banner('1D profiles creation -> Done')

    def save_individual_plots(self, histo, name, title, tcutg=None, draw_opt='', opt_stats=0, path='./', verbosity=False, opt_fit=0, addElem='', clone=False, doLogZ=False):
        if verbosity: self.print_banner('Saving {n}...'.format(n=name))
        gROOT.SetBatch(True)
        blabla = gROOT.ProcessLine("gErrorIgnoreLevel = {f};".format(f=kError))
        c0 = TCanvas('c_{n}'.format(n=name), title, 2100, 1500)
        c0.SetLeftMargin(0.1)
        c0.SetRightMargin(0.2)
        histo.SetStats(1)
        gStyle.SetOptFit(opt_fit)
        gStyle.SetOptStat(opt_stats)
        if addElem is not '':
            gStyle.SetStatX(0.4)
            gStyle.SetStatY(0.9)
            gStyle.SetStatW(0.15)
            gStyle.SetStatH(0.15)
        else:
            gStyle.SetStatX(0.8)
            gStyle.SetStatY(0.9)
            gStyle.SetStatW(0.15)
            gStyle.SetStatH(0.15)
        c0.cd()
        histo.Draw(draw_opt)
        if addElem is not '':
            c0.Update()
            st = c0.GetPrimitive('stats') if not clone else c0.GetPrimitive('stats')
            st.SetName('mystats') if not clone else st.SetName('mystats')
            lines = st.GetListOfLines()
            text = TLatex(0, 0, 'Correlation coef     {val}'.format(val=addElem))
            lines.Add(text)
            #st.AddText('Correlation coef     {val}'.format(val=addElem))
            histo.SetStats(0)
            st.Draw()
            c0.Modified()
        if tcutg is not None:
            tcutg.Draw('same')
        if not os.path.isdir('{dir}/Plots'.format(dir=path)):
            os.makedirs('{dir}/Plots'.format(dir=path))
        if not os.path.isdir('{dir}/Root'.format(dir=path)):
            os.makedirs('{dir}/Root'.format(dir=path))
        if doLogZ: c0.SetLogz()
        c0.SaveAs('{dir}/Root/c_{n}.root'.format(dir=path, n=name))
        c0.SaveAs('{dir}/Plots/c_{n}.png'.format(dir=path, n=name))
        c0.Close()
        gROOT.SetBatch(False)
        if verbosity: self.print_banner('{n} save -> Done'.format(n=name))
        del c0

    def save_cuts_distributions(self, histo1, histo2, name, title, draw_opt='', opt_stats=0, path='./', verbosity=False, histo3=''):
        if verbosity: self.print_banner('Saving {n}'.format(n=name))
        gROOT.SetBatch(True)
        blabla = gROOT.ProcessLine("gErrorIgnoreLevel = {f};".format(f=kError))
        c0 = TCanvas('c_{n}'.format(n=name), title, 2100, 1500)
        c0.SetLeftMargin(0.1)
        c0.SetRightMargin(0.1)
        histo1.SetStats(0)
        histo2.SetStats(1)
        gStyle.SetOptStat(opt_stats)
        gStyle.SetStatX(0.8)
        gStyle.SetStatY(0.9)
        gStyle.SetStatW(0.15)
        gStyle.SetStatH(0.15)
        c0.cd()
        histo1.Draw(draw_opt)
        histo2.Draw(draw_opt+'SAME')
        if histo3 != '':
            histo3.Draw(draw_opt+'SAME')
        c0.Update()
        c0.BuildLegend(0.65, 0.7, 0.9, 0.9)
        if not os.path.isdir('{dir}/Plots'.format(dir=path)):
            os.makedirs('{dir}/Plots'.format(dir=path))
        if not os.path.isdir('{dir}/Root'.format(dir=path)):
            os.makedirs('{dir}/Root'.format(dir=path))
        c0.SaveAs('{dir}/Root/c_{n}.root'.format(dir=path, n=name))
        c0.SaveAs('{dir}/Plots/c_{n}.png'.format(dir=path, n=name))
        c0.Close()
        gROOT.SetBatch(False)
        if verbosity: self.print_banner('{n} save -> Done'.format(n=name))
        del c0

    def save_cuts_overlay(self, histo0, histo1, histo2, histo3, histo4, histo5, histo6, histo7, name, title, draw_opt='', opt_stats=0, path='./', verbosity=False):
        if verbosity: self.print_banner('Saving {n}'.format(n=name))
        gROOT.SetBatch(True)
        blabla = gROOT.ProcessLine("gErrorIgnoreLevel = {f};".format(f=kError))
        c0 = TCanvas('c_{n}'.format(n=name), title, 2100, 1500)
        c0.SetLeftMargin(0.1)
        c0.SetRightMargin(0.1)
        histo1.SetStats(0)
        histo2.SetStats(1)
        gStyle.SetOptStat(opt_stats)
        gStyle.SetStatX(0.8)
        gStyle.SetStatY(0.9)
        gStyle.SetStatW(0.15)
        gStyle.SetStatH(0.15)
        c0.cd()
        s1 = THStack('s_{n}'.format(n=name), 's_{n}'.format(n=name))
        s1.Add(histo0)
        s1.Add(histo1)
        s1.Add(histo2)
        s1.Add(histo3)
        s1.Add(histo4)
        s1.Add(histo5)
        s1.Add(histo6)
        s1.Add(histo7)
        # histo0.Draw(draw_opt)
        # histo1.Draw(draw_opt+'SAME')
        # histo2.Draw(draw_opt+'SAME')
        # histo3.Draw(draw_opt+'SAME')
        # histo4.Draw(draw_opt+'SAME')
        # histo5.Draw(draw_opt+'SAME')
        # histo6.Draw(draw_opt+'SAME')
        # histo7.Draw(draw_opt+'SAME')
        s1.Draw('nostack')
        c0.Update()
        c0.SetLogy()
        c0.BuildLegend(0.65, 0.7, 0.9, 0.9)
        if not os.path.isdir('{dir}/Plots'.format(dir=path)):
            os.makedirs('{dir}/Plots'.format(dir=path))
        if not os.path.isdir('{dir}/Root'.format(dir=path)):
            os.makedirs('{dir}/Root'.format(dir=path))
        c0.SaveAs('{dir}/Root/c_{n}.root'.format(dir=path, n=name))
        c0.SaveAs('{dir}/Plots/c_{n}.png'.format(dir=path, n=name))
        c0.Close()
        gROOT.SetBatch(False)
        if verbosity: self.print_banner('{n} save -> Done'.format(n=name))
        del c0

    def clone_correlation_histograms(self, rocx, rocy, verbosity=False, optfit=0, extra1='', extra2='', extra3='', extra4=''):
        self.correl_col[rocy][rocx] = self.correl_col[rocx][rocy].Clone('corr_{ry}_{rx}_col'.format(ry=rocy, rx=rocx))
        self.correl_col[rocy][rocx].SetTitle('Col correlation between ROCs {ry} and {rx}'.format(ry=rocy, rx=rocx))
        self.save_individual_plots(self.correl_col[rocy][rocx], self.correl_col[rocy][rocx].GetName(), self.correl_col[rocy][rocx].GetTitle(), None, 'colz', 1000000011, self.save_dir, verbosity, optfit, extra1, True)
        self.correl_row[rocy][rocx] = self.correl_row[rocx][rocy].Clone('corr_{ry}_{rx}_row'.format(ry=rocy, rx=rocx))
        self.correl_row[rocy][rocx].SetTitle('Row correlation between ROCs {ry} and {rx}'.format(ry=rocy, rx=rocx))
        self.save_individual_plots(self.correl_row[rocy][rocx], self.correl_row[rocy][rocx].GetName(), self.correl_row[rocy][rocx].GetTitle(), None, 'colz', 1000000011, self.save_dir, verbosity, optfit, extra2, True)
        
        self.correl_x[rocy][rocx] = self.correl_x[rocx][rocy].Clone('corr_{ry}_{rx}_x'.format(ry=rocy, rx=rocx))
        self.correl_x[rocy][rocx].SetTitle('X correlation between ROCs {ry} and {rx}'.format(ry=rocy, rx=rocx))
        self.save_individual_plots(self.correl_x[rocy][rocx], self.correl_x[rocy][rocx].GetName(), self.correl_x[rocy][rocx].GetTitle(), None, 'colz', 1000000011, self.save_dir, verbosity, optfit, extra3, True)
        self.correl_y[rocy][rocx] = self.correl_y[rocx][rocy].Clone('corr_{ry}_{rx}_y'.format(ry=rocy, rx=rocx))
        self.correl_y[rocy][rocx].SetTitle('Y correlation between ROCs {ry} and {rx}'.format(ry=rocy, rx=rocx))
        self.save_individual_plots(self.correl_y[rocy][rocx], self.correl_y[rocy][rocx].GetName(), self.correl_y[rocy][rocx].GetTitle(), None, 'colz', 1000000011, self.save_dir, verbosity, optfit, extra4, True)

    def check_plot_existence(self, path, name):
        if os.path.isfile('{p}/Plots/{n}.png'.format(p=path, n=name)):
            return True
        else:
            return False
