# ==============================================
# IMPORTS
# ==============================================
from ROOT import TGraphErrors, TCanvas, TH1D, TH2D, gStyle, TH1F, gROOT, TLegend, TCut, TGraph, TProfile2D, TH2F, TProfile, TCutG, kRed, kBlack, kBlue, kMagenta, kGreen, kOrange, TF1, TPie, gPad, TLatex
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
        gStyle.SetPalette(53)
        gStyle.SetNumberContours(999)
        self.run = run
        self.runinfo = self.run.RunInfo
        self.binning = binning
        self.num_devices = num_devices
        self.num_entries = num_entries
        self.plot_settings = {
            'ph1Dbins': 200,
            'ph1Dmin': 0,
            'ph1Dmax': 100000,
            'nEventsAv': 1000,
            'event_bins': int(ceil(float(self.num_entries)/10)),
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
            'num_diff_cluster_sizes': 4
        }
        self.plot_settings['event_bins'] = int(ceil(float(self.num_entries)/10)) if self.num_entries <= 100000 else \
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

    def create_1D_histogram(self, type='landau', name='histo', title='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0):
        ph1Dbins = self.plot_settings['ph1Dbins']
        ph1Dmin = self.plot_settings['ph1Dmin']
        ph1Dmax = self.plot_settings['ph1Dmax']
        histo1D = TH1D(name, title, int(ph1Dbins + 1), ph1Dmin - float(ph1Dmax - ph1Dmin)/(2*ph1Dbins),
                       ph1Dmax + float(ph1Dmax - ph1Dmin)/(2*ph1Dbins))
        self.set_1D_options(type, histo1D, xTitle, yTitle, color, min_val)
        return (histo1D)

    def create_1D_profile(self, type='event', name='histo', title='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0):
        nbins = int(ceil(float(self.plot_settings['event_max'] - self.plot_settings['event_min'])/self.plot_settings['nEventsAv']))
        xmin = self.plot_settings['event_min']
        xmax = self.plot_settings['event_max']
        histo1D = TProfile(name, title, int(nbins + 1), xmin - float(xmax - xmin)/(2*nbins), xmax + float(xmax - xmin)/(2*nbins))
        self.set_1D_options(type, histo1D, xTitle, yTitle, color, min_val)
        return (histo1D)

    def set_1D_options(self, type='event', histo='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0):
        histo.GetXaxis().SetTitle(xTitle)
        histo.GetYaxis().SetTitle(yTitle)
        histo.GetYaxis().SetTitleOffset(1.3)
        if type is 'event': histo.SetMaximum(self.plot_settings['ph1Dmax'])
        histo.SetMinimum(min_val)
        histo.SetLineColor(color)
        histo.SetLineWidth(3*gStyle.GetLineWidth())

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

    def create_2D_histogram(self, type='spatial', name='histo', title='histo', xTitle='X', yTitle='Y', zTitle='Z', min_val=0, max_val=-1):
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
            ybins = self.plot_settings['ph1Dbins']
            ymin = self.plot_settings['ph1Dmin']
            ymax = self.plot_settings['ph1Dmax']
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
        devini = 0 if doTlscp else 4
        self.print_banner('Creating 1D histograms...')
        self.phROC_all = {i: self.create_1D_histogram('landau', 'phROC{n}_all'.format(n=i),
                                                      'Pulse Height ROC {n} all cluster sizes'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kBlack) for i in xrange(devini, self.num_devices)}
        self.phROC_all_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_all_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} all cluster sizes after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kBlack) for i in xrange(devini, self.num_devices)}
        self.phROC_1cl = {i: self.create_1D_histogram('landau', 'phROC{n}_1cl'.format(n=i),
                                                      'Pulse Height ROC {n} 1pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kBlue) for i in xrange(devini, self.num_devices)}
        self.phROC_1cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_1cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 1pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kBlue) for i in xrange(devini, self.num_devices)}
        self.phROC_2cl = {i: self.create_1D_histogram('landau', 'phROC{n}_2cl'.format(n=i),
                                                      'Pulse Height ROC {n} 2pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kGreen) for i in xrange(devini, self.num_devices)}
        self.phROC_2cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_2cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 2pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kGreen) for i in xrange(devini, self.num_devices)}
        self.phROC_3cl = {i: self.create_1D_histogram('landau', 'phROC{n}_3cl'.format(n=i),
                                                      'Pulse Height ROC {n} 3pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kRed) for i in xrange(devini, self.num_devices)}
        self.phROC_3cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_3cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 3pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kRed) for i in xrange(devini, self.num_devices)}
        self.phROC_M4cl = {i: self.create_1D_histogram('landau', 'phROC{n}_M4cl'.format(n=i),
                                                       'Pulse Height ROC {n} 4 or more pix cluster'.format(n=i), 'Charge (e)',
                                                       'Num Clusters', kMagenta) for i in xrange(devini, self.num_devices)}
        self.phROC_M4cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_M4cl_cuts'.format(n=i),
                                                            'Pulse Height ROC {n} 4 or more pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                            'Num Clusters', kMagenta) for i in xrange(devini, self.num_devices)}
        self.print_banner('1D histograms creation -> Done')
        # 2D Histograms
        self.print_banner('Creating 2D histograms...')
        self.hitMap = {i: self.create_2D_histogram('pixel', 'hitMapROC{n}'.format(n=i), 'Hit Map ROC {n}'.format(n=i),
                                                   'Column', 'Row', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.hitMap_cuts = {i: self.create_2D_histogram('pixel', 'hitMapROC{n}_cuts'.format(n=i), 'Hit Map ROC {n} after cuts'.format(n=i),
                                                        'Column', 'Row', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.ph1cl_vs_event = {i: self.create_2D_histogram('event', 'phCl1VsEventROC{n}'.format(n=i),
                                                           'PH 1 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.ph1cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl1VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 1 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.ph2cl_vs_event = {i: self.create_2D_histogram('event', 'phCl2VsEventROC{n}'.format(n=i),
                                                           'PH 2 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.ph2cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl2VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 2 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.ph3cl_vs_event = {i: self.create_2D_histogram('event', 'phCl3VsEventROC{n}'.format(n=i),
                                                           'PH 3 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.ph3cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl3VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 3 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.phM4cl_vs_event = {i: self.create_2D_histogram('event', 'phClM4VsEventROC{n}'.format(n=i),
                                                            'PH 4 or more pix cluster Vs Event ROC {n}'.format(n=i),
                                                            'Event', 'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.phM4cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phClM4VsEventROC{n}_cuts'.format(n=i),
                                                                 'PH 4 or more pix cluster Vs Event ROC {n} after cuts'.format(n=i),
                                                                 'Event', 'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.phAll_vs_event = {i: self.create_2D_histogram('event', 'phAllVsEventROC{n}'.format(n=i),
                                                           'PH all cluster sizes Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}
        self.phAll_vs_event_cuts = {i: self.create_2D_histogram('event', 'phAllVsEventROC{n}_cuts'.format(n=i),
                                                                'PH all cluster sizes Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(devini, self.num_devices)}

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
                                                        'Charge(e)', kBlack, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_all_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_all_cuts'.format(n=i),
                                                             'Mean PH ROC {n} all cluster sizes after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kBlack, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_1cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_1cl'.format(n=i),
                                                        'Mean PH ROC {n} 1 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kBlue, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_1cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_1cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 1 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kBlue, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_2cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_2cl'.format(n=i),
                                                        'Mean PH ROC {n} 2 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kGreen, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_2cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_2cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 2 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kGreen, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_3cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_3cl'.format(n=i),
                                                        'Mean PH ROC {n} 3 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kRed, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_3cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_3cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 3 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kRed, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_M4cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_M4cl'.format(n=i),
                                                         'Mean PH ROC {n} 4 or more pixs cluster'.format(n=i), 'Event',
                                                         'Charge(e)', kMagenta, 0) for i in xrange(devini, self.num_devices)}
        self.meanPhROC_M4cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_M4cl_cuts'.format(n=i),
                                                              'Mean PH ROC {n} 4 or more pixs cluster after cuts'.format(n=i), 'Event',
                                                              'Charge(e)', kMagenta, 0) for i in xrange(devini, self.num_devices)}
        self.print_banner('1D profiles creation -> Done')

    def save_individual_plots(self, histo, name, title, tcutg=None, draw_opt='', opt_stats=0, path='./', verbosity=False, opt_fit=0, addElem='', clone=False):
        if verbosity: self.print_banner('Saving {n}...'.format(n=name))
        gROOT.SetBatch(True)
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
