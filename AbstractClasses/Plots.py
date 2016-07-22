# ==============================================
# IMPORTS
# ==============================================
from ROOT import TGraphErrors, TCanvas, TH1D, TH2D, gStyle, TH1F, gROOT, TLegend, TCut, TGraph, TProfile2D, TH2F, TProfile, TCutG, kRed, kBlack, kBlue, kMagenta, kGreen, kOrange, TF1, TPie
# from TelescopeAnalysis import Analysis
# from CurrentInfo import Currents
# from numpy import array
from math import sqrt, ceil, log
from Elementary import Elementary
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
    def __init__(self, num_entries, run=None, num_devices=7, binning=-1):
        Elementary.__init__(self)
        gStyle.SetPalette(53)
        gStyle.SetNumberContours(999)
        self.run = run
        self.runinfo = self.run.RunInfo
        self.binning = binning
        self.num_devices = num_devices
        self.num_entries = num_entries
        self.plot_settings = {
            'ph1Dbins': 100,
            'ph1Dmin': 0,
            'ph1Dmax': 50000,
            'nEventsAv': 1000,
            'event_bins': int(ceil(float(self.num_entries)/10)),
            'event_min': 0,
            'event_max': self.num_entries,
            'maxphplots': int(ceil(8*self.num_entries/100)),  ## for landau histograms histograms
            'nBinsX': 80,  # 277
            'xmin': -6,
            'xmax': 6,
            'nBinsY': 120,  # 415
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
        self.plot_settings['deltaX'] = float(self.plot_settings['xmax']-self.plot_settings['xmin'])/self.plot_settings['nBinsX']
        self.plot_settings['deltaY'] = float(self.plot_settings['ymax']-self.plot_settings['ymin'])/self.plot_settings['nBinsY']

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
        if type is 'event': histo.SetMaximum(self.plot_settings['ph1Dmax'])
        histo.SetMinimum(min_val)
        histo.SetLineColor(color)

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
        else:
            xbins =  self.plot_settings['event_bins']
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
        histo.SetMinimum(min_val)
        if max_val is not -1: histo.SetMaximum(max_val)

    def create_histograms(self):
        # 1D Histograms
        self.print_banner('Creating 1D histograms...')
        self.phROC_all = {i: self.create_1D_histogram('landau', 'phROC{n}_all'.format(n=i),
                                                      'Pulse Height ROC {n} all cluster sizes'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kBlack) for i in xrange(self.num_devices)}
        self.phROC_all_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_all_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} all cluster sizes after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kBlack) for i in xrange(self.num_devices)}
        self.phROC_1cl = {i: self.create_1D_histogram('landau', 'phROC{n}_1cl'.format(n=i),
                                                      'Pulse Height ROC {n} 1pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kBlue) for i in xrange(self.num_devices)}
        self.phROC_1cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_1cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 1pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kBlue) for i in xrange(self.num_devices)}
        self.phROC_2cl = {i: self.create_1D_histogram('landau', 'phROC{n}_2cl'.format(n=i),
                                                      'Pulse Height ROC {n} 2pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kGreen) for i in xrange(self.num_devices)}
        self.phROC_2cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_2cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 2pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kGreen) for i in xrange(self.num_devices)}
        self.phROC_3cl = {i: self.create_1D_histogram('landau', 'phROC{n}_3cl'.format(n=i),
                                                      'Pulse Height ROC {n} 3pix cluster'.format(n=i), 'Charge (e)',
                                                      'Num Clusters', kRed) for i in xrange(self.num_devices)}
        self.phROC_3cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_3cl_cuts'.format(n=i),
                                                           'Pulse Height ROC {n} 3pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                           'Num Clusters', kRed) for i in xrange(self.num_devices)}
        self.phROC_M4cl = {i: self.create_1D_histogram('landau', 'phROC{n}_M4cl'.format(n=i),
                                                       'Pulse Height ROC {n} 4 or more pix cluster'.format(n=i), 'Charge (e)',
                                                       'Num Clusters', kMagenta) for i in xrange(self.num_devices)}
        self.phROC_M4cl_cuts = {i: self.create_1D_histogram('landau', 'phROC{n}_M4cl_cuts'.format(n=i),
                                                            'Pulse Height ROC {n} 4 or more pix cluster after cuts'.format(n=i), 'Charge (e)',
                                                            'Num Clusters', kMagenta) for i in xrange(self.num_devices)}
        self.print_banner('1D histograms creation -> Done')
        # 2D Histograms
        self.print_banner('Creating 2D histograms...')
        self.hitMap = {i: self.create_2D_histogram('pixel', 'hitMapROC{n}'.format(n=i), 'Hit Map ROC {n}'.format(n=i),
                                                   'Column', 'Row', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.hitMap_cuts = {i: self.create_2D_histogram('pixel', 'hitMapROC{n}_cuts'.format(n=i), 'Hit Map ROC {n} after cuts'.format(n=i),
                                                        'Column', 'Row', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.ph1cl_vs_event = {i: self.create_2D_histogram('event', 'phCl1VsEventROC{n}'.format(n=i),
                                                           'PH 1 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.ph1cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl1VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 1 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.ph2cl_vs_event = {i: self.create_2D_histogram('event', 'phCl2VsEventROC{n}'.format(n=i),
                                                           'PH 2 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.ph2cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl2VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 2 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.ph3cl_vs_event = {i: self.create_2D_histogram('event', 'phCl3VsEventROC{n}'.format(n=i),
                                                           'PH 3 pix cluster Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.ph3cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phCl3VsEventROC{n}_cuts'.format(n=i),
                                                                'PH 3 pix cluster Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.phM4cl_vs_event = {i: self.create_2D_histogram('event', 'phClM4VsEventROC{n}'.format(n=i),
                                                            'PH 4 or more pix cluster Vs Event ROC {n}'.format(n=i),
                                                            'Event', 'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.phM4cl_vs_event_cuts = {i: self.create_2D_histogram('event', 'phClM4VsEventROC{n}_cuts'.format(n=i),
                                                                 'PH 4 or more pix cluster Vs Event ROC {n} after cuts'.format(n=i),
                                                                 'Event', 'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.phAll_vs_event = {i: self.create_2D_histogram('event', 'phAllVsEventROC{n}'.format(n=i),
                                                           'PH all cluster sizes Vs Event ROC {n}'.format(n=i), 'Event',
                                                           'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.phAll_vs_event_cuts = {i: self.create_2D_histogram('event', 'phAllVsEventROC{n}_cuts'.format(n=i),
                                                                'PH all cluster sizes Vs Event ROC {n} after cuts'.format(n=i), 'Event',
                                                                'Charge (e)', 'Entries', 0, -1) for i in xrange(self.num_devices)}
        self.print_banner('2D histograms creation -> Done')
        # alternative 2D
        self.print_banner('Creating 2D profiles...')
        self.avPhROC_local_all = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_local_all'.format(n=i),
                                                            'Average Pulse Height ROC {n} Local Coord. all cluster sizes'.format(n=i),
                                                            'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(self.num_devices)}
        self.avPhROC_local_all_cuts = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_local_all_cuts'.format(n=i),
                                                                 'Average Pulse Height ROC {n} Local Coord. all cluster sizes after cuts'.format(n=i),
                                                                 'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(self.num_devices)}
        self.avPhROC_telescope_all = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_telescope_all'.format(n=i),
                                                                'Average Pulse Height ROC {n} telescope Coord. all cluster sizes'.format(n=i),
                                                                'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(self.num_devices)}
        self.avPhROC_telescope_all_cuts = {i: self.create_2D_profile('spatial', 'avPh_ROC{n}_telescope_all_cuts'.format(n=i),
                                                                     'Average Pulse Height ROC {n} telescope Coord. all cluster sizes after cuts'.format(n=i),
                                                                     'x (mm)', 'y (mm)', 'Charge (e)', 0, -1) for i in xrange(self.num_devices)}
        self.avPhROC_pixelated_all = {i: self.create_2D_profile('pixel', 'avPh_ROC{n}_pixelated_all'.format(n=i),
                                                                'Average Pulse Height ROC {n} pixelated Coord. all cluster sizes'.format(n=i),
                                                                'Column', 'Row', 'Charge (e)', 0, -1) for i in xrange(self.num_devices)}
        self.avPhROC_pixelated_all_cuts = {i: self.create_2D_profile('pixel', 'avPh_ROC{n}_pixelated_all_cuts'.format(n=i),
                                                                     'Average Pulse Height ROC {n} pixelated Coord. all cluster sizes after cuts'.format(n=i),
                                                                     'Column', 'Row', 'Charge (e)', 0, -1) for i in xrange(self.num_devices)}
        self.print_banner('2D profiles creation -> Done')

        # Alternative to TGraphErrors
        self.print_banner('Creating 1D profiles...')
        self.meanPhROC_all = {i: self.create_1D_profile('event', 'meanPHROC{n}_all'.format(n=i),
                                                        'Mean PH ROC {n} all cluster sizes'.format(n=i), 'Event',
                                                        'Charge(e)', kBlack, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_all_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_all_cuts'.format(n=i),
                                                             'Mean PH ROC {n} all cluster sizes after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kBlack, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_1cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_1cl'.format(n=i),
                                                        'Mean PH ROC {n} 1 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kBlue, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_1cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_1cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 1 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kBlue, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_2cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_2cl'.format(n=i),
                                                        'Mean PH ROC {n} 2 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kGreen, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_2cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_2cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 2 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kGreen, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_3cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_3cl'.format(n=i),
                                                        'Mean PH ROC {n} 3 pix cluster'.format(n=i), 'Event',
                                                        'Charge(e)', kRed, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_3cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_3cl_cuts'.format(n=i),
                                                             'Mean PH ROC {n} 3 pix cluster after cuts'.format(n=i), 'Event',
                                                             'Charge(e)', kRed, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_M4cl = {i: self.create_1D_profile('event', 'meanPHROC{n}_M4cl'.format(n=i),
                                                         'Mean PH ROC {n} 4 or more pixs cluster'.format(n=i), 'Event',
                                                         'Charge(e)', kMagenta, 0) for i in xrange(self.num_devices)}
        self.meanPhROC_M4cl_cuts = {i: self.create_1D_profile('event', 'meanPHROC{n}_M4cl_cuts'.format(n=i),
                                                              'Mean PH ROC {n} 4 or more pixs cluster after cuts'.format(n=i), 'Event',
                                                              'Charge(e)', kMagenta, 0) for i in xrange(self.num_devices)}
        self.print_banner('1D profiles creation -> Done')

    def save_individual_plots(self, histo, name, title, tcutg=None, draw_opt='', opt_stats=0, path='./'):
        self.print_banner('Saving {n}...'.format(n=name))
        c0 = TCanvas('c_{n}'.format(n=name), title)
        gStyle.SetOptStat(opt_stats)
        c0.cd()
        histo.Draw(draw_opt)
        if tcutg is not None:
            tcutg.Draw('same')
        c0.SaveAs('{dir}/c_{n}.root'.format(dir=path, n=name))
        c0.SaveAs('{dir}/c_{n}.png'.format(dir=path, n=name))
        c0.Close()
        self.print_banner('{n} save -> Done'.format(n=name))
        del c0
