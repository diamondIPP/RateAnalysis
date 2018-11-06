# ==============================================
# IMPORTS
# ==============================================
import os
from ROOT import TGraphErrors, TCanvas, TH1D, TH2D, gStyle, gROOT, kError, TProfile2D, TProfile, kBlack, TLatex, THStack
from math import ceil, sqrt
from numpy import arange, array
from Utils import increased_range, print_banner

from Elementary import Elementary

__author__ = 'DA'


# ==============================================
# MAIN CLASS
# ==============================================

class Plots(Elementary):
    def __init__(self, run=None, binning=-1, roc_tel=None):
        Elementary.__init__(self)
        # gStyle.SetPalette(53)  # kDarkBodyRadiator
        gStyle.SetPalette(55)  # kRainBow
        gStyle.SetNumberContours(255)
        self.run = run
        self.runinfo = self.run.RunInfo
        self.binning = binning
        self.num_devices = run.NPlanes
        self.NEntries = run.n_entries
        self.NCols = 52
        self.NRows = 80
        self.Settings = {'ph1Dmin': -5000,
                         'ph1Dmax': 60000,
                         'ph1Dbins': 65000 / 500,
                         'ph1DbinsSi': 160,
                         'ph1DminSi': 0,
                         'ph1DmaxSi': 90000,
                         'nEventsAv': 20000,
                         'event_bins': max(int(ceil(float(self.NEntries) / 10)), 200),
                         'event_min': 0,
                         'event_max': self.NEntries,
                         'maxphplots': int(ceil(8 * self.NEntries / 100)),
                         'nBinsX': 52,
                         'nBinsY': 80,
                         'xmin': -3900,
                         'xmax': 3900,
                         'ymin': -4000,
                         'ymax': 4000,
                         'nBinCol': 51,
                         'nCols': 52,
                         'nRows': 80,
                         'minCol': 0,
                         'maxCol': 51,
                         'nBinRow': 79,
                         'minRow': 0,
                         'maxRow': 79,
                         'num_diff_cluster_sizes': 4,
                         'vcalBins': [int((1350 * 47) / 500), -100, 1250],
                         'globalCoods': [-.5025, .5175, -.505, .515],
                         'xpix': .015,
                         'ypix': .01}
        self.Settings['phBins'] = [self.Settings['ph1Dbins'], self.Settings['ph1Dmin'], self.Settings['ph1Dmax']]
        self.Settings['2DBins'] = [self.Settings['nCols'], - .5, self.Settings['nCols'] - .5, self.Settings['nRows'], - .5, self.Settings['nRows'] - .5]
        self.Settings['2DBinsX'] = [self.Settings['nCols'], - .5, self.Settings['nCols'] - .5] * 2
        self.Settings['2DBinsY'] = [self.Settings['nRows'], - .5, self.Settings['nRows'] - .5] * 2
        self.Settings['event_bins'] = int(ceil(float(self.NEntries) / 5000)) if self.NEntries <= 100000 else \
            int(ceil(float(self.NEntries) / 100)) if self.NEntries <= 500000 else int(ceil(float(self.NEntries) / self.Settings['nEventsAv']))
        self.Settings['deltaX'] = float(self.Settings['xmax'] - self.Settings['xmin']) / self.Settings['nBinsX']
        self.Settings['deltaY'] = float(self.Settings['ymax'] - self.Settings['ymin']) / self.Settings['nBinsY']
        self.RocTel = range(4) if roc_tel is None else roc_tel
        self.save_dir = './'

    @staticmethod
    def get_arrays(lst):
        arrays = []
        for i in xrange(len(lst) / 3):
            step = (lst[3 * i + 2] - lst[3 * i + 1]) / float(lst[3 * i])
            arrays.append(arange(lst[3 * i + 1], lst[3 * i + 2] + step, step, 'd'))
        ret_val = []
        for arr in arrays:
            ret_val.append(len(arr) - 1)
            ret_val.append(arr)
        return ret_val

    def get_tcal_bins(self):
        return [int(round(sqrt(len(self.run.TCal))))] + increased_range([min(self.run.TCal), max(self.run.TCal)], .1, .1)

    def get_global_bins(self, res, mode=None, arrays=False):
        x, y = self.Settings['globalCoods'][:2], self.Settings['globalCoods'][2:]
        size_x = self.Settings['xpix'] if mode in ['x', None] else self.Settings['ypix']
        size_y = self.Settings['ypix'] if mode in ['y', None] else self.Settings['xpix']
        bins = [int(ceil((x[1] - x[0]) / size_x * sqrt(12) / res)), x[0], x[1], int(ceil((y[1] - y[0]) / size_y * sqrt(12) / res)), y[0], y[1]]
        x_arr, y_arr = arange(x[0], x[1], size_x, 'd'), arange(y[0], y[1], size_y, 'd')
        return bins if not arrays else [len(x_arr) - 1, x_arr, len(y_arr) - 1, y_arr]

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
            ph1Dbins = self.Settings['ph1DbinsD4'] if roc is self.roc_d1 else self.Settings['ph1DbinsD5'] if roc is self.roc_d2 else self.Settings['ph1DbinsSi']
            ph1Dmin = self.Settings['ph1DminD4'] if roc is self.roc_d1 else self.Settings['ph1DminD5'] if roc is self.roc_d2 else self.Settings['ph1DminSi']
            ph1Dmax = self.Settings['ph1DmaxD4'] if roc is self.roc_d1 else self.Settings['ph1DmaxD5'] if roc is self.roc_d2 else self.Settings['ph1DmaxSi']
        elif type is 'chi2':
            ph1Dbins = self.Settings['chi2_1Dbins']
            ph1Dmin = self.Settings['chi2_1Dmin']
            ph1Dmax = self.Settings['chi2_1Dmax']
        elif type is 'angle':
            ph1Dbins = self.Settings['angle_1Dbins']
            ph1Dmin = self.Settings['angle_1Dmin']
            ph1Dmax = self.Settings['angle_1Dmax']
        elif type is 'rhit':
            ph1Dbins = self.Settings['rhit_1Dbins']
            ph1Dmin = self.Settings['rhit_1Dmin']
            ph1Dmax = self.Settings['rhit_1Dmax']
        histo1D = TH1D(name, title, int(ph1Dbins + 1), ph1Dmin - float(ph1Dmax - ph1Dmin)/(2*ph1Dbins),
                       ph1Dmax + float(ph1Dmax - ph1Dmin)/(2*ph1Dbins))
        self.set_1D_options(type, histo1D, xTitle, yTitle, color, min_val)
        return (histo1D)

    def create_1D_profile(self, type='event', name='histo', title='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0, roc=4):
        nbins = int(ceil(float(self.Settings['event_max'] - self.Settings['event_min']) / self.Settings['nEventsAv']))
        xmin = self.Settings['event_min']
        xmax = self.Settings['event_max']
        histo1D = TProfile(name, title, int(nbins + 1), xmin - float(xmax - xmin)/(2*nbins), xmax + float(xmax - xmin)/(2*nbins))
        self.set_1D_options(type, histo1D, xTitle, yTitle, color, min_val, roc)
        return (histo1D)

    def set_1D_options(self, type='event', histo='histo', xTitle='X', yTitle='Y', color=kBlack, min_val=0, roc=4):
        histo.GetXaxis().SetTitle(xTitle)
        histo.GetYaxis().SetTitle(yTitle)
        histo.GetYaxis().SetTitleOffset(1.3)
        if type is 'event':
            if roc is self.roc_d1:
                histo.SetMaximum(self.Settings['ph1DmaxD4'])
            elif roc is self.roc_d2:
                histo.SetMaximum(self.Settings['ph1DmaxD5'])
            else:
                histo.SetMaximum(self.Settings['ph1DmaxSi'])
        histo.SetMinimum(min_val)
        histo.SetLineColor(color)
        histo.SetLineWidth(3*gStyle.GetLineWidth())
        if type is 'landaus':
            histo.SetFillColor(color)

    def create_2D_profile(self, type='spatial', name='histo', title='histo', xTitle='X', yTitle='Y', zTitle='Z', min_val=0, max_val=-1):
        xbins = self.Settings['nBinsX'] if type is 'spatial' else self.Settings['nBinCol']
        xmin = self.Settings['xmin'] if type is 'spatial' else self.Settings['minCol']
        xmax = self.Settings['xmax'] if type is 'spatial' else self.Settings['maxCol']
        ybins = self.Settings['nBinsY'] if type is 'spatial' else self.Settings['nBinRow']
        ymin = self.Settings['ymin'] if type is 'spatial' else self.Settings['minRow']
        ymax = self.Settings['ymax'] if type is 'spatial' else self.Settings['maxRow']
        histo2D = TProfile2D(name, title, int(xbins + 1), xmin - float(xmax-xmin)/(2*xbins),
                             xmax + float(xmax-xmin)/(2*xbins), int(ybins + 1), ymin - float(ymax-ymin)/(2*ybins),
                             ymax + float(ymax-ymin)/(2*ybins))
        self.set_2D_options(histo2D, xTitle, yTitle, zTitle, min_val, max_val)
        return histo2D

    def create_2D_histogram(self, type='spatial', name='histo', title='histo', xTitle='X', yTitle='Y', zTitle='Z', min_val=0, max_val=-1, roc=4):
        if type is 'spatial':
            xbins = self.Settings['nBinsX']
            xmin = self.Settings['xmin']
            xmax = self.Settings['xmax']
            ybins = self.Settings['nBinsY']
            ymin = self.Settings['ymin']
            ymax = self.Settings['ymax']
        elif type is 'pixel':
            xbins = self.Settings['nBinCol']
            xmin = self.Settings['minCol']
            xmax = self.Settings['maxCol']
            ybins = self.Settings['nBinRow']
            ymin = self.Settings['minRow']
            ymax = self.Settings['maxRow']
        elif type is 'correlpixcol':
            xbins = self.Settings['nBinCol']
            xmin = self.Settings['minCol']
            xmax = self.Settings['maxCol']
            ybins = self.Settings['nBinCol']
            ymin = self.Settings['minCol']
            ymax = self.Settings['maxCol']
        elif type is 'correlpixx':
            xbins = self.Settings['nBinsX']
            xmin = self.Settings['xmin']
            xmax = self.Settings['xmax']
            ybins = self.Settings['nBinsX']
            ymin = self.Settings['xmin']
            ymax = self.Settings['xmax']
        elif type is 'correlpixrow':
            xbins = self.Settings['nBinRow']
            xmin = self.Settings['minRow']
            xmax = self.Settings['maxRow']
            ybins = self.Settings['nBinRow']
            ymin = self.Settings['minRow']
            ymax = self.Settings['maxRow']
        elif type is 'correlpixy':
            xbins = self.Settings['nBinsY']
            xmin = self.Settings['ymin']
            xmax = self.Settings['ymax']
            ybins = self.Settings['nBinsY']
            ymin = self.Settings['ymin']
            ymax = self.Settings['ymax']
        else:
            xbins = self.Settings['event_bins']
            xmin = self.Settings['event_min']
            xmax = self.Settings['event_max']
            ybins = self.Settings['ph1DbinsD4'] if roc is self.roc_d1 else self.Settings['ph1DbinsD5'] if roc is self.roc_d2 else self.Settings['ph1DbinsSi']
            ymin = self.Settings['ph1DminD4'] if roc is self.roc_d1 else self.Settings['ph1DminD5'] if roc is self.roc_d2 else self.Settings['ph1DminSi']
            ymax = self.Settings['ph1DmaxD4'] if roc is self.roc_d1 else self.Settings['ph1DmaxD5'] if roc is self.roc_d2 else self.Settings['ph1DmaxSi']
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
        if min_val is not 'auto': histo.SetMinimum(min_val)
        if max_val is not -1: histo.SetMaximum(max_val)


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
        c0.Update()
        c0.Modified()
        c0.SaveAs('{dir}/Root/c_{n}.root'.format(dir=path, n=name))
        c0.SaveAs('{dir}/Plots/c_{n}.png'.format(dir=path, n=name))
        c0.Close()
        gROOT.SetBatch(False)
        if verbosity: self.print_banner('{n} save -> Done'.format(n=name))
        del c0

    def save_cuts_distributions(self, histo1, histo2, name, title, draw_opt='', opt_stats=0, path='./', verbosity=False, histo3='', doLogY=False):
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
        if doLogY: c0.SetLogy()
        c0.Update()
        c0.Modified()
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
        if verbosity: print_banner('{n} save -> Done'.format(n=name))
        del c0

    def check_plot_existence(self, path, name):
        if os.path.isfile('{p}/Plots/{n}.png'.format(p=path, n=name)):
            return True
        else:
            return False

    def get_binning(self, evts_per_bin, start_event=0, end_event=None, time_bins=False):
        binning = arange(start_event, self.NEntries if end_event is None else end_event, self.NEntries / (self.NEntries / evts_per_bin), 'd')
        return [len(binning) - 1, array([self.run.get_time_at_event(int(ev)) for ev in binning], 'd') if time_bins else binning]

    def get_time_binning(self, bin_width, rel_time=False, start_time=0, end_time=None):
        """ returns bins with fixed time width. bin_width in secons """
        end_time = self.run.EndTime if end_time is None else end_time
        bins = range(self.run.StartTime + start_time, end_time, bin_width)
        bins += [end_time] if bins[-1] != end_time else []
        return [len(bins) - 1, array(bins, 'd') - (self.run.StartTime if rel_time else 0)]
