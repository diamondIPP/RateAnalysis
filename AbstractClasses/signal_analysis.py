# ==============================================
# IMPORTS
# ==============================================
from ROOT import TGraphErrors, TCanvas, TH2D, gStyle, TH1F, gROOT, TLegend, TCut, TGraph, TProfile2D, TH2F, TProfile, TCutG
from newAnalysis import Analysis
from numpy import array
from math import sqrt, ceil, log
from argparse import ArgumentParser
from Extrema import Extrema2D
from time import time, sleep
from collections import OrderedDict
from sys import stdout

__author__ = 'micha'


# ==============================================
# MAIN CLASS
# ==============================================
class SignalAnalysis(Analysis):
    def __init__(self, run, channel, low_rate_run=None, binning=20000):
        Analysis.__init__(self, run, low_rate=low_rate_run)

        # main
        self.channel = channel
        self.run_number = self.run.run_number
        self.diamond_name = self.run.diamond_names[channel]
        self.bias = self.run.bias[channel]
        self.cut = self.cuts[channel]
        self.save_dir = '{tc}_{run}_{dia}'.format(tc=self.TESTCAMPAIGN[2:], run=self.run_number, dia=self.diamond_name)

        # stuff
        self.draw_option = 'COLZ'
        self.bin_size = binning
        self.binning = self.__get_binning()
        self.time_binning = self.get_time_binning()
        self.n_bins = len(self.binning)
        self.polarity = self.polarities[channel]

        # names
        self.signal_name = self.signal_names[channel]
        self.pedestal_name = self.pedestal_names[channel]

        # projection
        self.signal_projections = {}
        # graphs
        self.PulseHeight = None
        self.Pedestal = None
        self.signals = {}
        self.tmp_histos = {}
        # histograms
        self.signaltime = None
        self.SignalMapHisto = None
        self.MeanSignalHisto = None
        self.PeakValues = None

    def set_channel(self, ch):
        self.channel = ch
        self.diamond_name = self.run.diamondname[ch]
        self.bias = self.run.bias[ch]
        self.cut = self.cuts[ch]
        self.save_dir = '{tc}_{run}_{dia}'.format(tc=self.TESTCAMPAIGN[2:], run=self.run_number, dia=self.run.diamondname[ch])
        self.polarity = self.polarities[ch]
        self.signal_name = self.signal_names[ch]
        self.pedestal_name = self.pedestal_names[ch]

    def __set_bin_size(self, value):
        self.bin_size = value
        self.binning = self.__get_binning()
        self.time_binning = self.get_time_binning()
        self.n_bins = len(self.binning)
        return value

    # ==========================================================================
    # region 2D SIGNAL DISTRIBUTION
    def draw_signal_map(self, draw_option='surf2', show=True, factor=4):
        margins = self.find_diamond_margins(show_plot=False)
        x = [margins['x'][0], margins['x'][1]]
        y = [margins['y'][0], margins['y'][1]]
        nr = 1 if not self.channel else 2
        # get bin size via digital resolution of the telescope pixels
        x_bins = int(ceil(((x[1] - x[0]) / 0.015 * sqrt(12) / factor)))
        y_bins = int(ceil((y[1] - y[0]) / 0.01 * sqrt(12) / factor))
        h = TProfile2D('signal_map', 'Signal Map', x_bins, x[0], x[1], y_bins, y[0], y[1])
        if not show:
            gROOT.SetBatch(1)
        signal = '{sig}-{pol}*{ped}'.format(sig=self.signal_name, ped=self.pedestal_name, pol=self.polarity)
        print 'drawing signal map of {dia} for Run {run}...'.format(dia=self.diamond_name, run=self.run_number)
        self.tree.Draw('{z}:diam{nr}_track_x:diam{nr}_track_y>>signal_map'.format(z=signal, nr=nr), self.cut.all_cut, 'goff')
        c = TCanvas('c', 'Signal Map', 1000, 1000)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.12)
        gStyle.SetPalette(53)
        self.format_histo(h, x_tit='track_x [cm]', y_tit='track_y [cm]', y_off=1.6)
        if draw_option == 'surf2':
            self.format_histo(h, x_off=1.6, y_off=2.2, x_tit='track_x [cm]', y_tit='track_y [cm]')
        h.SetStats(0)
        h.Draw(draw_option)
        self.save_plots('SignalMap2D_' + draw_option, sub_dir=self.save_dir)
        gROOT.SetBatch(0)
        self.SignalMapHisto = h
        self.canvases[0] = c
        return h

    def make_region_cut(self):
        self.draw_mean_signal_distribution(show=False)
        return self.cut.generate_region(self.SignalMapHisto, self.MeanSignalHisto)

    def find_2d_regions(self):
        self.draw_mean_signal_distribution(show=False)
        extrema = Extrema2D(self.SignalMapHisto, self.MeanSignalHisto)
        extrema.clear_voting_histos()
        extrema.region_scan()
        extrema.show_voting_histos()
        self.save_plots('Regions2D', sub_dir=self.save_dir)
        return extrema

    def find_2d_extrema(self, size=1, histo=None, show=True):
        self.draw_mean_signal_distribution(show=False)
        extrema = Extrema2D(self.SignalMapHisto, self.MeanSignalHisto)
        extrema.clear_voting_histos()
        extrema.square_scan(size, histo)
        if show:
            extrema.show_voting_histos()
        self.save_plots('Extrema2D', sub_dir=self.save_dir)
        return extrema

    def draw_mean_signal_distribution(self, show=True):
        """
        Draws the distribution of the mean pulse height values of the bins from the signal map
        :param show: shows a plot of the canvas if True
        """
        sig_map = self.SignalMapHisto if self.SignalMapHisto is not None else self.draw_signal_map(show=False)
        x = [int(sig_map.GetMinimum()) / 10 * 10, int(sig_map.GetMaximum() + 10) / 10 * 10]
        h = TH1F('h', 'Mean Signal Distribution', 50, x[0], x[1])
        for bin_ in xrange((sig_map.GetNbinsX() + 2) * (sig_map.GetNbinsY() + 2)):
            h.Fill(sig_map.GetBinContent(bin_))
        if show:
            self.canvases[0] = TCanvas('c', 'Mean Signal Distribution', 1000, 1000)
            self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.2)
            h.Draw()
            self.save_plots('MeanSignalHisto', sub_dir=self.save_dir)
        self.MeanSignalHisto = h

    def draw_error_signal_map(self, show=False):
        self.draw_mean_signal_distribution(show=False)
        h = self.SignalMapHisto.ProjectionXY('', 'c=e')
        if show:
            c = TCanvas('c', 'Signal Map Errors', 1000, 1000)
            c.SetLeftMargin(0.12)
            c.SetRightMargin(0.11)
            self.format_histo(h, name='sig_map_errors', title='Signal Map Errors', x_tit='track_x [cm]', y_tit='track_y [cm]', y_off=1.6)
            h.SetStats(0)
            h.Draw('colz')
            self.save_plots('SignalMapErrors', sub_dir=self.save_dir, canvas=c)
            self.canvases[0] = c
            self.histos[0] = h
        return h

    def fit_mean_signal_distribution(self):
        pickle_path = self.get_program_dir() + self.pickle_dir + 'MeanSignalFit/{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            self.draw_mean_signal_distribution(show=False)
            return self.MeanSignalHisto.Fit('gaus', 'qs')

        fit = self.do_pickle(pickle_path, func)
        return fit

    def get_mean_fwhm(self):
        fit = self.fit_mean_signal_distribution()
        conversion_factor = 2 * sqrt(2 * log(2))  # sigma to FWHM
        return fit.Parameter(2) * conversion_factor

    def draw_diamond_hitmap(self):
        self.find_diamond_margins(show_frame=True)

    def find_diamond_margins(self, show_plot=True, show_frame=False):
        pickle_path = self.get_program_dir() + self.pickle_dir + 'Margins/{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            print 'getting margins for {dia} of run {run}...'.format(dia=self.diamond_name, run=self.run_number)
            if not show_plot:
                gROOT.SetBatch(1)
            h = TH2F('h', 'Diamond Margins', 80, -.3, .3, 52, -.3, .3)
            nr = 1 if not self.channel else 2
            self.tree.Draw('diam{nr}_track_x:diam{nr}_track_y>>h'.format(nr=nr), self.cut.all_cut, 'goff')
            projections = [h.ProjectionX(), h.ProjectionY()]
            efficient_bins = [[], []]
            zero_bins = [[], []]
            bin_low = [[], []]
            bin_high = [[], []]
            for i, proj in enumerate(projections):
                last_bin = None
                for bin_ in xrange(proj.GetNbinsX()):
                    efficiency = proj.GetBinContent(bin_) / float(proj.GetMaximum())
                    if efficiency > .3:
                        efficient_bins[i].append(proj.GetBinCenter(bin_))
                        bin_low[i].append(proj.GetBinLowEdge(bin_))
                        bin_high[i].append(proj.GetBinLowEdge(bin_ + 1))
                    if bin_:
                        if efficiency and not last_bin:
                            zero_bins[i].append(proj.GetBinCenter(bin_ - 1))
                        elif not efficiency and last_bin:
                            zero_bins[i].append((proj.GetBinCenter(bin_)))
                    last_bin = proj.GetBinContent(bin_)
            if show_plot:
                print zero_bins
                self.canvas = TCanvas('c', 'Diamond Hit Map', 1000, 1000)
                h.GetXaxis().SetRangeUser(zero_bins[0][0], zero_bins[0][1])
                h.GetYaxis().SetRangeUser(zero_bins[1][0], zero_bins[1][1])
                h.Draw('colz')
                if show_frame:
                    self.__show_frame(bin_low, bin_high)
                self.histos[0] = h
                self.save_plots('diamond_hitmap', sub_dir=self.save_dir, canvas=self.canvas)
            gROOT.SetBatch(0)
            return {name: [efficient_bins[i][0], efficient_bins[i][-1]] for i, name in enumerate(['x', 'y'])}

        margins = self.do_pickle(pickle_path, func) if not show_plot else func()
        return margins

    def __show_frame(self, bin_low, bin_high):
        frame = TCutG('frame', 4)
        frame.SetLineColor(2)
        frame.SetLineWidth(4)
        frame.SetVarX('x')
        frame.SetVarY('y')
        frame.SetPoint(0, bin_low[0][0], bin_low[1][0])
        frame.SetPoint(1, bin_high[0][-1], bin_low[1][0])
        frame.SetPoint(2, bin_high[0][-1], bin_high[1][-1])
        frame.SetPoint(3, bin_low[0][0], bin_high[1][-1])
        frame.SetPoint(4, bin_low[0][0], bin_low[1][0])
        frame.Draw('same')
        self.histos[1] = frame

    def calc_signal_spread(self, min_percent=5, max_percent=99):
        """
        Calculates the relative spread of mean signal response from the 2D signal response map.
        :param min_percent: min quantile
        :param max_percent: max quantile
        :return: relative spread [%]
        """
        if self.MeanSignalHisto is None:
            self.draw_mean_signal_distribution(show=False)
        q = array([min_percent / 100., max_percent / 100.])
        y = array([0., 0.])
        self.MeanSignalHisto.GetQuantiles(2, y, q)
        max_min_ratio = (y[1] / y[0] - 1) * 100
        delta_y = self.draw_error_signal_map(show=False).GetMinimum()
        # error propagation
        err = 100 * delta_y / y[0] * (1 + y[1] / y[0])
        print 'Relative Signal Spread is: {spr} +- {err}'.format(spr=max_min_ratio, err=err)
        return [max_min_ratio, err]

    # endregion

    # ==========================================================================
    # region PEAK VALUES
    def draw_peak_values(self, region='b', draw=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        num = self.get_signal_numbers(region, self.peak_integral)[self.channel]
        peak_val = 'IntegralPeaks[{num}]'.format(num=num)
        title = 'Peak Values {reg}{int}'.format(reg=region, int=self.peak_integral)
        x = self.run.signal_regions[region]
        h = TH1F('peakvalues', title, x[1] - x[0], x[0] / 2., x[1] / 2.)
        self.format_histo(h, x_tit='time [ns]', y_tit='Entries', y_off=2)
        cut = self.cut.all_cut
        self.tree.Draw(peak_val + '/2.>>peakvalues', cut, 'goff')
        if draw:
            c = TCanvas('c', 'Signal Peak Distribution', 1000, 1000)
            c.SetLeftMargin(0.14)
            h.Draw()
            self.save_plots('peak_values_{reg}{int}'.format(reg=region, int=self.peak_integral), 'png', canvas=c, sub_dir=self.save_dir)
            self.canvases[0] = c
        self.PeakValues = h
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')

    def fit_peak_values(self, draw=True):
        pickle_path = self.get_program_dir() + self.pickle_dir + 'PeakValues/Fit_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            print 'Getting peak value fit for {dia} of run {run}...'.format(run=self.run_number, dia=self.diamond_name)
            if self.PeakValues is None:
                self.draw_peak_values(draw=False)
            h = self.PeakValues
            max_bin = h.GetMaximumBin()
            x = [h.GetBinCenter(max_bin + i) for i in [-7, 1]]
            return h.Fit('gaus', 'qs0', '', x[0], x[1])

        fit = func() if draw else self.do_pickle(pickle_path, func)
        return fit

    def calc_peak_value_fwhm(self):
        pickle_path = self.get_program_dir() + self.pickle_dir + 'PeakValues/FWHM_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            print 'Getting peak value FWHM for {dia} of run {run}...'.format(run=self.run_number, dia=self.diamond_name)
            if self.PeakValues is None:
                self.draw_peak_values(draw=False)
            return self.calc_fwhm(self.PeakValues)

        fwhm = self.do_pickle(pickle_path, func)
        return fwhm

    # endregion

    # ==========================================================================
    # region SIGNAL/PEDESTAL
    def make_signal_time_histos(self, ped=False, signal=None, corr=False, show=True):
        gROOT.SetBatch(1)
        signal = self.signal_name if signal is None else signal
        signal = signal if not ped else self.pedestal_name
        signal = '{sig}-{pol}*{ped}'.format(sig=signal, ped=self.pedestal_name, pol=self.polarity) if corr else signal
        # 2D Histogram
        name = "signaltime_" + str(self.run_number)
        xbins = array(self.time_binning)
        x_min = -50 if not ped else -20
        x_max = 300 if not ped else 20
        bins = 1000 if not ped else 80
        h = TH2D(name, "signaltime", len(xbins) - 1, xbins, bins, x_min, x_max)
        self.tree.Draw("{name}:time>>{histo}".format(histo=name, name=signal), self.cut.all_cut, 'goff')
        if show:
            gROOT.SetBatch(0)
            c = TCanvas('c', 'Pulse Height vs Time', 1000, 1000)
            c.SetLeftMargin(.12)
            self.format_histo(h, x_tit='time [min]', y_tit='Pulse Height [au]', y_off=1.4)
            h.Draw('colz')
            self.save_plots('SignalTime', sub_dir=self.save_dir)
            self.signaltime = h
            self.canvases[0] = c
        gROOT.SetBatch(0)
        return h

    def draw_pedestal(self, binning=None, draw=True):
        bin_size = binning if binning is not None else self.bin_size
        picklepath = 'Configuration/Individual_Configs/Pedestal/{tc}_{run}_{ch}_{bins}_Ped_Means.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, ch=self.channel, bins=bin_size)
        gr = self.make_tgrapherrors('pedestal', 'Pedestal')

        def func():
            print 'calculating pedestal of ch', self.channel
            if binning is not None:
                self.__set_bin_size(binning)
            ped_time = self.make_signal_time_histos(ped=True, show=False)
            gROOT.SetBatch(1)
            means = []
            empty_bins = 0
            count = 0
            for i in xrange(self.n_bins):
                h_proj = ped_time.ProjectionY(str(i), i + 1, i + 1)
                if h_proj.GetEntries() > 0:
                    fit = self.fit_fwhm(h_proj)
                    gr.SetPoint(count, (self.time_binning[i] - self.run.startTime) / 60e3, fit.Parameter(1))
                    gr.SetPointError(count, 0, fit.ParError(1))
                    count += 1
                    means.append(fit.Parameter(1))
                else:
                    empty_bins += 1
            if draw:
                gROOT.SetBatch(0)
            print 'Empty proj. bins:\t', str(empty_bins) + '/' + str(self.n_bins)
            fit_pars = gr.Fit('pol0', 'qs')
            print 'mean:', fit_pars.Parameter(0), '+-', fit_pars.ParError(0)
            c = TCanvas('bla', 'blub', 1000, 1000)
            c.SetLeftMargin(.14)
            gStyle.SetOptFit(1)
            self.format_histo(gr, x_tit='time [min]', y_tit='Mean Pulse Height [au]', y_off=1.6)
            gr.Draw('alp')
            gr.Draw()
            self.save_plots('Pedestal', sub_dir=self.save_dir)
            self.Pedestal = gr
            self.canvases[0] = c
            gROOT.SetBatch(0)
            return means

        all_means = self.do_pickle(picklepath, func)
        if draw and not gROOT.FindObject('pedestal'):
            func()
        return all_means

    def draw_pulse_height(self, binning=None, draw=True, ped_corr=False, eventwise_corr=False, sig=None):
        signal = self.signal_name if sig is None else sig
        bin_size = binning if binning is not None else self.bin_size
        correction = ''
        if ped_corr:
            correction = 'binwise'
        elif eventwise_corr:
            correction = 'eventwise'
        suffix = '{bins}_{cor}_{sig}'.format(bins=bin_size, cor=correction, sig=self.get_all_signal_names()[signal])
        picklepath = 'Configuration/Individual_Configs/Ph_fit/{tc}_{run}_{ch}_{suf}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, ch=self.channel, suf=suffix)

        self.signaltime = None

        def func():
            print 'calculating pulse height fit of ch', self.channel
            tit_suffix = 'with {cor} Pedestal Correction'.format(cor=correction.title()) if ped_corr or eventwise_corr else ''
            gr = self.make_tgrapherrors('signal', 'Pulse Height Evolution ' + tit_suffix)
            if binning is not None:
                self.__set_bin_size(binning)
            sig_time = self.make_signal_time_histos(corr=eventwise_corr, signal=signal, show=False)
            mode = 'mean'
            empty_bins = 0
            count = 0
            means = self.draw_pedestal(bin_size, draw=False) if ped_corr else None
            gROOT.SetBatch(1)
            for i in xrange(self.n_bins):
                h_proj = sig_time.ProjectionY(str(i), i + 1, i + 1)
                if h_proj.GetEntries() > 0:
                    if mode in ["mean", "Mean"]:
                        mean = h_proj.GetMean()
                        mean -= self.polarity * means[count] if ped_corr else 0
                        gr.SetPoint(count, (self.time_binning[i] - self.run.startTime) / 60e3, mean)
                        gr.SetPointError(count, 0, h_proj.GetRMS() / sqrt(h_proj.GetEntries()))
                    elif mode in ["fit", "Fit"]:
                        h_proj.GetMaximum()
                        maxposition = h_proj.GetBinCenter(h_proj.GetMaximumBin())
                        h_proj.Fit("landau", "Q", "", maxposition - 50, maxposition + 50)
                        fitfun = h_proj.GetFunction("landau")
                        mpv = fitfun.GetParameter(1)
                        mpverr = fitfun.GetParError(1)
                        gr.SetPoint(count, (i + 0.5) * self.run.totalMinutes / self.n_bins, mpv)
                        gr.SetPointError(count, 0, mpverr)
                    count += 1
                else:
                    empty_bins += 1
            print 'Empty proj. bins:\t', str(empty_bins) + '/' + str(self.n_bins)
            if draw:
                gROOT.SetBatch(0)
            c = TCanvas('bla', 'blub', 1000, 1000)
            c.SetLeftMargin(.14)
            gStyle.SetOptFit(1)
            self.format_histo(gr, x_tit='time [min]', y_tit='Mean Pulse Height [au]', y_off=1.6)
            fit_par = gr.Fit('pol0', 'qs') if draw else gr.Fit('pol0', 'qs0')
            gr.Draw('alp')
            self.save_plots('PulseHeight', sub_dir=self.save_dir)
            self.PulseHeight = gr
            self.canvas = c
            gROOT.SetBatch(0)
            return fit_par

        fit = self.do_pickle(picklepath, func)
        if draw and not gROOT.FindObject('bla'):
            func()
        return fit

    def show_signal_histo(self, cut=None, corr=True):
        suffix = 'with Pedestal Correction' if corr else ''
        h = TH1F('signal b2', 'Pulse Height ' + suffix, 350, -50, 300)
        cut = self.cut.all_cut if cut is None else cut
        signal = '{sig}-{pol}*{ped}'.format(sig=self.signal_name, ped=self.pedestal_name, pol=self.polarity) if corr else self.signal_name
        self.tree.Draw('{name}>>signal b2'.format(name=signal), cut, 'goff')
        c = TCanvas('c', 'Signal Distribution', 1000, 1000)
        c.SetLeftMargin(.13)
        self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.8)
        h.Draw()
        self.save_plots('SignalDistribution', sub_dir=self.save_dir)
        self.histos[0] = h
        self.canvases[0] = c

    def show_pedestal_histo(self, region='ab', peak_int='2', cut=True, fwhm=True, draw=True):
        cut = self.cut.all_cut if cut else TCut()
        fw = 'fwhm' if fwhm else 'full'
        suffix = '{reg}_{fwhm}_{cut}'.format(reg=region + str(peak_int), cut=cut.GetName(), fwhm=fw)
        picklepath = 'Configuration/Individual_Configs/Pedestal/{tc}_{run}_{ch}_{suf}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, ch=self.channel, suf=suffix)

        def func():
            gROOT.SetBatch(1)
            print 'making pedestal histo for region {reg}{int}...'.format(reg=region, int=peak_int)
            h = TH1F('ped1', 'Pedestal Distribution', 100, -20, 20)
            name = self.get_pedestal_names(region, peak_int)[self.channel]
            self.tree.Draw('{name}>>ped1'.format(name=name), cut, 'goff')
            fit_pars = self.fit_fwhm(h, do_fwhm=fwhm, draw=draw)
            gStyle.SetOptFit(1)
            if draw:
                gROOT.SetBatch(0)
            c = TCanvas('c', 'Pedestal Distribution', 1000, 1000)
            c.SetLeftMargin(.13)
            self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.8)
            h.Draw()
            save_name = 'Pedestal_{reg}{cut}'.format(reg=region, cut=cut.GetName())
            self.save_plots(save_name, 'png', canvas=c, sub_dir=self.save_dir)
            self.tmp_histos[0] = h
            self.canvas = c
            gROOT.SetBatch(0)
            return fit_pars

        fit_par = self.do_pickle(picklepath, func)
        if draw and not gROOT.FindObject('ped1'):
            func()
        return fit_par

    def compare_pedestals(self):
        legend = TLegend(0.7, 0.7, 0.98, .9)
        gr1 = TGraph()
        gr1.SetTitle('pedestal comparison')
        gr1.SetMarkerStyle(20)
        gr2 = TGraph()
        gr2.SetTitle('pedestal comparison with cuts')
        gr2.SetMarkerStyle(20)
        gr2.SetMarkerColor(2)
        gr2.SetLineColor(2)
        gr3 = TGraph()
        gr3.SetTitle('pedestal comparison with cuts full fit')
        gr3.SetMarkerStyle(20)
        gr3.SetMarkerColor(3)
        gr3.SetLineColor(3)
        gROOT.SetBatch(1)
        gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
        for i, reg in enumerate(self.run.pedestal_regions):
            print 'calculation region', reg
            mean1 = self.show_pedestal_histo(reg).keys()[1]
            mean2 = self.show_pedestal_histo(reg, 'median').keys()[1]
            mean3 = self.show_pedestal_histo(reg, 'all').keys()[1]
            gr1.SetPoint(i, i, mean1)
            gr2.SetPoint(i, i, mean2)
            gr3.SetPoint(i, i, mean3)
        gROOT.SetBatch(0)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        for i, reg in enumerate(self.run.pedestal_regions):
            bin_x = gr1.GetXaxis().FindBin(i)
            gr1.GetXaxis().SetBinLabel(bin_x, reg)
        self.canvases[0] = TCanvas('bla', 'blub', 1000, 1000)
        self.tmp_histos[1] = gr1
        self.tmp_histos[0] = gr2
        self.tmp_histos[2] = gr3
        gr1.Draw('alp')
        gr2.Draw('lp')
        gr3.Draw('lp')
        legend.AddEntry(gr1, 'mean fit fwhm w/ cuts 2', 'lp')
        legend.AddEntry(gr2, 'mean fit fwhm w/ cuts median', 'lp')
        legend.AddEntry(gr3, 'mean fit fwhm w/ cuts all', 'lp')
        legend.Draw()
        self.tmp_histos[4] = legend

    # endregion

    def draw_pulser_rate(self, binning=200):
        """
        Shows the fraction of accepted pulser events as a function of event numbers. Peaks appearing in this graph are most likely beam interruptions.
        :param binning:
        """
        nbins = self.run.n_entries / binning
        h = TProfile('h', 'Pulser Rate', nbins, 0, z.run.n_entries)
        self.tree.Draw('(pulser!=0)*100:Entry$>>h', '', 'goff')
        c = TCanvas('c', 'Pulser Rate Canvas', 1000, 1000)
        self.format_histo(h, name='pulser_rate', title='Pulser Rate', x_tit='Event Number', y_tit='Pulser Fraction [%]', y_off=1.3)
        h.Draw('hist')
        self.run.draw_run_info(canvas=c, channel=self.channel)
        self.save_plots('pulser_rate', canvas=c, sub_dir=self.save_dir)
        self.canvases[0] = c
        self.histos[0] = h

    # ==========================================================================
    # region CUTS
    def compare_single_cuts(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        c1 = TCanvas('single', '', 1000, 1000)
        c2 = TCanvas('all', '', 1000, 1000)
        c2.SetLeftMargin(0.15)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        histos = []
        drawn_first = False
        for key, value in self.cut.cut_strings.iteritems():
            if str(value) or key == 'raw':
                print 'saving plot', key
                save_name = 'signal_distribution_{cut}'.format(cut=key)
                histo_name = 'signal {range}{peakint}'.format(range=self.signal_region, peakint=self.peak_integral)
                histo_title = 'signal with cut ' + key
                histo = TH1F(histo_name, histo_title, 350, -50, 300)
                # safe single plots
                c1.cd()
                self.tree.Draw("{name}>>{histo}".format(name=self.signal_name, histo=histo_name), value)
                self.save_plots(save_name, 'png', canvas=c1, sub_dir=self.save_dir)
                # draw all single plots into c2
                c2.cd()
                histo.SetLineColor(self.get_color())
                if not drawn_first:
                    self.format_histo(histo, title='Signal Distribution of Different Single Cuts', x_tit='Pulse Height [au]', y_tit='Entries', y_off=2)
                    histo.SetStats(0)
                    histo.Draw()
                    drawn_first = True
                else:
                    if key == 'all_cuts':
                        histo.SetLineWidth(2)
                    histo.Draw('same')
                histos.append(histo)
                legend.AddEntry(histo, key, 'l')
        # save c2
        legend.Draw()
        self.save_plots('all', 'png', canvas=c2, sub_dir=self.save_dir)
        self.save_plots('all', 'root', canvas=c2, sub_dir=self.save_dir)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        gROOT.SetBatch(0)

    def compare_normalised_cuts(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        c1 = TCanvas('single', '', 1000, 1000)
        c2 = TCanvas('normalised', '', 1000, 1000)
        c2.SetLeftMargin(0.15)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        histos = []
        drawn_first = False
        for key, value in self.cut.cut_strings.iteritems():
            if str(value) or key == 'raw':
                print 'saving plot', key
                save_name = 'signal_distribution_normalised_{cut}'.format(cut=key)
                histo_name = 'signal {range}{peakint}'.format(range=self.signal_region, peakint=self.peak_integral)
                histo_title = 'normalised signal with cut ' + key
                histo = TH1F(histo_name, histo_title, 350, -50, 300)
                # safe single plots
                c1.cd()
                self.tree.Draw("{name}>>{histo}".format(name=self.signal_name, histo=histo_name), value)
                histo = self.normalise_histo(histo)
                histo.Draw()
                self.save_plots(save_name, 'png', canvas=c1, sub_dir=self.save_dir)
                # draw all single plots into c2
                c2.cd()
                histo.SetLineColor(self.get_color())
                if not drawn_first:
                    self.format_histo(histo, title='Normalised Signal Distribution with Single Cuts', x_tit='Pulse Height [au]', y_tit='Normalised Integral', y_off=2)
                    histo.SetStats(0)
                    histo.Draw()
                    drawn_first = True
                else:
                    if key == 'all_cuts':
                        histo.SetLineWidth(2)
                    histo.Draw('same')
                histos.append(histo)
                legend.AddEntry(histo, key, 'l')
        # save c2
        legend.Draw()
        self.save_plots('normalised', 'png', canvas=c2, sub_dir=self.save_dir)
        self.save_plots('normalised', 'root', canvas=c2, sub_dir=self.save_dir)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        gROOT.SetBatch(0)

    def compare_consecutive_cuts(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        c1 = TCanvas('consecutive', '', 1000, 1000)
        c2 = TCanvas('all', '', 1000, 1000)
        c2.SetLeftMargin(0.15)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        histos = []
        drawn_first = False
        ind = 0
        cut = TCut('consecutive', '')
        for key, value in self.cut.cut_strings.iteritems():
            if (str(value) or key == 'raw') and key != 'all_cuts':
                cut += value
                print 'saving plot with {n} cuts'.format(n=ind)
                save_name = 'signal_distribution_{n}cuts'.format(n=ind)
                histo_name = 'signal {range}{peakint}'.format(range=self.signal_region, peakint=self.peak_integral)
                histo_title = 'signal with {n} cuts'.format(n=ind)
                histo = TH1F(histo_name, histo_title, 550, -50, 500)
                # safe single plots
                c1.cd()
                self.tree.Draw("{name}>>{histo}".format(name=self.signal_name, histo=histo_name), cut)
                self.save_plots(save_name, 'png', canvas=c1, sub_dir=self.save_dir)
                # draw all single plots into c2
                c2.cd()
                color = self.get_color()
                histo.SetLineColor(color)
                histo.SetFillColor(color)
                if not drawn_first:
                    self.format_histo(histo, title='Signal Distribution with Consecutive Cuts', x_tit='Pulse Height [au]', y_tit='Entries', y_off=2)
                    histo.SetStats(0)
                    histo.Draw()
                    drawn_first = True
                    legend.AddEntry(histo, key, 'f')
                else:
                    histo.Draw('same')
                    legend.AddEntry(histo, '+ ' + key, 'f')
                histos.append(histo)
                ind += 1
        # save c2
        legend.Draw()
        self.save_plots('consecutive', 'png', canvas=c2, sub_dir=self.save_dir)
        self.save_plots('consecutive', 'root', canvas=c2, sub_dir=self.save_dir)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        gROOT.SetBatch(0)

    # endregion

    # ==========================================================================
    # region SHOW
    def draw_bucket_pedestal(self, show=True):
        gROOT.SetBatch(1)
        reg_name = 'e2'
        three_bucket_num = self.get_signal_numbers(reg_name[0], reg_name[1])[self.channel]
        reg_margins = self.run.signal_regions[reg_name[0]]
        x_bins = (reg_margins[1] - reg_margins[0])
        h = TH2F('h', 'Bucket Pedestal', x_bins, reg_margins[0] / 2., reg_margins[1] / 2., 550, -50, 500)
        self.tree.Draw('{sig}:IntegralPeaks[{num}]/2>>h'.format(sig=self.signal_name, num=three_bucket_num), z.cut.cut_strings['tracks'] + z.cut.cut_strings['pulser'], 'goff')
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Bucket Pedestal', 1000, 1000)
        c.SetLogz()
        self.format_histo(h, x_tit='Highest Peak Timing [ns]', y_tit='Pulse Height [au]', y_off=1.25)
        h.SetStats(0)
        h.Draw('colz')
        self.save_plots('BucketPedestal', sub_dir=self.save_dir)
        gROOT.SetBatch(0)
        self.histos[0] = [c, h]

    def draw_waveforms(self, n=1000, start_event=None, cut_string=None, show=True, ret_event=False):
        """
        Draws stacked waveforms.
        :param n: number of waveforms
        :param cut_string:
        :param start_event: event to start
        :param show:
        :param ret_event: return number of valid events if True
        :return: histo with waveform
        """
        gROOT.SetBatch(1)
        t = time()
        start = self.start_event if start_event is None else start_event
        assert self.run.n_entries >= start >= 0, 'The start event is not within the range of tree events!'
        if not self.run.wf_exists(self.channel):
            return
        cut = self.cut.all_cut if cut_string is None else cut_string
        n_events = self.find_n_events(n, cut, start)
        h = TH2F('wf', 'Waveform', 1024, 0, 511, 1000, -500, 500)
        h.SetStats(0)
        gStyle.SetPalette(55)
        self.tree.Draw('wf0:Iteration$/2>>wf', cut, 'goff', n_events, start)
        h.GetYaxis().SetRangeUser(-500 + h.FindFirstBinAbove(0, 2) / 50 * 50, -450 + h.FindLastBinAbove(0, 2) / 50 * 50)
        self.print_elapsed_time(t)
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'WaveForm', 1000, 500)
        c.SetRightMargin(.045)
        self.format_histo(h, x_tit='Time [ns]', y_tit='Signal [au]')
        h.Draw('col')
        gROOT.SetBatch(0)
        if n > 1:
            self.save_plots('WaveForms{n}'.format(n=n), sub_dir=self.save_dir)
        self.histos[0] = [c, h]
        return h if not ret_event else n_events

    def show_single_waveforms(self, n, cut=None):
        ev = self.start_event
        for i in xrange(n):
            ev += self.draw_waveforms(n=1, cut_string=cut, ret_event=True, start_event=ev)
            sleep(1)

    # endregion

    def find_n_events(self, n_events, cut, start):
        """
        Finds the amount of events from the startevent that are not subject to the cut.
        :param n_events: number of wanted events
        :param cut:
        :param start:
        :return: actual number of events s.t. n_events are drawn
        """
        print 'Finding the correct number of events',
        n = self.tree.Draw('1', cut, 'goff', n_events, start)
        new_events = n_events
        while n != n_events:
            diff = n_events - n
            new_events += diff * 2 if abs(diff) > 1 else diff
            print '\b.',
            stdout.flush()
            n = self.tree.Draw('1', cut, 'goff', new_events, start)
        print
        return new_events

    @staticmethod
    def normalise_histo(histo):
        h = histo
        h.GetXaxis().SetRangeUser(0, 30)
        min_bin = h.GetMinimumBin()
        h.GetXaxis().UnZoom()
        max_bin = h.GetNbinsX() - 1
        h.Scale(1 / h.Integral(min_bin, max_bin))
        return h

    def analyise_signal_histograms(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        # gROOT.SetBatch(1)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        gr1 = TGraphErrors()
        gr1.SetTitle('mean values')
        gr1.SetMarkerStyle(20)
        gr2 = TGraph()
        gr2.SetTitle('median values')
        gr2.SetMarkerStyle(21)
        gr2.SetMarkerColor(2)
        gr3 = TGraph()
        gr3.SetMarkerStyle(22)
        gr3.SetMarkerColor(3)
        histos = []
        i = 0
        for key, value in self.cut.cut_strings.iteritems():
            if str(value) or key == 'raw':
                print 'process cut ' + key
                h = TH1F('h', '', 600, -100, 500)
                self.tree.Draw("{name}>>h".format(name=self.signal_name), value)
                mean = self.__get_mean(h)
                median = self.__get_median(h)
                mpv = self.__get_mpv(h)
                # print mean, median, mpv
                gr1.SetPoint(i, i, mean[0])
                gr1.SetPointError(i, 0, mean[1])
                gr2.SetPoint(i, i, median)
                gr3.SetPoint(i, i, mpv)
                histos.append(h)
                i += 1
        # rename bins
        legend.AddEntry(gr1, 'mean', 'lp')
        legend.AddEntry(gr2, 'median', 'lp')
        legend.AddEntry(gr3, 'mpv', 'lp')
        xaxis = gr1.GetXaxis()
        i = 0
        for key, value in self.cut.cut_strings.iteritems():
            if str(value) or key == 'raw':
                bin_x = xaxis.FindBin(i)
                gr1.GetXaxis().SetBinLabel(bin_x, key[:7])
                i += 1
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        # gROOT.SetBatch(0)
        c1 = TCanvas('c1', '', 1000, 1000)
        c1.cd()
        gr1.GetXaxis().SetRangeUser(-1, len(histos) + 1)
        gr1.Draw('alp')
        gr2.Draw('lp')
        gr3.Draw('lp')
        legend.Draw()
        self.histos['legend'] = legend
        return [gr1, gr2, gr3]

    @staticmethod
    def __get_histo_without_pedestal(histo):
        h = histo
        h.GetXaxis().SetRangeUser(0, 30)
        min_bin = h.GetMinimumBin()
        min_x = h.GetBinCenter(min_bin)
        h.GetXaxis().SetRangeUser(min_x, 500)
        return h

    def __get_mean(self, histo):
        h = self.__get_histo_without_pedestal(histo)
        h.GetXaxis().SetRangeUser(0, 30)
        min_bin = h.GetMinimumBin()
        min_x = h.GetBinCenter(min_bin)
        h.GetXaxis().SetRangeUser(min_x, 500)
        return [h.GetMean(), h.GetMeanError()]

    def __get_median(self, histo):
        h = self.__get_histo_without_pedestal(histo)
        integral = h.GetIntegral()
        median_i = 0
        for j in range(h.GetNbinsX() - 1):
            if integral[j] < 0.5:
                median_i = j
            else:
                break
        weight = (0.5 - integral[median_i]) / (integral[median_i + 1] - integral[median_i])
        median_x = h.GetBinCenter(median_i) + (h.GetBinCenter(median_i + 1) - h.GetBinCenter(median_i)) * weight
        return median_x

    def __get_mpv(self, histo):
        h = self.__get_histo_without_pedestal(histo)
        max_bin = h.GetMaximumBin()
        return h.GetBinCenter(max_bin)

    def draw_snrs(self):
        gr = self.make_tgrapherrors('gr', 'Signal to Noise Ratios')
        l1 = TLegend(.7, .68, .9, .9)
        l1.SetHeader('Regions')
        l2 = TLegend(.7, .47, .9, .67)
        l2.SetHeader('PeakIntegrals')
        for i, name in enumerate(self.get_all_signal_names().iterkeys()):
            snr = self.calc_snr(name)
            gr.SetPoint(i, i + 1, snr[0])
            gr.SetPointError(i, 0, snr[1])
        # rename bins
        for i, region in enumerate(self.get_all_signal_names().itervalues(), 1):
            bin_x = gr.GetXaxis().FindBin(i)
            gr.GetXaxis().SetBinLabel(bin_x, region)
        c = TCanvas('c', 'SNR', 1000, 1000)
        [l1.AddEntry(0, '{reg}:  {val}'.format(reg=reg, val=value), '') for reg, value in self.run.signal_regions.iteritems() if len(reg) < 2]
        [l2.AddEntry(0, '{reg}:  {val}'.format(reg=integ, val=value), '') for integ, value in self.run.peak_integrals.iteritems() if len(integ) < 2]
        self.format_histo(gr, y_tit='SNR', y_off=1.2, color=self.get_color())
        gr.SetLineColor(2)
        gr.Draw('bap')
        l1.Draw()
        l2.Draw()
        self.save_plots('SNR', sub_dir=self.save_dir)
        self.canvases[0] = c
        self.histos[0] = gr
        self.histos['legend'] = [l1, l2]

    def find_best_snr(self, show=True):
        gROOT.SetBatch(1)
        gr = self.make_tgrapherrors('gr', 'Signal to Noise Ratios')
        peak_integrals = OrderedDict(sorted({key:value for key, value in self.run.peak_integrals.iteritems() if len(key) < 3}.items()))
        for i, name in enumerate(peak_integrals.iterkeys()):
            signal = self.get_signal_names(self.get_signal_numbers('b', name))[self.channel]
            snr = self.calc_snr(signal)
            gr.SetPoint(i, i + 1, snr[0])
            gr.SetPointError(i, 0, snr[1])
        [gr.GetXaxis().SetBinLabel(gr.GetXaxis().FindBin(i), '{0}/{1}'.format(lst[0], lst[1])) for i, lst in enumerate(peak_integrals.itervalues(), 1)]
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        gr.Draw('ap')
        gROOT.SetBatch(0)
        self.save_plots('BestSNR', sub_dir=self.save_dir)
        self.histos[0] = [gr, c]

    def calc_snr(self, sig=None):
        signal = self.signal_name if sig is None else sig
        peak_int = self.get_all_signal_names()[signal][-1]
        ped_fit = self.show_pedestal_histo(draw=False, peak_int=peak_int)
        sig_fit = self.draw_pulse_height(eventwise_corr=True, draw=False, sig=signal)
        sig_mean = sig_fit.Parameter(0)
        ped_sigma = ped_fit.Parameter(2)

        snr = sig_mean / ped_sigma
        snr_err = snr * (sig_fit.ParError(0) / sig_mean + ped_fit.ParError(2) / ped_sigma)
        print 'SNR is: {snr} +- {err}'.format(snr=snr, err=snr_err)
        return [snr, snr_err]

    # ============================================
    # region MISCELLANEOUS
    def get_peak_position(self, event=None, region='b', peak_int='2'):
        num = self.get_signal_numbers(region=region, integral=peak_int)[self.channel]
        ev = self.start_event if event is None else event
        self.tree.GetEntry(ev)
        return self.tree.IntegralPeaks[num]

    def get_polarity(self):
        self.tree.GetEntry(0)
        return self.tree.polarities[self.channel]

    def get_all_signal_names(self):
        names = OrderedDict()
        regions = [reg for reg in self.run.signal_regions if len(reg) < 3]
        integrals = [integral for integral in self.run.peak_integrals if len(integral) < 3]
        for region in regions:
            for integral in integrals:
                if len(integral) > 2:
                    integral = '_' + integral
                name = 'ch{ch}_signal_{reg}_PeakIntegral{int}'.format(ch=self.channel, reg=region, int=integral)
                num = self.integral_names[name]
                reg = region + integral
                names['{pol}*IntegralValues[{num}]'.format(pol=self.polarity, num=num)] = reg
        return names

    @staticmethod
    def fit_fwhm(histo, fitfunc='gaus', do_fwhm=True, draw=False):
        h = histo
        if do_fwhm:
            peak_pos = h.GetBinCenter(h.GetMaximumBin())
            bin1 = h.FindFirstBinAbove(h.GetMaximum() / 2)
            bin2 = h.FindLastBinAbove(h.GetMaximum() / 2)
            fwhm = h.GetBinCenter(bin2) - h.GetBinCenter(bin1)
            option = 'qs' if draw else 'qs0'
            fit = h.Fit(fitfunc, option, '', peak_pos - fwhm / 2, peak_pos + fwhm / 2)
        else:
            fit = h.Fit(fitfunc, 'qs')
        return fit

    def __get_binning(self):
        jumps = self.cut.jump_ranges
        n_jumps = len(jumps['start'])
        bins = [0, self.GetMinEventCut()]
        ind = 0
        for start, stop in zip(jumps['start'], jumps['stop']):
            gap = stop - start
            # continue if first start and stop outside min event
            if stop < bins[-1]:
                ind += 1
                continue
            # if there is a jump from the start
            if start < bins[-1] < stop:
                bins.append(stop)
                ind += 1
                continue
            # add bins until hit interrupt
            while bins[-1] + self.bin_size < start:
                bins.append(bins[-1] + self.bin_size)
            # two jumps shortly after one another
            if ind < n_jumps - 2:
                next_start = jumps['start'][ind + 1]
                next_stop = jumps['stop'][ind + 1]
                if bins[-1] + self.bin_size + gap > next_start:
                    gap2 = next_stop - next_start
                    bins.append(bins[-1] + self.bin_size + gap + gap2)
                else:
                    bins.append(bins[-1] + self.bin_size + gap)
            else:
                bins.append(bins[-1] + self.bin_size + gap)
            ind += 1
        return bins

    def get_time_binning(self):
        time_bins = []
        for event in self.binning:
            time_bins.append(self.run.get_time_at_event(event))
        return time_bins

    def print_info_header(self):
        header = ['Run', 'Type', 'Diamond', 'HV [V]', 'Region']
        for info in header:
            print self.adj_length(info),
        print

    def print_information(self, header=True):
        if header:
            self.print_info_header()
        infos = [self.run_number, self.run.RunInfo['type'], self.diamond_name.ljust(4), self.bias, self.signal_region + self.peak_integral + '   ']
        for info in infos:
            print self.adj_length(info),
        print

    # endregion

    def __placeholder(self):
        pass


if __name__ == "__main__":
    st = time()
    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=392, type=int)
    parser.add_argument('ch', nargs='?', default=0, type=int)
    args = parser.parse_args()
    test_run = args.run
    print '\nAnalysing run', test_run, '\n'
    z = SignalAnalysis(test_run, args.ch)
    z.print_elapsed_time(st, 'Instantiation')
