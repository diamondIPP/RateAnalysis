# ==============================================
# IMPORTS
# ==============================================
from ROOT import TGraphErrors, TCanvas, TH2D, gStyle, TH1F, gROOT, TLegend, TCut, TGraph, TProfile2D, TH2F, TProfile, TCutG, kGreen, TF1, TPie, THStack
from TelescopeAnalysis import Analysis
from Elementary import Elementary
from CurrentInfo import Currents
from numpy import array, mean
from math import sqrt, ceil, log
from argparse import ArgumentParser
from Extrema import Extrema2D
from ChannelCut import ChannelCut
from time import time, sleep
from collections import OrderedDict
from sys import stdout
from copy import deepcopy

__author__ = 'micha'


# ==============================================
# MAIN CLASS
# ==============================================
class SignalAnalysis(Analysis):
    def __init__(self, run, channel, high_low_rate_run=None, binning=20000):

        self.channel = channel
        Analysis.__init__(self, run, high_low_rate=high_low_rate_run)

        # main
        self.diamond_name = self.run.diamond_names[channel]
        self.bias = self.run.bias[channel]
        self.save_dir = '{tc}_{run}_{dia}'.format(tc=self.TESTCAMPAIGN[2:], run=self.run_number, dia=self.diamond_name)

        # stuff
        self.BinSize = binning
        self.binning = self.__get_binning()
        self.time_binning = self.get_time_binning()
        self.n_bins = len(self.binning)
        self.Polarity = self.get_polarity()

        # names
        self.SignalName = self.get_signal_name(region=self.SignalRegion, peak_integral=self.PeakIntegral)
        self.SignalNumber = self.get_signal_number(self.SignalRegion, self.PeakIntegral)
        self.PedestalName = self.get_pedestal_name(region=self.PedestalRegion, peak_int=self.PeakIntegral)
        self.PulserName = self.get_pulser_name()

        # cuts
        self.Cut = ChannelCut(self, channel)

        # currents
        self.Currents = Currents(self)

        # graphs
        self.PulseHeight = None
        self.Pedestal = None
        # histograms
        self.SignalTime = None
        self.SignalMapHisto = None
        self.MeanSignalHisto = None
        self.PeakValues = None

    def __del__(self):
        for obj in [self.PulseHeight, self.Pedestal, self.SignalMapHisto, self.SignalTime, self.PeakValues, self.MeanSignalHisto]:
            self.del_rootobj(obj)
        for c in gROOT.GetListOfCanvases():
            c.Close()
        for lst in self.histos.itervalues():
            if not type(lst) is list:
                lst = [lst]
            for obj in lst:
                self.del_rootobj(obj)

    def show_current(self, relative_time=True):
        self.Currents.draw_graphs(relative_time=relative_time)

    # ==========================================================================
    # region INIT
    def get_polarity(self):
        self.tree.GetEntry(0)
        return self.tree.polarities[self.channel]

    def get_signal_number(self, region='b', peak_integral='2', sig_type='signal'):
        assert sig_type in ['signal', 'pedestal', 'pulser'], 'Invalid type of signal'
        if sig_type == 'signal':
            assert region in self.run.signal_regions, 'Invalid signal region {reg}!'.format(reg=region)
        elif sig_type == 'pedestal':
            assert region in self.run.pedestal_regions, 'Invalid pedestal region {reg}!'.format(reg=region)
        assert str(peak_integral) in self.run.peak_integrals, 'Invalid peak integral {reg}!'.format(reg=peak_integral)
        int_name = 'ch{ch}_{type}{reg}_PeakIntegral{int}'.format(ch=self.channel, reg='_' + region if region else '', int=peak_integral, type=sig_type)
        return self.IntegralNames[int_name]

    def get_signal_name(self, region='b', peak_integral='2', sig_type='signal'):
        num = self.get_signal_number(region, peak_integral, sig_type)
        return '{pol}*IntegralValues[{num}]'.format(pol=self.Polarity, num=num)

    def get_pedestal_name(self, region='ab', peak_int='2'):
        return self.get_signal_name(region=region, peak_integral=peak_int, sig_type='pedestal')

    def get_pulser_name(self):
        return self.get_signal_name(region='', sig_type='pulser')
    # endregion

    def set_channel(self, ch):
        self.channel = ch
        self.diamond_name = self.run.diamondname[ch]
        self.bias = self.run.bias[ch]
        self.Cut = ChannelCut(self, ch)
        self.save_dir = '{tc}_{run}_{dia}'.format(tc=self.TESTCAMPAIGN[2:], run=self.run_number, dia=self.run.diamondname[ch])
        self.Polarity = self.get_polarity()
        self.SignalName = self.get_signal_name()
        self.PedestalName = self.get_pedestal_name()

    def __set_bin_size(self, value):
        self.BinSize = value
        self.binning = self.__get_binning()
        self.time_binning = self.get_time_binning()
        self.n_bins = len(self.binning)
        return value

    # ==========================================================================
    # region BEAM PROFILE

    def draw_beam_profile(self, mode='x', show=True, fit=True, fit_margin=.6):
        assert mode.lower() in ['x', 'y'], 'Mode has to be either "x" or "y"!'
        margins = self.find_diamond_margins(show_plot=False, make_histo=True)
        h = deepcopy(self.histos[0])
        if not show:
            gROOT.SetBatch(1)
        prof = h.ProjectionX() if mode.lower() == 'x' else h.ProjectionY()
        margins[mode] = [prof.GetBinLowEdge(prof.FindBin(margins[mode][0])), prof.GetBinLowEdge(prof.FindBin(margins[mode][1]) + 1)]
        center = (margins[mode][1] + margins[mode][0]) / 2.
        width = (prof.FindBin(margins[mode][1]) - prof.FindBin(margins[mode][0])) / 2. * fit_margin * prof.GetBinWidth(1)
        fit_range = [center - width, center + width]
        c = TCanvas('c', 'Beam Profile', 1000, 1000)
        c.SetLeftMargin(.145)
        self.format_histo(prof, 'prof', 'Profile ' + mode.title(), y_tit='Entries', y_off=2, x_tit='Track Position {mod} [cm]'.format(mod=mode.title()))
        prof.GetXaxis().SetRangeUser(prof.GetBinCenter(prof.FindFirstBinAbove(0) - 1), prof.GetBinCenter(prof.FindLastBinAbove(0) + 1))
        prof.Draw()
        sleep(.1)
        lines = [self.make_tgaxis(x, c.GetUymin(), c.GetUymax(), '', 2, 2) for x in margins[mode]]
        fit_result = self.__fit_beam_profile(prof, fit_range, show) if fit else 0
        fits = None
        if fit:
            f1 = gROOT.GetFunction('gaus')
            f2 = deepcopy(f1)
            f2.SetLineColor(2)
            f2.SetLineStyle(1)
            f1.SetLineColor(kGreen + 1)
            f2.SetRange(fit_range[0], fit_range[1])
            f1.SetLineStyle(7)
            f1.Draw('same')
            f2.Draw('same')
            prof.GetXaxis().UnZoom()
            fits = [f1, f2]
        for line in lines:
            line.Draw()
        c.RedrawAxis()
        gROOT.SetBatch(0)
        self.save_plots('BeamProfile{mod}{fit}'.format(mod=mode.title(), fit='Fit' if fit else ''), sub_dir=self.save_dir)
        self.histos[1] = [prof, c, lines, fits]
        return fit_result if fit else prof

    @staticmethod
    def __fit_beam_profile(histo, fit_range, show=True):
        h = histo
        fit = h.Fit('gaus', 'qs{0}'.format('' if show else '0'), '', fit_range[0], fit_range[1])
        return fit

    def fit_beam_profile(self, mode='x', show=True, fit_margin=.6):
        pickle_path = self.PickleDir + 'BeamProfile/Fit{mod}_{tc}_{run}_{dia}_{mar}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name, mod=mode.title(), mar=fit_margin)

        def func():
            return self.draw_beam_profile(mode=mode, show=show, fit_margin=fit_margin)

        return self.do_pickle(pickle_path, func)

    def draw_beam_fit_properties(self, show=True, mode='x', sigma=True):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gr = self.make_tgrapherrors('gr', 'Beam Profile {0} {mod}'.format(mode.title(), mod='Fit #chi^{2}s / NDF' if not sigma else 'Sigma'))
        max_range = 11 if sigma else 10
        index = 0
        for i in xrange(1, max_range):
            perc = i / 10.
            fit = self.fit_beam_profile(mode=mode, show=False, fit_margin=perc)
            if fit.Ndf():
                y = fit.Parameter(2) if sigma else fit.Chi2() / fit.Ndf()
                gr.SetPoint(index, perc * 100, y)
                t = self.make_tlatex(perc * 100 - 2, y, str(fit.Ndf()), color=807, size=.04, align=32)
                gr.GetListOfFunctions().Add(t)
                index += 1
        c = TCanvas('c', 'Beam Chi2', 1000, 1000)
        self.format_histo(gr, x_tit='Range [%]', y_tit='#chi^{2} / NDF' if not sigma else 'Sigma', y_off=1.4)
        one = TF1('one', '1', 0, 100)
        t1 = self.make_tlatex(15, .95 * gr.GetYaxis().GetXmax(), 'NDF:', color=807, size=0.04, align=12)
        gr.GetListOfFunctions().Add(t1)
        gr.GetXaxis().SetRangeUser(-5, 105)
        gr.Draw('alp')
        one.Draw('same')

        self.histos[1] = [gr, c, t1]
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('BeamProf{mod}{dir}'.format(mod='Sigmas' if sigma else 'Chi2s', dir=mode.title()), sub_dir=self.save_dir)

    # endregion

    # ==========================================================================
    # region 2D SIGNAL DISTRIBUTION
    def draw_signal_map(self, draw_option='surf3z', show=True, factor=4):
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
        signal = '{sig}-{pol}*{ped}'.format(sig=self.SignalName, ped=self.PedestalName, pol=self.Polarity)
        print 'drawing signal map of {dia} for Run {run}...'.format(dia=self.diamond_name, run=self.run_number)
        self.tree.Draw('{z}:diam{nr}_track_x:diam{nr}_track_y>>signal_map'.format(z=signal, nr=nr), self.Cut.all_cut, 'goff')
        c = TCanvas('c', 'Signal Map', 1000, 1000)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.12)
        gStyle.SetPalette(53)
        self.format_histo(h, x_tit='track_x [cm]', y_tit='track_y [cm]', y_off=1.6)
        if draw_option.lower().startswith('surf'):
            self.format_histo(h, x_off=2, y_off=2.4, x_tit='track_x [cm]', y_tit='track_y [cm]')
        h.SetStats(0)
        h.SetContour(50)
        h.Draw(draw_option)
        self.save_plots('SignalMap2D_' + draw_option, sub_dir=self.save_dir)
        gROOT.SetBatch(0)
        self.SignalMapHisto = h
        self.canvases[0] = c
        return h

    def make_region_cut(self):
        self.draw_mean_signal_distribution(show=False)
        return self.Cut.generate_region(self.SignalMapHisto, self.MeanSignalHisto)

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
        gStyle.SetEndErrorSize(4)
        gr1 = self.make_tgrapherrors('gr', 'errors', width=3, marker_size=0, color=kGreen + 2)
        gr2 = self.make_tgrapherrors('gr', 'errors', width=3, marker_size=0, color=2)
        gr1.SetPoint(0, h.GetXaxis().GetXmin() + 5, h.GetMaximum() - 2)
        gr2.SetPoint(0, h.GetXaxis().GetXmin() + 5, h.GetMaximum() - 2)
        errors = self.SignalMapHisto.ProjectionXY('', 'c=e')
        gr1.SetPointError(0, errors.GetMinimum(), 0)
        gr2.SetPointError(0, errors.GetMaximum(), 0)
        l = self.make_tlatex(gr1.GetX()[0], gr1.GetY()[0] + 0.5, 'Errors', align=20, size=0.03)
        gr1.GetListOfFunctions().Add(l)
        if show:
            self.canvases[0] = TCanvas('c', 'Mean Signal Distribution', 1000, 1000)
            self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.2)
            h.Draw()
            gr2.Draw('[]')
            gr1.Draw('[]')
            gr2.Draw('p')
            gr1.Draw('p')
            self.save_plots('MeanSignalHisto', sub_dir=self.save_dir)
        self.MeanSignalHisto = h
        self.histos[0] = [gr1, gr2]

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
        pickle_path = self.PickleDir + 'MeanSignalFit/{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            self.draw_mean_signal_distribution(show=False)
            return self.MeanSignalHisto.Fit('gaus', 'qs')

        fit = self.do_pickle(pickle_path, func)
        return fit

    def get_mean_fwhm(self):
        fit = self.fit_mean_signal_distribution()
        conversion_factor = 2 * sqrt(2 * log(2))  # sigma to FWHM
        return fit.Parameter(2) * conversion_factor

    def draw_diamond_hitmap(self, cut=None, show_frame=True):
        self.find_diamond_margins(show_frame=show_frame, cut=cut)

    def find_diamond_margins(self, show_plot=True, show_frame=False, cut=None, make_histo=False):
        pickle_path = self.PickleDir + 'Margins/{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            print 'getting margins for {dia} of run {run}...'.format(dia=self.diamond_name, run=self.run_number)
            cut_string = self.Cut.all_cut if cut is None else cut
            if not show_plot:
                gROOT.SetBatch(1)
            h = TH2F('h', 'Diamond Margins', 52, -.3, .3, 80, -.3, .3)
            nr = 1 if not self.channel else 2
            self.tree.Draw('diam{nr}_track_x:diam{nr}_track_y>>h'.format(nr=nr), cut_string, 'goff')
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
                    if bin_ > 1:
                        if efficiency and not last_bin:
                            zero_bins[i].append(proj.GetBinCenter(bin_ - 1))
                        elif not efficiency and last_bin:
                            zero_bins[i].append((proj.GetBinCenter(bin_)))
                    last_bin = proj.GetBinContent(bin_)
            if show_plot:
                c = TCanvas('c', 'Diamond Hit Map', 1000, 1000)
                c.SetRightMargin(.14)
                c.SetBottomMargin(.15)
                h.GetXaxis().SetRangeUser(zero_bins[0][0], zero_bins[0][1])
                h.GetYaxis().SetRangeUser(zero_bins[1][0], zero_bins[1][1])
                h.SetStats(0)
                h.Draw('colz')
                if show_frame:
                    self.__show_frame(bin_low, bin_high)
                self.save_plots('DiamondHitmap', sub_dir=self.save_dir)
            self.histos[0] = h
            gROOT.SetBatch(0)
            return {name: [efficient_bins[i][0], efficient_bins[i][-1]] for i, name in enumerate(['x', 'y'])}

        margins = func() if show_plot or make_histo else 0
        return self.do_pickle(pickle_path, func, margins)

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
    # region SIGNAL PEAK POSITION
    def draw_peak_values(self, region=None, type_='signal', draw=True, ucut=None):
        num = self.get_signal_number('b', '2') if region is None else self.get_signal_number(region=region, sig_type=type_)
        region = 'b' if region is None else region
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        peak_val = 'IntegralPeaks[{num}]'.format(num=num)
        title = 'Peak Values {reg}'.format(reg=region)
        x = self.run.signal_regions[region] if type_ == 'signal' else self.run.get_regions('pulser')['pulser']
        h = TH1F('peakvalues', title, x[1] - x[0], x[0] / 2., x[1] / 2.)
        self.format_histo(h, x_tit='time [ns]', y_tit='Entries', y_off=2, fill_color=17)
        cut = self.Cut.all_cut if type_ == 'signal' else '!({0})'.format(self.Cut.CutStrings['pulser'])
        cut = cut if ucut is None else ucut
        self.tree.Draw(peak_val + '/2.>>peakvalues', cut, 'goff')
        if draw:
            c = TCanvas('c', 'Signal Peak Distribution', 1000, 1000)
            c.SetLeftMargin(0.14)
            h.Draw()
            self.save_plots('peak_values_{reg}{int}'.format(reg=region, int=self.PeakIntegral), 'png', canvas=c, sub_dir=self.save_dir)
            self.histos[0] = c
        self.PeakValues = h
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')

    def fit_peak_values(self, draw=True, pulser=False):
        pickle_path = self.PickleDir + 'PeakValues/Fit_{tc}_{run}_{dia}{pul}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name, pul='_pulser' if pulser else '')

        def func():
            print 'Getting peak value fit for {dia} of run {run}...'.format(run=self.run_number, dia=self.diamond_name)
            self.draw_peak_values(draw=draw) if not pulser else self.draw_pulser_peakvalues(draw=draw)
            h = self.PeakValues
            max_bin = h.GetMaximumBin()
            x = [h.GetBinCenter(max_bin + i) for i in [-7, 1]] if not pulser else [h.GetXaxis().GetXmin() + 1, h.GetXaxis().GetXmax() - 1]
            return h.Fit('gaus', 'qs{0}'.format('' if draw else '0'), '', x[0], x[1])

        mean_val = func() if draw else 0
        return self.do_pickle(pickle_path, func, mean_val)

    def calc_peak_value_fwhm(self):
        pickle_path = self.PickleDir + 'PeakValues/FWHM_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

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
    def __generate_signal_name(self, signal, evnt_corr, off_corr, bin_corr, cut=None):
        sig_name = signal
        if bin_corr:
            return sig_name
        elif off_corr:
            ped_fit = self.show_pedestal_histo(cut=cut, draw=False)
            sig_name += '-{0}'.format(ped_fit.Parameter(1))
        elif evnt_corr:
            sig_name += '-{ped}'.format(sig=sig_name, ped=self.PedestalName, pol=self.Polarity)
        return sig_name

    def make_signal_time_histos(self, ped=False, signal=None, evnt_corr=False, off_corr=False, show=True, bin_corr=False):
        gROOT.SetBatch(1)
        signal = self.SignalName if signal is None else signal
        signal = signal if not ped else self.PedestalName
        signal = self.__generate_signal_name(signal, evnt_corr, off_corr, bin_corr)
        # 2D Histogram
        name = "signaltime_" + str(self.run_number)
        xbins = array(self.time_binning)
        x_min = -50 if not ped else -20
        x_max = 300 if not ped else 20
        bins = 1000 if not ped else 80
        h = TH2D(name, "signaltime", len(xbins) - 1, xbins, bins, x_min, x_max)
        self.tree.Draw("{name}:time>>{histo}".format(histo=name, name=signal), self.Cut.all_cut, 'goff')
        if show:
            gROOT.SetBatch(0)
            c = TCanvas('c', 'Pulse Height vs Time', 1000, 1000)
            c.SetLeftMargin(.12)
            self.format_histo(h, x_tit='time [ms]', y_tit='Pulse Height [au]', y_off=1.4)
            h.Draw('colz')
            self.save_plots('SignalTime', sub_dir=self.save_dir)
            self.SignalTime = h
            self.canvases[0] = c
        gROOT.SetBatch(0)
        return h

    def draw_pedestal(self, binning=None, draw=True):
        bin_size = binning if binning is not None else self.BinSize
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

    def draw_pulse_height(self, binning=None, show=True, save_graph=False, evnt_corr=True, bin_corr=False, off_corr=False, sig=None):
        signal = self.SignalName if sig is None else sig
        bin_size = binning if binning is not None else self.BinSize
        correction = ''
        if bin_corr:
            correction = 'binwise'
        elif off_corr:
            correction = 'constant'
        elif evnt_corr:
            correction = 'eventwise'
        suffix = '{bins}_{cor}_{sig}'.format(bins=bin_size, cor=correction, sig=self.get_all_signal_names()[signal])
        picklepath = 'Configuration/Individual_Configs/Ph_fit/{tc}_{run}_{ch}_{suf}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, ch=self.channel, suf=suffix)

        self.SignalTime = None

        def func():
            print 'drawing pulse height fit for run {run} and {dia}...'.format(run=self.run_number, dia=self.diamond_name)
            if binning is not None:
                self.__set_bin_size(binning)
            tit_suffix = 'with {cor} Pedestal Correction'.format(cor=correction.title()) if bin_corr or evnt_corr or off_corr else ''
            gr = self.make_tgrapherrors('signal', 'Pulse Height Evolution Bin{0} '.format(self.BinSize) + tit_suffix)
            sig_time = self.make_signal_time_histos(evnt_corr=evnt_corr, signal=signal, show=False, off_corr=off_corr, bin_corr=bin_corr)
            mode = 'mean'
            empty_bins = 0
            count = 0
            means = self.draw_pedestal(bin_size, draw=False) if bin_corr else None
            gROOT.SetBatch(1)
            for i in xrange(self.n_bins - 1):
                h_proj = sig_time.ProjectionY(str(i), i + 1, i + 1)
                if h_proj.GetEntries() > 10:
                    if mode in ["mean", "Mean"]:
                        i_mean = h_proj.GetMean()
                        i_mean -= means[count] if bin_corr else 0
                        gr.SetPoint(count, (self.time_binning[i] - self.run.startTime) / 60e3, i_mean)
                        gr.SetPointError(count, 0, h_proj.GetRMS() / sqrt(h_proj.GetEntries()))
                        count += 1
                    elif mode in ["fit", "Fit"]:
                        h_proj.GetMaximum()
                        maxposition = h_proj.GetBinCenter(h_proj.GetMaximumBin())
                        h_proj.Fit("landau", "Q", "", maxposition - 50, maxposition + 50)
                        fitfun = h_proj.GetFunction("landau")
                        mpv = fitfun.GetParameter(1)
                        mpverr = fitfun.GetParError(1)
                        gr.SetPoint(count, (i + 0.5) * self.run.totalMinutes / self.n_bins, mpv)
                        gr.SetPointError(count, 0, mpverr)
                else:
                    empty_bins += 1
            if empty_bins:
                print 'Empty proj. bins:\t', str(empty_bins) + '/' + str(self.n_bins)
            if show:
                gROOT.SetBatch(0)
            c = TCanvas('bla', 'blub', 1000, 1000)
            c.SetLeftMargin(.14)
            gStyle.SetOptFit(1)
            self.format_histo(gr, x_tit='time [min]', y_tit='Mean Pulse Height [au]', y_off=1.6)
            fit_par = gr.Fit('pol0', 'qs')
            gr.Draw('apl')
            self.save_plots('PulseHeight{0}'.format(self.BinSize), sub_dir=self.save_dir)
            self.PulseHeight = gr
            self.canvas = c
            gROOT.SetBatch(0)
            return fit_par

        fit = func() if show or save_graph else 0
        return self.do_pickle(picklepath, func, fit)

    def draw_ph_distribution(self, binning=None, show=True, fit=True):
        if binning is not None:
            self.__set_bin_size(binning)
        sig_time = self.make_signal_time_histos(evnt_corr=True, show=False)
        if not show:
            gROOT.SetBatch(1)
        means = [h_proj.GetMean() for h_proj in [sig_time.ProjectionY(str(i), i + 1, i + 1) for i in xrange(self.n_bins - 1)] if h_proj.GetEntries() > 10]
        extrema = [int(min(means)), int(max(means))]
        h = TH1F('h', 'Signal Bin{0} Distribution'.format(self.BinSize), int(log(len(means), 2) * 2), extrema[0], extrema[1] + 2)
        for mean_ in means:
            h.Fill(mean_)
        c = TCanvas('c', 'Pulse Height Distribution', 1000, 1000)
        c.SetLeftMargin(.12)
        self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.5)
        if fit:
            h.Fit('gaus', 'q')
        h.SetFillColor(kGreen - 9)
        h.Draw()
        gROOT.SetBatch(0)
        self.save_plots('SignalBin{0}Disto'.format(self.BinSize), sub_dir=self.save_dir)
        self.histos[0] = [h, c]
        return h

    def show_ph_overview(self, binning=None):
        self.draw_pulse_height(binning=binning, show=False, save_graph=True)
        h1 = self.PulseHeight
        self.format_histo(h1, y_off=1.4)
        h2 = self.draw_ph_distribution(binning=binning, show=False)
        print h1, h2
        c = TCanvas('c', 'Pulse Height Distribution', 1500, 750)
        c.Divide(2, 1)
        for i, h in enumerate([h1, h2], 1):
            pad = c.cd(i)
            pad.SetBottomMargin(.15)
            h.Draw()
        self.save_plots('PHEvolutionOverview{0}'.format(self.BinSize), sub_dir=self.save_dir)
        self.histos[1] = [c]

    def show_signal_histo(self, cut=None, evnt_corr=True, off_corr=False, show=True, sig=None, binning=350, events=None, start=None):
        print 'drawing signal distribution for run {run} and {dia}...'.format(run=self.run_number, dia=self.diamond_name)
        suffix = 'with Pedestal Correction' if evnt_corr else ''
        h = TH1F('signal b2', 'Pulse Height ' + suffix, binning, -50, 300)
        cut = self.Cut.all_cut if cut is None else cut
        sig_name = self.SignalName if sig is None else sig
        sig_name = self.__generate_signal_name(sig_name, evnt_corr, off_corr, False, cut)
        start_event = int(float(start)) if start is not None else 0
        n_events = self.find_n_events(n_events=events, cut=str(cut), start=start_event) if events is not None else self.run.n_entries
        self.tree.Draw('{name}>>signal b2'.format(name=sig_name), str(cut), 'goff', n_events, start_event)
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Signal Distribution', 1000, 1000)
        c.SetLeftMargin(.13)
        self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.8)
        h.Draw()
        self.save_plots('SignalDistribution', sub_dir=self.save_dir)
        self.histos[0] = [c, h]
        gROOT.SetBatch(0)
        return h

    def draw_signal_vs_peakpos(self, show=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gr = self.make_tgrapherrors('gr', 'Signal vs Peak Position')
        i = 0
        for peak_pos in xrange(120, 151):
            print '\rcalculating peak pos: {0:03d}'.format(peak_pos),
            self.Cut.add_signal_peak_pos_cut([peak_pos, peak_pos + 1])
            events = self.tree.Draw('1', self.Cut.all_cut, 'goff')
            print '({0:5d})'.format(events),
            stdout.flush()
            if events > 500:
                ph_fit = self.draw_pulse_height(show=False, save_graph=True)
                print ph_fit.Parameter(0)
                gr.SetPoint(i, peak_pos, ph_fit.Parameter(0))
                gr.SetPointError(i, 0, ph_fit.ParError(0))
                i += 1
        gStyle.SetPalette(55)
        gROOT.SetBatch(1) if not show else self.do_nothing()
        c = TCanvas('c', 'Signal Distribution', 1000, 1000)
        c.SetLeftMargin(.11)
        self.format_histo(gr, x_tit='Signal Peak Position', y_tit='Pulse Height [au]', y_off=1.2)
        gr.Draw('alp')
        self.save_plots('SignalVsPeakPos', sub_dir=self.save_dir)
        self.histos[0] = [c, gr]
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')

    def draw_landau_vs_peakpos(self, show=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        hs = THStack('lpp', 'Landau vs. Signal Peak Postion;pulse height;#')
        for peak_pos in xrange(120, 151):
            print '\rcalculating peak pos: {0:03d}'.format(peak_pos),
            self.Cut.add_signal_peak_pos_cut([peak_pos, peak_pos + 1])
            events = self.tree.Draw('1', self.Cut.all_cut, 'goff')
            print '({0:05d})'.format(events),
            stdout.flush()
            if events > 500:
                h = self.show_signal_histo(show=False)
                h.SetLineColor(self.get_color())
                h.Scale(1/h.GetMaximum())
                hs.Add(h)
        gStyle.SetPalette(55)
        self.format_histo(hs, y_tit='Pulse Height [au]', y_off=1.2)
        self.histos.append(self.draw_histo(hs, 'SignalVsTriggerCell', show, self.save_dir, lm=.11, draw_opt='nostack'))
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
    def show_pedestal_histo(self, region='ab', peak_int='2', cut=None, fwhm=True, draw=True):
        cut = self.Cut.all_cut if cut is None else cut
        cut = TCut('', cut) if type(cut) is str else cut
        fw = 'fwhm' if fwhm else 'full'
        suffix = '{reg}_{fwhm}_{cut}'.format(reg=region + str(peak_int), cut=cut.GetName(), fwhm=fw)
        picklepath = 'Configuration/Individual_Configs/Pedestal/{tc}_{run}_{ch}_{suf}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, ch=self.channel, suf=suffix)

        def func():
            gROOT.SetBatch(1)
            print 'making pedestal histo for region {reg}{int}...'.format(reg=region, int=peak_int)
            h = TH1F('ped1', 'Pedestal Distribution', 100, -20, 20)
            name = self.get_pedestal_name(region, peak_int)
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
            self.histos[0] = [h, c]
            gROOT.SetBatch(0)
            return fit_pars

        fit_par = func() if draw else 0
        return self.do_pickle(picklepath, func, fit_par)

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
        c = TCanvas('bla', 'blub', 1000, 1000)
        gr1.Draw('alp')
        gr2.Draw('lp')
        gr3.Draw('lp')
        legend.AddEntry(gr1, 'mean fit fwhm w/ cuts 2', 'lp')
        legend.AddEntry(gr2, 'mean fit fwhm w/ cuts median', 'lp')
        legend.AddEntry(gr3, 'mean fit fwhm w/ cuts all', 'lp')
        legend.Draw()
        self.histos[0] = [gr1, gr2, gr3, c, legend]

    # endregion

    # ==========================================================================
    # region CUTS
    def show_cut_contributions(self):
        main_cut = [self.Cut.CutStrings['event_range'], self.Cut.CutStrings['beam_interruptions']]
        contributions = {}
        cutted_events = 0
        cuts = TCut('consecutive', '')
        for cut in main_cut + self.Cut.CutStrings.values():
            name = cut.GetName()
            if not name.startswith('old') and name != 'all_cuts' and name not in contributions and str(cut):
                cuts += cut
                events = int(z.tree.Draw('1', '!({0})'.format(cuts), 'goff'))
                events -= cutted_events
                print name, events
                contributions[cut.GetName()] = events
                cutted_events += events
        sorted_contr = OrderedDict()
        while contributions:
            for key, value in contributions.iteritems():
                if value == max(contributions.values()):
                    sorted_contr[key] = value
                    contributions.pop(key)
                    break
            for key, value in contributions.iteritems():
                if value == min(contributions.values()):
                    sorted_contr[key] = value
                    contributions.pop(key)
                    break
        contributions = sorted_contr
        values = contributions.values() + [self.run.n_entries - cutted_events]
        i = 0
        colors = [self.get_color() for i in xrange(1, len(values) + 1)]
        print values, i
        pie = TPie('pie', 'Cut Contributions', len(values), array(values, 'f'), array(colors, 'i'))
        for i, label in enumerate(contributions.iterkeys()):
            # if pie.GetEntryVal(i) < 2000:
            pie.SetEntryRadiusOffset(i, .05)
            pie.SetEntryLabel(i, label.title())
        pie.SetEntryRadiusOffset(i + 1, .05)
        pie.SetEntryLabel(i + 1, 'Good Events')
        pie.SetHeight(.04)
        pie.SetRadius(.2)
        pie.SetTextSize(.025)
        pie.SetAngle3D(70)
        pie.SetLabelFormat('#splitline{%txt}{%perc}')
        pie.SetAngularOffset(240)
        c = TCanvas('c', 'Cut Pie', 1000, 1000)
        pie.Draw('3drsc')
        self.save_plots('CutContributions', sub_dir=self.save_dir)
        self.histos[0] = [pie, c]
        return contributions

    def show_bucket_histos(self):
        h = TH1F('h', 'Bucket Cut Histograms', 250, -50, 300)
        self.tree.Draw('{name}>>h'.format(name=self.SignalName), '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['old_bucket'], pul=self.Cut.CutStrings['pulser']), 'goff')
        h1 = deepcopy(h)
        fit = self.Cut.triple_gauss_fit(h1, show=False)
        sig_fit = TF1('f1', 'gaus', -50, 300)
        sig_fit.SetParameters(fit.GetParameters())
        ped1_fit = TF1('f2', 'gaus', -50, 300)
        ped2_fit = TF1('f2', 'gaus', -50, 300)
        ped1_fit.SetParameters(*[fit.GetParameter(i) for i in xrange(3, 6)])
        ped2_fit.SetParameters(*[fit.GetParameter(i) for i in xrange(6, 9)])
        h_sig = deepcopy(h)
        h_ped1 = deepcopy(h)
        h_ped2 = deepcopy(h)
        h_sig.Add(ped1_fit, -1)
        h_sig.Add(ped2_fit, -1)
        h_ped1.Add(ped2_fit, -1)
        h_ped2.Add(ped1_fit, -1)
        h_ped1.Add(h_sig, -1)
        h_ped2.Add(h_sig, -1)
        c = TCanvas('c', 'Bucket Histos', 1000, 1000)
        for i, h in enumerate([h_ped1, h_ped2, h_sig]):
            h.SetStats(0)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw('same') if i else h.Draw()
        self.save_plots('BucketHistos', sub_dir=self.save_dir)
        self.histos[0] = [h, h_sig, h_ped1, h_ped2, c]

    def show_bucket_numbers(self, show=True):
        pickle_path = self.PickleDir + 'Cuts/BucketEvents_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            print 'getting number of bucket events for run {run} and {dia}...'.format(run=self.run_number, dia=self.diamond_name)
            n_new = self.tree.Draw('1', '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['bucket'], pul=self.Cut.CutStrings['pulser']), 'goff')
            n_old = self.tree.Draw('1', '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['old_bucket'], pul=self.Cut.CutStrings['pulser']), 'goff')
            if show:
                print 'New Bucket: {0} / {1} = {2:4.2f}%'.format(n_new, self.run.n_entries, n_new / float(self.run.n_entries) * 100)
                print 'Old Bucket: {0} / {1} = {2:4.2f}%'.format(n_old, self.run.n_entries, n_old / float(self.run.n_entries) * 100)
            return {'old': n_old, 'new': n_new, 'all': float(self.run.n_entries)}

        return self.do_pickle(pickle_path, func)

    def show_bucket_hits(self, show=True):
        # hit position
        h = TH2F('h', 'Diamond Margins', 80, -.3, .3, 52, -.3, .3)
        nr = 1 if not self.channel else 2
        cut = '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['old_bucket'], pul=self.Cut.CutStrings['pulser'])
        self.tree.Draw('diam{nr}_track_x:diam{nr}_track_y>>h'.format(nr=nr), cut, 'goff')
        projections = [h.ProjectionX(), h.ProjectionY()]
        zero_bins = [[], []]
        for i, proj in enumerate(projections):
            last_bin = None
            for bin_ in xrange(proj.GetNbinsX()):
                efficiency = proj.GetBinContent(bin_) / float(proj.GetMaximum())
                if bin_ > 1:
                    if efficiency > .05 and last_bin < 5:
                        zero_bins[i].append(proj.GetBinCenter(bin_ - 1))
                    elif efficiency < .05 and last_bin > 5:
                        zero_bins[i].append((proj.GetBinCenter(bin_)))
                last_bin = proj.GetBinContent(bin_)
        if show:
            print zero_bins
            c = TCanvas('c', 'Diamond Hit Map', 1000, 1000)
            h.GetXaxis().SetRangeUser(zero_bins[0][0], zero_bins[0][-1])
            h.GetYaxis().SetRangeUser(zero_bins[1][0], zero_bins[1][-1])
            h.Draw('colz')
            self.histos[0] = [c, h]
        return h

    def show_bucket_means(self, show=True, plot_histos=True):
        pickle_path = self.PickleDir + 'Cuts/BucketMeans_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name)

        def func():
            gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
            cuts_nobucket = TCut('no_bucket', '')
            cuts_oldbucket = TCut('old_bucket', '')
            for key, value in self.Cut.CutStrings.iteritems():
                if not key.startswith('old') and key not in ['all_cuts', 'bucket']:
                    cuts_nobucket += value
                if key not in ['all_cuts', 'bucket']:
                    cuts_oldbucket += value
            h1 = self.show_signal_histo(show=False, evnt_corr=True)
            h2 = self.show_signal_histo(show=False, evnt_corr=True, cut=cuts_nobucket)
            h3 = self.show_signal_histo(show=False, evnt_corr=True, cut=cuts_oldbucket)
            if plot_histos:
                c = TCanvas('c', 'Bucket Histos', 1000, 1000)
                self.format_histo(h1, color=self.get_color(), lw=1, x_tit='Pulse Height [au]', y_tit='Entries')
                h1.Draw()
                self.format_histo(h2, color=self.get_color(), lw=1)
                h2.Draw('same')
                self.format_histo(h3, color=self.get_color(), lw=1)
                h3.Draw('same')
                self.histos[0] = [h1, h2, h3, c]
            result = {name: [h.GetMean(), h.GetMeanError()] for name, h in zip(['new', 'no', 'old'], [h1, h2, h3])}
            gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
            if show:
                print result
            return result

        res = func() if plot_histos else 0
        return self.do_pickle(pickle_path, func, res)

    def compare_single_cuts(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        c1 = TCanvas('single', '', 1000, 1000)
        c2 = TCanvas('all', '', 1000, 1000)
        c2.SetLeftMargin(0.15)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        histos = []
        drawn_first = False
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                print 'saving plot', key
                save_name = 'signal_distribution_{cut}'.format(cut=key)
                histo_name = 'signal {range}{peakint}'.format(range=self.SignalRegion, peakint=self.PeakIntegral)
                histo_title = 'signal with cut ' + key
                histo = TH1F(histo_name, histo_title, 350, -50, 300)
                # safe single plots
                c1.cd()
                self.tree.Draw("{name}>>{histo}".format(name=self.SignalName, histo=histo_name), value)
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
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                print 'saving plot', key
                save_name = 'signal_distribution_normalised_{cut}'.format(cut=key)
                histo_name = 'signal {range}{peakint}'.format(range=self.SignalRegion, peakint=self.PeakIntegral)
                histo_title = 'normalised signal with cut ' + key
                histo = TH1F(histo_name, histo_title, 350, -50, 300)
                # safe single plots
                c1.cd()
                self.tree.Draw("{name}>>{histo}".format(name=self.SignalName, histo=histo_name), value)
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
        for key, value in self.Cut.CutStrings.iteritems():
            if (str(value) or key == 'raw') and key != 'all_cuts':
                cut += value
                print 'saving plot with {n} cuts'.format(n=ind)
                save_name = 'signal_distribution_{n}cuts'.format(n=ind)
                histo_name = 'signal {range}{peakint}'.format(range=self.SignalRegion, peakint=self.PeakIntegral)
                histo_title = 'signal with {n} cuts'.format(n=ind)
                histo = TH1F(histo_name, histo_title, 550, -50, 500)
                # safe single plots
                c1.cd()
                self.tree.Draw("{name}>>{histo}".format(name=self.SignalName, histo=histo_name), cut)
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
    # region PULSER
    def draw_pulser_rate(self, binning=200, cut=None, show=True):
        """
        Shows the fraction of accepted pulser events as a function of event numbers. Peaks appearing in this graph are most likely beam interruptions.
        :param binning:
        :param cut:
        :param show:
        """
        gROOT.SetBatch(1)
        cut = '' if cut is None else cut
        nbins = self.run.n_entries / binning
        h = TProfile('h', 'Pulser Rate', nbins, 0, self.run.n_entries)
        self.tree.Draw('(pulser!=0)*100:Entry$>>h', cut, 'goff')
        gROOT.SetBatch(0) if show else self.do_nothing()
        c = TCanvas('c', 'Pulser Rate Canvas', 1000, 1000)
        self.format_histo(h, name='pulser_rate', title='Pulser Rate', x_tit='Event Number', y_tit='Pulser Fraction [%]', y_off=1.3)
        h.Draw('hist')
        self.save_plots('pulser_rate', canvas=c, sub_dir=self.save_dir)
        self.canvases[0] = c
        self.histos[0] = h
        gROOT.SetBatch(0)
        return h

    def fit_pulser_rate(self, binning=2000, show=True):
        cut = self.Cut.CutStrings['event_range'] + self.Cut.CutStrings['beam_interruptions']
        h = self.draw_pulser_rate(show=show, cut=cut, binning=binning)
        gStyle.SetOptFit(1)
        fit = h.Fit('pol0', 'qs{0}'.format('' if show else 0))
        print 'Pulser Rate: {0} +- {1}'.format(fit.Parameter(0), fit.ParError(0))
        return fit

    def draw_pulser_pulseheight(self, binning=20000, draw_opt='hist'):
        """
        Shows the average pulse height as a function of event numbers.
        :param binning:
        :param draw_opt:
        """
        nbins = self.run.n_entries / binning
        h = TProfile('h', 'Pulser Pulse Height', nbins, 0, self.run.n_entries)
        cut = self.Cut.generate_pulser_cut()
        signal = self.__generate_signal_name(self.PulserName, False, True, False, cut)
        print signal
        self.tree.Draw('{0}:Entry$>>h'.format(signal), cut, 'goff')
        c = TCanvas('c', 'Pulser Rate Canvas', 1000, 1000)
        c.SetGridy()
        self.format_histo(h, name='h', title='Pulser Pulse Height', x_tit='Event Number', y_tit='Pulse Height [au]', y_off=1.3)
        gStyle.SetStatY(.65)
        gStyle.SetStatX(.65)
        h.Draw(draw_opt)
        self.run.draw_run_info(canvas=c, channel=self.channel)
        self.save_plots('PulserPulseHeight', sub_dir=self.save_dir)
        self.canvases[0] = c
        self.histos[0] = h
        return h

    def fit_pulser_pulseheight(self):
        h = self.draw_pulser_pulseheight(draw_opt='')
        gStyle.SetOptFit()
        fit = h.Fit('pol0', 'qs')
        return fit

    def show_pulser_histo(self, show=True, corr=True, beam_on=True, binning=700, events=None, start=None):
        cut = self.Cut.generate_pulser_cut(beam_on)
        h = self.show_signal_histo(cut=cut, sig=self.PulserName, show=show, off_corr=corr, evnt_corr=False, binning=binning, events=events, start=start)
        c = gROOT.GetListOfCanvases()[-1]
        c.SetLogy()
        c.Update()
        return h

    def calc_pulser_fit(self, show=True, corr=True, beam_on=True, events=None, start=None, binning=350):
        start_string = '_{0}'.format(start) if start is not None else ''
        events_string = '_{0}'.format(events) if events is not None else ''
        suffix = '{corr}_{beam}{st}{ev}'.format(corr='_ped_corr' if corr else '', beam='BeamOff' if not beam_on else 'BeamOn', st=start_string, ev=events_string)
        pickle_path = self.PickleDir + 'Pulser/HistoFit_{tc}_{run}_{dia}{suf}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, dia=self.diamond_name, suf=suffix)

        def func():
            gStyle.SetOptFit(1)
            h = self.show_pulser_histo(show=show, corr=corr, beam_on=beam_on, binning=binning, events=events, start=start)
            fit_func = h.Fit('gaus', 'qs{0}'.format('' if show else '0'), '', 0, h.GetBinCenter(h.GetMaximumBin() + 2))
            f = gROOT.GetFunction('gaus')
            f.SetLineStyle(7)
            f.SetRange(0, 500)
            f.Draw('same')
            self.histos[1] = [f]
            self.save_plots('PulserHisto', sub_dir=self.save_dir)
            return fit_func

        fit = func() if show else 0
        return self.do_pickle(pickle_path, func, fit)

    def draw_pulser_peakvalues(self, draw=True):
        self.draw_peak_values('', 'pulser', ucut=self.Cut.generate_pulser_cut(), draw=draw)

    def show_pulser_pedestal(self, show=True):
        return z.show_pedestal_histo(cut=self.Cut.generate_pulser_cut(), draw=show)

    def compare_pulser_pedestal(self, show=True):
        self.show_pedestal_histo()
        h1 = deepcopy(self.histos[0][0])
        self.show_pulser_pedestal()
        h2 = self.histos[0][0]
        print h1, h2
        gROOT.SetBatch(0) if show else gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Pedestal Comparison', 1000, 1000)
        legend = TLegend(.7, .78, .88, .88)
        names = ['Signal', 'Pulser']
        for i, h in enumerate([h1, h2]):
            self.format_histo(h, color=self.get_color(), lw=2)
            h.SetStats(0)
            h.Scale(1 / h.GetMaximum())
            h.GetYaxis().SetRangeUser(0, 1.1)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, names[i])
        legend.Draw()
        gROOT.SetBatch(0)
        self.save_plots('PulserPedestalComparison', sub_dir=self.save_dir)
        self.histos[1] = [c, h1, h2, legend]

    def draw_pulser_waveform(self, n=1, start_event=None, add_buckets=False, cut=None, fixed_range=None):
        cut = self.Cut.generate_pulser_cut() if cut is None else cut
        start = self.StartEvent + self.count if start_event is None else start_event + self.count
        print 'Start at event number:', start
        cnt = self.draw_waveforms(n=n, start_event=start, add_buckets=add_buckets, cut_string=cut, ret_event=True, fixed_range=fixed_range)
        print cnt
        if cnt is None:
            return
        self.count += cnt

    def save_pulser_shapes(self, n_pics=5, fixed_range=None):
        frange = [-100, 50] if fixed_range is None else fixed_range
        events_spacing = (self.EndEvent - self.StartEvent) / n_pics
        start_events = [self.StartEvent + events_spacing * i for i in xrange(n_pics)]
        for i, start in enumerate(start_events):
            self.count = 0
            self.draw_pulser_waveform(n=1000, start_event=start, fixed_range=frange)
            self.save_plots('WaveForms{0}'.format(i), sub_dir='{0}/WaveForms'.format(self.save_dir))

    def draw_pulser_vs_time(self, n_points=5, _mean=True, show=True, corr=True, events=5000):
        events_spacing = (self.EndEvent - self.StartEvent) / n_points
        start_events = [self.StartEvent + events_spacing * i for i in xrange(n_points)]
        mode = 'Mean' if _mean else 'Sigma'
        gr = self.make_tgrapherrors('gr', '{0} of Pulser vs Time'.format(mode))
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        for i, start in enumerate(start_events):
            fit = self.calc_pulser_fit(show=False, start=start, events=events, binning=200, corr=corr)
            par = 1 if _mean else 2
            gr.SetPoint(i, (self.run.get_time_at_event(start) - self.run.startTime) / 60e3, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
        gROOT.SetBatch(0) if show else gROOT.SetBatch(1)
        c = TCanvas('c', '{0} of Pulser vs Time'.format(mode), 1000, 1000)
        gr.Draw('alp')
        self.save_plots('Pulser{0}VsTime'.format(mode), sub_dir=self.save_dir)
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.histos[2] = [c, gr]

    # endregion

    # ==========================================================================
    # region SHOW
    def draw_bucket_pedestal(self, show=True):
        gROOT.SetBatch(1)
        reg_name = 'e2'
        three_bucket_num = self.get_signal_number(reg_name[0], reg_name[1])
        reg_margins = self.run.signal_regions[reg_name[0]]
        x_bins = (reg_margins[1] - reg_margins[0])
        h = TH2F('h', 'Bucket Pedestal', x_bins, reg_margins[0] / 2., reg_margins[1] / 2., 550, -50, 500)
        self.tree.Draw('{sig}:IntegralPeaks[{num}]/2>>h'.format(sig=self.SignalName, num=three_bucket_num), self.Cut.CutStrings['tracks'] + self.Cut.CutStrings['pulser'], 'goff')
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

    def draw_single_wf(self, event=None):
        cut = '!({0})&&!pulser'.format(self.Cut.CutStrings['old_bucket'])
        return self.draw_waveforms(n=1, cut_string=cut, add_buckets=True, ret_event=True, start_event=event)

    def draw_waveforms(self, n=1000, start_event=None, cut_string=None, show=True, ret_event=False, add_buckets=False, fixed_range=None):
        """
        Draws stacked waveforms.
        :param n: number of waveforms
        :param cut_string:
        :param start_event: event to start
        :param show:
        :param ret_event: return number of valid events if True
        :param add_buckets: draw buckets and most probable peak values if True
        :param fixed_range: fixes x-range to given value if set
        :return: histo with waveform
        """
        gROOT.SetBatch(1)
        t = time()
        start = self.StartEvent if start_event is None else start_event
        assert self.run.n_entries >= start >= 0, 'The start event is not within the range of tree events!'
        if not self.run.wf_exists(self.channel):
            return
        cut = self.Cut.all_cut if cut_string is None else cut_string
        n_events = self.find_n_events(n, cut, start)
        h = TH2F('wf', 'Waveform', 1024, 0, 511, 1000, -500, 500)
        h.SetStats(0)
        gStyle.SetPalette(55)
        self.tree.Draw('wf0:Iteration$/2>>wf', cut, 'goff', n_events, start)
        if fixed_range is None:
            h.GetYaxis().SetRangeUser(-500 + h.FindFirstBinAbove(0, 2) / 50 * 50, -450 + h.FindLastBinAbove(0, 2) / 50 * 50)
        else:
            assert type(fixed_range) is list, 'Range has to be a list!'
            h.GetYaxis().SetRangeUser(fixed_range[0], fixed_range[1])
        self.print_elapsed_time(t)
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'WaveForm', 1000, 500)
        c.SetRightMargin(.045)
        self.format_histo(h, x_tit='Time [ns]', y_tit='Signal [au]', markersize=.4)
        draw_option = 'scat' if n == 1 else 'col'
        h.Draw(draw_option)
        if add_buckets:
            sleep(.2)
            h.GetXaxis().SetNdivisions(26)
            c.SetGrid()
            c.SetBottomMargin(.186)
            self._add_buckets(c)
        gROOT.SetBatch(0)
        if n > 1:
            self.save_plots('WaveForms{n}'.format(n=n), sub_dir=self.save_dir)
        else:
            self.save_plots('SingleWaveForm', sub_dir=self.save_dir)
        self.histos[0] = [c, h]
        return h if not ret_event else n_events

    def show_single_waveforms(self, n=1, cut=None):
        ev = self.StartEvent
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
        n = mean([self.tree.Draw('1', cut, 'goff', n_events, start + i * n_events) for i in xrange(4)])
        new_events = n_events
        ratio = n_events / n if n else 5
        i = 0
        while n != n_events:
            diff = n_events - n
            # print n, diff, new_events
            if i < 3:
                new_events += int(diff * ratio)
            else:
                new_events += int(diff * ratio / i)
            print '\b.',
            stdout.flush()
            n = self.tree.Draw('1', cut, 'goff', new_events, start)
            i += 1
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

    def analyse_signal_histograms(self):
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
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                print 'process cut ' + key
                # h = TH1F('h', '', 600, -100, 500)
                # self.tree.Draw("{name}>>h".format(name=self.signal_name), value)
                h = self.show_signal_histo(evnt_corr=True, cut=value, show=False)
                i_mean = self.__get_mean(h)
                median = self.__get_median(h)
                mpv = self.__get_mpv(h)
                # print mean, median, mpv
                gr1.SetPoint(i, i, i_mean[0])
                gr1.SetPointError(i, 0, i_mean[1])
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
        for key, value in self.Cut.CutStrings.iteritems():
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
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
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
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('SNR', sub_dir=self.save_dir)
        self.canvases[0] = c
        self.histos[0] = gr
        self.histos['legend'] = [l1, l2]

    def calc_snr(self, sig=None):
        signal = self.SignalName if sig is None else sig
        peak_int = self.get_all_signal_names()[signal][-2:] if self.get_all_signal_names()[signal][-2].isdigit() else self.get_all_signal_names()[signal][-1]
        ped_fit = self.show_pedestal_histo(draw=False, peak_int=peak_int)
        sig_fit = self.draw_pulse_height(evnt_corr=True, show=False, sig=signal)
        sig_mean = sig_fit.Parameter(0)
        ped_sigma = ped_fit.Parameter(2)

        snr = sig_mean / ped_sigma
        snr_err = snr * (sig_fit.ParError(0) / sig_mean + ped_fit.ParError(2) / ped_sigma)
        print 'SNR is: {snr} +- {err}'.format(snr=snr, err=snr_err)
        return [snr, snr_err]

    # ============================================
    # region PEAK INTEGRAL

    def find_best_snr(self, show=True, same_width=False):
        gROOT.SetBatch(1)
        gr = self.make_tgrapherrors('gr', 'Signal to Noise Ratios')
        peak_integrals = OrderedDict(sorted({key: value for key, value in self.run.peak_integrals.iteritems() if len(key) < 3}.items()))
        i = 0
        for name, value in peak_integrals.iteritems():
            signal = self.get_signal_name('b', name)
            snr = self.calc_snr(signal)
            print value
            x = (value[1] + value[0]) / 2. if not same_width else value[0] / 2.
            gr.SetPoint(i, x, snr[0])
            gr.SetPointError(i, 0, snr[1])
            i += 1
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        self.format_histo(gr, x_tit='Integralwidth [ns]', y_tit='SNR')
        gr.Draw('ap')
        gROOT.SetBatch(0)
        self.save_plots('BestSNR', sub_dir=self.save_dir)
        self.save_plots('BestSNR', file_type='root', sub_dir=self.save_dir)
        self.histos[0] = [gr, c]

    def signal_vs_peakintegral(self, show=True, ped=False):
        gROOT.SetBatch(1)
        gr = self.make_tgrapherrors('gr', '{sig} vs Peak Integral'.format(sig='Signal' if not ped else 'Pedestal'))
        peak_integrals = OrderedDict(sorted({key: value for key, value in self.run.peak_integrals.iteritems() if len(key) < 3}.items()))
        i = 0
        ratio = '{0}{1}'.format(self.run.peak_integrals.values()[0][0], self.run.peak_integrals.values()[0][1])
        for name, value in peak_integrals.iteritems():
            sig_name = self.get_signal_name(region='b', peak_integral=name)
            signal = self.draw_pulse_height(evnt_corr=True, show=False, sig=sig_name) if not ped else self.show_pedestal_histo(draw=False, peak_int=name)
            par = 2 if ped else 0
            gr.SetPoint(i, (value[1] + value[0]) / 2., signal.Parameter(par))
            gr.SetPointError(i, 0, signal.ParError(par))
            i += 1
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Signal vs Peak Integral', 1000, 1000)
        self.format_histo(gr, x_tit='Integralwidth [ns]', y_tit='Signal [au]', y_off=1.3)
        gr.Draw('ap')
        gROOT.SetBatch(0)
        self.save_plots('{sig}PeakInt_{rat}'.format(rat=ratio, sig='Ped' if ped else 'Sig'), sub_dir=self.save_dir)
        self.save_plots('{sig}PeakInt_{rat}'.format(rat=ratio, sig='Ped' if ped else 'Sig'), sub_dir=self.save_dir, file_type='root')
        self.histos[0] = [gr, c]

    # endregion

    # ============================================
    # region MISCELLANEOUS

    def get_cut(self):
        """ :return: full cut_string """
        return self.Cut.all_cut

    def get_peak_position(self, event=None, region='b', peak_int='2'):
        num = self.get_signal_number(region, peak_int)
        ev = self.StartEvent if event is None else event
        self.tree.GetEntry(ev)
        return self.tree.IntegralPeaks[num]

    def get_all_signal_names(self):
        names = OrderedDict()
        regions = [reg for reg in self.run.signal_regions if len(reg) < 3]
        integrals = [integral for integral in self.run.peak_integrals if len(integral) < 3]
        for region in regions:
            for integral in integrals:
                if len(integral) > 2:
                    integral = '_' + integral
                name = 'ch{ch}_signal_{reg}_PeakIntegral{int}'.format(ch=self.channel, reg=region, int=integral)
                num = self.IntegralNames[name]
                reg = region + integral
                names['{pol}*IntegralValues[{num}]'.format(pol=self.Polarity, num=num)] = reg
        return names

    def __get_binning(self):
        jumps = self.Cut.jump_ranges
        n_jumps = len(jumps['start'])
        bins = [self.Cut.get_min_event()]
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
            while bins[-1] + self.BinSize < start:
                bins.append(bins[-1] + self.BinSize)
            # two jumps shortly after one another
            if ind < n_jumps - 2:
                next_start = jumps['start'][ind + 1]
                next_stop = jumps['stop'][ind + 1]
                if bins[-1] + self.BinSize + gap > next_start:
                    gap2 = next_stop - next_start
                    bins.append(bins[-1] + self.BinSize + gap + gap2)
                else:
                    bins.append(bins[-1] + self.BinSize + gap)
            else:
                bins.append(bins[-1] + self.BinSize + gap)
            # fill up the end
            if ind == n_jumps - 1 and bins[-1] >= stop:
                while bins[-1] + self.BinSize < self.run.n_entries:
                    bins.append(bins[-1] + self.BinSize)
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
        infos = [self.run_number, self.run.RunInfo['type'], self.diamond_name.ljust(4), self.bias, self.SignalRegion + self.PeakIntegral + '   ']
        for info in infos:
            print self.adj_length(info),
        print

    def print_integral_names(self):
        for key, value in self.IntegralNames.iteritems():
            if key.startswith('ch{ch}'.format(ch=self.channel)):
                print str(value).zfill(3), key
        return

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
    message = 'STARTING PAD-ANALYSIS OF RUN {0}'.format(test_run)
    print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(message)) * '=', msg=message)
    a = Elementary()
    a.print_testcampaign()
    z = SignalAnalysis(test_run, args.ch)
    z.print_elapsed_time(st, 'Instantiation')
