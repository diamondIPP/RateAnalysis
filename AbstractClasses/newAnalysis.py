import copy
import json
import numpy as np
from numpy import array, zeros
import types as t

import ROOT
from ROOT import TCanvas, TH2F, gROOT, TProfile, TH1F, TLegend, gStyle, kGreen, TArrow, kOrange, kViolet, kCyan
from AbstractClasses.ConfigClass import *
from AbstractClasses.RunClass import Run
from AbstractClasses.Cut import Cut
from BinCollection import BinCollection
from collections import OrderedDict
from argparse import ArgumentParser
from ConfigParser import ConfigParser
from copy import deepcopy
from time import sleep


class Analysis(Elementary):
    """ Class for the analysis of the non-channel specific stuff of a single run. """

    def __init__(self, run, diamonds=3, verbose=False, high_low_rate=None):
        """
        An Analysis Object collects all the information and data for the analysis of one single run.
        The most important fields are:
            .run        - instance of Run class, containing run info and data
            .cut        - dict containing two instances of Cut class
            .RunInfo    - dict containing run infos
            .Pads       - dict containing two intances of
        The most important methods are:
            .CreateMeanSignalHistogram(channel)
            .Draw(varexp)
            .DrawRunInfo()
            .GetCut(channel)
            .GetEventAtTime(dt)
            .MakePreAnalysis()
            .ShowPulserRate()
        :param run:         run object of type "Run" or integer run number
        :param diamonds:    An integer number defining the diamonds activated for analysis: 0x1=ch0 (diamond 1) 0x2=ch3 (diamond 2)
        :param verbose:     if True, verbose printing is activated
        """
        Elementary.__init__(self, verbose=verbose)

        # basics
        self.diamonds = diamonds
        self.run = self.init_run(run)
        self.run.analysis = self
        self.run_number = self.run.run_number
        self.RunInfo = deepcopy(self.run.RunInfo)
        self.lowest_rate_run = high_low_rate['min'] if high_low_rate is not None else self.run.run_number
        self.highest_rate_run = high_low_rate['max'] if high_low_rate is not None else self.run.run_number
        self.parser = self.load_parser()
        self.pickle_dir = self.get_program_dir() + self.parser.get('SAVE', 'pickle_dir')
        # self.saveMCData = parser.getboolean("SAVE", "SaveMCData")
        self.ana_save_dir = '{tc}_{run}'.format(tc=self.TESTCAMPAIGN[2:], run=self.run.run_number)

        # tree
        self.tree = self.run.tree

        # miscellaneous
        self.Polarity = self.get_polarity() if hasattr(self, 'get_polarity') else 1
        self.channel = self.channel if hasattr(self, 'channel') else None

        # names
        self.IntegralNames = self.get_integral_names()
        self.SignalRegion = self.parser.get('BASIC', 'signal_region')
        self.PedestalRegion = self.parser.get('BASIC', 'pedestal_region')
        self.PeakIntegral = self.parser.get('BASIC', 'peak_integral')
        # todo FIX THAT
        self.SignalName = self.get_signal_name(region=self.SignalRegion, peak_integral=self.PeakIntegral) if hasattr(self, 'get_signal_name') else None
        self.PedestalName = self.get_pedestal_name(region=self.PedestalRegion, peak_int=self.PeakIntegral) if hasattr(self, 'get_pedestal_name') else None
        self.PulserName = self.get_pulser_name() if hasattr(self, 'get_pulser_name') else None

        # self.signal_num = self.get_signal_numbers(self.signal_region, self.peak_integral)
        # self.signal_names = self.get_signal_names()
        # self.pedestal_num = self.get_pedestal_numbers(self.pedestal_region, self.peak_integral)
        # self.pedestal_names = self.get_pedestal_names()

        self.Checklist = self.init_checklist()
        self.ExtremaResults = self.init_extrema_results()
        self.SignalHistoFitResults = self.init_signal_hist_fit_results()

        self.Cut = Cut(self)
        self.start_event = self.Cut.CutConfig['EventRange'][0]
        self.pedestalFitMean = {}

        # save histograms // canvases
        self.signal_canvas = None
        self.histos = {}
        self.canvases = {}
        self.lines = {}
        self.pedestalSigma = {}

        # alignment
        self.IsAligned = self.check_alignment(draw=False)

    # ============================================================================================
    # region INIT

    def load_parser(self):
        parser = ConfigParser()
        parser.read("Configuration/AnalysisConfig_" + self.TESTCAMPAIGN + ".cfg")
        return parser

    def get_integral_names(self):
        names = OrderedDict()
        self.tree.GetEntry(0)
        for i, name in enumerate(self.tree.IntegralNames):
            names[name] = i
        return names

    def init_extrema_results(self):
        names = ['TrueNPeaks', 'FoundNMaxima', 'FoundMaxima', 'FoundNMinima', 'FoundMinima', 'Ninjas', 'Ghosts']
        sub_dic = {name: False for name in names}
        sub_dic['SignalHeight'] = 0.
        return {ch: sub_dic for ch in self.run.channels}

    def init_signal_hist_fit_results(self):
        names = ['FitFunction', 'Peak', 'FWHM', 'Chi2', 'NDF']
        sub_dic = {name: False for name in names}
        return {ch: sub_dic for ch in self.run.channels}

    def init_checklist(self):
        dic = {}
        ch_dic = {ch: False for ch in self.run.channels}
        dic['LoadTrackData'] = False
        dic['HitsDistribution'] = False
        dic['GlobalPedestalCorrection'] = ch_dic
        dic['MeanSignalHisto'] = ch_dic
        dic['FindExtrema'] = {'Minima': ch_dic, 'Maxima': ch_dic}
        return dic

    def init_run(self, run):
        if not isinstance(run, Run):
            assert type(run) is int, 'run has to be either a Run instance or an integer run number'
            return Run(run, self.diamonds)
        else:
            assert run.run_number is not None, 'No run selected, choose run.SetRun(run_nr) before you pass the run object'
            return run

    def load_bincollection(self):
        dic = {}
        for ch in self.run.channels:
            dic[ch] = BinCollectionConfig(run=self.run, channel=ch)
        return dic
    # endregion

    # ============================================================================================
    # region REGIONS AND PEAK INTEGRAL

    def __draw_single_wf(self, event=None, show=True):
        start = self.start_event if event is None else event
        if hasattr(self, 'draw_waveforms') and self.run.wf_exists(self.channel):
            h = self.draw_waveforms(n=1, show=show, start_event=start)
        else:
            h = TH2F('regions', '', 1024, 0, 511, 1000, -200, 50)
            if self.run.wf_exists(0):
                self.tree.Draw('wf0:Iteration$/2>>regions', self.Cut.all_cut, 'goff', 1, start)
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c2', 'Regions', 1000, 500)
        c.SetMargin(.075, .045, .1, .1)
        c.SetGrid()
        h.SetStats(0)
        h.GetXaxis().SetNdivisions(26)
        self.format_histo(h, markersize=0.3, x_tit='Time [ns]', y_tit='Signal [au]')
        h.Draw()
        gROOT.SetBatch(0)
        self.histos[0] = [c, h]
        return h

    def draw_regions(self, ped=True, event=None):
        h = self.__draw_single_wf(event=event, show=False)
        c = TCanvas('c1', 'Regions', 1000, 500)
        c.SetMargin(.075, .045, .2, .1)
        c.SetGrid()
        h.Draw()
        tit = 'Pedestal Regions' if ped else 'Signal Regions'
        h.SetTitle(tit)
        lines = []
        starts = []
        titles = []
        regions = self.run.pedestal_regions if ped else self.run.signal_regions
        gr = self.make_tgrapherrors('gr', '', color=2, marker_size=0, width=3)
        i = 0
        gStyle.SetEndErrorSize(4)
        sleep(.5)
        for reg, lst in regions.iteritems():
            if len(reg) < 3:
                offset = 40 if not lst[0] in starts else 20
                if lst[1] - lst[0] > 1:
                    gr.SetPoint(i, (lst[1] + lst[0]) / 4., c.GetUymax() - offset)
                    gr.SetPointError(i, (lst[1] - lst[0]) / 4., 0)
                    l = self.make_tlatex(gr.GetX()[i], gr.GetY()[i] + 3, '{sig}{reg}'.format(reg=reg, sig='p' if ped else 's'), color=2, size=.04)
                    gr.GetListOfFunctions().Add(l)
                    i += 1
                l1 = self.make_tgaxis(lst[0] / 2, c.GetUymin(), c.GetUymax() - offset, '', 2)
                l2 = self.make_tgaxis(lst[1] / 2, c.GetUymin(), c.GetUymax() - offset, '', 2) if lst[1] - lst[0] > 1 else 0
                if not lst[1] - lst[0] > 1:
                    l1.SetLineColor(4)
                    l1.SetLineWidth(2)
                    l1.SetTitleColor(4)
                    l1.SetY2(c.GetUymax() - 100)
                    tit = self.make_tlatex(lst[0] / 2, c.GetUymax() - 97, '{sig}{reg}'.format(reg=reg, sig='p' if ped else 's'), size=.04, color=4)
                    tit.Draw()
                    titles.append(tit)
                l1.Draw()
                l2.Draw() if l2 else self.do_nothing()
                lines.append([l1, l2])
                starts.append(lst[0])
        gr.Draw('[]')
        gr.Draw('p')
        self._add_buckets()
        save_name = 'PedestalRegions' if ped else 'SignalRegions'
        self.save_plots(save_name, sub_dir=self.ana_save_dir, ch=None)
        self.histos[0] = [h, c, gr, lines, titles]

    def _add_buckets(self, canvas=None):
        c = gROOT.GetSelectedPad() if canvas is None else canvas
        axis = []
        labels = []
        arrows = []
        start = self.run.signal_regions['b'][0] % 40
        stop = int(.8 * c.GetUxmax()) if c.GetUxmax() > 500 else int(c.GetUxmax())
        print 'sta-sto:', start, stop
        bucket0 = self.run.signal_regions['b'][0] / 40
        x_range = c.GetUxmax() - c.GetUxmin()
        print 'xrange', x_range
        l = self.make_tlatex(c.GetUxmin() - .015 * x_range, c.GetUymin() - 30, 'Bucket:', align=30, color=kGreen + 2, size=0.03)
        l.Draw()
        labels.append(l)
        l1 = self.make_tlatex(c.GetUxmin() - .015 * x_range, c.GetUymin() + 8, 'Peak:', align=30, color=kOrange + 7, size=0.03)
        l1.Draw()
        labels.append(l1)
        peak_fit = self.fit_peak_values(draw=False) if hasattr(self, 'fit_peak_values') else 0
        for i, x in enumerate(xrange(start, stop, 20), -bucket0):
            a = self.make_tgaxis(x, c.GetUymin() - 30, c.GetUymin() - 12, '', kGreen + 2)
            if x <= stop - 20:
                l = self.make_tlatex(x + 10, c.GetUymin() - 30, str(i), align=20, color=kGreen + 2, size=0.03)
                labels.append(l)
                l.Draw()
                if peak_fit:
                    pos = peak_fit.Parameter(1) % 20
                    ar = TArrow(x + pos, c.GetUymin() + 1, x + pos, c.GetUymin() + 10, .005, '<|')
                    ar.SetLineWidth(2)
                    ar.SetFillColor(kOrange + 7)
                    ar.SetLineColor(kOrange + 7)
                    ar.Draw()
                    arrows.append(ar)
            a.Draw()
            axis.append(a)
        self.histos[1] = [axis, labels, arrows]

    def draw_peak_integrals(self, event=None):
        h = self.__draw_single_wf(event=event, show=False)
        c = TCanvas('c1', 'Regions', 1000, 500)
        c.SetMargin(.075, .045, .2, .1)
        c.SetGrid()
        self.format_histo(h, title='Peak Integrals', markersize=.5)
        h.GetXaxis().SetRangeUser(self.run.signal_regions['e'][0] / 2 - 20, self.run.signal_regions['e'][1] / 2)
        h.Draw()
        sleep(.5)
        # draw line at found peak and pedestal 'ab'
        peak_pos = self.get_peak_position(event) / 2. if hasattr(self, 'get_peak_position') else self.run.signal_regions['a'][0] / 2.
        ped_pos = self.run.pedestal_regions['ab'][1] / 2.
        l = self.make_tgaxis(peak_pos, c.GetUymin(), c.GetUymax() - 100, '', 4, 2)
        l2 = self.make_tgaxis(ped_pos, c.GetUymin(), c.GetUymax() - 100, '', kViolet + 3, 2)
        l2.Draw()
        l.Draw()
        t1 = self.make_tlatex(peak_pos, c.GetUymax() - 97, 'found peak', color=4)
        t2 = self.make_tlatex(ped_pos, c.GetUymax() - 97, 'ab', color=kViolet + 3)
        t1.Draw()
        t2.Draw()
        # draw error bars
        gr1 = self.make_tgrapherrors('gr1', '', color=kGreen + 2, marker_size=0, asym_err=True, width=3)
        gr2 = self.make_tgrapherrors('gr2', '', color=kCyan - 3, marker_size=0, asym_err=True, width=3)
        gStyle.SetEndErrorSize(4)
        i = 0
        for int_, lst in self.run.peak_integrals.iteritems():
            if len(int_) < 3:
                gr1.SetPoint(i, peak_pos, c.GetUymax() - 30 * (i + 1) - 100)
                gr2.SetPoint(i, ped_pos, c.GetUymax() - 33 * (i + 1) - 100)
                gr1.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr1.SetPointError(i, .5, .5, 0, 0)
                gr2.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr2.SetPointError(i, .5, .5, 0, 0)
                l1 = self.make_tlatex(gr1.GetX()[i], gr1.GetY()[i] + 5, ' ' + int_, color=kGreen + 2, align=10)
                gr1.GetListOfFunctions().Add(l1)
                i += 1
        for gr in [gr1, gr2]:
            gr.Draw('[]')
            gr.Draw('p')
        self._add_buckets()
        self.save_plots('IntegralPeaks', sub_dir=self.ana_save_dir, ch=None)
        self.histos[0] = [gr1, gr2, c, l, t1, h, l2, t2]
    # endregion

    # ============================================================================================
    # region TRACKS
    def show_chi2(self, mode=None, show=True):
        gROOT.SetBatch(1)
        assert mode in ['x', 'y', None], 'mode has to be in {lst}!'.format(lst=['x', 'y', None])
        n_bins = 500 if mode is None else 1000
        mode = 'tracks' if mode is None else mode
        h = TH1F('h', '#chi^{2} in ' + mode, n_bins, 0, 100)
        self.tree.Draw('chi2_{mod}>>h'.format(mod=mode), '', 'goff')
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Chi2 in ' + mode, 1000, 1000)
        c.SetLeftMargin(.13)
        if show or mode == 'tracks':
            yq = zeros(1)
            h.GetQuantiles(1, yq, array([.9]))
            h.GetXaxis().SetRangeUser(0, yq[0])
        self.format_histo(h, x_tit='#chi^{2}', y_tit='Entries', y_off=1.8)
        h.Draw()
        self.histos[0] = h
        self.canvases[0] = c
        gROOT.SetBatch(0)
        return h

    def show_all_chi2(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [self.show_chi2(mode, show=False) for mode in [None, 'x', 'y']]
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        max_chi2 = int(max([h.GetMaximum() for h in histos])) / 1000 * 1000 + 1000
        histos[0].GetYaxis().SetRangeUser(0, max_chi2)
        histos[0].SetTitle('All #chi^{2}')
        legend = TLegend(.7, .7, .9, .9)
        leg_names = ['#chi^{2} ' + mode for mode in ['of both', 'in x', 'in y']]
        for i, h in enumerate(histos):
            h.SetStats(0)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, leg_names[i], 'l')
            self.histos[i] = h
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.histos['legend'] = legend
        self.save_plots('Chi2', canvas=c, sub_dir=self.ana_save_dir, ch=None)
        self.canvases[0] = c

    def show_angle(self, mode='x', show=True):
        """
        Displays the angle distribution of the tracks.
        :param mode: has to be eiher 'x' or 'y'
        :param show:
        :return: histogram
        """
        assert mode in ['x', 'y']
        gROOT.SetBatch(1)
        h = TH1F('h', 'Track Angle Distribution in ' + mode, 320, -4, 4)
        self.tree.Draw('slope_{mod}>>h'.format(mod=mode), '', 'goff')
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Angle in ' + mode, 1000, 1000)
        c.SetLeftMargin(.13)
        self.format_histo(h, x_tit='Track Angle [deg]', y_tit='Entries', y_off=1.8, lw=2)
        h.Draw()
        self.histos[0] = h
        self.canvases[0] = c
        gROOT.SetBatch(0)
        # a = gROOT.GetListOfCanvases()
        # print a[0]
        self.save_plots('TrackAngle{mod}'.format(mod=mode.upper()), sub_dir=self.ana_save_dir, ch=None)
        return h

    def calc_angle_fit(self, mode='x', show=True):
        pickle_path = self.pickle_dir + 'Tracks/AngleFit_{tc}_{run}_{mod}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, mod=mode)

        def func():
            h = self.show_angle(mode, show=show)
            return self.fit_fwhm(h, draw=show)

        fit = self.do_pickle(pickle_path, func)
        return fit

    def show_both_angles(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        self.get_color()
        histos = [self.show_angle(mode, show=False) for mode in ['x', 'y']]
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        max_angle = int(max([h.GetMaximum() for h in histos])) / 1000 * 1000 + 1000
        histos[0].GetYaxis().SetRangeUser(0, max_angle)
        legend = TLegend(.7, .7, .9, .9)
        leg_names = ['Angle in ' + mode for mode in ['x', 'y']]
        for i, h in enumerate(histos):
            h.SetStats(0)
            h.SetTitle('Track Angle Distributions')
            h.SetLineColor(self.get_color())
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, leg_names[i], 'l')
            self.histos[i] = h
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.canvases[0] = c
        self.histos['legend'] = legend
        self.save_plots('TrackAngles', sub_dir=self.ana_save_dir, ch=None)
    # endregion

    # ==============================================
    # region ALIGNMENT
    def check_alignment(self, binning=5000, draw=True):
        pickle_path = 'Configuration/Individual_Configs/Alignment/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run.run_number)

        def func():
            gROOT.SetBatch(1)
            nbins = self.run.n_entries / binning
            h = TProfile('h', 'Pulser Rate', nbins, 0, self.run.n_entries)
            self.tree.Draw('(@col.size()>1)*100:Entry$>>h', 'pulser', 'goff')
            self.format_histo(h, name='align', title='Event Alignment', x_tit='Event Number', y_tit='Hits per Event @ Pulser Events [%]', y_off=1.3)
            h.GetYaxis().SetRangeUser(0, 100)
            if draw:
                gROOT.SetBatch(0)
            c = TCanvas('c', 'Pulser Rate Canvas', 1000, 1000)
            h.Draw('hist')
            self.save_plots('EventAlignment', sub_dir=self.ana_save_dir, ch=None)
            gROOT.SetBatch(0)
            self.histos[0] = [h, c]
            align = self.__check_alignment_histo(h)
            return align

        aligned = func() if draw else self.do_pickle(pickle_path, func)
        if not aligned:
            msg = 'The events are not aligned!'
            print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(msg)) * '!', msg=msg)
        return aligned

    def find_alignment_offset(self):
        offsets = [i for i in xrange(-3, 4) if i]
        h = TH1F('h', 'Pixel Hits @ Pulser Events', 20, 0, 20)
        right_offset = None
        for offset in offsets:
            pulser_events = 0
            for event in xrange(self.start_event, self.run.n_entries):
                print '\rpulser events: {0:04d}'.format(pulser_events),
                if pulser_events >= 1000:
                    break
                self.tree.GetEntry(event)
                if self.tree.pulser:
                    pulser_events += 1
                    self.tree.GetEntry(event + offset)
                    hits = len(self.tree.col)
                    h.Fill(hits)
            if self.__check_alignment_histo(h):
                right_offset = offset
                break
            h.Reset()
        h.Draw()
        self.histos[0] = h
        print '\nThe event offset is {off}'.format(off=right_offset)

    def __check_alignment_histo(self, histo):
        h = histo
        for bin_ in xrange(h.FindBin(self.start_event), h.GetNbinsX()):
            if h.GetBinContent(bin_) > 40:
                return False
        return True
    # endregion

    # ==============================================
    # region SHOW & PRINT

    def draw_preliminary(self):
        c = gROOT.GetListOfCanvases()[-1]
        text = ROOT.TText((c.GetUxmax() - c.GetUxmin()) / 2., (c.GetUymax() - c.GetUymin()) / 2., "Preliminary")
        text.SetTextColor(19)
        text.SetTextSize(.2)
        text.SetTextAngle(30)
        text.SetTextAlign(20)
        h = None
        for obj in c.GetListOfPrimitives():
            print obj.IsA().GetName()
            if obj.IsA().GetName() in ['TH1F', 'TH2F', 'TGraph', 'TGraphErrors']:
                h = obj
        text.Draw()
        h.Draw('same')
        c.RedrawAxis()
        self.histos[1] = text

    def ShowFFT(self, drawoption="", cut=None, channel=None, startevent=0, endevent=10000000, savePlots=True, canvas=None):
        '''
        Draws the fft_mean VS 1/fft_max scatter plot.
        For the drawoption ROOT TH2D Drawoption can be used, as well as
        "mc" for a special multicolor plot, which differentiates between
        different waveform types.
        The cut can be None, for using a default cut, True for using the
        cut from the Analysis instance or can be set customly.
        :param drawoption:
            "" draws all events with the given cut
            "mc" for multicolor: draw different colors for different
                waveform types
            "col" gives colored TH2D plot (also possible: "colz")
        :param cut:
            if True: apply cut from Analysis instance
            if None: "!pulser" set as default (except for "mc"
                drawoption)
            custom: cut="!pulser&&!is_saturated[{channel}]"
                (--> syntax "[{channel}]" important for vectors!)
        :param channel:
            if None: takes the activated channels (diamonds) from
                Analysis instance
            else: 0 or 3 to select one particular channel
        :param startevent:
            first event number to draw
        :param endevent:
            last event number to draw
        :param savePlots:
            True to save the generated scatter plots
        :return:
        '''
        channels = self.GetChannels(channel=channel)

        if drawoption in ["MC", "mc", "Mc", "mC"]:
            drawoption = ""
            multicolor = True
            if cut == None:
                cut = ""
        else:
            multicolor = False

        if cut == None:
            cut = "!pulser"

        events = endevent - startevent
        if events < 10000000:
            if endevent >= 10000000: endevent = self.run.tree.GetEntries()
            comment = "Events: {start}-{end}".format(start=startevent, end=endevent)

        else:
            comment = None
        if canvas == None:
            self.fftcanvas = ROOT.TCanvas("fftcanvas", "fftcanvas", len(channels) * 500, 500)
        else:
            # assert(isinstance(canvas, ROOT.TCanvas))
            self.fftcanvas = canvas

        self.fftcanvas.Divide(len(channels), 1)

        if len(channels) == 2:
            namesuffix = ""
        else:
            namesuffix = "_Ch{ch}".format(ch=channels[0])

        self.fftcanvas.cd()
        i = 1
        self.fftHistos = {}
        self.fftlegend = {}
        ROOT.gStyle.SetOptStat(0)
        for ch in channels:
            if cut == True:
                thiscut = self.GetCut(channel=ch)
                thisusercut = self.GetUserCutString(channel=ch)
                events = self.GetNEventsCut(channel=ch)
                startevent = self.GetMinEventCut(channel=ch)
            else:
                thiscut = cut.format(channel=ch)
                thisusercut = thiscut
            self.fftHistos[ch] = ROOT.TH2D("fft_ch{channel}".format(channel=ch), "FFT {" + thisusercut + "}", 5000, 2e-6, 0.0025, 5000, 1e1, 1e4)
            self.fftcanvas.cd(i)
            ROOT.gPad.SetLogy()
            ROOT.gPad.SetLogx()
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}".format(channel=ch), thiscut, "", events, startevent)
            self.fftHistos[ch].SetTitle("{diamond} ".format(diamond=self.run.diamond_names[ch]) + self.fftHistos[ch].GetTitle())
            self.fftHistos[ch].GetXaxis().SetTitle("1/fft_max")
            self.fftHistos[ch].GetYaxis().SetTitle("fft_mean")
            # self.fftHistos[ch].Draw(drawoption)
            if multicolor:
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_sat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch),
                                       "is_saturated[{channel}]".format(channel=ch), "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_sat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch),
                                       (thiscut + "&&is_saturated[{channel}]").format(channel=ch), "", events, startevent)
                saturated_histo = ROOT.gROOT.FindObject("fft_ch{channel}_sat".format(channel=ch))
                saturated_histo.SetMarkerStyle(1)
                saturated_histo.SetMarkerColor(6)
                saturated_histo.SetFillColor(6)
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_med(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch),
                                       "abs(median[{channel}])>8".format(channel=ch), "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_med(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch),
                                       (thiscut + "&&abs(median[{channel}])>8").format(channel=ch), "", events, startevent)
                median_histo = ROOT.gROOT.FindObject("fft_ch{channel}_med".format(channel=ch))
                median_histo.SetMarkerStyle(1)
                median_histo.SetMarkerColor(8)
                median_histo.SetFillColor(8)
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_flat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch),
                                       "sig_spread[{channel}]<10".format(channel=ch), "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_flat(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch),
                                       (thiscut + "&&sig_spread[{channel}]<10").format(channel=ch), "", events, startevent)
                flat_histo = ROOT.gROOT.FindObject("fft_ch{channel}_flat".format(channel=ch))
                flat_histo.SetMarkerStyle(1)
                flat_histo.SetMarkerColor(4)
                flat_histo.SetFillColor(4)
                if cut == "":
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_pulser(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch), "pulser", "", events, startevent)
                else:
                    self.run.tree.Draw("fft_mean[{channel}]:1./fft_max[{channel}]>>fft_ch{channel}_pulser(5000, 2e-6, 0.0025, 5000, 1e1, 1e4)".format(channel=ch),
                                       (thiscut + "&&pulser").format(channel=ch), "", events, startevent)
                pulser_histo = ROOT.gROOT.FindObject("fft_ch{channel}_pulser".format(channel=ch))
                pulser_histo.SetMarkerStyle(1)
                pulser_histo.SetMarkerColor(2)
                pulser_histo.SetFillColor(2)
                self.fftHistos[ch].Draw("")
                saturated_histo.Draw("same")
                median_histo.Draw("same")
                pulser_histo.Draw("same")
                flat_histo.Draw("same")
                self.fftlegend[ch] = ROOT.TLegend(0.1, 0.1, 0.3, 0.3)
                self.fftlegend[ch].AddEntry(pulser_histo, "pulser", "f")
                self.fftlegend[ch].AddEntry(saturated_histo, "saturated", "f")
                self.fftlegend[ch].AddEntry(median_histo, "|median|>8", "f")
                self.fftlegend[ch].AddEntry(flat_histo, "sig_spread<10", "f")
                self.fftlegend[ch].Draw("same")
                self.fftcanvas.Update()
                self.DrawRunInfo(channel=ch, comment=comment, canvas=self.fftcanvas)
            else:
                self.fftHistos[ch].Draw(drawoption)
                self.fftcanvas.Update()
                self.DrawRunInfo(channel=ch, comment=comment, canvas=self.fftcanvas)
            i += 1
        self.fftcanvas.Update()
        if savePlots:
            savename = "Run{run}_FFT" + namesuffix
            self.save_plots(savename.format(run=self.run.run_number), file_type="png", canvas=self.fftcanvas, sub_dir="png/")
            # self.SavePlots(savename.format(run=self.run.run_number), ending="root", canvas=self.fftcanvas, subDir="root/")

        self.if_wait("FFT shown..")

    def GetIncludedEvents(self, maxevent, channel=None):
        '''
        Returns list of all event numbers, which are neither excluded by
        beaminerruptions, by the excludeFirst cut, or are outside of a
        event range window.
        :return: list of included event numbers
        '''
        if channel == None:
            minevent0 = self.Cut[0].get_min_event()
            minevent3 = self.Cut[3].get_min_event()
            maxevent0 = self.Cut[0].get_max_event()
            maxevent3 = self.Cut[3].get_max_event()
            minevent = min(minevent0, minevent3)
            if maxevent == None:
                maxevent = max(maxevent0, maxevent3)
            else:
                maxevent = min(max(maxevent0, maxevent3), maxevent)

            excluded = [i for i in np.arange(0, minevent)]  # first events
            # if self.cuts[0].cut_types["noBeamInter"] and self.cuts[3].cut_types["noBeamInter"]:
            #     self.cuts[0].get_beam_interruptions()
            for i in xrange(len(self.Cut[0].jump_ranges["start"])):
                excluded += [i for i in np.arange(self.Cut[0].jump_ranges["start"][i], self.Cut[0].jump_ranges["stop"][i] + 1)]  # events around jumps
            excluded.sort()
            all_events = np.arange(0, maxevent)
            included = np.delete(all_events, excluded)
            return included
            # cut0 = self.cut[0].GetIncludedEvents(maxevent=maxevent)
            # cut3 = self.cut[3].GetIncludedEvents(maxevent=maxevent)
            # assert(cut0 == cut3), "Included Events not The same for both channels, select a particular channel"
            # return cut0
        else:
            assert (channel in [0, 3])
            return self.Cut[channel].get_included_events(maxevent=maxevent)

    # def GetMinEventCut(self, channel=None):
    #     if channel == None:
    #         opt0 = self.Cut[0].get_min_event()
    #         opt3 = self.Cut[3].get_min_event()
    #         assert (opt0 == opt3), "GetMinEvent Not the same for both cut channels. Choose a specific channel."
    #         return opt0
    #     else:
    #         assert (channel in [0, 3])
    #         return self.Cut[channel].get_min_event()

    def GetMaxEventCut(self, channel=None):
        '''
        Returns the highest event number, which satisfies the cut
        conditions. It is recommended to handle the channel ID. If no
        channel will be passed and the event numbers differ for both
        channels, an assertion error will be raised.
        :param channel:
        :return:
        '''
        if channel == None:
            opt0 = self.Cut[0].get_max_event()
            opt3 = self.Cut[3].get_max_event()
            assert (opt0 == opt3), "GetMaxEvent Not the same for both cut channels. Choose a specific channel."
            return opt0
        else:
            assert (channel in [0, 3])
            return self.Cut[channel].get_max_event()

    def GetNEventsCut(self, channel=None):
        '''
        Returns the number of events, which satisfies the cut
        conditions. It is recommended to handle the channel ID. If no
        channel will be passed and the event numbers differ for both
        channels, an assertion error will be raised.
        :param channel:
        :return:
        '''
        if channel == None:
            opt0 = self.Cut[0].get_n_events()
            opt3 = self.Cut[3].get_n_events()
            assert (opt0 == opt3), "GetNEvents Not the same for both cut channels. Choose a specific channel."
            return opt0
        else:
            assert (channel in [0, 3])
            return self.Cut[channel].get_n_events()

    def ShowSignalHisto(self, channel=None, canvas=None, drawoption="", cut="", color=None, normalized=True, drawruninfo=True, binning=600, xmin=None, xmax=None, savePlots=False, logy=False,
                        gridx=False):
        '''
        Shows a histogram of the signal responses for the selected
        channels.
        :param channel:
        :param canvas:
        :param drawoption:
        :param cut:
        :param color:
        :param normalized:
        :param drawruninfo:
        :param binning:
        :param xmin:
        :param xmax:
        :param savePlots:
        :param logy:
        :param gridx:
        :return:
        '''
        self._ShowHisto(self.GetSignalDefinition(channel=channel), channel=channel, canvas=canvas, drawoption=drawoption, cut=cut, color=color, normalized=normalized, infoid="SignalHisto",
                        drawruninfo=drawruninfo, binning=binning, xmin=xmin, xmax=xmax, savePlots=savePlots, logy=logy, gridx=gridx)

    def ShowPedestalHisto(self, channel=None, canvas=None, drawoption="", color=None, normalized=True, drawruninfo=False, binning=600, xmin=None, xmax=None, savePlots=False, logy=False, cut=""):
        '''
        Shows a histogram of the pedestal signal responses for the
        selected channels.
        :param channel:
        :param canvas:
        :param drawoption:
        :param color:
        :param normalized:
        :param drawruninfo:
        :param binning:
        :param xmin:
        :param xmax:
        :param savePlots:
        :param logy:
        :return:
        '''
        if channel != None:
            pedestalname = self.pedestaldefinition[channel]
        else:
            pedestalname = self.pedestalname + "[{channel}]"
        self._ShowHisto(pedestalname, channel=channel, canvas=canvas, drawoption=drawoption, color=color, cut=cut, normalized=normalized, infoid="PedestalHisto", drawruninfo=drawruninfo,
                        binning=binning, xmin=xmin, xmax=xmax, savePlots=savePlots, logy=logy)

    def ShowMedianHisto(self, channel=None, canvas=None, drawoption="", color=None, normalized=True, drawruninfo=True, binning=100, xmin=-30, xmax=30, savePlots=False, logy=True):
        '''
        Shows a histogram of the median for the selected
        channels.
        :param channel:
        :param canvas:
        :param drawoption:
        :param color:
        :param normalized:
        :param drawruninfo:
        :param binning:
        :param xmin:
        :param xmax:
        :param savePlots:
        :param logy:
        :return:
        '''
        self.medianhisto = self._ShowHisto("median[{channel}]", logy=logy, infoid="median", drawruninfo=drawruninfo, xmin=xmin, xmax=xmax, binning=binning, channel=channel, canvas=canvas,
                                           drawoption=drawoption, color=color, normalized=normalized, savePlots=False)
        if canvas == None:
            canvas = ROOT.gROOT.FindObject("mediancanvas")

        canvas.SetGridx()

        if self.medianhisto:
            self.medianhisto.GetXaxis().SetNdivisions(20)

        canvas.Update()
        if savePlots:
            if channel != None:
                dia = "_" + self.run.diamond_names[channel]
            else:
                dia = ""
            self.save_plots("Run{run}_MedianHisto{dia}.png".format(run=self.run.run_number, dia=dia))

    def ShowSignalPedestalHisto(self, channel, canvas=None, savePlots=True, cut="", normalized=True, drawruninfo=True, binning=600, xmin=None, xmax=None, logy=False, gridx=True):
        if canvas == None:
            self.signalpedestalcanvas = ROOT.TCanvas("signalpedestalcanvas{run}{channel}".format(run=self.run.run_number, channel=channel), "signalpedestalcanvas")
        else:
            self.signalpedestalcanvas = canvas

        self.ResetColorPalette()
        self.ShowPedestalHisto(channel=channel, canvas=self.signalpedestalcanvas, savePlots=False, cut=cut, normalized=normalized, drawruninfo=False, binning=binning, xmin=xmin, xmax=xmax, logy=logy)
        self.ShowSignalHisto(channel=channel, canvas=self.signalpedestalcanvas, drawoption="sames", savePlots=False, cut=cut, normalized=normalized, drawruninfo=False, binning=binning, xmin=xmin,
                             xmax=xmax, logy=logy, gridx=gridx)

        if drawruninfo:
            self.DrawRunInfo(channel=channel, canvas=self.signalpedestalcanvas, infoid="signalpedestal" + str(self.run.run_number) + str(channel))

        if savePlots:
            for ending in ["png", "root"]:
                self.save_plots(savename="Run{run}_SignalPedestal_Ch{channel}.{end}".format(run=self.run.run_number, channel=channel, end=ending), sub_dir=ending, canvas=self.signalpedestalcanvas)

        self.if_wait("ShowSignalPedeslHisto")

    def _ShowHisto(self, signaldef, channel=None, canvas=None, drawoption="", cut="", color=None, normalized=True, infoid="histo", drawruninfo=False, binning=600, xmin=None, xmax=None,
                   savePlots=False, logy=False, gridx=False):
        '''

        :param channel:
        :param canvas:
        :param drawoption:
        :param normalized:
        :return:
        '''
        channels = self.GetChannels(channel=channel)

        if canvas == None:
            canvas = ROOT.TCanvas(infoid + "canvas")
            self.ResetColorPalette()
        else:
            pass
            # drawoption = "sames"
        canvas.cd()

        if xmin == None:
            xmin = -100
        if xmax == None:
            xmax = 500

        ROOT.gStyle.SetOptStat(1)

        if color == None: color = self.get_new_color()
        for ch in channels:
            if len(channels) > 1 and drawoption == "" and ch == 3:
                drawoption = "sames"
                color = self.get_new_color()
            if cut == "":
                thiscut = self.GetCut(ch)
                thisusercut = self.GetUserCutString(channel=ch)
            else:
                thiscut = cut.format(channel=ch)
                thisusercut = thiscut
            print "making " + infoid + " using\nSignal def:\n\t{signal}\nCut:\n\t({usercut})\n\t{cut}".format(signal=signaldef, usercut=thisusercut, cut=thiscut)
            self.run.tree.Draw(
                (signaldef + ">>{infoid}{run}({binning}, {min}, {max})").format(infoid=(self.run.diamond_names[ch] + "_" + infoid), channel=ch, run=self.run.run_number, binning=binning, min=xmin,
                                                                                max=xmax), thiscut, drawoption, self.GetNEventsCut(channel=ch), self.GetMinEventCut(channel=ch))
            canvas.Update()
            if logy: canvas.SetLogy()
            if gridx: canvas.SetGridx()
            histoname = "{infoid}{run}".format(infoid=(self.run.diamond_names[ch] + "_" + infoid), run=self.run.run_number)
            histo = ROOT.gROOT.FindObject(histoname)

            if histo:
                histo.GetXaxis().SetTitle(signaldef.format(channel=""))
                histo.SetLineColor(color)
                histo.Draw(drawoption)
                print "\n"
                print "__" + infoid + "__"
                print "\tNEvents: ", histo.GetEntries()
                print "\tMean:    ", histo.GetMean()
                print "\n"
                histo.SetTitle("{signal} {cut}".format(signal=infoid, cut="{" + thisusercut + "}"))
                stats = histo.FindObject("stats")
                if stats:
                    stats.SetTextColor(color)
                    if channels.index(ch) == 1:  # shift second statbox down
                        point = stats.GetBBoxCenter()
                        point.SetY(150)
                        stats.SetBBoxCenter(point)
            canvas.Modified()

            if normalized:
                histo.Scale(1. / histo.GetMaximum())
                histo.Draw(drawoption)

        if drawruninfo:
            if len(channels) == 2:
                self.DrawRunInfo(canvas=canvas)
            else:
                self.DrawRunInfo(channel=channels[0], canvas=canvas)
        canvas.Update()
        self.if_wait(infoid + " shown")
        if savePlots:
            self.save_plots("Run{run}_{signal}.png".format(run=self.run.run_number, signal=infoid), canvas=canvas, sub_dir=infoid + "/png/")
            self.save_plots("Run{run}_{signal}.root".format(run=self.run.run_number, signal=infoid), canvas=canvas, sub_dir=infoid + "/root/")
        return histo

    def ShowDiamondCurrents(self):
        # todo
        pass

    def LoadTrackData(self, minimum_bincontent=None):  # min_bincontent in config file
        '''
        Create a bin collection object as self.Pads and load data from ROOT TTree
        into the Pad object. Then get the 2-dim signal distribution from self.Pads
        :param minimum_bincontent: Bins with less hits are ignored
        :return: -
        '''
        print "Loading Track information with \n\tmin_bincontent: {mbc}\n\tfirst: {first}\n\tmaxevent: {me}".format(mbc=self.minimum_bincontent, first=self.GetMinEventCut(channel=0),
                                                                                                                    me=self.loadMaxEvent)
        if minimum_bincontent != None: self.minimum_bincontent = minimum_bincontent
        assert (self.minimum_bincontent > 0), "minimum_bincontent has to be a positive integer"  # bins with less hits are ignored

        # create a bin collection object:
        self.Pads = {}
        for ch in [0, 3]:  # self.run.GetChannels():
            self.Pads[ch] = BinCollection(self, ch, *self.config_object[ch].Get2DAttributes())

        # fill two 2-dim histograms to collect the hits and signal strength
        x_ = {}
        y_ = {}
        channels = [0, 3]  # self.run.GetChannels()
        includedCh = {}
        if self.loadMaxEvent > 0:
            included = self.GetIncludedEvents(maxevent=self.loadMaxEvent)  # all event numbers without jump events and initial cut
            includedCh[0] = self.GetIncludedEvents(maxevent=self.loadMaxEvent, channel=0)
            includedCh[3] = self.GetIncludedEvents(maxevent=self.loadMaxEvent, channel=3)
        else:
            included = self.GetIncludedEvents()  # all event numbers without jump events and initial cut
            includedCh[0] = self.GetIncludedEvents(channel=0)
            includedCh[3] = self.GetIncludedEvents(channel=3)

        # if not self.pedestal_correction:
        #     signaldef = "self.run.tree.{signal}[{channel}]"
        # else:
        signaldef = "self.run.tree.{signal}[{channel}] - self.run.tree.{pedestal}[{channel}]"

        print "Loading tracking informations..."
        print "Signal def:\n\t{signal}".format(signal=signaldef)
        try:
            self.run.tree.diam1_track_x
        except AttributeError:
            print "\n\n--------------- ERROR --------------"
            print "NO TRACKING INFORMATION IN ROOT FILE"
            print "------------------------------------\n\n"

        cutfunctions = {}
        # cutfunctions = lambda pulser, is_saturated, n_tracks, fft_mean, INVfft_max: 1
        exec ("cutfunctions[0] = {cutf}".format(cutf=self.Cut[0].get_function_def()))
        exec ("cutfunctions[3] = {cutf}".format(cutf=self.Cut[3].get_function_def()))

        for i in included:
            if i % 10000 == 0: print "processing Event ", i
            # read the ROOT TTree
            self.run.tree.GetEntry(i)
            x_[0] = self.run.tree.diam1_track_x
            y_[0] = self.run.tree.diam1_track_y
            x_[3] = self.run.tree.diam2_track_x
            y_[3] = self.run.tree.diam2_track_y
            for ch in channels:
                signal_ = eval(signaldef.format(signal=self.signal_names[ch], channel=ch, pedestal=self.pedestal_names[ch]))
                pulser = self.run.tree.pulser
                is_saturated = self.run.tree.is_saturated[ch]
                fft_mean = self.run.tree.fft_mean[ch]
                INVfft_max = 1. / self.run.tree.fft_max[ch]
                n_tracks = self.run.tree.n_tracks
                sig_time = self.run.tree.sig_time[ch]
                sig_spread = self.run.tree.sig_spread[ch]
                median = self.run.tree.median[ch]
                if cutfunctions[ch](pulser, is_saturated, n_tracks, fft_mean, INVfft_max, sig_time, sig_spread,
                                    median):  # (not pulser and not is_saturated and fft_mean>50 and fft_mean<500 and INVfft_max>1e-4):
                    self.Pads[ch].Fill(x_[ch], y_[ch], signal_)
        print "Tracking information of ", len(included), " events loaded"


        # self.Pads[channel].MakeFits()
        self.Signal2DDistribution = {}
        for ch in channels:
            self.Signal2DDistribution[ch] = self.Pads[ch].GetMeanSignalDistribution(self.minimum_bincontent)
            self.Signal2DDistribution[ch].SetDirectory(0)
            self.Signal2DDistribution[ch].SetStats(False)
            self.Signal2DDistribution[ch].GetXaxis().SetTitle("pos x / cm")
            self.Signal2DDistribution[ch].GetYaxis().SetTitle("pos y / cm")
            self.Signal2DDistribution[ch].GetYaxis().SetTitleOffset(1.4)

        self.Checklist["LoadTrackData"] = True

    def set_channels(self, diamonds):
        self.run.set_channels(diamonds=diamonds)

    def get_event_at_time(self, time_sec):
        return self.run.get_event_at_time(time_sec)

    def get_flux(self):
        return self.run.get_flux()

    def ShowSignalMaps(self, draw_minmax=True, saveplots=False, savename="Run{run}_SignalMaps", ending="png", saveDir="Results/", show3d=False):
        """
        Shows a 2-dimensional mean signal response map for each channel which is activated for analysis.
        :param saveplots: if True, save the plot
        :param savename: filename if saveplots = True
        :param ending: datatype of file if saveplots = True
        :param saveDir: directory to save the plot - has to end with '/'
        :return: -
        """
        if not self.Checklist['LoadTrackData']:
            self.LoadTrackData()

        channels = self.run.get_active_channels()
        self.signal_canvas = ROOT.TCanvas('signal_canvas{run}', 'Mean Signal Maps', len(channels) * 500, 500)
        self.signal_canvas.Divide(len(channels), 1)
        if len(channels) == 2:
            namesuffix = ''
        else:
            namesuffix = '_Ch{ch}'.format(ch=channels[0])
        self.signal_canvas.cd(1)
        ROOT.SetOwnership(self.signal_canvas, False)

        ROOT.gStyle.SetPalette(53)  # Dark Body Radiator palette
        ROOT.gStyle.SetNumberContours(999)

        for channel in channels:
            pad = self.signal_canvas.cd(channels.index(channel) + 1)
            # Plot the Signal2D TH2D histogram

            if show3d:
                self.Signal2DDistribution[channel].Draw('SPEC dm(2,10) pa(1,1,1) ci(1,1,1) a(15,45,0) s(1,1)')
                savename += '3D'
            else:
                self.Signal2DDistribution[channel].Draw('colz')

            if draw_minmax and not show3d:
                self._DrawMinMax(pad=self.signal_canvas.cd(channels.index(channel) + 1), channel=channel)

            if not show3d: 
                self.DrawRunInfo(canvas=pad, infoid='signalmap{run}{ch}'.format(run=self.run.run_number, ch=channel))

        self.signal_canvas.Update()
        self.if_wait('2d drawn')
        if saveplots:
            savename = savename.format(run=self.run.run_number) + namesuffix
            self.save_plots(savename, ending, canvas=self.signal_canvas, sub_dir=ending, save_dir=saveDir)
            self.save_plots(savename, 'root', canvas=self.signal_canvas, sub_dir="root", save_dir=saveDir)

    def _DrawMinMax(self, pad, channel, theseMaximas=None, theseMinimas=None):
        pad.cd()

        if theseMaximas is None:
            if self.ExtremaResults[channel]['FoundMaxima'] == None: self.FindMaxima(channel=channel)
        else:
            assert (type(theseMaximas) is t.ListType)
        if theseMinimas is None:
            if self.ExtremaResults[channel]['FoundMinima'] == None: self.FindMinima(channel=channel)
        else:
            assert (type(theseMinimas) is t.ListType)

        if self.run.IsMonteCarlo:
            print "Run is MONTE CARLO"
            if self.run.SignalParameters[0] > 0:
                height = self.run.SignalParameters[4]
                self.Pads[channel].GetSignalInRow(height, show=True)
            # self.combined_canvas.cd(1)
            self.Pads[channel].MaximaSearch.real_peaks.SetMarkerColor(ROOT.kBlue)
            self.Pads[channel].MaximaSearch.real_peaks.SetLineColor(ROOT.kBlue)
            self.Pads[channel].MaximaSearch.real_peaks.Draw('SAME P0')
        # self.Pads[channel].MaximaSearch.found_extrema.SetMarkerColor(ROOT.kGreen+2)
        # self.Pads[channel].MaximaSearch.found_extrema.Draw('SAME P0')
        if theseMinimas == None:
            minima = self.ExtremaResults[channel]['FoundMinima']
        else:
            minima = theseMinimas

        for i in xrange(len(minima)):
            text = ROOT.TText()
            text.SetTextColor(ROOT.kBlue - 4)
            text.DrawText(minima[i][0] - 0.01, minima[i][1] - 0.005, 'low')

        if theseMaximas == None:
            maxima = self.ExtremaResults[channel]['FoundMaxima']
        else:
            maxima = theseMaximas

        for i in xrange(len(maxima)):
            text = ROOT.TText()
            text.SetTextColor(ROOT.kRed)
            text.DrawText(maxima[i][0] - 0.02, maxima[i][1] - 0.005, 'high')
        if len(maxima) * len(minima) > 0:
            maxbin = self.Pads[channel].GetBinByCoordinates(*(maxima[0]))
            maxbin.FitLandau()
            minbin = self.Pads[channel].GetBinByCoordinates(*(minima[0]))
            minbin.FitLandau()
            print '\nApproximated Signal Amplitude: {0:0.0f}% - (high/low approximation)\n'.format(100. * (maxbin.Fit['MPV'] / minbin.Fit['MPV'] - 1.))

    def ShowExtendedSignalMaps(self, draw_minmax=True, saveplots=False, savename="SignalDistribution", ending="png", saveDir=None, PS=False, test=""):
        '''
        Shows the 2-dimensional mean signal map, as well as the
        corresponding histogram of the mean signals inside this map.
        (i.e. the mean of each bin corresponds to one entry in the
        histogram)
        :param draw_minmax:
        :param saveplots:
        :param savename:
        :param ending:
        :param saveDir:
        :param PS:
        :param test:
        :return:
        '''
        self.combined_canvas_name = "combined_canvas" + test
        self.combined_canvas = ROOT.gROOT.GetListOfCanvases().FindObject(self.combined_canvas_name)

        channels = self.run.get_active_channels()

        if not self.combined_canvas:
            self.combined_canvas = ROOT.TCanvas(self.combined_canvas_name, "Combined Canvas", 1000, len(channels) * 500)
            ROOT.SetOwnership(self.combined_canvas, False)
            self.combined_canvas.Divide(2, len(channels))

        for channel in channels:
            # self.CreatePlots(False)
            self.CreateMeanSignalHistogram(channel=channel, saveplots=False, show=False)

            self.combined_canvas.cd(1 + channels.index(channel) * 2)  # 2D Signal Map Pad
            # ROOT.gStyle.SetPalette(55)
            # ROOT.gStyle.SetNumberContours(999)
            ROOT.gStyle.SetPalette(53)
            ROOT.gStyle.SetNumberContours(999)
            self.Signal2DDistribution[channel].Draw("colz")  # "CONT1Z")#)'colz')
            if draw_minmax:
                self._DrawMinMax(pad=self.combined_canvas.cd(1 + channels.index(channel) * 2), channel=channel)

            self.combined_canvas.cd(2 + channels.index(channel) * 2)  # Histo Pad

            if PS:  # if photoshop mode, fill histogram pink
                ROOT.gStyle.SetHistFillColor(6)
                ROOT.gStyle.SetHistFillStyle(1001)
            else:
                ROOT.gStyle.SetHistFillColor(7)
                ROOT.gStyle.SetHistFillStyle(3003)

            self.MeanSignalHisto[channel].UseCurrentStyle()
            self.MeanSignalHisto[channel].GetXaxis().SetTitle("Signal response")
            self.MeanSignalHisto[channel].Draw()
        self.combined_canvas.cd()
        self.combined_canvas.Update()

        savename = self.run.diamond_names[channel] + "_" + savename + "_" + str(self.run.run_number)  # diamond_irradiation_savename_runnr
        if saveplots:
            self.save_plots(savename, ending, saveDir, canvas=self.combined_canvas)
            self.save_plots(savename, "root", saveDir, canvas=self.combined_canvas)
        if PS:
            ROOT.gStyle.SetHistFillColor(7)
            ROOT.gStyle.SetHistFillStyle(3003)
        self.if_wait("Combined 2D Signal DistributionsShown")

    def FindMaxima(self, channel, show=False):
        '''
        Conducts a maxima search using the class FindMaxima. The Results
        are stored in the dict:
            self.extremaResults[channel]
        :param channel:
        :param show:
        :return:
        '''
        if not hasattr(self, "Pads"):
            self.LoadTrackData()
        self.Pads[channel].FindMaxima(minimum_bincount=self.minimum_bincontent, show=show)

    def FindMinima(self, channel, show=False):
        if not hasattr(self, "Pads"):
            self.LoadTrackData()
        self.Pads[channel].FindMinima(minimum_bincount=self.minimum_bincontent, show=show)

    def GetMPVSigmas(self, channel, show=False):
        '''
        Returns MPVs and sigmas form all Bins of current run as Lists.
        If shown=True then a scatter plot is shown
        :return:
        '''
        if not hasattr(self, "Pads"):
            self.LoadTrackData()

        MPVs = []
        MPVErrs = []
        Sigmas = []
        SigmaErrs = []
        for bin in self.Pads[channel].listOfBins:
            entries = bin.GetEntries()
            if entries >= self.minimum_bincontent:
                if bin.Fit["MPV"] == None:
                    bin.FitLandau()
                MPVs.append(bin.Fit["MPV"])
                MPVErrs.append(bin.Fit["MPVErr"])
                Sigmas.append(bin.Fit["Sigma"])
                SigmaErrs.append(bin.Fit["SigmaErr"])

        if show:
            canvas = ROOT.TCanvas("MPVSigmaCanvas", "MPVSigmaCanvas")
            canvas.cd()
            graph = ROOT.TGraphErrors()
            graph.SetNameTitle("mpvsigmascatterplot", "Landau Fit Parameters of all Bins")
            count = 0
            for i in xrange(len(MPVs)):
                graph.SetPoint(count, MPVs[i], Sigmas[i])
                graph.SetPointError(count, MPVErrs[i], SigmaErrs[i])
                count += 1
            graph.Draw("AP")
            graph.GetYaxis().SetTitle("Sigma of Landau fit")
            graph.GetXaxis().SetTitle("MPV of Landau fit")
            self.DrawRunInfo(channel=channel)
            canvas.Update()
            self.if_wait("MPV vs Sigma shown")

        return MPVs, Sigmas, MPVErrs, SigmaErrs

    def ShowTotalSignalHistogram(self, channel, save=False, scale=False, showfit=False):
        '''
        The same as self.ShowSignalHisto(), but concerning only data
        with tracking information (i.e. the data that is stored in
        self.Pads[channel]).
        In addition,  this method also has the Langaus fit implemented.
        :param channel: 
        :param save: 
        :param scale: 
        :param showfit: produce a 'Langaus' fit
        :return:
        '''
        if hasattr(self, "Pads"):
            self.Pads[channel].CreateTotalSignalHistogram(saveplot=save, scale=scale, showfit=showfit)
            # if showfit:
            #     self.GetSignalHistoFitResults() # Get the Results from self.signalHistoFitResults to self.signalHistoFitResults
        else:
            self.LoadTrackData()
            self.Pads[channel].CreateTotalSignalHistogram(saveplot=save, scale=scale, showfit=showfit)

    def GetSignalHistoFitResults(self, channel, show=False):
        '''
        Returns fit results of the signal histogram fit ('Langau'). If
        the fits of the signal histogram does not exist, it will create
        it and show the histogram containing the fit.
        :param channel: 
        :return: FitFunction, Peak, FWHM, Chi2, NDF
        '''
        if self.SignalHistoFitResults[channel]["Peak"] != None:
            FitFunction = self.SignalHistoFitResults[channel]["FitFunction"]
            Peak = self.SignalHistoFitResults[channel]["Peak"]
            FWHM = self.SignalHistoFitResults[channel]["FWHM"]
            Chi2 = self.SignalHistoFitResults[channel]["Chi2"]
            NDF = self.SignalHistoFitResults[channel]["NDF"]
            if show: self.ShowTotalSignalHistogram(channel=channel, save=False, showfit=True)
            return FitFunction, Peak, FWHM, Chi2, NDF
        else:
            self.ShowTotalSignalHistogram(channel=channel, save=False, showfit=True)
            if (self.SignalHistoFitResults[channel]["Peak"] != None) or (hasattr(self, "ExtremeAnalysis") and self.SignalHistoFitResults[channel]["Peak"] != None):
                self.GetSignalHistoFitResults()
            else:
                assert (False), "BAD SignalHistogram Fit, Stop program due to possible infinity loop"

    def SignalTimeEvolution(self, channel, Mode="Mean", show=True, time_spacing=3, save=True, binnumber=None, RateTimeEvolution=False,
                            nameExtension=None):  # CUTS!! not show: save evolution data / Comment more
        '''
        Creates Signal vs time plot. The time is bunched into time buckets of width time_spaceing,
        from all Signals inside one time bucket the Mean or the MPV is evaluated and plotted in a TGraph.
        The Mean is taken from the histogramm of the particular time bucket. The MPV is simply the
        Position (Center) of the Bin containing the most entries in the time bucket histogram.
        :param Mode: What should be plotted as y: either "Mean" or "MPV"
        :param show:
        :param time_spacing: time bucket width in minutes
        :param save:
        :return:
        '''
        False  # MAKE FASTER USING self.GetEventAtTime()
        # assert(Mode in ["Mean", "MPV"]), "Wrong Mode, Mode has to be `Mean` or `MPV`"
        # if not hasattr(self, "Pads"):
        #     self.LoadTrackData()
        # if True: #not self.run.IsMonteCarlo:
        #
        #     # Names
        #     if binnumber != None:
        #         binCenter = self.Pads[channel].GetBinCenter(binnumber)
        #         if nameExtension is None:
        #             nameExtension = "_Bin{0:.3f}_{1:.3f}".format(*binCenter)
        #     else:
        #         binCenter = 0,0
        #         if nameExtension is None:
        #             nameExtension = "OverAll"
        #     if RateTimeEvolution:
        #         type_ = "Rate"
        #     else:
        #         type_ = "Signal"
        #
        #     results = {}
        #     self.SignalEvolution = ROOT.TGraphErrors()
        #     self.SignalEvolution.SetNameTitle(type_+"Evolution", type_+" Time Evolution "+nameExtension)
        #
        #     self.run.tree.GetEntry(1)
        #     starttime = self.run.tree.time
        #     self.run.tree.GetEntry(self.run.tree.GetEntries())
        #     # endtime = self.run.tree.time
        #     ok_deltatime = 0
        #     TimeERROR = False
        #     for i in xrange(self.run.tree.GetEntries()-1):
        #         i += 1
        #         # read the ROOT TTree
        #         self.run.tree.GetEntry(i)
        #         signal_ = self.run.tree.sig_spread[0]
        #         time_ = self.run.tree.time
        #         deltatime_min = (time_-starttime)/60000.
        #         if deltatime_min < -1: # possible time_stamp reset
        #             deltatime_min = ok_deltatime + time_/60000.
        #             TimeERROR = True
        #         else:
        #             ok_deltatime = deltatime_min
        #         time_bucket = int(deltatime_min)/int(time_spacing)*int(time_spacing)
        #         pulser = self.run.tree.pulser
        #
        #         # check if hit is inside bin (ONLY FOR BIN SIGNAL TIME EVOLUTION)
        #         if binnumber != None:
        #             x_ = self.run.tree.diam1_track_x
        #             y_ = self.run.tree.diam1_track_y
        #             binnumber_ = self.Pads[channel].GetBinNumber(x_, y_)
        #         else:
        #             binnumber_ = None
        #
        #         if not pulser and (binnumber == None or binnumber == binnumber_ ):
        #             try:
        #                 results[time_bucket].append(signal_)
        #             except KeyError:
        #                 results[time_bucket] = [signal_]
        #     if TimeERROR:
        #         print "\n\n\n\n\nWARNING: Error occured in time flow of Run "+str(self.run.run_number)+". Possible reset of timestamp during data taking!\n\n\n\n\n"
        #
        #     time_buckets = results.keys()
        #     time_buckets.sort()
        #
        #     # drop last bucket if Rate Time Scan (last bucket may be overshooting the data time window)
        #     if RateTimeEvolution:
        #         time_buckets = time_buckets[:-1]
        #
        #     count = 0
        #     _, xmin, xmax, _, ymin, ymax =  self.Pads[channel].Get2DAttributes()
        #     area = (xmax-xmin)*(ymax-ymin)
        #     for t in time_buckets:
        #         histo = ROOT.TH1D(type_+"Evolution_time_"+str(t), type_+" Time Histogram for t = {0:0.0f}-{1:0.0f}min".format(t, t+int(time_spacing)), 500, 0, 500)
        #         for i in xrange(len(results[t])):
        #             histo.Fill(results[t][i])
        #         if Mode == "Mean" or RateTimeEvolution:
        #             if RateTimeEvolution:
        #                 if binnumber != None:
        #                     c = 1./(60.*time_spacing*(self.config_object[channel].config["2DHist"]["binsize"])**2)
        #                     N = histo.GetEntries()
        #                     self.SignalEvolution.SetPoint(count, t, c*N) # Rate/cm^2 = Entries/(seconds*(binsize)**2)
        #                 else:
        #                     c = 1./(60.*time_spacing*area)
        #                     N = histo.GetEntries()
        #                     self.SignalEvolution.SetPoint(count, t, c*N) # Rate/cm^2 = Entries/(seconds*Area)
        #                 self.SignalEvolution.SetPointError(count,0,c*np.sqrt(N))
        #                 signalname = "Rate / Hz/cm^2"
        #             else:
        #                 self.SignalEvolution.SetPoint(count, t, histo.GetMean())
        #                 self.SignalEvolution.SetPointError(count, 0, histo.GetRMS()/np.sqrt(histo.GetEntries()))
        #                 signalname = "Mean of Signal Response"
        #         elif Mode == "MPV":
        #             self.SignalEvolution.SetPoint(count, t, histo.GetBinCenter(histo.GetMaximumBin()))
        #             self.SignalEvolution.SetPointError(count, 0, 0)
        #             signalname = "MPV of Signal Response"
        #         count += 1
        #         del histo
        #         ROOT.gROOT.Delete(type_+"Evolution_time_"+str(t))
        #
        #     self.SignalEvolution.GetXaxis().SetTitle("Time / min")
        #     self.SignalEvolution.GetYaxis().SetTitle(signalname)
        #     self.SignalEvolution.GetYaxis().SetTitleOffset(1.1)
        #     if RateTimeEvolution: # Start y axis from 0 when rate time evolution
        #         if binnumber != None: # just one bin
        #             self.SignalEvolution.GetYaxis().SetRangeUser(0, 4000)#self.SignalEvolution.GetYaxis().GetXmax())
        #         else: #overall
        #             self.SignalEvolution.GetYaxis().SetRangeUser(0, 1700)#self.SignalEvolution.GetYaxis().GetXmax())
        #     canvas = ROOT.TCanvas("signaltimeevolution", "signaltimeevolution")
        #
        #     self.SignalEvolution.Draw("ALP*")
        #     # Draw a second axis if it is a RateTimeEvolution:
        #     if RateTimeEvolution:
        #         canvas.Update()
        #         rightaxis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUymin()/c, ROOT.gPad.GetUymax()/c, 20210, '+L')
        #         rightaxis.SetTitle('# Hits')
        #         rightaxis.SetTitleOffset(1.2)
        #         ROOT.SetOwnership(rightaxis, False)
        #         rightaxis.Draw('SAME')
        #     if save:
        #         self.SavePlots(type_+"TimeEvolution"+Mode+nameExtension+".png", canvas=canvas)
        #     if TimeERROR:
        #         self.SignalEvolution.SaveAs(self.save_directory+"ERROR_"+type_+"TimeEvolution"+Mode+nameExtension+".root")
        #     self.IfWait("Showing "+type_+" Time Evolution..")
        #     canvas.Close()
        # else:
        #     print "Run is Monte Carlo. Signal- and Rate Time Evolution cannot be created."

    def _ShowPreAnalysisOverview(self, channel=None, savePlot=True):

        channels = self.GetChannels(channel=channel)

        for ch in channels:
            self.pAOverviewCanv = TCanvas("PAOverviewCanvas", "PAOverviewCanvas", 1500, 900)
            self.pAOverviewCanv.Divide(2, 1)

            PApad = self.pAOverviewCanv.cd(1)
            rightPad = self.pAOverviewCanv.cd(2)
            rightPad.Divide(1, 3)

            PApad.cd()
            self.MakePreAnalysis(channel=ch, savePlot=True, canvas=PApad)

            upperpad = rightPad.cd(1)
            upperpad.Divide(2, 1)
            pulserPad = upperpad.cd(1)
            self.ShowPulserRate(canvas=pulserPad, savePlot=False)
            fftPad = upperpad.cd(2)
            self.ShowFFT("mc", cut=True, channel=ch, savePlots=False, canvas=fftPad)

            self.ResetColorPalette()
            middlepad = rightPad.cd(2)
            middlepad.Divide(2, 1)
            spreadPad = middlepad.cd(1)
            self._ShowHisto(signaldef="sig_spread[{channel}]", channel=ch, canvas=spreadPad, infoid="Spread", drawruninfo=True, savePlots=False, logy=True, gridx=True, binning=150, xmin=0, xmax=150)
            medianPad = middlepad.cd(2)
            self.ShowMedianHisto(channel=ch, canvas=medianPad)

            lowerPad = rightPad.cd(3)
            lowerPad.Divide(2, 1)
            peakPosPad = lowerPad.cd(1)
            self.ShowPeakPosition(channel=ch, canvas=peakPosPad, savePlot=False)
            sigpedPad = lowerPad.cd(2)
            self.CalculateSNR(channel=ch, name="", savePlots=False, canvas=sigpedPad)

            if savePlot:
                self.save_plots(savename="Run{run}_PreAnalysisOverview_{dia}.png".format(run=self.run.run_number, dia=self.run.diamond_names[ch]), sub_dir=self.run.diamond_names[ch],
                                canvas=self.pAOverviewCanv)

    def SetIndividualCuts(self, showOverview=True, savePlot=False):
        '''
        With this method the run-specific cuts can be configured.
        The cuts can be set in the terminal window, where a canvas pops-
        up in order to decide on the cut configurations.
        The individual cut configuration is stored in
            Configuration/Individual_Configs/
        and will be applied the next time an Analysis object with the
        same run number is created.
        In order to ignore the individual cuts, the file has to be
        removed from the directory.
        :param showOverview:
        :param savePlot:
        :return:
        '''
        default_dict_ = {  #
                           "EventRange": None,  # [1234, 123456]
                           "ExcludeFirst": None,  # 50000 events
                           "noPulser": None,  # 1: nopulser, 0: pulser, -1: no cut
                           "notSaturated": None,
                           "noBeamInter": None,
                           "FFT": None,
                           "Tracks": None,
                           "peakPos_high": None,
                           "spread_low": None,
                           "absMedian_high": None
                           }
        self.individualCuts = {
            0: copy.deepcopy(default_dict_),
            3: copy.deepcopy(default_dict_)
        }
        try:
            for ch in [0, 3]:
                if showOverview:
                    self._ShowPreAnalysisOverview(channel=ch, savePlot=savePlot)
                range_min = raw_input("{dia} - Event Range Cut. Enter LOWER Event Number: ".format(dia=self.run.diamond_names[ch]))
                range_max = raw_input("{dia} - Event Range Cut. Enter UPPER Event Number: ".format(dia=self.run.diamond_names[ch]))
                peakPos_high = raw_input("{dia} - Peak Position Cut. Enter maximum Peak Position Sample Point: ".format(dia=self.run.diamond_names[ch]))
                spread_low = raw_input("{dia} - Spread Cut. Enter minimum Spread (max-min): ".format(dia=self.run.diamond_names[ch]))
                absMedian_high = raw_input("{dia} - Median Cut. Enter maximum abs(median) value: ".format(dia=self.run.diamond_names[ch]))

                if range_max != "" and range_min != "":
                    self.individualCuts[ch]["EventRange"] = [int(range_min), int(range_max)]
                elif range_max != "":
                    self.individualCuts[ch]["EventRange"] = [0, int(range_max)]
                elif range_min != "":
                    self.individualCuts[ch]["ExcludeFirst"] = int(range_min)

                if peakPos_high != "":
                    self.individualCuts[ch]["peakPos_high"] = int(peakPos_high)

                if spread_low != "":
                    self.individualCuts[ch]["spread_low"] = int(spread_low)

                if absMedian_high != "":
                    self.individualCuts[ch]["absMedian_high"] = int(absMedian_high)
        except:
            print "Aborted."
        else:
            print "Creating run-specific config file:"
            path = "Configuration/Individual_Configs/"
            filename = "{testcp}_Run{run}.json".format(testcp=self.TESTCAMPAIGN, run=self.run.run_number)
            print "\t" + path + filename

            f = open(path + filename, "w")
            json.dump(self.individualCuts, f, indent=2, sort_keys=True)
            f.close()
            print "done."

    def ShowWaveFormsPulser(self, nevents=1000, startevent=None, channels=None):
        nWFChannels, draw_waveforms = self.__check_wv_channels()
        self.pulserWaveformCanvas = ROOT.TCanvas("pulserWaveformCanvas", "Pulser Waveform Canvas", 1500, nWFChannels * 300)
        self.pulserWaveformCanvas.cd()
        self.pulserWaveformCanvas.Divide(2, 1)
        pad1 = self.pulserWaveformCanvas.cd(1)
        pad2 = self.pulserWaveformCanvas.cd(2)
        self.ShowWaveForms(nevents=nevents, cut="!pulser", startevent=startevent, channels=channels, canvas=pad1, infoid="notpulser")
        # self.pulserWaveformCanvas.cd(2)
        # pad2 = self.pulserWaveformCanvas.cd(2)
        # pad2.Clear()
        self.ShowWaveForms(nevents=nevents, cut="pulser", startevent=startevent, channels=channels, canvas=pad2, infoid="pulser")

    def GetUserCutString(self, channel=None):
        '''
        Returns a short, more user-friendly cut string, which can be
        used to display the cut configuration as terminal prompt or to
        mention in canvases.
        :param channel:
        :return:
        '''
        if channel == None:
            opt0 = self.Cut[0].get_easy_cutstring()
            opt3 = self.Cut[3].get_easy_cutstring()
            assert (opt0 == opt3), "GetUserCutString Not the same for both cut channels. Choose a specific channel."
            return opt0
        else:
            assert (channel in [0, 3])
            return self.Cut[channel].get_easy_cutstring()

    def ShowPeakPosition(self, channel=None, cut="", canvas=None, savePlot=True):
        '''
        Generates a plot, which shows the locations of the peak
        positions in the waveforms. In the 2-dimensional plot the pulse-
        height vs sample point of the peak is drawn.
        :param channel:
        :param cut:
        :param canvas:
        :return:
        '''
        if channel == None:
            channels = self.run.get_active_channels()
            namesuffix = ""
        else:
            channels = [channel]
            namesuffix = "_ch{ch}".format(ch=channel)

        min_ = self.run.signalregion_low - 10
        max_ = self.run.signalregion_high + 10
        binning = max_ - min_

        if canvas == None:
            canvas = ROOT.TCanvas("{run}signalpositioncanvas".format(run=self.run.run_number), "{run}signalpositioncanvas".format(run=self.run.run_number), 800, len(channels) * 300)
        canvas.Divide(1, len(channels))

        ROOT.gStyle.SetPalette(55)  # rainbow palette
        ROOT.gStyle.SetNumberContours(200)

        for ch in channels:
            pad = canvas.cd(channels.index(ch) + 1)
            if cut == "":
                thiscut = self.GetCut(ch)
            else:
                thiscut = cut.format(channel=ch)

            self.Draw(("(" + self.signaldefinition[ch] + "):sig_time[{channel}]>>signalposition{run}{channel}({bins}, {low}, {high}, 600, -100, 500)").format(channel=ch, run=self.run.run_number,
                                                                                                                                                              bins=binning, low=min_, high=max_),
                      thiscut, "colz")
            hist = ROOT.gROOT.FindObject("signalposition{run}{channel}".format(channel=ch, run=self.run.run_number))
            pad.SetLogz()
            pad.SetGridx()
            if hist:
                hist.SetStats(0)
                hist.SetTitle("Peak Position {" + self.Cut[ch].get_easy_cutstring() + "}")
                hist.GetXaxis().SetTitle("Sample Point of Peak")
                hist.GetXaxis().SetTitleSize(0.05)
                hist.GetXaxis().SetLabelSize(0.05)
                hist.GetXaxis().SetNdivisions(20)
                hist.GetYaxis().SetTitle("Pulse Height ({sigdef})".format(sigdef=self.signaldefinition[ch].format(channel=ch)))
                hist.GetYaxis().SetTitleSize(0.05)
                hist.GetYaxis().SetLabelSize(0.05)
            self.DrawRunInfo(channel=ch, canvas=pad, infoid="peakpos{run}{ch}".format(run=self.run.run_number, ch=ch), userWidth=0.15, userHeight=0.15)

        canvas.Update()
        if savePlot:
            self.save_plots("Run{run}_PeakPosition{ns}.png".format(run=self.run.run_number, ns=namesuffix), canvas=canvas, sub_dir="PeakPosition")
        self.if_wait("Peak Position shown")

    def ShowSignalSpread(self, channel=None, cut=""):
        if channel == None:
            channels = self.run.get_active_channels()
            namesuffix = ""
        else:
            channels = [channel]
            namesuffix = "_ch{ch}".format(ch=channel)

        canvas = ROOT.TCanvas("{run}signalspreadcanvas".format(run=self.run.run_number), "{run}signalspreadcanvas".format(run=self.run.run_number), 800, len(channels) * 300)
        canvas.Divide(1, len(channels))

        for ch in channels:
            canvas.cd(channels.index(ch) + 1)
            if cut == "":
                thiscut = self.GetCut(ch)
            else:
                thiscut = cut.format(channel=ch)

            self.Draw(("sig_spread[{channel}]>>signalspread{run}{channel}(400, 0, 400)").format(channel=ch, run=self.run.run_number), thiscut)
            hist = ROOT.gROOT.FindObject("signalspread{run}{channel}".format(channel=ch, run=self.run.run_number))
            if hist: hist.SetStats(0)

        canvas.Update()
        self.save_plots("Run{run}_SignalSpread{ns}.png".format(run=self.run.run_number, ns=namesuffix), canvas=canvas, sub_dir="Cuts")
        self.if_wait("Peak Position shown")

    def Draw(self, varexp, selection="", drawoption="", nentries=1000000000, firstentry=0):
        '''
        The ROOT TTree Draw method as a shortcut for
        self.run.tree.Draw(varexp, selection, drawoption, nentries,
        firstentry)
        :param varexp:
        :param selection:
        :param drawoption:
        :param nentries:
        :param firstentry:
        :return:
        '''
        self.run.tree.Draw(varexp, selection, drawoption, nentries, firstentry)

    def MakeSNRAnalyis(self, channel, name="", binning=1000, xmin=-500, xmax=500):
        name = name + str(self.run.run_number) + str(channel)
        cut0 = self.GetCut(channel=channel)
        cut = cut0
        SNRs = {}

        channelstring = "[{channel}]"
        windownames = ["a", "b", "c"]
        integralnames = ["1", "2", "3"]
        signaldefs = {
            "a1": "sig_a1[{channel}]-ped_integral1[{channel}]",
            "a2": "sig_a2[{channel}]-ped_integral2[{channel}]",
            "a3": "sig_a3[{channel}]-ped_integral3[{channel}]",
            "b1": "sig_b1[{channel}]-ped_integral1[{channel}]",
            "b2": "sig_b2[{channel}]-ped_integral2[{channel}]",
            "b3": "sig_b3[{channel}]-ped_integral3[{channel}]",
            "c1": "sig_integral1[{channel}]-ped_integral1[{channel}]",
            "c2": "sig_integral2[{channel}]-ped_integral2[{channel}]",
            "c3": "sig_integral3[{channel}]-ped_integral3[{channel}]",
        }

        for windowname in windownames:

            for integralname in integralnames:
                self.CalculateSNR(signaldefinition=signaldefs[windowname + integralname], pedestaldefinition="ped_integral" + integralname + channelstring, name=name + "_" + windowname + integralname,
                                  binning=binning, xmin=xmin, xmax=xmax, channel=channel)
                pedestal_5_sigma_range = [self.pedestalFitMean[channel] - 5 * self.pedestalSigma[channel], self.pedestalFitMean[channel] + 5 * self.pedestalSigma[channel]]
                cut = cut0 + "&&" + self.pedestalname + "[{channel}]>" + str(pedestal_5_sigma_range[0]) + "&&" + self.pedestalname + "[{channel}]<" + str(pedestal_5_sigma_range[1])
                SNRs[windowname + integralname] = self.CalculateSNR(signaldefinition=signaldefs[windowname + integralname], pedestaldefinition="ped_integral" + integralname + channelstring, cut=cut,
                                                                    name=name + "_" + windowname + integralname, binning=binning, xmin=xmin, xmax=xmax, channel=channel)

                if windowname == "c":
                    SNRs[windowname + integralname + "b"] = self.CalculateSNR(signaldefinition=signaldefs[windowname + integralname], pedestaldefinition="ped_integral" + integralname + channelstring,
                                                                              cut=cut + "&&sig_time[{channel}]<250", name=name + "_" + windowname + integralname + "b", binning=binning, xmin=xmin,
                                                                              xmax=xmax, channel=channel)

        signaldefs2 = {
            "spread": "sig_spread[{channel}]",
            "int": "sig_int[{channel}]-ped_int[{channel}]"
        }

        # for key in ["int"]:#signaldefs2:
        #     self.CalculateSNR(signaldefinition=signaldefs2[key], pedestaldefinition="ped_"+key+channelstring, name=name+"_"+key, binning=binning, xmin=xmin, xmax=xmax, channel=channel, fitwindow=40)
        #     pedestal_5_sigma_range = [self.pedestalFitMean[channel]-5*self.pedestalSigma[channel], self.pedestalFitMean[channel]+5*self.pedestalSigma[channel]]
        #     cut = cut0+"&&"+self.pedestalname+"[{channel}]>"+str(pedestal_5_sigma_range[0])+"&&"+self.pedestalname+"[{channel}]<"+str(pedestal_5_sigma_range[1])
        #     SNRs[key] = self.CalculateSNR(signaldefinition=signaldefs2[key], cut=cut, pedestaldefinition="ped_"+key+channelstring, name=name+"_"+key, binning=binning, xmin=xmin, xmax=xmax, channel=channel, fitwindow=40)
        #     SNRs[key+"-b"] = self.CalculateSNR(signaldefinition=signaldefs2[key], pedestaldefinition="ped_"+key+channelstring, name=name+"_"+key+"-b", cut=cut+"&&sig_time[{channel}]<250", binning=binning, xmin=xmin, xmax=xmax, channel=channel, fitwindow=40)

        # key = "spread"
        # self.MakeGlobalPedestalCorrection(channel=channel)
        # cut=cut0
        # SNRs[key] = self.CalculateSNR(name=name+"_"+key, binning=binning, xmin=xmin, xmax=xmax, channel=channel, fitwindow=40)
        # SNRs[key+"-b"] = self.CalculateSNR(name=name+"_"+key+"-b", cut=cut+"&&sig_time[{channel}]<250", binning=binning, xmin=xmin, xmax=xmax, channel=channel, fitwindow=40)

        # for key in SNRs.keys():
        #     print key, " - ", SNRs[key]

        print "OUTPUT:"
        print "a1"
        print "a2"
        print "a3"
        print "b1"
        print "b2"
        print "b3"
        print "d1"
        print "d2"
        print "d3"
        print "d1b"
        print "d2b"
        print "d3b"
        print "spread"
        print "spread-b"
        print "int"
        print "int-b"

        print SNRs["a1"]
        print SNRs["a2"]
        print SNRs["a3"]
        print SNRs["b1"]
        print SNRs["b2"]
        print SNRs["b3"]
        print SNRs["c1"]
        print SNRs["c2"]
        print SNRs["c3"]
        print SNRs["c1b"]
        print SNRs["c2b"]
        print SNRs["c3b"]
        print "spread"  # SNRs["spread"]
        print "spread-b"  # SNRs["spread-b"]
        print "int"  # SNRs["int"]
        print "int-b"  # SNRs["int-b"]

    def MakeGlobalPedestalCorrection(self, channel=None):

        channels = self.GetChannels(channel=channel)

        for ch in channels:
            if not self.Checklist["GlobalPedestalCorrection"][ch]:
                self.CalculateSNR(channel=ch, signaldefinition=None, pedestaldefinition=None, fitwindow=20, savePlots=True, name="GlobalPedestalCorrection")

                self.signaldefinition[ch] = self.signaldefinition[ch] + "-" + str(self.pedestalFitMean[ch])
                self.pedestaldefinition[ch] = self.pedestaldefinition[ch] + "-" + str(self.pedestalFitMean[ch])

                self.Checklist["GlobalPedestalCorrection"][ch] = True

    def GetChannels(self, channel=None):

        if channel == None:
            channels = self.run.get_active_channels()
        else:
            channels = [channel]

        return channels

    def AnalyzePedestalContribution(self, channel, refactor=5):
        self.pedestal_analysis_canvas = ROOT.TCanvas("pedestal_analysis_canvas", "pedestal_analysis_canvas")

        cut = self.GetCut(channel=channel)
        if not hasattr(self, "pedestalFitMean"):
            self.CalculateSNR(channel=channel, signaldefinition=None, pedestaldefinition=None, logy=True, cut="", name="", fitwindow=40, binning=1000, xmin=-500, xmax=500, savePlots=False,
                              canvas=None)
        else:
            if not self.pedestalFitMean.has_key(channel):
                self.CalculateSNR(channel=channel, signaldefinition=None, pedestaldefinition=None, logy=True, cut="", name="", fitwindow=40, binning=1000, xmin=-500, xmax=500, savePlots=False,
                                  canvas=None)

        sigma = self.pedestalSigma[channel]
        pedestalmean = self.pedestalFitMean[channel]

        self.pedestal_n_sigma_range = [pedestalmean - refactor * sigma, pedestalmean + refactor * sigma]

        cut = cut + "&&" + self.pedestalname + "[{channel}]>" + str(self.pedestal_n_sigma_range[0]) + "&&" + self.pedestalname + "[{channel}]<" + str(self.pedestal_n_sigma_range[1])
        self.ShowSignalPedestalHisto(channel, canvas=self.pedestal_analysis_canvas, savePlots=False, cut=cut, normalized=False, drawruninfo=True, binning=1000, xmin=-500, xmax=500, logy=True,
                                     gridx=True)

        histo = ROOT.gROOT.FindObject("{dia}_SignalHisto{run}".format(dia=self.run.diamond_names[channel], run=self.run.run_number))
        self.h_nopedestal = copy.deepcopy(histo)

        landauMax = histo.GetMaximum()
        landauMaxPos = histo.GetBinCenter(histo.GetMaximumBin())

        xmin = pedestalmean - 5 * sigma - 20
        xmax = landauMaxPos + 20
        name = "{dia}_LandauGaus{run}".format(dia=self.run.diamond_names[channel], run=self.run.run_number)

        flandaupedestal = ROOT.TF1('f_%s' % name, 'gaus(0)+landau(3)', xmin, xmax)
        flandaupedestal.SetParLimits(0, -.001, landauMax * 0.1)  # gaus: height
        flandaupedestal.SetParLimits(1, -5, 0)  # gaus: mean
        flandaupedestal.SetParLimits(2, sigma * 0.5, sigma * 1.5)  # gaus: sigma
        flandaupedestal.SetParLimits(3, landauMax * 0.5, landauMax * 10)  # landau: height
        flandaupedestal.SetParLimits(4, landauMaxPos * 0.5, landauMaxPos * 1.5)  # landau: MPV
        flandaupedestal.SetParLimits(5, 5, 50)  # landau: sigma
        # flandaupedestal.FixParameter(1,0)

        histo.Fit(flandaupedestal, "", "", xmin, xmax)
        landaugausfit = histo.GetFunction("f_" + name)

        fpedestal = ROOT.TF1('fped_%s' % name, 'gaus', -500, 500)
        fpedestal.SetParLimits(1, -.001, 0.001)
        fpedestal.FixParameter(1, 0)

        for i in range(3):
            fpedestal.SetParameter(i, landaugausfit.GetParameter(i))

        fpedestal.SetLineColor(ROOT.kGreen)
        fpedestal.Draw("same")

        self.h_nopedestal.Add(fpedestal, -1)
        self.h_nopedestal.SetLineColor(ROOT.kGray + 3)
        self.h_nopedestal.FindObject("stats").SetTextColor(ROOT.kGray + 3)
        self.h_nopedestal.SetName("{dia}_PedCorrected{run}".format(dia=self.run.diamond_names[channel], run=self.run.run_number))
        self.h_nopedestal.Draw("sames")
        self.pedestal_analysis_canvas.Update()

        mean = histo.GetMean()
        mean_nopedestal = self.h_nopedestal.GetMean()

        return mean, mean_nopedestal

    def AnalyzeTrackContribution(self, channel=3):
        '''
        -->config:

        [CUT]
        EventRange_min = 0
        EventRange_max = 0
        cut0 = -1
        cut3 = -1
        notPulser = 1
        notSaturated = True
        noBeamInter = True
        excludeBeforeJump = 2000
        excludeAfterJump = 20000
        FFT = False
        hasTracks = False
        excludefirst = 0
        peakPos_high = -1
        spread_low = -1
        absMedian_high = -1
        pedestalsigma = 5
        beaminterruptions_folder = beaminterruptions
        '''
        c = ROOT.TCanvas("trackanalysis", "trackanalysis")
        self.ResetColorPalette()

        def makeit(name, cut):
            h = self._ShowHisto(self.signaldefinition[channel], channel=channel, canvas=c, drawoption="", cut=self.GetCut(channel) + cut, color=None, normalized=True, infoid=name, drawruninfo=False,
                                binning=600, xmin=None, xmax=None, savePlots=False, logy=False, gridx=False)
            entries = h.GetEntries()
            mean = h.GetMean()
            print "\n"
            print name + ":"
            print "\tNEvents: ", entries
            print "\tMean:    ", mean
            return entries, mean

        # d2:
        d2 = makeit("d2", "")
        # d2 | track:
        d2_t = makeit("d2_track", "&&n_tracks")
        # d2B:
        d2b = makeit("d2b", "&&sig_time[{channel}]<250")
        # d2B:
        d2b_t = makeit("d2b_track", "&&sig_time[{channel}]<250&&n_tracks")

        results_mean = str(self.run.get_flux()) + "\t" + str(d2b[1]) + "\t" + str(d2b_t[1]) + "\t" + str(d2[1]) + "\t" + str(d2_t[1])
        results_events = str(self.run.get_flux()) + "\t" + str(d2b[0]) + "\t" + str(d2b_t[0]) + "\t" + str(d2[0]) + "\t" + str(d2_t[0])

        print "\n"
        print "Rate, ", "mean d2b ", "mean d2b_t, ", "mean d2, ", "mean d2_t"
        print self.run.get_flux(), d2b[1], d2b_t[1], d2[1], d2_t[1]
        print "\n"
        print "\n"
        print "Rate, ", "# d2b ", "# d2b_t, ", "# d2, ", "# d2_t"
        print self.run.get_flux(), d2b[0], d2b_t[0], d2[0], d2_t[0]

        return results_mean, results_events

    def AnalyzeMultiHitContribution(self):

        single_particle = self.run.tree.Draw("1",
                                             "!pulser&&clusters_per_plane[0]>=1&&clusters_per_plane[1]>=1&&clusters_per_plane[2]>=1&&clusters_per_plane[3]>=1&&(clusters_per_plane[0]>1||clusters_per_plane[1]>1||clusters_per_plane[2]>1||clusters_per_plane[3]>1)")
        multiparticle = self.run.tree.Draw("1", "!pulser&&clusters_per_plane[0]<=1&&clusters_per_plane[1]<=1&&clusters_per_plane[2]<=1&&clusters_per_plane[3]<=1")

        total = single_particle + multiparticle

        return 1. * multiparticle / total

if __name__ == "__main__":
    ana_parser = ArgumentParser()
    ana_parser.add_argument('run', nargs='?', default=392, type=int)
    args = ana_parser.parse_args()
    this_run = args.run
    print '\nAnalysing run', this_run, '\n'
    z = Analysis(this_run)


'''
htemp = copy.deepcopy(histo)
h_nopedestal = copy.deepcopy(histo)
xmin = -.5
xmax = 4
fpedestal = ROOT.TF1('fped_%s'%name,'gaus',xmin,xmax)
# print 'set par limits 1'
fpedestal.SetParLimits(1,-.001,0.001)
# print 'fix par limits 1'
fpedestal.FixParameter(1,0)

sigma = pedestal_width*refactor
# print 'SIGMA is: ',pedestal_width,refactor,sigma
# print 'set Par limits 2'
fpedestal.SetParLimits(2,sigma*.01,sigma*1.001)
# print 'fix Par 2'
fpedestal.SetParameter(2,sigma)
# print 'fit'
bin = find_local_minimum_left_to(0,htemp)
max_fit_range = htemp.GetXaxis().GetBinCenter(bin)
htemp.Fit(fpedestal,'QB','',-.2,min(.1,max_fit_range))
htemp.Fit(fpedestal,'QB','',-.2,min(.1,max_fit_range))
fpedestal_full = fpedestal.Clone()
fpedestal_full.SetRange(htemp.GetXaxis().GetXmin(),htemp.GetXaxis().GetXmax())
fpedestal_full.SetLineStyle(2)

h_nopedestal.Add(fpedestal_full,-1,'')
h_nopedestal.SetLineColor(ROOT.kBlue)
h_nopedestal.Draw('')

c = ROOT.TCanvas()
fit = htemp.FindObject("PrevFitTMP")
fit.GetXmin()
fit.GetXmax()
fit2 = ROOT.TF1("fbkg","gaus",-500,500)
for i in range(3):
     fit2.SetParameter(i,fit.GetParameter(i+3))
fit2.SetLineColor(ROOT.kGreen)
fit2.Draw("same")
htemp2 = htemp.Clone()
htemp2.Add(fit2,-1)
htemp2.SetLineColor(ROOT.kBlue)
htemp2.Draw("same")
htemp.GetMean()
htemp2.GetMean()
m = htemp.GetMean()
m2 = htemp2.GetMean()
d = m2-m
d/m*100
'''
