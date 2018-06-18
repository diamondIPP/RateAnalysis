#!/usr/bin/env python
# ==============================================
# IMPORTS
# ==============================================
from argparse import ArgumentParser
from math import log
from sys import stdout
from json import loads

from ROOT import TGraphErrors, TCanvas, TH2D, gStyle, TH1F, gROOT, TLegend, TCut, TGraph, TProfile2D, TH2F, TProfile, TCutG, kGreen, TF1, \
    THStack, TArrow, kOrange, TSpectrum, TMultiGraph, Long, TH2I, gRandom

from CutPad import CutPad
from CurrentInfo import Currents
from Elementary import Elementary
from Extrema import Extrema2D
from TelescopeAnalysis import Analysis
from Pulser import PulserAnalysis
from Pedestal import PedestalAnalysis
from Peaks import PeakAnalysis
from Timing import TimingAnalysis
from Run import Run
from Utils import *

__author__ = 'micha'


# ==============================================
# MAIN CLASS
# ==============================================
class PadAnalysis(Analysis):
    def __init__(self, run, dia, high_low_rate_run=None, binning=10000):

        self.RunNumber = run.RunNumber
        Analysis.__init__(self, run, high_low_rate=high_low_rate_run, binning=binning)
        self.channel = self.load_channel(dia)

        # main
        self.DiamondName = self.run.DiamondNames[dia - 1]
        self.DiamondNumber = dia
        self.Bias = self.run.Bias[dia - 1]
        self.save_dir = '{dia}/{run}/'.format(run=str(self.RunNumber).zfill(3), dia=self.DiamondName)

        # stuff
        if run.tree:
            # polarities
            self.Polarity = self.get_polarity()
            self.PulserPolarity = self.get_pulser_polarity()

            # regions // ranges
            self.IntegralNames = self.get_integral_names()
            self.IntegralRegions = self.load_regions()
            self.SignalRegion = self.IntegralRegions['signal']
            self.PedestalRegion = self.IntegralRegions['pedestal']
            self.PeakIntegral = self.load_peak_integral()

            # names
            self.SignalDefinition = '({pol}*TimeIntegralValues[{num}])'
            self.SignalNumber = self.get_signal_number()
            self.SignalName = self.get_signal_name()
            self.PedestalName = self.get_pedestal_name()

            # cuts
            self.Cut = CutPad(self, self.channel)
            self.AllCuts = self.Cut.all_cut

            # alignment
            self.IsAligned = self.check_alignment(show=False)

            # subclasses
            self.Pulser = PulserAnalysis(self)
            self.Pedestal = PedestalAnalysis(self)
            self.Peaks = PeakAnalysis(self)
            self.Timing = TimingAnalysis(self)

        # currents
        self.Currents = Currents(self)

        # histograms
        self.PedestalHisto = None
        self.SignalTime = None
        self.SignalMapHisto = None
        self.MeanSignalHisto = None
        self.PeakValues = None

    def __del__(self):
        for obj in [self.Pedestal, self.SignalMapHisto, self.SignalTime, self.PeakValues, self.MeanSignalHisto]:
            self.del_rootobj(obj)
        for c in gROOT.GetListOfCanvases():
            c.Close()
        for lst in self.histos + self.RootObjects:
            if not type(lst) is list:
                lst = [lst]
            for obj in lst:
                self.del_rootobj(obj)

    def draw_current(self, relative_time=True, averaging=1):
        self.Currents.draw_indep_graphs(rel_time=relative_time, averaging=averaging)

    # ==========================================================================
    # region INIT

    def load_channel(self, dia):
        assert dia in [1, 2], 'You have to choose either diamond 1 or 2'
        return self.run.Channels[dia - 1]

    def get_integral_names(self):
        if self.run.TreeConfig.has_section('Integral Names'):
            return [str(name) for name in loads(self.run.TreeConfig.get('Integral Names', 'Names'))]
        self.tree.GetEntry(0)
        return [str(name) for name in self.tree.IntegralNames]

    def get_polarity(self):
        self.tree.GetEntry(0)
        return self.tree.polarities[self.channel]

    def get_pulser_polarity(self):
        self.tree.GetEntry(0)
        return self.tree.pulser_polarities[self.channel]

    def load_regions(self):
        all_regions = {}
        for name in ['signal', 'pedestal', 'pulser']:
            option = '{}_region'.format(name)
            region = '{name}_{region}'.format(name=name, region=self.ana_config_parser.get('BASIC', option)) if option in self.ana_config_parser.options('BASIC') else ''
            regions = [reg for reg in self.run.IntegralRegions[self.DiamondNumber - 1] if reg.startswith(name)]
            all_regions[name] = region if region in regions else regions[0]
        return all_regions

    def load_peak_integral(self):
        peak_int = 'PeakIntegral{}'.format(self.ana_config_parser.get('BASIC', 'peak_integral'))
        return peak_int if peak_int in self.run.PeakIntegrals[self.DiamondNumber - 1] else self.run.PeakIntegrals[self.DiamondNumber - 1].keys()[0]

    def get_signal_number(self, region=None, peak_integral=None, sig_type='signal'):
        region = self.IntegralRegions[sig_type] if region is None else self.make_region(sig_type, region)
        peak_integral = self.PeakIntegral if peak_integral is None else 'PeakIntegral{}'.format(peak_integral)
        int_name = 'ch{ch}_{reg}_{int}'.format(ch=self.channel, reg=region, int=peak_integral)
        return self.IntegralNames.index(int_name)

    def get_signal_name(self, region=None, peak_integral=None, sig_type='signal'):
        num = self.get_signal_number(region, peak_integral, sig_type)
        return self.SignalDefinition.format(pol=self.Polarity, num=num)

    def set_signal_definitions(self, use_time=True, sig_region=None, peak_int=None):
        signal = 'TimeIntegralValues' if use_time else 'IntegralValues'
        signal = '({{pol}}*{sig}[{{num}}])'.format(sig=signal)
        print 'changed SignalDefinition to:', signal
        self.SignalDefinition = signal
        self.update_signal_definitions(sig_region, peak_int)

    def update_signal_definitions(self, sig_region=None, peak_int=None):
        self.SignalNumber = self.get_signal_number(sig_region, peak_int)
        self.SignalName = self.get_signal_name(sig_region, peak_int)

    def get_pedestal_name(self, region=None, peak_int=None):
        return self.get_signal_name(region=region, peak_integral=peak_int, sig_type='pedestal')

    def get_peak_name(self, region=None, type_='signal', t_corr=True):
        peak_name = 'IntegralPeakTime' if t_corr else 'IntegralPeaks'
        return '{name}[{num}]'.format(name=peak_name, num=self.get_signal_number(region, None, type_))
    # endregion

    def set_channel(self, ch):
        self.channel = ch
        self.DiamondName = self.run.diamondname[ch]
        self.Bias = self.run.bias[ch]
        self.Cut = CutPad(self, ch)
        self.save_dir = '{tc}_{run}_{dia}'.format(tc=self.TESTCAMPAIGN[2:], run=self.RunNumber, dia=self.run.diamondname[ch])
        self.Polarity = self.get_polarity()
        self.SignalName = self.get_signal_name()

    # ==========================================================================
    # region BEAM PROFILE

    def draw_beam_profile(self, mode='x', show=True, fit=True, fit_margin=.6):
        assert mode.lower() in ['x', 'y'], 'Mode has to be either "x" or "y"!'
        h = deepcopy(self.histos[-1])
        if not show:
            gROOT.SetBatch(1)
        prof = h.ProjectionX() if mode.lower() == 'x' else h.ProjectionY()
        margins = [prof.GetBinLowEdge(prof.FindBin(-.4)), prof.GetBinLowEdge(prof.FindBin(.4) + 1)]
        center = (margins[1] + margins[0]) / 2.
        width = (prof.FindBin(margins[1]) - prof.FindBin(margins[0])) / 2. * fit_margin * prof.GetBinWidth(1)
        fit_range = [center - width, center + width]
        c = TCanvas('c', 'Beam Profile', 1000, 1000)
        c.SetLeftMargin(.145)
        self.format_histo(prof, 'prof', 'Profile ' + mode.title(), y_tit='Entries', y_off=2, x_tit='Track Position {mod} [cm]'.format(mod=mode.title()))
        prof.GetXaxis().SetRangeUser(prof.GetBinCenter(prof.FindFirstBinAbove(0) - 1), prof.GetBinCenter(prof.FindLastBinAbove(0) + 1))
        prof.Draw()
        sleep(.1)
        lines = [self.draw_axis(x, c.GetUymin(), c.GetUymax(), '', 2, 2) for x in margins]
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
        self.histos.append([prof, c, lines, fits])
        return fit_result if fit else prof

    @staticmethod
    def __fit_beam_profile(histo, fit_range, show=True):
        h = histo
        fit = h.Fit('gaus', 'qs{0}'.format('' if show else '0'), '', fit_range[0], fit_range[1])
        return fit

    def fit_beam_profile(self, mode='x', show=True, fit_margin=.6):
        pickle_path = self.PickleDir + 'BeamProfile/Fit{mod}_{tc}_{run}_{dia}_{mar}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.RunNumber, dia=self.DiamondName, mod=mode.title(), mar=fit_margin)

        def func():
            return self.draw_beam_profile(mode=mode, show=show, fit_margin=fit_margin)

        return do_pickle(pickle_path, func)

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
                t = self.draw_tlatex(perc * 100 - 2, y, str(fit.Ndf()), color=807, size=.04, align=32)
                gr.GetListOfFunctions().Add(t)
                index += 1
        c = TCanvas('c', 'Beam Chi2', 1000, 1000)
        self.format_histo(gr, x_tit='Range [%]', y_tit='#chi^{2} / NDF' if not sigma else 'Sigma', y_off=1.4)
        one = TF1('one', '1', 0, 100)
        t1 = self.draw_tlatex(15, .95 * gr.GetYaxis().GetXmax(), 'NDF:', color=807, size=0.04, align=12)
        gr.GetListOfFunctions().Add(t1)
        gr.GetXaxis().SetRangeUser(-5, 105)
        gr.Draw('alp')
        one.Draw('same')

        self.histos.append([gr, c, t1])
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('BeamProf{mod}{dir}'.format(mod='Sigmas' if sigma else 'Chi2s', dir=mode.title()), sub_dir=self.save_dir)

    # endregion

    # ==========================================================================
    # region 2D SIGNAL DISTRIBUTION

    def draw_efficiency_map(self, res=1.5, cut='all', show=True):
        cut_string = TCut(cut) + self.Cut.CutStrings['tracks']
        cut_string = self.Cut.generate_special_cut(excluded=['fiducial']) if cut == 'all' else cut_string
        p = TProfile2D('p_em', 'Efficiency Map {d}'.format(d=self.DiamondName), *self.Plots.get_global_bins(res))
        self.tree.Draw('({s}>10)*100:dia_track_y[{r1}]:dia_track_x[{r1}]>>p_em'.format(s=self.generate_signal_name(), r1=self.DiamondNumber - 1), cut_string, 'goff')
        set_statbox(entries=4, opt=1000000010, x=.81)
        self.format_histo(p, x_tit='Track x [cm]', y_tit='Track y [cm]', z_tit='Efficiency [%]', y_off=1.4, z_off=1.5, ncont=100)
        self.save_histo(p, 'Efficiency Map', show, lm=.13, rm=.17, draw_opt='colz')

    def draw_signal_map(self, res=1.5, cut=None, fid=False, hitmap=False, redo=False, show=True, prnt=True):
        cut = self.Cut.generate_special_cut(excluded=['fiducial'], prnt=prnt) if not fid and cut is None else cut
        cut = self.Cut.all_cut if cut is None else TCut(cut)
        suf = '{c}_{ch}'.format(c=cut.GetName(), ch=self.Cut.CutConfig['chi2X'])
        pickle_path = self.make_pickle_path('SignalMaps', 'Hit' if hitmap else 'Signal', run=self.RunNumber, ch=self.DiamondNumber, suf=suf)

        def func():
            set_root_output(0)
            name = 'h_hm' if hitmap else 'h_sm'
            h1 = TH2I(name, 'Diamond Hit Map', *self.Plots.get_global_bins(res)) if hitmap else TProfile2D(name, 'Signal Map', *self.Plots.get_global_bins(res))
            self.log_info('drawing {mode}map of {dia} for Run {run}...'.format(dia=self.DiamondName, run=self.RunNumber, mode='hit' if hitmap else 'signal '))
            sig = self.generate_signal_name()
            x_var, y_var = (self.Cut.get_track_var(self.DiamondNumber - 1, v) for v in ['x', 'y'])
            self.tree.Draw('{z}{y}:{x}>>{h}'.format(z=sig + ':' if not hitmap else '', x=x_var, y=y_var, h=name), cut, 'goff')
            self.set_dia_margins(h1)
            self.set_ph_range(h1)
            z_tit = 'Number of Entries' if hitmap else 'Pulse Height [au]'
            self.format_histo(h1, x_tit='track_x [cm]', y_tit='track_y [cm]', y_off=1.4, z_off=1.3, z_tit=z_tit, ncont=50, ndivy=510, ndivx=510)
            self.SignalMapHisto = h1
            return h1

        set_statbox(only_entries=True, x=0.82)
        gStyle.SetPalette(1 if hitmap else 53)
        h = func() if redo else None
        h = do_pickle(pickle_path, func, h)
        self.draw_histo(h, '', show, lm=.12, rm=.16, draw_opt='colzsame')
        self.draw_fiducial_cut()
        self.save_canvas(canvas=get_last_canvas(), name='HitMap' if hitmap else 'SignalMap2D', print_names=prnt)
        return h

    def draw_dia_hitmap(self, show=True, res=1.5, cut=None, fid=False, redo=False, prnt=True):
        return self.draw_signal_map(show=show, res=res, cut=cut, fid=fid, hitmap=True, redo=redo, prnt=prnt)

    def set_dia_margins(self, h, size=.3):
        # find centers in x and y
        xmid, ymid = [(p.GetBinCenter(p.FindFirstBinAbove(0)) + p.GetBinCenter(p.FindLastBinAbove(0))) / 2 for p in [h.ProjectionX(), h.ProjectionY()]]
        self.format_histo(h, x_range=[xmid - size, xmid + size], y_range=[ymid - size, ymid + size])

    def set_ph_range(self, h):
        values = [h.GetBinContent(bin_) for bin_ in xrange(h.GetNbinsX() * h.GetNbinsY()) if h.GetBinContent(bin_)]
        m, s = calc_mean(values)
        n_sig = 3
        if 2 * s > m:
            self.format_histo(h, z_range=[min(values), 0.8 * max(values)])
        else:
            self.format_histo(h, z_range=[m - n_sig * s, m + n_sig * s])

    def draw_sig_map_disto(self, show=True, factor=1.5, cut=None, fid=True, redo=False):
        source = self.draw_signal_map(factor, cut, fid, hitmap=False, redo=redo, show=False)
        h = TH1F('h_smd', 'Signal Map Distribution', 100, -50, 350)
        [h.SetBinContent(source.GetBinContent(ibin)) for ibin in xrange(source.GetNbinsX() * source.GetNbinsY()) if source.GetBinContent(ibin)]
        self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=2, fill_color=self.FillColor)
        self.save_histo(h, 'SignalMapDistribution', lm=.15, show=show)

    def draw_sig_map_profiles(self, mode='x', factor=1.5, cut=None, fid=False, hitmap=False, redo=False, show=True):
        s = self.draw_signal_map(factor, cut, fid, hitmap=hitmap, redo=redo, show=False)
        g = self.make_tgrapherrors('g_smp', 'Signal Map Profile')
        values = [[] for _ in xrange(s.GetNbinsX() if mode == 'x' else s.GetNbinsY())]
        for xbin in xrange(s.GetNbinsX()):
            for ybin in xrange(s.GetNbinsY()):
                value = s.GetBinContent(xbin, ybin)
                if value:
                    values[(xbin if mode == 'x' else ybin)].append(value)
        for i, lst in enumerate(values):
            m, sigma = calc_mean(lst) if lst else (0., 0.)
            xval = s.GetXaxis().GetBinCenter(i) if mode == 'x' else s.GetYaxis().GetBinCenter(i)
            g.SetPoint(i, xval, m)
            g.SetPointError(i, 0, sigma)

        self.format_histo(g, x_tit='Track in {m} [cm]'.format(m=mode), y_tit='Pulse Height [au]', y_off=1.5, ndivx=515)
        self.save_histo(g, 'SignalMapProfile', draw_opt='ap', lm=.14, show=show, gridx=True)

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
        # todo: save mean
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
        l = self.draw_tlatex(gr1.GetX()[0], gr1.GetY()[0] + 0.5, 'Errors', align=20, size=0.03)
        gr1.GetListOfFunctions().Add(l)
        if show:
            c = TCanvas('c', 'Mean Signal Distribution', 1000, 1000)
            self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.2)
            h.Draw()
            gr2.Draw('[]')
            gr1.Draw('[]')
            gr2.Draw('p')
            gr1.Draw('p')
            self.save_plots('MeanSignalHisto', sub_dir=self.save_dir)
            self.histos.append([gr1, gr2, c])
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
            self.histos.append([h, c])
        return h

    def fit_mean_signal_distribution(self):
        pickle_path = self.PickleDir + 'MeanSignalFit/{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.RunNumber, dia=self.DiamondName)

        def func():
            self.draw_mean_signal_distribution(show=False)
            return self.MeanSignalHisto.Fit('gaus', 'qs')

        fit = do_pickle(pickle_path, func)
        return fit

    def get_mean_fwhm(self):
        fit = self.fit_mean_signal_distribution()
        conversion_factor = 2 * sqrt(2 * log(2))  # sigma to FWHM
        return fit.Parameter(2) * conversion_factor

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
        self.histos.append(frame)

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
    def draw_peak_timing(self, region=None, sig_type='signal', show=True, cut=None, corr=True, draw_cut=True):
        xmin, xmax = self.run.IntegralRegions[self.DiamondNumber - 1][self.SignalRegion if region is None else self.make_region(sig_type, region)]
        # increase range for timing correction and convert to ns
        fac = 2. if self.run.Digitiser == 'drs4' else 2.5
        xmin = xmin / fac - (10 if corr else 10)
        xmax = xmax / fac + (10 if corr else 10)
        print int((xmax - xmin) * (4 if corr else 1)), xmin, xmax
        h = TH1F('h_pv', '{typ} Peak Positions'.format(typ=sig_type.title()), int((xmax - xmin) * (4 if corr else 4)), xmin, xmax)
        self.format_histo(h, x_tit='Signal Peak Timing [ns]', y_tit='Number of Entries', y_off=1.3, stats=0)
        cut = self.Cut.generate_special_cut(excluded=['timing']) if cut is None else TCut(cut)
        dic = self.Cut.calc_timing_range(show=False)
        t_correction = '({p1}* trigger_cell + {p2} * trigger_cell*trigger_cell)'.format(p1=dic['t_corr'].GetParameter(1), p2=dic['t_corr'].GetParameter(2))
        draw_string = '{peaks}{op}>>h_pv'.format(peaks=self.get_peak_name(region, sig_type, corr), op='/{f}'.format(f=fac) if not corr else '-' + t_correction)
        print draw_string
        self.tree.Draw(draw_string, cut, 'goff')
        self.draw_histo(h, show=show, sub_dir=self.save_dir, lm=.12, logy=True)
        f, fit, fit1 = self.fit_peak_timing(h)
        self.__draw_timing_legends(draw_cut, f, fit, fit1)
        h.Draw('same')
        h.GetXaxis().SetRangeUser(f.Parameter(1) - 5 * f.Parameter(2), f.Parameter(1) + 10 * f.Parameter(2))
        self.save_plots('{typ}PeakPositions'.format(typ=sig_type.title()))
        self.PeakValues = h
        return f

    def __draw_timing_legends(self, draw_cut, f, fit, fit1):
        l1 = None
        if draw_cut:
            l1 = self.make_legend(.66, .7, nentries=3, name='l1')
            g = self.__draw_timing_cut()
            l1.AddEntry(g, 'Timing Cut', 'fl')
            l1.AddEntry(fit1, 'Fitting Range', 'l')
            l1.AddEntry(fit, 'Fit Function', 'l')
            l1.Draw()
        l2 = self.make_legend(.66, nentries=3, name='fr', margin=.05, clean=False)
        l2.SetHeader('Fit Results')
        l2.AddEntry(0, 'Mean:', '')
        l2.AddEntry(0, '{0:5.2f} #pm {1:5.2f} ns'.format(f.Parameter(1), f.ParError(1)), '').SetTextAlign(32)
        l2.AddEntry(0, 'Sigma:', '')
        l2.AddEntry(0, '{0:5.2f} #pm {1:5.2f} ns'.format(f.Parameter(2), f.ParError(2)), '').SetTextAlign(32)
        l2.SetNColumns(2)
        l2.Draw()
        l2.GetListOfPrimitives().First().SetTextAlign(22)
        self.RootObjects.append([l1, l2])

    def __draw_timing_cut(self):
        timing_fit = self.Cut.calc_timing_range(show=False)['timing_corr']
        xmin, xmax = timing_fit.GetParameter(1) - 3 * timing_fit.GetParameter(2), timing_fit.GetParameter(1) + 3 * timing_fit.GetParameter(2)
        g = TCutG('timing', 5)
        g.SetVarX('y')
        g.SetVarY('x')
        ymin, ymax = -10, 1e7
        for i, (x, y) in enumerate([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)]):
            g.SetPoint(i, x, y)
        g.SetLineColor(827)
        g.SetLineWidth(2)
        g.SetFillColor(827)
        g.SetFillStyle(3001)
        g.Draw('f')
        g.Draw('l')
        self.RootObjects.append(g)
        return g

    def fit_peak_timing(self, histo):
        h = histo
        fit1 = h.Fit('gaus', 'qs0')
        mean_, sigma = fit1.Parameter(1), fit1.Parameter(2)
        fit = h.Fit('gaus', 'qs', '', mean_ - sigma, mean_ + sigma)
        fit2 = TF1('f1', 'gaus', mean_ - 5 * sigma, mean_ + 5 * sigma)
        fit3 = TF1('f2', 'gaus', mean_ - sigma, mean_ + sigma)
        pars = [fit.Parameter(i) for i in xrange(3)]
        fit2.SetParameters(*pars)
        fit3.SetParameters(*pars)
        fit3.SetLineWidth(2)
        fit3.SetLineColor(2)
        fit2.SetLineStyle(2)
        fit2.Draw('same')
        self.RootObjects.append([fit2, fit3])
        return fit, fit2, fit3

    def calc_peak_value_fwhm(self):
        pickle_path = self.PickleDir + 'PeakValues/FWHM_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.RunNumber, dia=self.DiamondName)

        def func():
            print 'Getting peak value FWHM for {dia} of run {run}...'.format(run=self.RunNumber, dia=self.DiamondName)
            if self.PeakValues is None:
                self.draw_peak_timing(show=False)
            return self.calc_fwhm(self.PeakValues)

        fwhm = do_pickle(pickle_path, func)
        return fwhm

    def draw_forc_times(self, show=True, corr=False):
        self.tree.Draw('forc_pos', 'forc_pos[0]>20', 'goff')
        htemp = gROOT.FindObject('htemp')
        x = [int(htemp.GetBinCenter(htemp.FindFirstBinAbove(5000))) - 10, int(htemp.GetBinCenter(htemp.FindLastBinAbove(5000))) + 10]
        h = TH1F('ft', 'FORC Timing', x[1] - x[0], x[0] / 2., x[1] / 2.)
        forc = 'forc_pos/2.' if not corr else 'forc_time'
        self.tree.Draw('{forc}>>ft'.format(forc=forc), self.Cut.all_cut, 'goff')
        self.format_histo(h, x_tit='time [ns]', y_tit='Entries', y_off=2, fill_color=17)
        self.histos.append(self.save_histo(h, 'FORCTiming', show, sub_dir=self.save_dir, lm=.14))

    # endregion

    # ==========================================================================
    # region TRIGGER CELL
    def draw_trigger_cell(self, show=True, cut=None):
        h = TH1F('tc', 'Trigger Cell', 1024, 0, 1024)
        cut = self.Cut.all_cut if cut is None else cut
        self.tree.Draw('trigger_cell>>tc', cut, 'goff')
        self.format_histo(h, x_tit='trigger cell', y_tit='Entries', y_off=1.7, fill_color=17)
        h.SetStats(0)
        h.GetYaxis().SetRangeUser(0, h.GetMaximum() * 1.05)
        h.Fit('pol0', 'qs')
        self.histos.append(self.save_histo(h, 'TriggerCell', show, sub_dir=self.save_dir, lm=.11))

    def draw_trigger_cell_vs_peakpos(self, show=True, cut=None, tprofile=False, corr=True, t_corr=False):
        x = self.run.signal_regions[self.SignalRegion]
        if not tprofile:
            ybins = (x[1] - x[0]) if not corr else 4 * (x[1] - x[0])
            h = TH2D('tcpp', 'Trigger Cell vs. Signal Peak Position', 1024, 0, 1024, ybins, x[0] / 2., x[1] / 2.)
        else:
            h = TProfile2D('tcpp', 'Trigger Cell vs. Signal Peak Position', 1024, 0, 1024, x[1] - x[0], x[0] / 2., x[1] / 2.)
        h1 = TProfile('hpr', 'hpr', 100, 0, 1024)

        cut = self.Cut.generate_special_cut(excluded=['timing']) if cut is None else cut
        # cut = self.Cut.all_cut if cut is None else cut
        prof = '' if not tprofile else ':'
        sig = '' if not tprofile else '{sig}-{ped}'.format(sig=self.SignalName, ped=self.Pedestal.SignalName)
        gStyle.SetPalette(55)
        peaks = 'IntegralPeaks[{num}]/2.' if not corr else 'IntegralPeakTime[{num}]'
        peaks = peaks.format(num=self.SignalNumber)
        dic = self.Cut.calc_timing_range(show=False)
        t_correction = '-({p1}* trigger_cell + {p2} * trigger_cell*trigger_cell)'.format(p1=dic['t_corr'].GetParameter(1), p2=dic['t_corr'].GetParameter(2)) if t_corr else ''
        self.tree.Draw('{z}{prof}{peaks}{tc}:trigger_cell>>tcpp'.format(z=sig, prof=prof, peaks=peaks, tc=t_correction), cut, 'goff')
        self.tree.Draw('{peaks}{tc}:trigger_cell>>hpr'.format(peaks=peaks, tc=t_correction), self.AllCuts, 'goff')
        self.format_histo(h, x_tit='trigger cell', y_tit='Signal Peak Timing [ns]', y_off=1.25, z_tit='Pulse Height [au]' if tprofile else 'Number of Entries', z_off=1.2, stats=0)
        self.format_histo(h1, color=1, lw=3)
        h.GetZaxis().SetRangeUser(60, 120) if tprofile else do_nothing()
        fit = h.ProjectionY().Fit('gaus', 'qs0')
        h.GetYaxis().SetRangeUser(fit.Parameter(1) - 4 * fit.Parameter(2), fit.Parameter(1) + 5 * fit.Parameter(2))
        self.histos.append(self.draw_histo(h, 'TriggerCellVsPeakPos{0}'.format('Signal' if tprofile else ''), show, self.save_dir, lm=.11, draw_opt='colz', rm=.15, logz=True))
        h1.Draw('hist same')
        self.save_plots('TriggerCellVsPeakPos{0}{1}{2}'.format('Signal' if tprofile else '', 'BothCorr' if t_corr else '', 'Corr' if corr else ''), self.save_dir)
        self.RootObjects.append(h1)

    def draw_trigger_cell_vs_forc(self, show=True, cut=None, full_range=False, corr=False):
        if not full_range:
            self.tree.Draw('forc_pos', 'forc_pos[0]>20', 'goff')
            htemp = gROOT.FindObject('htemp')
            x = [int(htemp.GetBinCenter(htemp.FindFirstBinAbove(5000))) - 10, int(htemp.GetBinCenter(htemp.FindLastBinAbove(5000))) + 10]
        else:
            x = [0, 1024]
        h = TH2D('tcf', 'Trigger Cell vs. FORC Timing', 1024, 0, 1024, x[1] - x[0], x[0] / 2., x[1] / 2.)
        cut = self.AllCuts if cut is None else cut
        gStyle.SetPalette(55)
        forc = 'forc_pos/2.' if not corr else 'forc_time'
        self.tree.Draw('{forc}:trigger_cell>>tcf'.format(forc=forc), cut, 'goff')
        self.format_histo(h, x_tit='trigger cell', y_tit='forc timing [ns]', y_off=1.4)
        h.SetStats(0)
        self.histos.append(self.save_histo(h, 'TriggerCellVsFORC{0}'.format('FullRange' if full_range else ''), show, self.save_dir, lm=.11, draw_opt='colz', rm=.15))

    def draw_intlength_vs_triggercell(self, show=True, bin_size=2, prof=False):
        if prof:
            h = TProfile('hltc', 'Integral Length vs. Triggercell', 1024 / bin_size, 0, 1024)
        else:
            y_expect = (self.run.peak_integrals[self.PeakIntegral][0] + self.run.peak_integrals[self.PeakIntegral][1]) * .5
            h = TH2F('hltc', 'Integral Length vs. Triggercell', 1024 / bin_size, 0, 1024, 100, y_expect - 2, y_expect + 2)
        self.tree.Draw('IntegralLength[{num}]:trigger_cell>>hltc'.format(num=self.SignalNumber), self.Cut.all_cut, 'goff')
        self.format_histo(h, x_tit='Triggercell', y_tit='Integral Length [ns]', y_off=1.4, z_tit='Number of Entries', z_off=1.2)
        self.RootObjects.append(self.draw_histo(h, 'IntLengthVsTriggerCell', show, draw_opt='' if prof else 'colz', lm=.12, rm=.16 if not prof else .1))
        if not prof:
            gStyle.SetOptFit(1)
            gStyle.SetOptStat(0)
            gStyle.SetPalette(53)
            set_statbox(.82, .88, .15, 5)
            h_y = h.ProjectionY()
            fit = h_y.Fit('gaus', 'qs0')
            h.GetYaxis().SetRangeUser(fit.Parameter(1) - 5 * fit.Parameter(2), fit.Parameter(1) + 5 * fit.Parameter(2))
            f = TF1('f', '[0]*sin([1]*x - [2]) + [3]')
            f.SetLineColor(600)
            for i, name in enumerate(['y_sc', 'x_sc', 'x_off', 'y_off']):
                f.SetParName(i, name)
            f.SetParLimits(0, .1, 3)
            f.SetParLimits(1, 1e-4, 1e-2)
            h.Fit(f, 'q')
        self.save_plots('IntLengthVsTriggerCell', self.save_dir)
        gStyle.SetPalette(1)

    def draw_intdiff_vs_triggercell(self, show=True):
        h = TH2F('hdtc', 'Difference of the Integral Definitions vs Triggercell', 1024 / 2, 0, 1024, 200, 0, 25)
        hprof = TProfile('hdtc_p', 'Difference of the Integral Definitions vs Triggercell', 1024 / 8, 0, 1024)
        self.tree.Draw('(TimeIntegralValues[{num}]-IntegralValues[{num}]):trigger_cell>>hdtc'.format(num=self.SignalNumber), self.Cut.all_cut, 'goff')
        self.tree.Draw('(TimeIntegralValues[{num}]-IntegralValues[{num}]):trigger_cell>>hdtc_p'.format(num=self.SignalNumber), self.Cut.all_cut, 'goff')
        gStyle.SetPalette(53)
        self.format_histo(h, x_tit='Triggercell', y_tit='Integral2 - Integral1 [au]', z_tit='Number of Entries', stats=0, y_off=1.4, z_off=1.1)
        self.RootObjects.append(self.draw_histo(h, '', show, draw_opt='colz', lm=.12, rm=.15))
        self.format_histo(hprof, lw=3, color=600)
        hprof.Draw('hist same')
        p = h.ProjectionY()
        h.GetYaxis().SetRangeUser(0, p.GetBinCenter(p.FindLastBinAbove(p.GetMaximum() / 15.)))
        self.RootObjects.append(hprof)
        self.save_plots('IntDiffVsTriggerCell', self.save_dir)
        gStyle.SetPalette(1)

    # endregion

    # ==========================================================================
    # region SIGNAL/PEDESTAL
    def print_off_results(self, prnt=True):
        ph, ped, pul = self.draw_pulse_height(show=False)[1], self.Pedestal.draw_disto_fit(save=False), self.Pulser.draw_distribution_fit(show=False)
        string = '{0:3.2f} {1:3.2f} {2:3.2f}'.format(ph.Parameter(0), ped.Parameter(1), pul.Parameter(1))
        if prnt:
            print 'Signal\tPedest.\tPulser'
            print string
        else:
            return string

    def generate_signal_name(self, signal=None, evnt_corr=True, off_corr=False, bin_corr=False, cut=None):
        sig_name = signal if signal is not None else self.SignalName
        # pedestal polarity is always the same as signal polarity
        ped_pol = '1'
        # change polarity if pulser has opposite polarity to signal
        if hasattr(self, 'Pulser') and signal == self.Pulser.SignalName:
            ped_pol = '-1' if self.PulserPolarity != self.Polarity else ped_pol
        if bin_corr:
            return sig_name
        elif off_corr:
            ped_fit = self.Pedestal.draw_disto_fit(cut=cut, save=False)
            sig_name += '-{pol}*{ped}'.format(ped=ped_fit.Parameter(1), pol=ped_pol)
        elif evnt_corr:
            sig_name += '-{pol}*{ped}'.format(ped=self.PedestalName, pol=ped_pol)
        return sig_name

    def make_signal_time_histos(self, signal_name=None, evnt_corr=False, off_corr=False, bin_corr=False, rel_t=False, show=True):
        signal_name = self.generate_signal_name(self.SignalName if signal_name is None else signal_name, evnt_corr, off_corr, bin_corr)
        h = TH2D('h_st', 'Signal vs. Time', *(self.get_time_bins() + [225, -50, 500]))
        set_statbox(only_entries=True, x=.83)
        gStyle.SetPalette(53)
        self.tree.Draw('{name}:time/1000>>h_st'.format(name=signal_name), self.Cut.all_cut, 'goff')
        self.format_histo(h, x_tit='Time [min]', y_tit='Pulse Height [au]', y_off=1.4, t_ax_off=self.run.StartTime if rel_t else 0)
        self.save_histo(h, 'SignalTime', show, lm=.12, draw_opt='colz', rm=.15)
        return h

    def draw_pulse_height(self, binning=10000, redo=False, corr=True, sig=None, rel_t=True, show=True, save=True, prnt=True):

        sig = self.SignalName if sig is None else sig
        bin_size = binning if binning is not None else self.BinSize
        correction = '' if not corr else '_eventwise'
        suffix = '{bins}{cor}_{reg}'.format(bins=bin_size, cor=correction, reg=self.get_all_signal_names()[sig])
        picklepath = self.make_pickle_path('Ph_fit', None, self.RunNumber, self.channel, suf=suffix)

        def func():
            signal = self.generate_signal_name(self.SignalName if sig is None else sig, corr)
            prof = TProfile('pph', 'Pulse Height Evolution', *self.get_time_bins(binning))
            self.tree.Draw('{sig}:time/1000.>>pph'.format(sig=signal), self.Cut.all_cut, 'goff')
            self.PulseHeight = prof
            return prof

        p = func() if redo else None
        p = do_pickle(picklepath, func, p)
        set_statbox(entries=2, only_fit=True, w=.3)
        y_vals = [p.GetBinContent(i) for i in xrange(2, p.GetNbinsX() + 1)]
        self.format_histo(p, name='Fit Result', x_tit='Time [min]', y_tit='Mean Pulse Height [au]', y_off=1.6, x_range=[self.run.StartTime, self.get_time_bins()[1][-1]],
                          t_ax_off=self.run.StartTime if rel_t else 0, y_range=increased_range([min(y_vals), max(y_vals)], .5, .5), ndivx=505)
        self.draw_histo(p, show=show, lm=.14, prnt=save)
        fit = self.fit_pulse_height(p, picklepath)
        self.save_plots('PulseHeight{0}'.format(self.BinSize), show=show, save=save, prnt=prnt)
        return p, fit

    def fit_pulse_height(self, p, picklepath):
        fit = p.Fit('pol0', 'qs', '', 0, self.__get_max_fit_pos(p))
        kinder_pickle(picklepath, fit)
        return FitRes(fit)

    @staticmethod
    def __get_max_fit_pos(h):
        """ look for huge fluctiations in ph graph and return last stable point"""
        sum_ph = h.GetBinContent(1)
        for i in xrange(2, h.GetNbinsX() + 1):
            sum_ph += h.GetBinContent(i)
            if h.GetBinContent(i) < .7 * sum_ph / (i + 1):
                log_warning('Found huge ph fluctiation! Stopping Fit! y value = {y}, mean_y = {m}'.format(y=h.GetBinContent(i), m=sum_ph / (i + 1)))
                return h.GetBinCenter(i - 1)
        return h.GetBinCenter(h.GetNbinsX()) + 1000

    def draw_ph_distribution(self, binning=None, show=True, fit=True, xmin=0, xmax=270., bin_size=.5, save=True):
        if binning is not None:
            self.set_bin_size(binning)
        sig_time = self.make_signal_time_histos(evnt_corr=True, show=False)
        if not show:
            gROOT.SetBatch(1)
        means = [h_proj.GetMean() for h_proj in [sig_time.ProjectionY(str(i), i + 1, i + 1) for i in xrange(self.n_bins - 1)] if h_proj.GetEntries() > 10]
        nbins = int((xmax - xmin) / bin_size)
        h = TH1F('h', 'Signal Bin{0} Distribution'.format(self.BinSize), nbins, xmin, xmax)  # int(log(len(means), 2) * 2), extrema[0], extrema[1] + 2)
        for mean_ in means:
            h.Fill(mean_)
        self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.5, fill_color=407)
        h.Fit('gaus', 'q') if fit else do_nothing()
        if save:
            self.save_histo(h, 'SignalBin{0}Disto'.format(self.BinSize), lm=.12)
        return h

    def show_ph_overview(self, binning=None):
        self.draw_pulse_height(binning=binning, show=False)
        h1 = self.draw_pulse_height(show=False)
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
        self.histos.append([h2, c])

    def draw_signal_distribution(self, cut=None, evnt_corr=True, off_corr=False, show=True, sig=None, binning=1, events=None,
                                 start=None, x_range=None, redo=False, stats=True, prnt=True):
        cut = self.AllCuts if cut is None else TCut(cut)
        suffix = '{b}_{c}_{cut}'.format(b=binning, c=int(evnt_corr), cut=cut.GetName())
        pickle_path = self.make_pickle_path('PulseHeight', 'Histo', run=self.RunNumber, ch=self.DiamondNumber, suf=suffix)
        x_range = [-50, 300] if x_range is None else x_range

        def func():
            self.log_info('Drawing signal distribution for run {run} and {dia}...'.format(run=self.RunNumber, dia=self.DiamondName), prnt=prnt)
            set_root_output(False)
            h1 = TH1F('h_sd', 'Pulse Height {s}'.format(s='with Pedestal Correction' if evnt_corr else ''), int((x_range[1] - x_range[0]) * binning), *x_range)
            sig_name = self.generate_signal_name(sig, evnt_corr, off_corr, False, cut)
            start_event = int(float(start)) if start is not None else 0
            n_events = self.find_n_events(n=events, cut=str(cut), start=start_event) if events is not None else self.run.n_entries
            self.tree.Draw('{name}>>h_sd'.format(name=sig_name), str(cut), 'goff', n_events, start_event)
            self.format_histo(h1, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=2, fill_color=self.FillColor)
            return h1

        set_statbox(only_entries=True) if stats else do_nothing()
        h = func() if redo else None
        h = do_pickle(pickle_path, func, h)
        self.save_histo(h, 'SignalDistribution', lm=.15, show=show, prnt=prnt, save=cut)
        return h

    def draw_signal_vs_peakpos(self, show=True, corr=False):
        gr = self.make_tgrapherrors('gr', 'Signal vs Peak Position')
        i = 0
        x = self.run.signal_regions[self.SignalRegion]
        self.draw_peak_timing(show=False, corr=corr)
        h = self.PeakValues
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        for peak_pos in xrange(x[0] + 2, x[1] - 2):
            print '\rcalculating peak pos: {0:03d}'.format(peak_pos),
            self.Cut.set_signal_peak_pos(peak_pos, peak_pos + 1) if not corr else self.Cut.set_signal_peak_time(peak_pos / 2., (peak_pos + 1) / 2.)
            print peak_pos / 2., (peak_pos + 1) / 2.
            events = int(h.GetBinContent(h.FindBin(peak_pos / 2.)))
            print '({0:05d})'.format(events),
            stdout.flush()
            if events > 500:
                ph_fit = self.draw_pulse_height(show=False)[1]
                gr.SetPoint(i, peak_pos / 2., ph_fit.Parameter(0))
                gr.SetPointError(i, 0, ph_fit.ParError(0))
                i += 1
        gr.GetXaxis().SetLimits(x[0] / 2., x[1] / 2.)
        self.format_histo(gr, x_tit='Signal Peak Position [ns]', y_tit='Pulse Height [au]', y_off=1.4)
        self.histos.append(self.save_histo(gr, 'SignalVsPeakPos', show, self.save_dir, lm=.11, draw_opt='alp'))
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')

    def draw_sig_vs_corr_peaktiming(self, show=True, prof=False):
        x = self.run.signal_regions[self.SignalRegion]
        h = TProfile('hspt', 'Signal vs. Corrected Peak Timing', (x[1] - x[0]), x[0] / 2, x[1] / 2)
        if not prof:
            h = TH2F('hspt', 'Signal vs. Corrected Peak Timing', (x[1] - x[0]), x[0] / 2, x[1] / 2, 350, -50, 300)
        dic = self.Cut.calc_timing_range(show=False)
        t_correction = '({p1}* trigger_cell + {p2} * trigger_cell*trigger_cell)'.format(p1=dic['t_corr'].GetParameter(1), p2=dic['t_corr'].GetParameter(2))
        draw_string = '{sig}:IntegralPeakTime[{num}]-{tc}>>hspt'.format(sig=self.SignalName, num=self.SignalNumber, tc=t_correction)
        exluded_cuts = ['timing', 'bucket', 'tracks', 'chi2X', 'chi2Y', 'track_angle']
        cut = self.Cut.generate_special_cut(excluded=exluded_cuts)
        self.tree.Draw(draw_string, cut, 'goff')
        self.format_histo(h, fill_color=1)
        self.RootObjects.append(self.draw_histo(h, show=show, draw_opt='colz'))
        self.__draw_timing_cut()

    def draw_landau_vs_peakpos(self, show=True, bins=2):
        hs = THStack('lpp', 'Landau vs. Signal Peak Postion;pulse height;entries')
        x = self.run.signal_regions[self.SignalRegion]
        self.Cut.reset_cut('signal_peak_pos')
        self.draw_peak_timing(show=False)
        h_pv = self.PeakValues
        l = TLegend(.7, .38, .90, .88)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        for peak_pos in xrange(x[0] + 2, x[1] - 2, bins):
            print '\rcalculating peak pos: {0:03d}'.format(peak_pos),
            self.Cut.set_signal_peak_pos(peak_pos, peak_pos + bins)
            events = 0
            for pp in xrange(peak_pos, peak_pos + bins):
                events += int(h_pv.GetBinContent(h_pv.FindBin(peak_pos / 2.)))
            print '({0:05d})'.format(events),
            stdout.flush()
            if events > 10000:
                h = self.draw_signal_distribution(show=False, binning=100)
                h.SetLineColor(self.get_color())
                h.Scale(1 / h.GetMaximum())
                l.AddEntry(h, '[{0},{1}] ns'.format(int(peak_pos / 2.), int(peak_pos / 2. + bins / 2.)), 'l')
                hs.Add(h)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.reset_colors()
        self.format_histo(hs, y_tit='Pulse Height [au]', y_off=1.2)
        self.histos.append(self.save_histo(hs, 'LandauVsPeakPos', show, self.save_dir, lm=.11, draw_opt='nostack', l=l))

    def draw_signal_vs_triggercell(self, show=True, bins=10):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gr = self.make_tgrapherrors('stc', 'Signal vs Trigger Cell')
        i = 0
        for tcell in xrange(0, 1024 - bins, bins):
            if tcell:
                print '\033[F',
            print '\rcalculating pulse height for trigger cell: {0:03d}'.format(tcell),
            self.Cut.set_trigger_cell(tcell, tcell + bins)
            stdout.flush()
            ph_fit = self.draw_pulse_height(show=False)[1]
            gr.SetPoint(i, tcell, ph_fit.Parameter(0))
            gr.SetPointError(i, 0, ph_fit.ParError(0))
            i += 1
        print
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.format_histo(gr, x_tit='trigger cell', y_tit='pulse height [au]', y_off=1.2)
        self.histos.append(self.save_histo(gr, 'SignalVsTriggerCell', show, self.save_dir, lm=.11, draw_opt='alp'))

    # endregion

    # ==========================================================================
    # region CUTS

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
        self.histos.append([h, h_sig, h_ped1, h_ped2, c])

    def show_bucket_numbers(self, show=True):
        pickle_path = self.PickleDir + 'Cuts/BucketEvents_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.RunNumber, dia=self.DiamondName)

        def func():
            print 'getting number of bucket events for run {run} and {dia}...'.format(run=self.RunNumber, dia=self.DiamondName)
            n_new = self.tree.Draw('1', '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['bucket'], pul=self.Cut.CutStrings['pulser']), 'goff')
            n_old = self.tree.Draw('1', '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['old_bucket'], pul=self.Cut.CutStrings['pulser']), 'goff')
            if show:
                print 'New Bucket: {0} / {1} = {2:4.2f}%'.format(n_new, self.run.n_entries, n_new / float(self.run.n_entries) * 100)
                print 'Old Bucket: {0} / {1} = {2:4.2f}%'.format(n_old, self.run.n_entries, n_old / float(self.run.n_entries) * 100)
            return {'old': n_old, 'new': n_new, 'all': float(self.run.n_entries)}

        return do_pickle(pickle_path, func)

    def show_bucket_hits(self, show=True):
        # hit position
        h = TH2F('h', 'Diamond Margins', 80, -.3, .3, 52, -.3, .3)
        nr = 1 if not self.channel else 2
        cut = '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['old_bucket'], pul=self.Cut.CutStrings['pulser'])
        self.tree.Draw('dia_track_x[{nr}]:dia_track_y[{nr}]>>h'.format(nr=nr), cut, 'goff')
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
            self.histos.append([h, c])
        return h

    def draw_bucket_pedestal(self, show=True, corr=True, additional_cut='', draw_option='colz'):
        gStyle.SetPalette(55)
        cut_string = self.Cut.generate_special_cut(included=['tracks', 'pulser', 'saturated'])
        cut_string += additional_cut
        h = self.draw_signal_vs_peak_position('e', '2', show, corr, cut_string, draw_option, 1, save=False)
        self.format_histo(h, x_range=[self.run.signal_regions[self.SignalRegion][0] / 2, self.run.signal_regions['e'][1] / 2], stats=0)
        self.save_plots('BucketPedestal')

    def draw_bucket_waveforms(self, show=True, t_corr=True, start=100000):
        good = self.draw_waveforms(1, show=False, start_event=None, t_corr=t_corr)[0]
        cut = self.Cut.generate_special_cut(excluded=['bucket', 'timing']) + TCut('!({0})'.format(self.Cut.CutStrings['bucket']))
        bucket = self.draw_waveforms(1, cut=cut, show=False, start_event=start, t_corr=t_corr)[0]
        cut = self.Cut.generate_special_cut(excluded=['bucket', 'timing']) + TCut('{buc}&&!({old})'.format(buc=self.Cut.CutStrings['bucket'], old=self.Cut.CutStrings['old_bucket']))
        bad_bucket = self.draw_waveforms(1, cut=cut, show=False, t_corr=t_corr, start_event=None)[0]
        self.reset_colors()
        mg = TMultiGraph('mg_bw', 'Bucket Waveforms')
        l = self.make_legend(.85, .4, nentries=3)
        names = ['good wf', 'bucket wf', 'both wf']
        for i, gr in enumerate([good, bucket, bad_bucket]):
            self.format_histo(gr, color=self.get_color(), markersize=.5)
            mg.Add(gr, 'lp')
            l.AddEntry(gr, names[i], 'lp')
        self.format_histo(mg, draw_first=True, x_tit='Time [ns]', y_tit='Signal [mV]')
        x = [self.run.signal_regions['e'][0] / 2, self.run.signal_regions['e'][1] / 2 + 20]
        self.format_histo(mg, x_range=x, y_off=.7)
        y = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
        self.draw_histo(mg, show=show, draw_opt='A', x=1.5, y=0.75, lm=.07, rm=.045, bm=.2, l=l)
        self._add_buckets(y[0], y[1], x[0], x[1], avr_pos=-1, full_line=True)
        self.save_plots('BucketWaveforms')
        self.reset_colors()

    def show_bucket_means(self, show=True, plot_histos=True):
        pickle_path = self.PickleDir + 'Cuts/BucketMeans_{tc}_{run}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.RunNumber, dia=self.DiamondName)

        def func():
            gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
            cuts_nobucket = TCut('no_bucket', '')
            cuts_oldbucket = TCut('old_bucket', '')
            for key, value in self.Cut.CutStrings.iteritems():
                if not key.startswith('old') and key not in ['all_cuts', 'bucket']:
                    cuts_nobucket += value
                if key not in ['all_cuts', 'bucket']:
                    cuts_oldbucket += value
            h1 = self.draw_signal_distribution(show=False, evnt_corr=True)
            h2 = self.draw_signal_distribution(show=False, evnt_corr=True, cut=cuts_nobucket)
            h3 = self.draw_signal_distribution(show=False, evnt_corr=True, cut=cuts_oldbucket)
            if plot_histos:
                c = TCanvas('c', 'Bucket Histos', 1000, 1000)
                self.format_histo(h1, color=self.get_color(), lw=1, x_tit='Pulse Height [au]', y_tit='Entries')
                h1.Draw()
                self.format_histo(h2, color=self.get_color(), lw=1)
                h2.Draw('same')
                self.format_histo(h3, color=self.get_color(), lw=1)
                h3.Draw('same')
                self.histos.append([h1, h2, h3, c])
            result = {name: [h.GetMean(), h.GetMeanError()] for name, h in zip(['new', 'no', 'old'], [h1, h2, h3])}
            gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
            if show:
                print result
            return result

        res = func() if plot_histos else None
        return do_pickle(pickle_path, func, res)

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
                self.save_plots(save_name, canvas=c1, sub_dir=self.save_dir)
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
        self.save_plots('all', canvas=c2, sub_dir=self.save_dir)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        gROOT.SetBatch(0)

    def compare_normalised_cuts(self, scale=False, show=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        self.reset_colors()
        c1 = TCanvas('single', '', 1000, 1000)
        name = 'sCutComparison'
        if scale:
            name += "_scaled"
        else:
            name += "_noarmalized"
        if scale:
            title = 'Scaled Signal Distribution with Single Cuts'
        else:
            title = 'Normalised Signal Distribution with Single Cuts'
        title += ';Pulse Height [au];Normalised Entries'

        stack = THStack(name, title)

        entries = 0
        for value in self.Cut.CutStrings.itervalues():
            if str(value):
                entries += 1
        legend = self.make_legend(x1=.57, nentries=entries - 2)
        histos = []
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                save_name = 'signal_distribution_normalised_{cut}'.format(cut=key)
                histo_name = 'signal {range}{peakint}'.format(range=self.SignalRegion, peakint=self.PeakIntegral)
                histo_title = 'normalized' if not scale else 'scaled'
                histo_title += ' signal with cut ' + key
                histo = TH1F(histo_name, histo_title, 350, -50, 300)
                # safe single plots
                c1.cd()
                self.tree.Draw("{name}>>{histo}".format(name=self.SignalName, histo=histo_name), value)
                if scale:
                    histo = self.scale_histo(histo)
                else:
                    histo = self.normalise_histo(histo)
                histo.Draw()
                c1.Update()
                self.save_plots(save_name, canvas=c1, sub_dir=self.save_dir)
                # draw all single plots into c2
                histo.SetLineColor(self.get_color())

                if key == 'all_cuts':
                    histo.SetLineWidth(2)
                stack.Add(histo)
                histos.append(histo)
                legend.AddEntry(histo, key, 'l')
        stack.Draw()
        gROOT.SetBatch(0)

        for h in histos:
            h.SetStats(False)
        name = '{0}Cuts'.format('Normalised' if not scale else 'Scaled')
        self.format_histo(stack, y_off=1.4, x_off=1.1)
        self.RootObjects.append(self.save_histo(stack, name, show, self.save_dir, lm=.15, l=legend, draw_opt='nostack'))
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        gROOT.SetBatch(0)

    def compare_consecutive_cuts(self, scale=False, show=True, save_single=True, short=False, x_range=None):
        x_range = [-50, 500] if x_range is None else x_range
        self.reset_colors()
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        short_cuts = ['raw', 'saturated', 'timing', 'pulser', 'tracks', 'bucket', 'fiducial']
        legend = self.make_legend(.75 if short else .71, .88, nentries=len(self.Cut.ConsecutiveCuts) + 1 if not short else len(short_cuts) + 1)
        cut = TCut('consecutive', '')
        stack = THStack('scc', 'Signal Distribution with Consecutive Cuts')
        i = 0
        leg_style = 'l' if scale else 'f'
        for key, value in self.Cut.ConsecutiveCuts.iteritems():
            if short:
                if key not in short_cuts:
                    continue
            self.log_info('adding cut {0}'.format(key))
            key = 'beam_stops' if key.startswith('beam') else key
            cut += value
            save_name = 'signal_distribution_{n}cuts'.format(n=i)
            h = TH1F('h_{0}'.format(i), 'signal with {n} cuts'.format(n=i), 550, *x_range)
            self.tree.Draw('{name}>>h_{i}'.format(name=self.SignalName, i=i), cut, 'goff')
            if scale:
                self.scale_histo(h)
            self.save_histo(h, save_name, show=False, save=save_single)
            color = self.get_color()
            self.format_histo(h, color=color, stats=0)
            if not scale:
                h.SetFillColor(color)
            stack.Add(h)
            leg_entry = '+ {0}'.format(key) if i else key
            legend.AddEntry(h, leg_entry, leg_style)
            i += 1
        if short:
            h = self.draw_signal_distribution(show=False, binning=550, x_range=x_range)
            color = self.get_color()
            self.format_histo(h, color=color, stats=0)
            h.SetFillColor(color) if not scale else do_nothing()
            stack.Add(h)
            legend.AddEntry(h, '+ other', leg_style)
        self.format_histo(stack, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.9, draw_first=True)
        save_name = 'Consecutive{1}{0}'.format('Scaled' if scale else '', 'Short' if short else '')
        self.save_histo(stack, save_name, show, self.save_dir, l=legend, draw_opt='nostack', lm=0.14)
        stack.SetName(stack.GetName() + 'logy')
        # stack.SetMaximum(stack.GetMaximum() * 1.2)
        self.save_histo(stack, '{name}LogY'.format(name=save_name), show, self.save_dir, logy=True, draw_opt='nostack', lm=0.14)
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")

    def draw_fiducial_cut(self):
        self.Cut.draw_fid_cut()

    def draw_cut_means(self, show=True, short=False):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gr = self.make_tgrapherrors('gr_cm', 'Mean of Pulse Height for Consecutive Cuts')
        cut = TCut('consecutive', '')
        names = []
        i = 1
        gr.SetPoint(0, 0, 0)
        for key, value in self.Cut.CutStrings.iteritems():
            if (str(value) or key == 'raw') and key not in ['all_cuts', 'old_bucket']:
                if short:
                    self.log_info('adding cut {0}'.format(key))
                    if key not in ['raw', 'saturated', 'timing', 'bucket', 'pulser', 'tracks', 'fiducial']:
                        continue
                key = 'beam_stops' if key.startswith('beam') else key
                cut += value
                h = self.draw_signal_distribution(cut=cut, show=False)
                self.log_info('{0}, {1}, {2}'.format(key, h.GetMean(), h.GetMeanError()))
                gr.SetPoint(i, i, h.GetMean())
                gr.SetPointError(i, 0, h.GetMeanError())
                names.append(key)
                i += 1
        if short:
            h = self.draw_signal_distribution(show=False)
            gr.SetPoint(i, i, h.GetMean())
            gr.SetPointError(i, 0, h.GetMeanError())
            names.append('other')
        self.format_histo(gr, markersize=.2, fill_color=17, y_tit='Mean Pulse Height [au]', y_off=1.4)
        y = [gr.GetY()[i] for i in xrange(1, gr.GetN())]
        gr.GetYaxis().SetRangeUser(min(y) - 1, max(y) + 1)
        gr.GetXaxis().SetLabelSize(.05)
        for i in xrange(1, gr.GetN()):
            bin_x = gr.GetXaxis().FindBin(i)
            gr.GetXaxis().SetBinLabel(bin_x, names[i - 1])
        self.RootObjects.append(self.save_histo(gr, 'CutMeans{s}'.format(s='Short' if short else ''), show, self.save_dir, bm=.20, draw_opt='bap', lm=.12, x_fac=1.5))
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')

    def draw_distance_vs_ph(self, show=True, steps=10):
        h = self.draw_track_length(show=False, save=False)
        xmin, xmax = [h.GetBinCenter(i) for i in [h.FindFirstBinAbove(5), h.FindLastBinAbove(5)]]
        xvals = [xmin + i * (xmax - xmin) / steps for i in xrange(steps + 1)]
        gr = self.make_tgrapherrors('gr_aph', 'Pulse Height Vs Distance in Diamond')
        j = 0
        for i in xrange(len(xvals) - 1):
            cut = self.Cut.generate_distance(xvals[i], xvals[i + 1])
            self.Cut.all_cut += cut
            fit = self.draw_pulse_height(show=False)[1]
            if fit.Parameter(0):
                gr.SetPoint(j, xvals[i], fit.Parameter(0))
                gr.SetPointError(j, 0, fit.ParError(0))
                j += 1
            self.Cut.update_all_cut()
        self.draw_histo(gr, show)

    def test_landau_stats(self):
        gr = self.make_tgrapherrors('gr_ls', 'Landau Statistics')
        set_root_output(False)
        self.start_pbar(sum(int(pow(2, i / 2.)) for i in xrange(1, 40)))
        k = 0
        for j, i in enumerate(xrange(1, 40)):
            h = TH1F('h', 'h', 500, 0, 1000)
            for _ in xrange(int(pow(2, i / 2.))):
                k += 1
                h.Fill(gRandom.Landau(80, 5))
            self.ProgressBar.update(k)
            gr.SetPoint(j, pow(2, i), h.GetMean())
            gr.SetPointError(j, 0, h.GetMeanError())
        self.ProgressBar.finish()
        self.draw_histo(gr, draw_opt='alp', logx=True)

    def find_conv(self):
        gr = self.make_tgrapherrors('gr_c', 'chi2 vs nconv')
        for i, j in enumerate(xrange(10, 70, 5)):
            print j
            f = self.fit_langau(j, False)
            gr.SetPoint(i, j, f.Chi2 / f.NDF)
        self.draw_histo(gr)

    # endregion

    # ==========================================================================
    # region SHOW
    def draw_signal_vs_peak_position(self, region=None, peak_int=None, show=True, corr=True, cut=None, draw_opt='colz', nbins=4, save=True):
        region = self.SignalRegion if region is None else region
        peak_int = self.PeakIntegral if peak_int is None else peak_int
        cut = self.Cut.generate_special_cut(excluded=[self.Cut.CutStrings['timing']]) if cut is None else cut
        num = self.get_signal_number(region, peak_int)
        reg_margins = self.run.signal_regions[region]
        x_bins = (reg_margins[1] - reg_margins[0]) * nbins
        h = TH2F('h_spp', 'Signal Vs Peak Positions', x_bins, reg_margins[0] / 2., reg_margins[1] / 2., 550, -50, 500)
        peak_string = 'IntegralPeaks' if not corr else 'IntegralPeakTime'
        draw_string = '{sig}:{peaks}[{num}]{scale}>>h_spp'.format(sig=self.SignalName, num=num, peaks=peak_string, scale='/2.' if not corr else '')
        self.tree.Draw(draw_string, cut, 'goff')
        self.format_histo(h, x_tit='Peak Timing [ns]', y_tit='Pulse Height [au]', y_off=1.35, z_off=1.2, stats=0, z_tit='Number of Entries')
        self.save_histo(h, 'SignalVsPeakPos', show, draw_opt=draw_opt, logz=True, rm=.15, lm=.12, save=save)
        return h

    def draw_signal_vs_signale(self, show=True):
        gStyle.SetPalette(53)
        cut = self.Cut.generate_special_cut(excluded=['bucket'])
        num = self.get_signal_number(region='e')
        cut += TCut('IntegralPeakTime[{0}]<94&&IntegralPeakTime[{0}]>84'.format(num))
        h = TH2F('hsse', 'Signal b vs Signal e', 62, -50, 200, 50, 0, 200)
        self.tree.Draw('{sige}:{sigb}>>hsse'.format(sigb=self.SignalName, sige=self.get_signal_name(region='e')), cut, 'goff')
        self.format_histo(h, x_tit='Signal s_b [au]', y_tit='Signal s_e [au]', z_tit='Number of Entries', z_off=1.1, y_off=1.5, stats=0)
        self.RootObjects.append(self.save_histo(h, 'SignalEvsSignalB', show, rm=.15, lm=.13, draw_opt='colz'))
        gStyle.SetPalette(1)

    def draw_waveforms(self, n=1, cut=None, start_event=None, t_corr=True, channel=None, show=True):
        """ Draws a stack of n waveforms. """
        channel = self.channel if channel is None else channel
        if not self.run.wf_exists(self.channel):
            return
        start_event = self.count + self.StartEvent if start_event is None else start_event
        self.log_info('Drawing {n} waveform, startint at event: {s}'.format(n=n, s=start_event))
        cut = self.Cut.all_cut if cut is None else TCut(cut)
        n_events = self.find_n_events(n, cut, start_event)
        self.tree.SetEstimate(n * 1024)
        n_entries = self.tree.Draw('wf{ch}:trigger_cell'.format(ch=channel), cut, 'goff', n_events, start_event)
        title = '{n}{tc} Waveform{p}'.format(n=n, tc=' Time Corrected' if t_corr else '', p='s' if n > 1 else '')
        h = TH2F('h_wf', title, 1024, 0, 512, 2048, -512, 512)
        values = [self.tree.GetV1()[i] for i in xrange(n_entries)]
        if t_corr:
            times = [self.run.get_calibrated_times(self.tree.GetV2()[1024 * i]) for i in xrange(n_entries / 1024)]
            times = [v for lst in times for v in lst]
        else:
            times = [(.4 if self.run.Digitiser == 'caen' else .5) * i for i in xrange(1024)] * n
        for v, t in zip(values, times):
            h.Fill(t, v)
        self.tree.SetEstimate()
        y_range = increased_range([min(values), max(values)], .1, .2)
        h = self.make_tgrapherrors('g_cw', title, x=times, y=values) if n == 1 else h
        self.format_histo(h, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, y_range=y_range, markersize=.5)
        self.save_histo(h, 'WaveForms{n}'.format(n=n), show=show, draw_opt='col' if n > 1 else 'apl', lm=.073, rm=.045, bm=.18, x_fac=1.5, y_fac=.5)
        self.count += n_events
        return h, n_events

    def draw_single_waveform(self, cut='', event=None, show=True):
        h, n = self.draw_waveforms(n=1, start_event=event, cut=cut, t_corr=True, show=False)
        self.draw_histo(h, show=show, gridy=1, gridx=1, lm=.073, rm=.045, bm=.18, x=1.5, y=.5)

    def show_single_waveforms(self, n=1, cut='', start_event=None):
        start = self.StartEvent + self.count if start_event is None else start_event + self.count
        activated_wfs = [wf for wf in xrange(4) if self.run.wf_exists(wf)]
        print 'activated wafeforms:', activated_wfs
        print 'Start at event number:', start
        wfs = [self.draw_waveforms(n=n, start_event=start, cut=cut, show=False, channel=wf) for wf in activated_wfs]
        n_wfs = len(activated_wfs)
        if not gROOT.GetListOfCanvases()[-1].GetName() == 'c_wfs':
            c = TCanvas('c_wfs', 'Waveforms', 2000, n_wfs * 500)
            c.Divide(1, n_wfs)
        else:
            c = gROOT.GetListOfCanvases()[-1]
        for i, wf in enumerate(wfs, 1):
            wf[0].SetTitle('{nam} WaveForm'.format(nam=self.run.DRS4Channels[activated_wfs[i - 1]]))
            c.cd(i)
            wf[0].Draw('aclp')
        self.RootObjects.append([c, wfs])
        cnt = wfs[0][1]
        # if cnt is None:
        #     return
        self.count += cnt

    # endregion

    def find_n_events(self, n, cut, start):
        total_events = self.tree.Draw('event_number', cut, 'goff', self.run.n_entries, start)
        evt_numbers = [self.tree.GetV1()[i] for i in xrange(total_events)]
        return int(evt_numbers[:n][-1] + 1 - start)

    def check_alignment(self, binning=5000, show=True):
        """ just check the number of pixel hits at pulser events for no offset """
        pickle_path = 'Configuration/Individual_Configs/Alignment/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run.RunNumber)

        def func():
            nbins = self.run.n_entries / binning
            h = TProfile('h', 'Pulser Rate', nbins, 0, self.run.n_entries)
            self.tree.Draw('(@col.size()>1)*100:Entry$>>h', 'pulser', 'goff')
            self.format_histo(h, title='Event Alignment', x_tit='Event Number', y_tit='Hit Efficiency @ Pulser Events [%]', y_off=1.3, stats=0, y_range=[0, 105], fill_color=self.FillColor)
            self.save_histo(h, 'EventAlignment', show, self.TelSaveDir, draw_opt='hist', prnt=show, rm=.08)
            return all(h.GetBinContent(bin_) < 40 for bin_ in xrange(5, h.GetNbinsX()))

        aligned = func() if show else None
        aligned = do_pickle(pickle_path, func, aligned)
        if not aligned:
            log_warning('Run {r} is misaligned :-('.format(r=self.RunNumber))
        return aligned

    def find_event_offsets(self, binning=5000, show=True):
        nbins = self.run.n_entries / binning
        histos = [TProfile('h{i}'.format(i=i), 'Pulser Rate', nbins, 0, self.run.n_entries) for i in xrange(5)]
        self.tree.SetEstimate(self.run.n_entries)
        self.tree.Draw('(@col.size()>1)*100', '', 'goff')
        cols = [self.tree.GetV1()[i] for i in xrange(self.run.n_entries)]
        n = self.tree.Draw('Entry$', 'pulser', 'goff')
        pulser_events = [int(self.tree.GetV1()[i]) for i in xrange(n)]
        for ev in pulser_events[:-1]:
            histos[0].Fill(ev, cols[ev])
            histos[1].Fill(ev, cols[ev - 1])
            histos[2].Fill(ev, cols[ev + 1])
            histos[3].Fill(ev, cols[ev - 2])
            histos[4].Fill(ev, cols[ev + 2])
        for h in histos:
            self.format_histo(h, title='Event Alignment', x_tit='Event Number', y_tit='Hits per Event @ Pulser Events [%]', y_off=1.3, stats=0, color=self.get_color(),
                              y_range=[0, 105], fill_color=self.FillColor)
        self.save_histo(histos[0], 'EventAlignment', show, self.TelSaveDir, draw_opt='hist', prnt=show, rm=.08)
        for h in histos[1:]:
            h.Draw('same')
        self.RootObjects.append([histos])
        self.reset_colors()

    @staticmethod
    def normalise_histo(histo, to100=False):
        h = histo
        h.GetXaxis().SetRangeUser(0, 30)
        min_bin = h.GetMinimumBin()
        h.GetXaxis().UnZoom()
        max_bin = h.GetNbinsX() - 1
        integral = h.Integral(min_bin, max_bin)
        if integral:
            fac = 100 if to100 else 1
            h.Scale(fac / integral)
        return h

    @staticmethod
    def scale_histo(histo):
        h = histo
        h.GetXaxis().SetRangeUser(30, 500)
        maximum = h.GetBinContent(h.GetMaximumBin())
        h.GetXaxis().UnZoom()
        if maximum:
            h.Scale(1. / maximum)
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
                h = self.draw_signal_distribution(evnt_corr=True, cut=value, show=False)
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
        self.histos.append(legend)
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

    def draw_snrs(self, show=True, lego=True, proj=False, draw_opt='lego2'):
        self.verbose = False
        gStyle.SetPaintTextFormat('5.4g')
        lego = False if proj else lego
        gr = self.make_tgrapherrors('gr', 'Signal to Noise Ratios')
        h = TProfile2D('h_snr', 'Signal to Noise Ratios', 70, 0, 70, 70, 0, 140)
        i = 0
        for name, region in self.get_all_signal_names().iteritems():
            if self.SignalRegion.split('_')[-1] in region:
                peak_integral = self.get_peak_integral(remove_letters(region))
                snr = self.calc_snr(name=name, reg=self.get_all_signal_names()[name])
                h.Fill(peak_integral[0] / 2., peak_integral[1] / 2., snr.n)
                gr.SetPoint(i, i + 1, snr.n)
                gr.SetPointError(i, 0, snr.s)
                gr.GetListOfFunctions().Add(self.draw_tlatex(i + 1, snr.n + snr.s * 1.5, str(peak_integral), align=22, size=.02))
                i += 1
        self.format_histo(gr, y_tit='SNR', y_off=1.2, color=self.get_color(), fill_color=1)
        vals = sorted([h.GetBinContent(i) for i in xrange(h.GetNbinsX() * h.GetNbinsY()) if h.GetBinContent(i)])
        x, y, z1 = Long(0), Long(0), Long(0)
        xmin, ymin = h.GetXaxis().GetXmin(), h.GetYaxis().GetXmin()
        h.GetBinXYZ(h.GetMaximumBin(), x, y, z1)
        x1, y1 = (x - 1) / 2. + xmin, (y - 1) / 2. + ymin
        self.__draw_profiles(h, x, y, proj)
        self.format_histo(h, x_tit='Left Length [ns]', x_off=1.45, y_tit='Right Length [ns]', y_off=1.6, z_tit='snr', z_off=1.6, stats=0, z_range=[vals[2], max(vals)])
        h.SetContour(50)
        gStyle.SetPalette(53)
        if draw_opt == 'coltext':
            self.show_best_snr(h, x1, y1, show)
        else:
            self.save_histo(h, 'SNRLego', show and lego, draw_opt=draw_opt, bm=.2, rm=.1, lm=.13, phi=-30, theta=40)
        gStyle.SetPalette(1)
        self.save_histo(gr, 'SNR', not (lego or proj) and show, draw_opt='bap')

    def show_best_snr(self, histo, x, y, show):
        h = histo
        self.format_histo(h, x_off=1, y_off=1.15, stats=0, z_tit='snr [au]', z_off=1.35)
        self.draw_histo(h, '', show, draw_opt='colztext', rm=.16)
        self.draw_vertical_line(x, -1e5, 1e5, color=2, style=2, name='a', w=2)
        self.draw_vertical_line(x + .5, -1e5, 1e5, color=2, style=2, name='b', w=2)
        self.draw_horizontal_line(y, 0, 10, color=418, style=2, w=2, name='c')
        self.draw_horizontal_line(y + .5, 0, 100, color=418, style=2, w=2, name='d')
        self.save_plots('SNRColText')

    def __draw_profiles(self, histo, x, y, show=True):
        h = histo
        py = h.ProfileY('Right Length', x, x)
        px = h.ProfileX('Left Length', y, y)
        vals = [py.GetBinContent(i) for i in xrange(py.GetNbinsX()) if py.GetBinContent(i)] + [px.GetBinContent(i) for i in xrange(px.GetNbinsX()) if px.GetBinContent(i)]
        self.format_histo(py, stats=0, lw=2)
        self.format_histo(px, stats=0, lw=2)
        py.SetLineColor(2)
        px.SetLineColor(418)
        l = self.make_legend(.68, .95)
        [l.AddEntry(p, p.GetName(), 'fp') for p in [py, px]]
        stack = THStack('s_sp', 'SNR Profiles')
        stack.Add(py, 'histe')
        stack.Add(px, 'histe')
        self.format_histo(stack, draw_first=True, x_tit='Integral Length [ns]', y_tit='snr [au]', y_off=1.35)
        stack.SetMinimum(increased_range([min(vals), max(vals)], .5, .5)[0])
        stack.SetMaximum(increased_range([min(vals), max(vals)], .5, .5)[1])
        self.save_histo(stack, 'SNRProfiles', show, draw_opt='nostack', l=l, lm=.13)

    def calc_snr(self, name=None, reg=''):
        signal_name = self.SignalName if name is None else name
        peak_int = remove_letters(self.get_all_signal_names()[signal_name])
        ped_sigma = make_ufloat(self.Pedestal.draw_disto_fit(save=False, name=self.Pedestal.get_signal_name(peak_int=peak_int), show=False), par=2)
        signal = make_ufloat(self.draw_pulse_height(corr=True, show=False, sig=signal_name)[1])
        snr = signal / ped_sigma
        print '{name} {0}\t| SNR is: {snr}\t {1} {2}'.format(self.get_peak_integral(peak_int), signal.n, ped_sigma.n, name=reg, snr=snr)
        return snr

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
        self.histos.append([gr, c])

    def signal_vs_peakintegral(self, show=True, ped=False):
        gROOT.SetBatch(1)
        gr = self.make_tgrapherrors('gr', '{sig} vs Peak Integral'.format(sig='Signal' if not ped else 'Pedestal'))
        peak_integrals = OrderedDict(sorted({key: value for key, value in self.run.peak_integrals.iteritems() if len(key) < 3}.items()))
        i = 0
        ratio = '{0}{1}'.format(self.run.peak_integrals.values()[0][0], self.run.peak_integrals.values()[0][1])
        for name, value in peak_integrals.iteritems():
            sig_name = self.get_signal_name(region='b', peak_integral=name)
            signal = self.draw_pulse_height(corr=True, show=False, sig=sig_name) if not ped else self.Pedestal.draw_disto(save=False, name=self.Pedestal.get_signal_name(peak_int=name))
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
        self.histos.append([gr, c])

    # endregion

    # ============================================
    # region MISCELLANEOUS
    def get_cut(self):
        """ :return: full cut_string """
        return self.Cut.all_cut

    def get_all_signal_names(self, sig_type='signal'):
        names = OrderedDict()
        for region in self.run.IntegralRegions[self.DiamondNumber - 1]:
            if sig_type in region:
                for integral in self.run.PeakIntegrals[self.DiamondNumber - 1]:
                    name = 'ch{ch}_{reg}_{int}'.format(ch=self.channel, reg=region, int=integral)
                    num = self.IntegralNames.index(name)
                    reg = region.replace(sig_type, '').strip('_') + integral.replace('PeakIntegral', '')
                    names[self.SignalDefinition.format(pol=self.Polarity, num=num)] = reg
        return names

    def print_info_header(self):
        header = ['Run', 'Type', 'Diamond', 'HV [V]', 'Region']
        for info in header:
            print self.adj_length(info),
        print

    def print_information(self, header=True):
        if header:
            self.print_info_header()
        infos = [self.RunNumber, self.run.RunInfo['type'], self.DiamondName.ljust(4), self.Bias, self.SignalRegion + self.PeakIntegral + '   ']
        for info in infos:
            print self.adj_length(info),
        print

    def show_integral_names(self):
        for i, name in enumerate(self.IntegralNames):
            if name.startswith('ch{}'.format(self.channel)):
                print str(i).zfill(3), name

    # endregion

    def spectrum(self, it=20, noise=20):
        decon = array(1024 * [0], 'f')
        s = TSpectrum(25)
        peaks = []
        for i in xrange(it):
            self.tree.GetEntry(300000 + i)
            data = array([-1 * self.tree.wf0[j] for j in xrange(1024)], 'f')
            thr = 100 * 2 * noise / max(data)
            print thr
            p = s.SearchHighRes(data, decon, 1024, 5, thr, True, 3, True, 5)
            xpos = [s.GetPositionX()[i] for i in xrange(p)]
            peaks.append(xpos)
        return decon, s, peaks

    def fixed_integrals(self):
        tcals = [0.4813, 0.5666, 0.3698, 0.6393, 0.3862, 0.5886, 0.5101, 0.5675, 0.4033, 0.6211, 0.4563, 0.5919, 0.4781, 0.5947, 0.417, 0.5269,
                 0.5022, 0.5984, 0.4463, 0.622, 0.4326, 0.5603, 0.3712, 0.6168, 0.5238, 0.5515, 0.514, 0.5949, 0.4198, 0.5711, 0.5344, 0.5856,
                 0.3917, 0.6125, 0.4335, 0.5817, 0.4658, 0.5338, 0.4442, 0.5865, 0.4482, 0.5778, 0.4755, 0.6118, 0.4113, 0.5609, 0.465, 0.6188,
                 0.3908, 0.5736, 0.5223, 0.5222, 0.5109, 0.493, 0.4421, 0.5908, 0.4555, 0.6737, 0.371, 0.5172, 0.5362, 0.5982, 0.5017, 0.4976,
                 0.5568, 0.5519, 0.416, 0.5788, 0.476, 0.5636, 0.4424, 0.5773, 0.4472, 0.6109, 0.4123, 0.616]
        sum_time = 0
        times = []
        for i in range(40):
            times.append(sum_time)
            sum_time += tcals[i]
        h = TH1F('h', 'Integral Length', len(times) - 1, array(times, 'f'))
        self.tree.GetEntry(200002)
        peak_pos = self.tree.IntegralPeaks[self.SignalNumber]
        wf = list(self.tree.wf0)
        mid = times[15] + tcals[15] / 2.
        for i in range(40):
            h.SetBinContent(i, abs((wf[peak_pos - 16 + i])))

        points_x1 = [mid - 4, mid - 4, h.GetBinLowEdge(9), h.GetBinLowEdge(9), mid - 4]
        points_x2 = [mid + 6, mid + 6, h.GetBinLowEdge(27), h.GetBinLowEdge(27), mid + 6]
        points_y1 = [0, -1 * wf[peak_pos - 8] - .3, -1 * wf[peak_pos - 8] - .3, 0, 0]
        points_y2 = [0, -1 * wf[peak_pos + 11] - .3, -1 * wf[peak_pos + 11] - .3, 0, 0]
        gr1 = TGraph(5, array(points_x1, 'd'), array(points_y1, 'd'))
        gr2 = TGraph(5, array(points_x2, 'd'), array(points_y2, 'd'))
        gr1.SetFillColor(kOrange + 7)
        gr2.SetFillColor(kOrange + 7)
        gr3 = TGraph(2, array([mid, mid], 'd'), array([0, -1 * wf[peak_pos]], 'd'))
        gr3.SetLineWidth(2)
        ar = TArrow(mid - 4, 50, mid + 6, 50, .015, '<|>')
        ar.SetLineWidth(2)
        ar.SetFillColor(1)
        ar.SetLineColor(1)

        c = TCanvas('c', 'c', 2500, 1500)
        self.format_histo(h, x_tit='Time [ns]', y_tit='Pulse Height [au]')
        h.SetStats(0)
        h1 = h.Clone()
        h1.GetXaxis().SetRangeUser(mid - 4 + .5, mid + 6 - .7)
        h1.SetFillColor(2)
        h.SetLineWidth(2)
        h.Draw()
        h1.Draw('same][')
        gr1.Draw('f')
        gr2.Draw('f')
        gr3.Draw('l')
        print mid - 4, mid + 6
        ar.Draw()
        self.histos.append([h, h1, gr1, gr2, gr3, ar, c])

    def draw_tcal(self, show=True):
        gr = self.make_tgrapherrors('gr_tcal', 'DRS4 Bin Sizes', marker_size=.5, x=range(len(self.run.TCal)), y=self.run.TCal)
        gr.Fit('pol0', 'qs') if show else do_nothing()
        self.format_histo(gr, x_tit='Bin number', y_tit='Length [ns]', y_off=1.5, y_range=[0, 1])
        set_statbox(only_fit=True)
        self.save_histo(gr, 'DRSBinSizes', x_fac=1.5, y_fac=.75, show=show)

    def draw_tcal_disto(self, show=True):
        h = TH1F('h_tcal', 'Bin Size Distribution', *self.Plots.get_tcal_bins())
        for value in self.run.TCal:
            h.Fill(value)
        self.format_histo(h, x_tit='time [ns]', y_tit='number of entries', y_off=1.5)
        self.save_histo(h, 'DRSBinSizeDisto', show, self.save_dir)

    def save_felix(self):
        self.save_dir = '{info}_{rp}'.format(info=self.make_info_string().strip('_'), rp=self.RunNumber)
        self.TelSaveDir = '{info}_{rp}'.format(info=self.make_info_string().strip('_'), rp=self.RunNumber)
        self.set_save_directory('PlotsFelix')
        self.Pulser.save_felix()

    def get_mpv_fwhm(self, show=True):
        h = self.draw_signal_distribution(show=show)
        return find_mpv_fwhm(h)

    def get_peak_integral(self, name):
        return  self.run.PeakIntegrals[self.DiamondNumber - 1]['PeakIntegral{}'.format(name) if not 'Peak' in str(name) else name]

    @staticmethod
    def make_region(signal, region=''):
        return '{}{}'.format(signal, '_' + region if region else '')

    def __placeholder(self):
        pass


if __name__ == "__main__":
    st = time()
    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=392, type=int)
    parser.add_argument('dia', nargs='?', default=1, type=int)
    parser.add_argument('-tc', '--testcampaign', nargs='?', default='')
    parser.add_argument('-v', '--verbose', action='store_false')
    parser.add_argument('-t', '--tree', action='store_false')
    args = parser.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    a = Elementary(tc)
    a.print_testcampaign()
    print_banner('STARTING PAD-ANALYSIS OF RUN {0}'.format(args.run))
    print
    run_class = Run(args.run, verbose=args.verbose, tree=None if args.tree else args.tree)
    z = PadAnalysis(run_class, args.dia)
    print_elapsed_time(st, 'Instantiation')
