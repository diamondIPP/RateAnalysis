#!/usr/bin/env python
# --------------------------------------------------------
#       Pedestal analysis of the pad waveforms
# created on May 23rd 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import TH2F, gStyle, TH1F
from Utils import set_statbox, set_drawing_range, FitRes
from copy import deepcopy
from numpy import array, mean
from collections import OrderedDict


class PedestalAnalysis(Elementary):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.Run = self.Ana.run
        self.Channel = self.Ana.channel
        self.Tree = self.Ana.tree
        self.Cut = self.Ana.Cut
        self.save_dir = self.Ana.save_dir
        self.Polarity = self.Ana.Polarity
        self.SignalName = self.get_signal_name()
        self.DiamondName = self.Ana.DiamondName
        self.DiamondNumber = self.Ana.DiamondNumber
        self.RunNumber = self.Ana.RunNumber

        self.ROOTObjects = []

    def get_signal_name(self, region=None, peak_int=None):
        return self.Ana.get_signal_name(region=region, peak_integral=peak_int, sig_type='pedestal')

    def get_all_signal_names(self):
        names = OrderedDict()
        regions = [reg for reg in self.Run.pedestal_regions if len(reg) < 3]
        integrals = [integral for integral in self.Run.peak_integrals if len(integral) < 3]
        for region in regions:
            for integral in integrals:
                name = 'ch{ch}_pedestal_{reg}_PeakIntegral{int}'.format(ch=self.Channel, reg=region, int=integral)
                num = self.Ana.IntegralNames[name]
                reg = region + integral
                names[self.Ana.SignalDefinition.format(pol=self.Polarity, num=num)] = reg
        return names

    def draw_disto(self, name=None, cut=None, logy=False, show=True, save=True, redo=False):
        show = False if not save else show
        cut = self.Cut.all_cut if cut is None else TCut(cut)
        picklepath = self.make_pickle_path('Pedestal', 'Disto', run=self.RunNumber, ch=self.DiamondNumber, suf=cut.GetName())

        def func():
            self.log_info('Drawing pedestal distribution for {d} of run {r}'.format(d=self.DiamondName, r=self.RunNumber))
            h1 = TH1F('h_pd', 'Pedestal Distribution', 400, -100, 100)
            self.Tree.Draw('{name}>>h_pd'.format(name=self.SignalName if name is None else name), cut, 'goff')
            self.format_histo(h1, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.8, fill_color=self.FillColor)
            return h1
        if show or save:
            set_statbox(entries=8, opt=1000000010)
        h = self.do_pickle(picklepath, func, redo=redo)
        set_drawing_range(h, rfac=.2)
        self.save_histo(h, 'PedestalDistribution', show, save=save, logy=logy, lm=.13)
        return h

    def draw_disto_fit(self, name=None, cut=None, logy=False, show=True, save=True, redo=False):
        show = False if not save else show
        set_statbox(only_fit=True, entries=4, w=.3)
        h = self.draw_disto(name, cut, logy, show=False, save=False, redo=redo)
        set_drawing_range(h)
        h.SetName('Fit Results')
        fit_pars = self.fit_fwhm(h, do_fwhm=True, draw=show)
        f = deepcopy(h.GetFunction('gaus'))
        f.SetNpx(1000)
        f.SetRange(h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        f.SetLineStyle(2)
        h.GetListOfFunctions().Add(f)
        self.save_histo(h, 'PedestalDistributionFit', show, save=save, logy=logy, lm=.13)
        return FitRes(fit_pars)

    def draw_signal_time(self, signal_name=None, show=True):
        signal_name = self.Ana.generate_signal_name(self.SignalName if signal_name is None else signal_name, evnt_corr=False)
        h = TH2F('h_st', 'Pedestal vs. Time', len(self.Ana.time_binning) - 1, array(self.Ana.get_minute_time_binning()), 160, -80, 80)
        set_statbox(only_entries=True, x=.83)
        gStyle.SetPalette(53)
        self.Tree.Draw('{name}:(time - {s}) / 60000>>h_st'.format(name=signal_name, s=self.Ana.run.startTime), self.Cut.all_cut, 'goff')
        self.format_histo(h, x_tit='Time [min]', y_tit='Pulse Height [au]', y_off=1.4)
        self.save_histo(h, 'PedestalTime', show, lm=.12, draw_opt='colz', rm=.15)
        return h

    def draw_pulse_height(self, signal_name=None, show=True, sigma=False, redo=False):
        signal_name = self.Ana.generate_signal_name(self.SignalName if signal_name is None else signal_name, evnt_corr=False)
        picklepath = self.make_pickle_path('Pedestal', 'Evolution', run=self.RunNumber, ch=self.DiamondNumber, suf=self.get_all_signal_names()[signal_name])

        def func():

            h = self.draw_signal_time(signal_name, False)
            g1 = self.make_tgrapherrors('g_pph', 'Pedestal{s} Pulse Height Evolution'.format(s=' Sigma' if sigma else ''))

            for ibin in xrange(1, h.GetNbinsX() - 1):
                proj_y = h.ProjectionY('_py{i}'.format(i=ibin), ibin, ibin)
                fit = proj_y.Fit('gaus', 'qs0')
                t_mean = mean(self.Ana.get_minute_time_binning()[(ibin - 1):ibin])
                g1.SetPoint(ibin - 1, t_mean, fit.Parameter(2 if sigma else 1))
                g1.SetPointError(ibin - 1, 0, fit.ParError(2 if sigma else 1))
            return g1

        g = self.do_pickle(picklepath, func, redo=redo)
        self.format_histo(g, x_tit='Time [min]', y_tit='Mean Pulse Height [au]', y_off=1.6)
        self.save_histo(g, 'Pedestal{s}PulseHeight'.format(s='Sigma' if sigma else ''), show=show, lm=.14, draw_opt='apl')
        return g

    def fit_pulse_height(self, signal_name=None, show=True, sigma=False):
        set_statbox(entries=4, only_fit=True)
        g = self.draw_pulse_height(signal_name, False, sigma)
        fit = g.Fit('pol0', 'qs')
        self.save_histo(g, '{n}Fit'.format(n=g.GetTitle()), show, lm=.12, draw_opt='apl')
        return FitRes(fit)

    def draw_sigma(self, signal_name=None):
        return self.draw_pulse_height(signal_name, sigma=True)
