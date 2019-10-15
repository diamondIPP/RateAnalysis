#!/usr/bin/env python
# --------------------------------------------------------
#       Pedestal analysis of the pad waveforms
# created on May 23rd 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from ROOT import TH2F, gStyle, TH1F, TCut, gROOT, gPad
from Utils import set_statbox, set_drawing_range, FitRes, server_pickle, do_pickle, make_ufloat
from copy import deepcopy
from numpy import mean, sqrt, log
from InfoLegend import InfoLegend


class PedestalAnalysis(Elementary):
    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.channel
        self.Tree = self.Ana.tree
        self.Cut = self.Ana.Cut
        self.save_dir = self.Ana.save_dir
        self.Polarity = self.Ana.Polarity
        self.SignalName = self.get_signal_name()
        self.DiamondName = self.Ana.DiamondName
        self.DiamondNumber = self.Ana.DiamondNumber
        self.RunNumber = self.Ana.RunNumber
        self.InfoLegend = InfoLegend(pad_analysis)

        self.Histogram = None

        self.ROOTObjects = []

    def get_signal_name(self, region=None, peak_int=None):
        return self.Ana.get_signal_name(region=region, peak_integral=peak_int, sig_type='pedestal')

    def get_all_signal_names(self):
        return self.Ana.get_all_signal_names('pedestal')

    def get_mean(self):
        return make_ufloat(self.draw_disto_fit(show=False, prnt=False), par=1)

    def get_noise(self):
        return make_ufloat(self.draw_disto_fit(show=False, prnt=False), par=2)

    def get_fwhm(self):
        return self.get_noise() * 2 * sqrt(2 * log(2))

    def draw_disto(self, name=None, cut=None, logy=False, show=True, save=True, redo=False, prnt=True, normalise=None):
        show = False if not save else show
        cut = self.Cut.all_cut if cut is None else TCut(cut)
        signal_name = self.SignalName if name is None else name
        picklepath = self.make_pickle_path('Pedestal', 'Disto', run=self.RunNumber, ch=self.DiamondNumber, suf='{c}_{r}'.format(c=cut.GetName(), r=self.get_all_signal_names()[signal_name]))

        def func():
            self.log_info('Drawing pedestal distribution for {d} of run {r}'.format(d=self.DiamondName, r=self.RunNumber), prnt=prnt)
            h1 = TH1F('h_pd', 'Pedestal Distribution', 2400, -150, 150)
            self.Tree.Draw('{name}>>h_pd'.format(name=signal_name), cut, 'goff')
            return h1

        if show:
            self.set_statbox(all_stat=True, n_entries=6, w=.3)
        h = do_pickle(picklepath, func, redo=redo)
        self.format_histo(h, 'Pedestal', x_tit='Pedestal [mV]', y_tit='Number of Entries', y_off=1.8, fill_color=self.FillColor, normalise=normalise)
        set_drawing_range(h, rfac=.2)
        self.save_histo(h, 'PedestalDistribution{}'.format(cut.GetName()), show, save=save, logy=logy, lm=.13, prnt=prnt)
        self.Histogram = h
        return h

    def draw_disto_fit(self, name=None, cut=None, logy=False, show=True, save=True, redo=False, prnt=True, draw_cut=False, normalise=None):
        cut = self.Cut.all_cut if cut is None else TCut(cut)
        cut = self.Cut.generate_special_cut(excluded='ped_sigma') if draw_cut else cut
        suffix = '{r}_fwhm_{c}'.format(c=cut.GetName(), r=self.get_all_signal_names()[self.SignalName if name is None else name])
        picklepath = self.make_pickle_path('Pedestal', run=self.RunNumber, ch=self.DiamondNumber, suf=suffix)
        show = False if not save else show
        self.set_statbox(fit=True, n_entries=6, w=.35)
        h = self.draw_disto(name, cut, logy, show=False, save=save, redo=redo, prnt=prnt, normalise=normalise)
        set_drawing_range(h, thresh=h.GetMaximum() * .02)
        fit_pars = self.fit_fwhm(h, do_fwhm=True, draw=True)
        f = deepcopy(h.GetFunction('gaus'))
        f.SetNpx(1000)
        f.SetRange(h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        f.SetLineStyle(2)
        h.GetListOfFunctions().Add(f)
        self.draw_histo(h, '', show, logy=logy, lm=.13, prnt=prnt)
        if draw_cut:
            b = self.__draw_cut(h)
            h.Draw('same')
            b.Draw('l')
            gPad.RedrawAxis()
        h.Sumw2(False)
        self.save_plots('PedestalDistributionFit{}'.format(cut.GetName()), save=save, prnt=prnt, show=show)
        self.Histogram = h
        server_pickle(picklepath, fit_pars)
        return fit_pars

    def __draw_cut(self, h):
        fit = h.GetListOfFunctions()[2]
        xmin, xmax = fit.GetParameter(1) - 3 * fit.GetParameter(2), fit.GetParameter(1) + 3 * fit.GetParameter(2)
        b = self.draw_box(xmin, -10, xmax, 1e7, color=2, width=2, fillstyle=3001, name='ped', style=7)
        legend = self.make_legend(.65, y2=.64, nentries=1, margin=.45, name='la', scale=1)
        legend.AddEntry(b, 'cut (3 sigma)', 'lf')
        legend.Draw()
        return b

    def draw_sigma_selection(self, show=True, redo=False):
        self.draw_disto_fit(cut=self.Cut.generate_special_cut(excluded=['ped_sigma']), logy=True, show=show, redo=redo)
        h = self.Histogram
        x = self.Cut.ped_range
        g = self.draw_box(x[0], -1e9, x[1], 1e9, color=self.FillColor, width=2, fillstyle=3001, show=False)
        l1 = self.make_legend(.7, .63, nentries=3)
        l1.AddEntry(g, 'Pedestal Cut', 'fl')
        l1.AddEntry(h.GetListOfFunctions()[1], 'Fitting Range', 'l')
        l1.AddEntry(h.GetListOfFunctions()[2], 'Fit Function', 'l')
        self.format_histo(h)
        set_drawing_range(h, rfac=1)
        g.Draw('f')
        g.Draw('l')
        self.save_histo(h, 'PedSigmaSelection', show=show, logy=True, l=l1, canvas=gROOT.GetListOfCanvases()[-1], draw_opt='same')
        # self.save_plots('PedSigmaSelection', show=show)

    def draw_pulse_height(self, bin_size=None, cut=None, y_range=None, redo=False, sig=None, rel_t=True, show=True, save=True, prnt=True):
        self.Ana.draw_pulse_height(bin_size=bin_size, cut=cut, y_range=y_range, redo=redo, corr=False, sig=self.SignalName if sig is None else sig, rel_t=rel_t, show=show, save=save, prnt=prnt)

    def draw_signal_time(self, signal_name=None, rel_t=False, show=True):
        signal_name = self.Ana.generate_signal_name(self.SignalName if signal_name is None else signal_name, evnt_corr=False)
        h = TH2F('h_st', 'Pedestal vs. Time', *(self.Ana.get_time_bins() + [160, -80, 80]))
        set_statbox(only_entries=True, x=.83)
        gStyle.SetPalette(53)
        self.Tree.Draw('{name}:time/1000 >> h_st'.format(name=signal_name), self.Cut.all_cut, 'goff')
        self.format_histo(h, x_tit='Time [min]', y_tit='Pulse Height [au]', y_off=1.4, t_ax_off=self.Ana.run.StartTime if rel_t else 0)
        self.save_histo(h, 'PedestalTime', show, lm=.12, draw_opt='colz', rm=.15)
        return h

    def draw_pulse_height_histo(self, signal_name=None, show=True, sigma=False, redo=False):
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

        g = do_pickle(picklepath, func, redo=redo)
        self.format_histo(g, x_tit='Time [min]', y_tit='Mean Pulse Height [au]', y_off=1.6)
        self.save_histo(g, 'Pedestal{s}PulseHeight'.format(s='Sigma' if sigma else ''), show=show, lm=.14, draw_opt='apl')
        return g

    def fit_pulse_height(self, signal_name=None, show=True, sigma=False):
        set_statbox(entries=4, only_fit=True)
        g = self.draw_pulse_height_histo(signal_name, False, sigma)
        fit = g.Fit('pol0', 'qs')
        self.save_histo(g, '{n}Fit'.format(n=g.GetTitle()), show, lm=.12, draw_opt='apl')
        return FitRes(fit)

    def draw_sigma(self, signal_name=None):
        return self.draw_pulse_height_histo(signal_name, sigma=True)

    def compare(self):
        g = self.make_tgrapherrors('g_cp', 'Pedestal Comparison')
        for i, reg in enumerate(self.Run.pedestal_regions):
            name = self.get_signal_name(region=reg)
            fit = self.draw_disto_fit(name, show=False)
            g.SetPoint(i, i, fit.Parameter(1))
            g.SetPointError(i, 0, fit.ParError(1))
        for i, reg in enumerate(self.Run.pedestal_regions):
            g.GetXaxis().SetBinLabel(g.GetXaxis().FindBin(i), reg)
        self.format_histo(g, x_tit='Region Integral', y_tit='Mean Pedestal', y_off=1.4)
        self.save_histo(g, 'PedestalComparison', lm=.12)