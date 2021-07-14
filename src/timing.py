#!/usr/bin/env python
# --------------------------------------------------------
#       timing analysis of the pad waveforms
# created on March 1st 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TF1, TCut, gPad
from src.sub_analysis import PadSubAnalysis
from helpers.draw import *
from string import ascii_lowercase


class TimingAnalysis(PadSubAnalysis):
    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='Timing')
        self.TimingCut = self.Cut.get_timing()

    def draw_all(self):
        self.draw_peaks(show=False, prnt=False, show_cut=True)
        self.draw_peaks_tc(cut=self.get_raw_cut(), show=False, prnt=False)
        self.draw_comparison(show=False)
        self.draw_fine_correction(show=False, prnt=False)

    # --------------------------
    # region GET
    def get(self, par=1, redo=False):
        def f():
            return FitRes(self.draw_peaks(show=False, prnt=False, redo=redo).GetListOfFunctions()[1])
        return fit2u(do_pickle(self.make_simple_pickle_path('PeakVals'), f, redo=redo), par=par)

    def get_peak_var(self, corr=True, fine_corr=False, cut=None, region=None, redo=False):
        fine_corr = ' - {}'.format(self.get_fine_correction(cut, redo=redo)) if fine_corr else ''
        return '{}{}'.format(self.Ana.get_peak_name(t_corr=corr, region=region), fine_corr)
    
    def get_cft_var(self, corr=False):
        return 'cft[{}]{}'.format(self.DUT.Number - 1, '- {}'.format(self.get_tc_correction()) if corr else '')

    def get_raw_var(self, corr=True, channel=None):
        return 'max_peak_{}[{}]'.format('time' if corr else 'position', choose(channel, self.Channel))

    def get_forc_var(self, corr=True):
        return 'forc_time' if corr else 'forc_pos * {}'.format(self.DigitiserBinWidth)

    def get_raw_cut(self, cut=None, nsigma=3):
        f = self.draw_peaks(show=False, fine_corr=False, prnt=False, save=False).GetListOfFunctions()[1]
        return TCut('RawTiming', '({} - {}) / {} < {}'.format(self.Ana.PeakName, f.GetParameter(1), f.GetParameter(2), nsigma)) + self.TimingCut(cut)

    def get_fine_correction(self, cut=None, redo=False):
        def f():
            fit = self.draw_peaks_tc(show=False, prnt=False, cut=self.get_raw_cut(cut), redo=redo).GetListOfFunctions()[1]
            good_pars = fit.GetChisquare() / fit.GetNDF() < 100 and abs(fit.GetParameter(0)) < 10  # require decent chi2 and a meaningful scaling of the sin(x)
            return '({0} * TMath::Sin({1} * (trigger_cell - {2})))'.format(*get_buf(fit.GetParameters(), 3)) if good_pars else '0'
        return do_pickle(self.make_simple_pickle_path('FineCorr', self.Cut(cut).GetName()), f, redo=redo)

    def get_channel_nr(self, name):
        i = where(array([ch_name.lower() for ch_name in self.Run.DigitizerChannels]) == name.lower())[0]
        return i[0] if i.size else warning('There is no Digitiser Channel with {n} in it'.format(n=name))

    def get_rf_channel(self):
        return self.get_channel_nr('rf')

    def get_scint_channel(self):
        return self.get_channel_nr('scint')
    # endregion GET
    # --------------------------

    # --------------------------
    # region RAW
    def draw_raw_peaks(self, xmin=100, xmax=400, bin_width=1., ch=None, corr=False, cut='', show=True):
        x = self.get_tree_vec(self.get_raw_var(corr, ch), self.Cut(cut))
        return self.Draw.distribution(x, make_bins(xmin, xmax, bin_width), x_tit='Time [ns]' if corr else 'Digitiser Bin', show=show, x_range=ax_range(x[(x > xmin) & (x < xmax)], thresh=4))

    @save_pickle('RawMean', suf_args=0)
    def get_raw(self, cut=None, _redo=False):
        x = self.get_tree_vec(self.get_raw_var(corr=True), self.Cut(cut))
        rx = ax_range(self.Ana.SignalRegion * self.DigitiserBinWidth, 0, .5, .5)
        x = x[(x > rx[0]) & (x < rx[1])]
        return fit_fwhm(self.Draw.distribution(x, make_bins(*rx, max(self.DigitiserBinWidth / 2, self.Bins.find_width(x))), show=False))[1]

    def draw_max_peaks(self, cut=None):
        h = self.draw_raw_peaks(0, 512, bin_width=.5, corr=True, show=False, cut=self.TimingCut(cut))
        format_histo(h, y_off=.5, tit_size=.07, lab_size=.06)
        self.Draw(h, 'MaxPeakTimings', logy=True, lm=.073, rm=.045, bm=.18, w=1.5, h=.5, stats=set_entries())

    def create_run_config(self, ch=None, off=0):
        h0, h1 = [self.draw_raw_peaks(max(5, xmin - off), xmax - off, ch=ch, show=False) for xmin, xmax in [(50, 450), (400, 900)]]
        m, bw = int(round(h0.GetMean())), int(round(self.BunchSpacing / self.DigitiserBinWidth))
        print('signal_b_region = [{}, {}]'.format(m - bw // 2, m + bw // 2))
        print('pedestal_ab_region = [{0}, {0}]'.format(m - bw))
        pmin, pmax = array(ax_range(h=h1, thresh=.1 * h1.GetMaximum(), fl=.3, fh=.3)).round().astype('i')
        print('pulser_region = [{}, {}]'.format(pmin, pmax))
        print('pedestal_ac_region = [{0}, {0}]'.format(pmin - bw // 2))
        format_histo(h0, x_range=[m - bw // 2, m + bw // 2], y_off=1.4)
        format_histo(h1, x_range=[pmin, pmax], y_off=1.4)
        c = Draw.canvas(w=2, divide=2)
        for i, h in enumerate((h0, h1), 1):
            c.cd(i)
            h.Draw()

    def print_regions(self, ch=None):
        h = self.draw_raw_peaks(50, 450, ch=ch, show=False)
        m, bw = int(round(h.GetMean())), self.BunchSpacing / self.DigitiserBinWidth
        letters = list(filter(lambda x: x != 'e', ascii_lowercase))
        for i in range(int((self.Run.NSamples - m - bw / 2) / bw) + 1):
            border1 = m + bw / 2 + i * bw
            print('signal_{}_region = [{}, {}]'.format(letters[i + 1], round(border1 - bw), round(border1)))

    # endregion RAW
    # --------------------------

    # --------------------------
    # region CONSTANT FRACTION
    def draw_cft(self, bin_size=.1, corr=False, cut=None, show=True):
        values = self.Run.get_tree_vec(var=self.get_cft_var(corr), cut=self.Cut(cut))
        return self.Draw.distribution(values, self.Ana.get_t_bins(bin_size), 'Constant Fraction Time', x_tit='Constant Fraction Time [ns]', show=show)

    def draw_cft_vs_triggercell(self, bin_size=10, show=True):
        x, y = self.Ana.get_tree_vec(var=['trigger_cell', self.get_cft_var()], cut=self.Cut())
        return self.Draw.profile(x, y, make_bins(0, self.Run.NSamples, bin_size), 'CFT vs. Trigger Cell', x_tit='Trigger Cell', y_tit='Constant Fraction Time [ns]', show=show)

    def get_tc_correction(self, redo=False, show=False):
        def f():
            p = self.draw_cft_vs_triggercell(show=show)
            fit = TF1('f0', '[0]*TMath::Sin([1]*(x - [2])) + [3]', -50, 1024)
            fit.SetParameters(1, 1, 100, mean(get_hist_vec(p, err=False)))
            fit.FixParameter(1, 1 / 1024. * 2 * pi)
            p.Fit(fit, 'q{}'.format('' if show else 0))
            return '({0} * TMath::Sin({1} * (trigger_cell - {2})))'.format(*[fit.GetParameter(i) for i in range(3)])
        return do_pickle(self.make_simple_pickle_path('CFTTC'), f, redo=redo or show)

    def draw_signal_vs_cft(self, bin_size=.1, show=True):
        x, y = self.get_tree_vec(var=[self.get_cft_var(), self.Ana.get_signal_var()], cut=self.Cut())
        self.Draw.profile(x, y, self.Ana.get_t_bins(bin_size), 'Signal vs. Constant Fraction Time', x_tit='Constant Fraction Time [ns]', y_tit='Pulse Height [mV]', show=show)
    # endregion CONSTANT FRACTION
    # --------------------------

    # --------------------------
    # region PEAK TIMING
    def draw_peaks(self, fit=True, cut=None, corr=True, fine_corr=True, show=True, prnt=True, redo=False, save=True, show_cut=False, normalise=None):
        def f():
            x = self.get_tree_vec(var=self.get_peak_var(corr, fine_corr, self.TimingCut(cut), redo=redo), cut=self.TimingCut(cut))
            x *= self.DigitiserBinWidth if not corr else 1
            bin_width = max(.1, min(.4, mean_sigma(x[(x > mean(x) - 5) & (x < mean(x) + 5)])[1].n / 8)) if corr else self.DigitiserBinWidth  # adjust bin width according to width of values
            return self.Draw.distribution(x, make_bins(*self.Ana.get_signal_range(.5, .5), bin_width), '{}Peak Positions'.format('Time Corrected ' if corr else ''), show=False)
        h = do_pickle(self.make_simple_pickle_path('Peak', '{}_{}{}'.format(self.TimingCut(cut).GetName(), int(corr), int(fine_corr))), f, redo=redo)
        if not h.GetEntries():
            return warning('histogram has no entries')
        format_histo(h, x_tit='Time [ns]', y_off=2.0, normalise=normalise, x_range=ax_range(h=h, thresh=.05 * h.GetMaximum(), fl=.2, fh=.7))
        self.fit_peaks(h, fit)
        self.Draw(h, lm=.13, show=show, prnt=prnt, stats=set_statbox(fit=fit, all_stat=not fit, entries=True))
        self.__draw_cut(h, show_cut)
        self.Draw.save_plots('{}PeakPos{}'.format('Raw' if not corr else 'Fine' if fine_corr else '', 'Fit' if fit else ''), show=show, prnt=prnt, save=save)
        return h

    @staticmethod
    def __draw_cut(h, show=True):
        if show:
            fit = h.GetListOfFunctions()[1]
            xmin, xmax = fit.GetParameter(1) - 3 * fit.GetParameter(2), fit.GetParameter(1) + 3 * fit.GetParameter(2)
            b = Draw.box(xmin, -10, xmax, 1e7, line_color=2, width=2, fillstyle=3001, style=7)
            h.Draw('same')
            b.Draw('l')
            Draw.legend([b], ['cut (3 sigma)'], 'lf', .59, y2=.656, margin=.45, scale=1.25)
            gPad.RedrawAxis()

    @staticmethod
    def fit_peaks(h, show=True):
        if show:
            max_val = h.GetBinCenter(h.GetMaximumBin())
            fit1 = h.Fit('gaus', 'qs0', '', max_val - 2, max_val + 2)
            mean_, sigma = fit1.Parameter(1), fit1.Parameter(2)
            fit = TF1('f', 'gaus', mean_ - sigma, mean_ + sigma)
            fit.SetParameters(*fit1.Parameters())
            h.Fit('f', 'q0', '', mean_ - sigma, mean_ + sigma)
            fit2 = TF1('f1', 'gaus', mean_ - 5 * sigma, mean_ + 5 * sigma)
            fit2.SetParameters(fit.GetParameters())
            fit2.SetLineStyle(2)
            fit.SetLineColor(1)
            fit2.SetLineColor(1)
            h.GetListOfFunctions().Add(fit2)

    def draw_comparison(self, show=True):
        histos = [self.draw_peaks(show=False, corr=c, fine_corr=fc, fit=f, prnt=False) for c, fc, f in array([(0, 0, 0), (1, 0, 1), (1, 1, 1)], dtype=bool)]
        titles = ['raw'] + ['sigma {} corr. {:1.2f}'.format(n, s) for n, s in zip(['time', 'fine'], [f.GetParameter(2) for f in [h.GetListOfFunctions()[1] for h in histos[1:]]])]
        s = self.Draw.stack(histos, 'Time Comparison', titles, scale=True, w=.3, show=show)
        format_histo(s, x_range=ax_range(h=histos[0], thresh=.02 * histos[0].GetMaximum(), fl=.2, fh=.2))

    def draw_peaks_tc(self, bin_size=4, corr=True, fine_corr=False, fit=True, cut=None, show=True, prnt=True, save=True, redo=False):
        def f():
            x, y = self.get_tree_vec(var=['trigger_cell', self.get_peak_var(corr, fine_corr)], cut=self.TimingCut(cut))
            return self.Draw.profile(x, y, make_bins(0, self.Run.NSamples, bin_size), '{}Peak Position vs Trigger Cell'.format('Corrected ' if corr else ''), show=False)
        h = do_pickle(self.make_simple_pickle_path('PeakTC', '{}{}{}'.format(self.TimingCut(cut).GetName(), bin_size, int(fine_corr))), f, redo=redo)
        self.fit_sin(h, fit)
        format_histo(h, x_tit='Trigger Cell', y_tit='Signal Peak Time [ns]', y_off=1.8, stats=fit)
        self.Draw(h, 'OriPeakPosVsTriggerCell', show, lm=.13, prnt=prnt, save=save, stats=set_statbox(fit=True, center_x=True))
        return h

    def fit_sin(self, h, fit):
        fit_func = TF1('fcor', '[0]*TMath::Sin([1]*(x - [2])) + [3]', -50, self.Run.NSamples)
        fit_func.SetParNames('Scale', 'Period', 'Phase', 'Offset')
        fit_func.SetParameters(1, 1 / self.Run.NSamples * 2 * pi, 100, mean([h.GetBinContent(i) for i in range(h.GetNbinsX())]))
        fit_func.FixParameter(1, 1 / self.Run.NSamples * 2 * pi)
        h.Fit(fit_func, 'q0')
        h.GetListOfFunctions().Add(fit_func) if fit else do_nothing()

    def draw_fine_correction(self, bin_size=12, show=True, prnt=True):
        profiles = [self.draw_peaks_tc(fine_corr=i, bin_size=bin_size, show=False, fit=False) for i in [0, 1]]
        for p in profiles:
            format_histo(p, y_range=ax_range(min(p.GetMinimum() for p in profiles), max(p.GetMaximum() for p in profiles), .1, .3), markersize=.5, color=self.Draw.get_color(2))
        self.Draw(profiles[0], show=show)
        profiles[1].Draw('same')
        Draw.legend(profiles, ['Uncorredted', 'Fine Correction'], 'pl')
        self.Draw.save_plots('FineCorrection', prnt=prnt)
    # endregion PEAK TIMING
    # --------------------------

    # --------------------------
    # region TRIGGER CELL
    def draw_trigger_cell(self, bin_width=1, show=True, cut=None):
        h = self.Draw.distribution(self.get_tree_vec(var='trigger_cell', cut=self.Cut(cut)), make_bins(0, self.Run.NSamples, bin_width), 'Trigger Cell', show=show)
        format_histo(h, x_tit='Trigger Cell', y_range=[0, 1.2 * h.GetMaximum()])
        h.Fit('pol0', 'qs')
        format_statbox(h, fit=True, entries=True)
        self.Draw.save_plots('TriggerCell')

    def draw_intlength_tc(self, bin_width=4, fit=True, show=True):
        x, y = self.get_tree_vec(var=['trigger_cell', 'IntegralLength[{}]'.format(self.Ana.get_signal_number())], cut=self.Cut())
        p = self.Draw.profile(x, y, self.Bins.get_wf(bin_width), 'Integral Length vs. Triggercell', x_tit='Triggercell', y_tit='Integral Length [ns]', show=show, stats=fit)
        self.fit_sin(p, fit)
        format_statbox(p, fit=fit)

    def draw_forc_tc(self, cut=None, corr=True, show=True):
        x, y = self.get_tree_vec(var=['trigger_cell', self.get_forc_var(corr)], cut=self.Cut(cut))
        self.Draw.histo_2d(x, y, self.Bins.get_wf() + self.Bins.get_wf(.5), 'FORC Timing vs. Trigger Cell', pal=55, x_tit='Trigger cell', y_tit='FORC Timing [ns]', stats=0,
                           y_range=ax_range(y, 0, .5, .5, thresh=2), show=show)

    def draw_tcal(self, bin_size=1, fit=False, show=True):
        x, y = arange(self.Run.TCal.size), self.Run.TCal
        g = self.Draw.profile(x, y, self.Bins.get_wf(bin_size), 'DRS4 Bin Sizes', markersize=.5, x_tit='Bin Number', y_tit='Size [ns]', y_range=ax_range(y, 0, .2, .2), w=1.5, h=.75, show=show)
        g.Fit('pol1', 'q') if fit else do_nothing()
        format_statbox(g, fit=fit, entries=True)

    def draw_tcal_disto(self, bin_size=.01, show=True):
        self.Draw.distribution(self.Run.TCal, make_bins(0, self.DigitiserBinWidth * 2, bin_size), 'TCal Distribution', x_tit='Bin Size [ns]', show=show)
    # endregion TRIGGER CELL
    # --------------------------

    # --------------------------
    # region RF
    def fit_rf(self, evnt=0, t_corr=True, show=True):
        gr, n = self.Ana.Waveform.draw(start_event=evnt, t_corr=t_corr, show=show, channel=self.Run.get_rf_channel(), cut='')
        fit = TF1('f', '[0] * TMath::Sin((x+[1])*2*pi/[2])+[3]', 0, 500)
        fit.SetParNames('Scale', 'Period', 'Phase', 'Offset')
        fit.SetNpx(1000)
        fit.SetParameters(100, 20, 3, -40)
        fit.SetParLimits(2, -20, 20)
        fit_res = gr.Fit('f', 'qs{}'.format('' if show else 0))
        format_statbox(gr, fit=True)
        print(evnt, fit.GetChisquare() / fit.GetNDF())
        return FitRes(fit_res)

    def draw_rf_phases(self, n=100, start_event=0, show=True):
        y = [fit2u(self.fit_rf(i + start_event, show=False), par=2) for i in range(n)]
        self.Draw.graph(arange(len(y)), y, title='RF Phases', x_tit='event', y_tit='RF Phase', show=show)

    def draw_rf_period(self, cut=None, show=True):
        x = self.get_tree_vec(var='rf_period', cut=self.Cut(cut))
        self.Draw.distribution(x, make_bins(19.7, 19.8, n=500), 'RF Period', x_tit='RF Period [ns]', ndivx=507, show=show)

    def draw_rf_phase(self, cut=None, show=True):
        x = self.get_tree_vec(var='rf_phase', cut=self.Cut(cut) + TCut('rf_chi2 < 100'))
        self.Draw.distribution(x, make_bins(-30, 30, n=500), 'RF Phase', x_tit='RF Phase [ns]', show=show)

    def draw_rf_vs_peak(self, show=True):
        x, y = self.get_tree_vec(var=['signal_peak_time[{}]'.format(self.Channel), 'rf_phase'], cut=self.Cut() + TCut('rf_chi2 < 100'))
        h = self.Draw.histo_2d(x, y, make_bins(0, 500, .25) + make_bins(-12, 12, .1), 'RF Phase vs. Peak Timings', x_tit='Signal Peak Timing [ns]', y_tit='RF Phase [ns]', stats=0, show=show)
        format_histo(h, **{n: v for n, v in zip(['x_range', 'y_range'], ax_range(5, 5, .2, .2, h))})
    # endregion RF
    # --------------------------

    # --------------------------
    # region MISCELLANEOUS
    def draw_fit_peak_timing(self, show=True):
        x = self.get_tree_vec(var='fit_peak_time[{}]'.format(self.Channel), cut=self.Cut())
        self.Draw.distribution(x, make_bins(*self.Ana.get_signal_range(), n=sqrt(x.size)), 'Fitted Peak Positions', x_tit='Time [ns]', show=show)

    def draw_peaking_time(self, show=True):
        x = self.get_tree_vec(var='peaking_time[{}]'.format(self.Channel), cut=self.Cut())
        self.Draw.distribution(x, make_bins(0, 20, n=sqrt(x.size)), 'Peaking Time', x_tit='Time [ns]', show=show)

    def draw_forc_times(self, bin_size=.5, corr=False, show=True):
        x = self.get_tree_vec(var=self.get_forc_var(corr), cut=self.Cut())
        self.Draw.distribution(x, self.Bins.get_wf(bin_size), 'FORC Timing', x_tit='Time [ns]', x_range=ax_range(x, 0, .2, .2, thresh=2), show=show)

    def draw_signal_peak(self, cut=None, corr=False, show=True):
        x = self.get_tree_vec(var='signal_peak_time[{}]{}'.format(self.Channel, '+rf_phase' if corr else ''), cut=self.Cut(cut))
        h = self.Draw.distribution(x, make_bins(0, 500, .1 if corr else .2), 'Signal Peak Timings', x_tit='Signal Peak Timing [ns]', show=show)
        format_histo(h, x_range=ax_range(5, 5, .1, .3, h))

    def draw_scint_inflexion(self, corr=False, show=True):
        return self.draw_inflexion_time(corr=corr, channel=self.get_scint_channel(), show=show)

    def draw_peak_width(self, show=True):
        x = abs(self.get_tree_vec(var='rise_width[{}]'.format(self.Channel), cut=self.Cut()))
        return self.Draw.distribution(x, make_bins(0, 1, n=500), 'Peak Width', x_tit='Rise Width [ns]', show=show, x_range=ax_range(x, 0, .1, .3, thresh=3))

    def draw_threshold(self, corr=False, channel=None, show=True):
        ch = choose(channel, self.Channel)
        x = self.get_tree_vec(var='signal_peak_time[{}]{}'.format(ch, '+rf_phase' if corr else ''), cut=self.Cut() + TCut('t_thresh[{}] > 10'.format(ch)))
        return self.Draw.distribution(x, make_bins(0, 500, .1 if corr else .2), 'Time over Threshold', x_tit='ToT [ns]', show=show, x_range=ax_range(x, 0, .1, .3, thresh=3))

    def draw_scint_threshold(self, corr=False, show=True):
        return self.draw_threshold(corr=corr, channel=self.get_scint_channel(), show=show)

    def draw_inflexion_time(self, corr=False, channel=None, show=True):
        x = self.get_tree_vec(var='rise_time[{}]{}'.format(choose(channel, self.Channel), '+rf_phase' if corr else ''), cut=self.Cut() + TCut('rise_time[{}]'.format(choose(channel, self.Channel))))
        h = self.Draw.distribution(x, make_bins(0, 20, .05), 'Inflexion Time', x_tit='Rise Time [ns]', show=show)
        format_histo(h, x_range=ax_range(5, h.GetMaximum() * .01, .1, .3, h))

    def draw_inter_dia_corr(self, show=True):
        x, y = self.get_tree_vec(var=['rise_time[{}]'.format(ch) for ch in self.Run.Channels], cut=self.Cut() + TCut('rise_time[{}]'.format(self.Channel)))
        h = self.Draw.histo_2d(x, y, make_bins(0, 20, .05) * 2, 'Inflextion Times of Diamond Signals', x_tit='Inflexion Time1 [ns]', y_tit='Inflexion Time2 [ns]', show=show)
        format_histo(h, **{n: v for n, v in zip(['x_range', 'y_range'], ax_range(5, 5, .2, .2, h))})

    def draw_inter_dia(self, show=True):
        x = self.get_tree_vec(var='rise_time[{}] - rise_time[{}]'.format(*self.Run.Channels), cut=self.Cut() + TCut('rise_time[{}]'.format(self.Channel)))
        h = self.Draw.distribution(x, make_bins(-10, 10, .05), 'Peak Width', x_tit='Rise Time Difference [ns]', show=show)
        format_histo(h, x_range=ax_range(5, 5, .1, .3, h))
    # endregion MISCELLANEOUS
    # --------------------------
