#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from warnings import simplefilter
from ROOT import TH1F, TCut, TProfile, THStack, TH2F, TMath
from ROOT.gRandom import Landau
from numpy import polyfit, RankWarning, vectorize, size, repeat, argmax, insert, array_split, invert
from numpy.random import normal, rand
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import poisson

from src.sub_analysis import PadSubAnalysis
from helpers.draw import *


class PeakAnalysis(PadSubAnalysis):

    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='Peaks')
        self.WF = self.Ana.Waveform
        self.NoiseThreshold = self.calc_noise_threshold()
        self.Threshold = max(self.NoiseThreshold, self.Ana.get_min_signal(self.Ana.get_signal_name(peak_int=1)))
        self.Prominence = 5
        self.Distance = 35
        self.BinWidth = self.DigitiserBinWidth
        self.StartAdditional = self.get_start_additional()
        self.NBunches = self.calc_n_bunches()
        self.BunchSpacing = self.Ana.BunchSpacing

    # ----------------------------------------
    # region INIT
    def calc_noise_threshold(self):
        """ return peak threshold, 5 times the raw noise + mean of the noise. """
        return abs(self.Ana.Pedestal.get_raw_mean() + 5 * self.Ana.Pedestal.get_raw_noise()).n

    def get_start_additional(self, b0=None):
        """Set the start of the additional peaks 2.5 bunches after the signal peak to avoid the biased bunches after the signal. Signalbunch = bunch 1 """
        return int(self.Run.IntegralRegions[self.DUT.Number - 1]['signal_a'][0] + self.Ana.BunchSpacing * (choose(b0, 4) - 1.5) / self.BinWidth)

    def get_bunch_range(self, b0=None, b_end=None, bins=False):
        a = array([self.get_start_additional(b0), self.Run.NSamples if b_end is None else self.get_start_additional(b_end)])
        return [int(i) for i in a] if bins else [float(i) for i in a * self.BinWidth]

    def calc_n_bunches(self):
        return arange(self.StartAdditional, self.Run.NSamples, self.BunchSpacing / self.BinWidth).size
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_binning(self, bin_size=None):
        return self.WF.get_binning(bin_size)

    def get_from_tree(self, redo=False, cut=...):
        return array(do_hdf5(self.make_simple_hdf5_path('Peaks'), self.Run.get_tree_vec, var=self.Ana.PeakName, dtype='f4', redo=redo))[cut]

    def get_cut(self, cut=None):
        return self.WF.get_cut(cut)

    def get_signal_values(self, f, ind=None, default=-1, *args, **kwargs):
        signal_ind, noind = self.get_signal_indices()
        values = insert(array(f(*args, **kwargs))[signal_ind], array(noind), default)
        return values if ind is None else values[ind]

    def get_all_cft(self, thresh=.5):
        return self.get_signal_values(self.find_all_cft, thresh=thresh)

    def get_cft(self, flat=False, ind=None):
        cft, n_peaks = self.find_all_cft(), self.find_all()[-1]
        if flat:
            return array(cft)
        cft = array(split(cft, cumsum(n_peaks)[:-1]))
        return cft if ind is None else cft[ind]

    def get_all_tot(self, thresh=None, fixed=True):
        return self.get_signal_values(self.calc_all_tot, thresh, fixed)

    def get_signal_heights(self):
        return self.get_signal_values(self.get_heights, default=-999, flat=True)

    def get_signal_times(self, fit=False, ind=None):
        return self.get_signal_values(self.get, ind, default=-999, flat=True, fit=fit)

    def get_signal_indices(self):
        """ :returns: all indices of the times in the signal region and the indices with no too many signals """
        values, n, m = self.get(flat=True), self.find_all()[2], mean(self.get_from_tree())
        w = self.BunchSpacing / 2
        signals = split((m - w < values) & (values < m + w), cumsum(n)[:-1])  # check if the times are in the signal region
        n_signals = array([count_nonzero(i) for i in signals])  # get the number of signals in the signal region
        no_signals = where((n_signals == 0) | (n_signals > 1))[0]
        n_signals = repeat(n_signals, n)
        return where((m - w < values) & (values < m + w) & (n_signals == 1))[0], no_signals - arange(no_signals.size)

    def get_all(self, cut=None, thresh=None, fit=False, redo=False):
        t, h, n = self.find_all(thresh, fit, redo)
        cut = self.get_cut(cut)
        lcut = cut.repeat(n)  # make cut for the flat arrays
        return array(t)[lcut], array(h)[lcut], array(n)[cut]

    def get(self, flat=False, thresh=None, fit=False, cut=None, i=0):
        cut = self.get_cut(cut)
        t, h, n = self.get_all(cut, thresh, fit)
        data = [t, h]
        return array(data[i]) if flat else array(split(data[i], cumsum(n)[:-1]), dtype=object)

    def get_heights(self, flat=False, thresh=None, fit=False, cut=None):
        return self.get(flat, thresh, fit, cut, i=1)

    def get_bunch_height(self, n=1, thresh=80, cut=-1, redo=False):
        def f():
            x = self.find_bunch_heights(n, cut=self.make_bunch_cut(n) if cut is None else self.Ana.get_pulser_cut() & (self.get_bunch_cut(n + cut) if isint(cut) else cut))
            return mean_sigma(x[(x > thresh) & (x < 500)])[0]
        return f() if is_iter(cut) else do_pickle(self.make_simple_pickle_path('BHeight', '{}_{}_{}'.format(n, thresh, cut)), f, redo=redo)

    def get_tbin_cut(self, ibin, bin_size=None, corr=True, fine_corr=False):
        cut_name = self.Ana.Timing.get_peak_var(corr, fine_corr)
        vmin, vmax = self.Ana.get_t_bins(bin_size)[1][ibin:ibin + 2]
        return TCut('tbin{}'.format(ibin), '{0} < {n} && {n} < {1}'.format(vmin, vmax, n=cut_name))

    def get_t_indices(self, ibin, bin_size=None):
        bins = self.Ana.get_t_bins(bin_size)[1]
        values = array(self.get_from_tree())
        return where((bins[ibin] <= values) & (values < bins[ibin + 1]))[0]

    def get_ph_indices(self, pmin, pmax):
        heights = self.get_signal_heights()
        return where((pmin <= heights) & (heights <= pmax))[0]

    def get_cft_indices(self, ibin, thresh=None, bin_size=None):
        bins = self.Ana.get_t_bins(bin_size)[1]
        values = self.get_all_cft(thresh)
        return where((bins[ibin] <= values) & (values < bins[ibin + 1]))[0]

    def get_indices(self, t_min, t_max):
        s1_indices = self.get_n()
        times = concatenate(self.get()[s1_indices])
        x = self.get_corrected_times(times[times > self.StartAdditional * self.BinWidth] % self.BunchSpacing, events=s1_indices)
        return s1_indices[where((t_min <= x) & (x <= t_max))]

    def get_event(self, i):
        return self.WF.get_all_times()[i], self.WF.get_all()[i], self.get_heights()[i], self.get()[i]

    def get_n_additional(self, start_bunch=None, end_bunch=None, thresh=None):
        def f():
            values = self.find_n_additional(start_bunch, end_bunch, thresh)
            m = mean(values)
            return ufloat(m, sqrt(m / values.size))
        suffix = '{}_{}'.format(start_bunch, end_bunch) if start_bunch is not None else ''
        return do_pickle(self.make_simple_pickle_path('NAdd', '{}{}'.format(suffix, '' if thresh is None else '{:1.0f}'.format(thresh))), f)

    def get_corrected_times(self, times, n_peaks=None, events=None, cut=None):
        signal_peak_times = repeat(self.get_from_tree(cut=self.get_cut(cut)), n_peaks) if events is None else array(self.get_from_tree())[events]
        return times - (signal_peak_times - self.get_from_tree(cut=self.get_cut())[0])  # subtract the first from the signal events

    def get_n_times(self, n=2, ret_indices=False):
        """ :returns all times with exactly [n] additional peaks. """
        times = self.get()
        for i in range(times.size):
            times[i] = times[i][(times[i] > self.StartAdditional * self.BinWidth)]
        return where(array([lst.size for lst in times]) == n)[0] if ret_indices else times[array([lst.size for lst in times]) == n]  # select all events with two additional peaks

    def get_n(self, n=1):
        return self.get_n_times(n, ret_indices=True)

    def get_n_total(self, cut=None, b0=None, b_end=None):
        return self.draw(show=False, cut=cut).Integral(*self.get_bunch_range(b0, b_end, bins=True))

    def get_lambda(self, n_peaks=None, lam=None, cut=None):
        cut = choose(cut, self.Ana.get_pulser_cut())
        lam = choose(lam, choose(self.get_n_total(cut), n_peaks) / cut[cut].size)
        return ufloat(lam, sqrt(lam / (cut[cut].size - 1)))

    def get_flux(self, n_peaks=None, lam=None, prnt=True):
        flux = self.get_lambda(n_peaks, lam) / (self.Ana.BunchSpacing * self.NBunches * self.get_area()) * 1e6
        self.info('Estimated Flux by number of peaks: {}'.format(make_flux_string(flux)), prnt=prnt)
        return flux

    def get_n_consecutive(self, n=2):
        cut = self.Ana.get_pulser_cut()
        lam = self.get_lambda(cut=cut).n / self.NBunches  # lambda for a single bunch
        n_events = cut[cut].size
        p = 1 - poisson.cdf(0, lam)  # probability to have at least one peak in a single bunch
        return sum(p ** (n - 1) * poisson.pmf(i, lam) * i for i in range(50)) * n_events

    def get_peak_fac(self):
        h = self.draw(cut=self.Ana.get_pulser_cut(), corr=False, show=False)
        facs = [self.get_n_bunch(i, h) / h.Integral(*self.get_bunch_range(i, i + 1, bins=True)) for i in range(4, self.NBunches + 3)]
        return mean_sigma(facs)[0]

    def get_bunch_eff(self, n=0):
        cut = self.Ana.get_pulser_cut() & self.Ana.get_event_cut('!bucket1[{}]'.format(self.Channel))
        return self.get_n_bunch(n, cut=cut) / self.get_peak_fac() / cut[cut].size

    def get_n_bunch(self, n=0, h=None, cut=None, show=False):
        h = self.draw(corr=False, show=show, cut=choose(cut, self.Ana.get_pulser_cut())) if h is None else h
        f = Draw.make_f('f', 'gaus', *self.get_bunch_range(n, n + 1))
        h.Fit(f, 'q', '', *self.get_bunch_range(n, n + 1))
        if show:
            f.Draw('same')
        return f.Integral(0, self.Run.NSamples) / self.BinWidth

    def get_area(self, bcm=False):
        """ :returns area of the DUT in cm^2"""
        return self.get_bcm_area() if bcm else self.DUT.ActiveArea * .01

    def get_bcm_area(self):
        """ return the total area of the BCM' pad sizes """
        i = int(self.Ana.DUT.Name.split('-')[-1]) - 1
        base_length = 0.0928125  # [cm]
        spacing = 0.0025
        radius = 0.0049568
        rounded_edge = radius ** 2 * (4 - pi)
        area = base_length ** 2
        return 2 ** i * area + self.get_spacings(i, spacing, base_length) - rounded_edge

    def get_spacings(self, i, spacing, length):
        """ return the additional spacings for the BCM' pad sizes """
        if i > 0:
            j = 2 ** (i / 2)
            return spacing * (j * length + (j - 1) * spacing) + 2 * self.get_spacings(i - 1, spacing, length)
        else:
            return 0

    def get_mean_sigma(self):
        v = self.get_signal_times()
        return mean_sigma(v[v > 0])

    def get_mpv_sigma_heights(self, redo=False):
        def f():
            h = self.draw_height_disto(show=False)
            max_bin = h.GetMaximumBin()
            fit = h.Fit('landau', 'qs0', '', h.GetBinCenter(max_bin - 10), h.GetBinCenter(max_bin + 30))
            return fit.Parameter(1), fit.Parameter(2)
        return do_pickle(self.make_simple_pickle_path('H'), f, redo=redo)

    def get_signal_ph(self):
        values = self.get_signal_heights()
        weights, bins = histogram(values, self.Ana.Bins.get_pad_ph()[1])
        bin_centers = (bins + (bins[1] - bins[0]) / 2)[:-1]
        m, s = mean_sigma(bin_centers, weights)
        return ufloat(m, s / sqrt(values.size)) - self.Ana.Pedestal.get_raw_mean()

    def get_bunch_cut(self, n, cut=...):
        n_peaks = self.find_n_additional(n, n + 1, cut=cut)
        return self.Ana.make_event_cut(self.Ana.get_events(cut)[n_peaks > 0])

    def make_bunch_cut(self, n=1):
        return self.Ana.get_pulser_cut() & invert(self.get_bunch_cut(n - 1)) & invert(self.get_bunch_cut(n + 1))
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_thresh(self, thresh=None):
        self.Draw.info('Threshold: {:.1f}'.format(choose(thresh, self.Threshold)))

    def draw(self, corr=False, split_=1, thresh=None, bin_size=None, y_range=None, cut=None, show=True, redo=False):
        times, heights, n_peaks = self.get_all(cut=cut, thresh=thresh, redo=redo)
        times = self.get_corrected_times(times, n_peaks, cut=cut) if corr else times
        times = array_split(times, split_)
        h = [TH1F(Draw.get_name('h'), 'Peak Times', *self.get_binning(bin_size)) for _ in range(split_)]
        for i in range(split_):
            v = times[i]
            h[i].FillN(v.size, array(v, 'd'), ones(v.size))
        h = h[0] if split_ == 1 else h
        if show and split_ == 1:
            set_statbox(entries=True)
            format_histo(h, x_tit='Time [ns]', y_tit='Number of Peaks', y_off=1.3, fill_color=Draw.FillColor, y_range=y_range)
            self.Draw(h, lm=.12, show=show, w=1.5, h=0.75, logy=True, stats=True)
        return h

    def draw_events(self, events=None, show=True):
        t, h, n = self.get_all(cut=events)
        h = self.Draw.distribution(t, self.get_binning(.5), 'Peak Times', x_tit='Time [ns]', y_tit='Number of Peaks', show=show, w=1.5, h=0.75)
        format_statbox(h, entries=True)

    def draw_signal(self, bin_size=.5, ind=None, fit=False, y=None, x=None, x_range=None, y_range=None, show=True, draw_ph=False, smear=None):
        times = choose(x, self.get_signal_times, fit=fit, ind=ind)
        self.smear_times(times, smear)
        h = self.Draw.distribution(times, self.Ana.get_t_bins(bin_size), 'Signal Peak Times', lm=.13, rm=.12 if draw_ph else None, show=show, x_tit='Signal Peak Time [ns]', y_off=1.8, x_range=x_range)
        c = self.draw_ph(get_last_canvas(), bin_size, times, y, x_range, y_range, draw_ph)
        format_statbox(h, all_stat=True, c=c)
        return h

    def draw_ph(self, c, bin_size, x, y, x_range, y_range, show):
        if show:
            p = self.Ana.draw_ph_peaktime(show=False, bin_size=bin_size, x=x, y=y)
            values = get_hist_vec(p, err=False)
            format_histo(p, title=' ', stats=0, x_tit='', l_off_x=1, y_range=choose(y_range, ax_range(min(values[values > 0]), max(values), .3, .3)), x_range=x_range)
            c.cd()
            Draw.tpad('psph', transparent=True, lm=.13, rm=.12)
            p.Draw('y+')
            update_canvas(c)
            return c

    def draw_heights_vs_time(self, bin_size=.5, corr=True, cft=False, show=True):
        times, heights, n_peaks = self.find_all()
        times = self.get_corrected_times(times, n_peaks) if corr else self.find_all_cft() if cft else times
        p = TProfile('pph', 'Peak Heights', *self.get_binning(bin_size))
        p.FillN(times.size, array(times).astype('d'), array(heights).astype('d'), ones(times.size))
        format_histo(p, x_tit='Time [ns]', y_tit='Peak Height [mV]', y_off=1.3, stats=0, fill_color=Draw.FillColor)
        self.Draw(p, lm=.12, show=show, w=1.5, h=0.75)

    def draw_signal_height_vs_time(self, bin_size=None, show=True):
        p = TProfile('pspt', 'Peak Height vs. Constant Fraction Time', *self.Ana.get_t_bins(bin_size))
        fill_hist(p, x=self.get(flat=True), y=self.get_heights(flat=True))
        format_histo(p, x_tit='Constrant Fraction Time [ns]', y_tit='Peak Height [mV]', y_off=1.4)
        self.Draw(p, show=show, lm=.12, stats=set_statbox(entries=True))
        return p

    def draw_combined_heights(self, hist=False, show=True):
        s1_indices = self.get_n()
        times, heights = concatenate(self.get()[s1_indices]), concatenate(self.get_heights()[s1_indices])
        x = self.get_corrected_times(times[times > self.StartAdditional * self.BinWidth] % self.BunchSpacing, events=s1_indices)
        y = heights[times > self.StartAdditional * self.BinWidth]
        p = TProfile('pcph', 'Combined Bunch Pulse Heights', int(ceil(self.BunchSpacing)) * 2, 0, ceil(self.BunchSpacing))
        p = TH1F('hcph', 'Combined Number of Peaks for all Bunches', int(ceil(self.BunchSpacing)) * 2, 0, ceil(self.BunchSpacing)) if hist else p
        p.FillN(x.size, x.astype('d'), ones(x.size)) if hist else p.FillN(x.size, x.astype('d'), y.astype('d'), ones(x.size))
        format_histo(p, x_tit='Time [ns]', y_tit='Number of Peaks' if hist else 'Peak Height [mV]', y_off=1.6, stats=0, fill_color=Draw.FillColor)
        self.Draw(p, lm=.12, show=show)

    def draw_n(self, do_fit=False, show=True):
        n_peaks = self.find_n_additional()
        h = self.Draw.distribution(n_peaks, make_bins(11), 'Number of Peaks', show=False)
        if do_fit:
            fit_poissoni(h, show=show)
        format_histo(h, x_tit='Number of Peaks', y_tit='Number of Entries', y_off=1.4, fill_color=Draw.FillColor, lw=2)
        self.Draw(h, 'PeakNumbers{}'.format('Fit' if do_fit else ''), show, logy=True, lm=.11, stats=set_statbox(fit=do_fit, entries=True, w=None if do_fit else .2))
        self.get_flux(n_peaks)
        return h

    def draw_n_bunch(self, b=0, y_range=None, show=True):
        bunches, times = self.find_bunches(), self.get_n_times(n=2)
        times = concatenate(times[[any((bunches[b] < lst) & (lst < bunches[b + 1])) for lst in times]])  # select events with a peak in bunch b
        set_statbox(entries=True, w=.2)
        self.Draw.distribution(times, self.get_binning(), '2 Additional Peaks for Bunch {}', x_tit='Time [ns]', y_range=y_range, show=show, w=1.5, h=0.75, logy=True, stats=True)
        return times

    def draw_bunch_systematics(self, n=None, show=True):
        n = self.NBunches - 1 if n is None else n
        bunches, times = self.find_bunches(), self.get_n_times(n=2)
        peaks_per_bunch = zeros(n + 1, dtype='u4')
        self.PBar.start(self.NBunches)
        for b in range(self.NBunches):
            t = times[[any((bunches[b] < lst) & (lst < bunches[b + 1])) for lst in times]]  # select events with a peak in bunch b
            for lst in t:
                for i in arange(n + 1):
                    if b + i < self.NBunches:
                        if any((bunches[i + b] < lst) & (lst < bunches[i + b + 1])):
                            peaks_per_bunch[i] += 1
            self.PBar.update()
        peaks_per_bunch = array([ufloat(v, sqrt(v)) for v in peaks_per_bunch]) / arange(self.NBunches, self.NBunches - n - 1, -1)  # add errors and normalise
        n_peaks = peaks_per_bunch[1:]  # exclude the signal peak
        self.Draw.graph(arange(1, n_peaks.size + 1), n_peaks, title='Bunch Systematics', x_tit='Bunch after a Signal', y_tit='Average Number of Peaks', show=show)

    def draw_height_disto(self, show=True):
        return self.Draw.distribution(self.get_signal_heights(), self.Bins.get_pad_ph(), 'Signal Heights', lm=.12, show=show, x_tit='Peak Height [mV]', y_off=1.5)

    def draw_flux_vs_threshold(self, tmin=None, tmax=None, steps=20):
        x = linspace(choose(tmin, self.NoiseThreshold), choose(tmax, self.Threshold), steps)
        y = array([self.get_flux(lam=self.get_n_additional(thresh=ix)) for ix in x]) / 1000.
        return self.Draw.graph(x, y, title='Flux vs. Peak Threshold', x_tit='Peak Finding Threshold [mV]', y_tit='Flux [MHz/cm^{2}]')

    def draw_spacing(self, overlay=False, bin_size=.5, show=True):
        values = concatenate(self.get() - self.get_signal_times())
        values = values[values != 0]
        values = ((values + self.BunchSpacing / 2) % self.BunchSpacing + self.BunchSpacing / 2) if overlay else values
        m, w = mean(values), self.Ana.get_t_bins()[1][-1] - self.Ana.get_t_bins()[1][0]
        bins = make_bins(m - w / 2, m + w / 2, bin_size) if overlay else self.get_binning(bin_size)
        x, y = [1, 1] if overlay else [1.5, .75]
        self.Draw.distribution(values, bins, 'Peak Spacing', x_tit='Peak Distance [ns]', w=x, h=y, show=show, stats=set_statbox(entries=True, w=.2))

    def draw_additional_disto(self, n_splits=4, show=True):
        x = [v.n for h in self.draw(split_=n_splits) for v in self.draw_additional(h, show=False)]
        self.Draw.distribution(x, make_bins(*ax_range(x, 0, .3, .3), n=sqrt(len(x))), 'Peak Heights', x_tit='Peak Height', show=show, stats=0)

    def draw_additional(self, h=None, show=True):
        values = get_hist_vec(self.draw(show=False) if h is None else h)[self.StartAdditional:]
        peaks = find_peaks([v.n for v in values], height=max(values).n / 2., distance=self.Ana.BunchSpacing)
        g = self.Draw.graph((peaks[0] + self.StartAdditional) / 2, values[peaks[0]], title='Additional Peak Heights', show=show, w=1.5, h=.75, gridy=True)
        g.Fit('pol0', 'qs')
        format_statbox(g, fit=True)
        return values[peaks[0]]

    def draw_bunch_height(self, n=1, bin_width=5, cut=None, show=True):
        x = self.find_bunch_heights(n, cut=self.make_bunch_cut(n) if cut is None else self.Ana.get_pulser_cut() & cut)
        return self.Draw.distribution(x, self.Bins.get_pad_ph(bin_width), 'Peak Height B{}'.format(n), x_tit='Peak Height [mV]', show=show)

    def compare_bunch_heights(self, n, bin_width=5, fill=None, show=True):
        h = [self.draw_bunch_height(i, bin_width, show=False) for i in make_list(n)]
        self.Draw.stack(h, 'Peak Height Comparison', ['B{}'.format(i) for i in make_list(n)], scale=True, show=show, fill=fill, opacity=.3)

    def draw_bunch_heights(self, thresh=80, cut=None, show=True):
        x, heights = arange(1, self.NBunches), []
        self.PBar.start(x.size)
        for i in x:
            heights.append(self.get_bunch_height(i, thresh=thresh, cut=cut))
            self.PBar.update()
        g = Draw.make_tgrapherrors(x, heights, x_tit='Bunch', y_tit='Peak Height [mV]')
        ga = Draw.make_tgrapherrors([ufloat(mean(x[3:]), (x[-1] - x[3]) / 2)], [mean_sigma(heights[3:])[0]], color=2)
        return self.Draw.multigraph([g, ga], 'Bunch Heights', ['single', 'average'], w=2, show=show, color=False, gridy=True)

    def draw_pre_bunch_heights(self, b=-1, thresh=80, y_range=None):
        g = [ig for cut in [None, b] for ig in self.draw_bunch_heights(thresh=thresh, show=False, cut=cut).GetListOfGraphs()]
        tits = ['single peak', 'single average', 'peak in bunch {}'.format(b), 'B{} average'.format(b)]
        self.Draw.multigraph(g, 'Bunch Peak Heights', tits, w=2, y_range=y_range, gridy=True)

    def draw_height_evolution(self, n=1, bin_size=10000, thresh=60):
        x = self.find_bunch_heights(n)
        x = x[(x > thresh) & (x < 500)]
        x = [mean_sigma(ix)[0] for ix in array_split(x, int(x.size / bin_size))]
        self.Draw.graph(arange(len(x)), x)

    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region FIND
    def find_bunches(self, center=False, bin_size=.5, show=False):
        h = self.draw(corr=False, show=show, bin_size=bin_size)
        values = get_hist_vec(self.draw(corr=False, show=False, bin_size=bin_size))[self.StartAdditional:]
        peaks, d = find_peaks([v.n for v in values], height=max(values).n / 2., distance=self.BunchSpacing)
        peaks += self.StartAdditional
        fit_peaks = []
        for p, ph in zip(peaks, d['peak_heights']):
            low, hi = int(p - 4 / bin_size), int(p + 4 / bin_size)
            f = h.Fit('gaus', 'qs', '', *[h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(ph / 3, 1, low, hi), h.FindLastBinAbove(ph / 2, 1, low, hi)]])
            fit_peaks.append(f.Parameter(1))
        fit_peaks = array(fit_peaks) - (self.BunchSpacing / 2. if not center else 0)
        return concatenate([fit_peaks, [fit_peaks[-1] + self.BunchSpacing]])

    def find_n_additional(self, start_bunch=None, end_bunch=None, thresh=None, cut=None):
        times = self.get(cut=cut, thresh=thresh)
        start, end = self.get_bunch_range(start_bunch, end_bunch)
        for i in range(times.size):
            times[i] = times[i][(times[i] > start) & (times[i] < end)]
        return array([lst.size for lst in times], dtype='u2')

    def find_bunch_heights(self, n=1, thresh=None, cut=None, excl=.05):
        t, h = [self.get(False, thresh, False, cut, i) for i in range(2)]
        t0 = self.get_start_additional(n) * self.BinWidth
        for i in range(t.size):
            h[i] = h[i][(t[i] > t0) & (t[i] < t0 + self.BunchSpacing)]
        values = concatenate(h).astype('d')
        return values[int(excl * values.size):]

    def find_all(self, thresh=None, fit=False, redo=False):
        suf = '{:1.0f}'.format(self.Threshold) if thresh is None and not fit else '{:1.0f}_{}'.format(thresh, int(fit)) if thresh is not None else int(fit)
        hdf5_path = self.make_simple_hdf5_path(suf=suf, dut=self.Channel)
        if file_exists(hdf5_path) and not redo:
            f = h5py.File(hdf5_path, 'r')
            return f['times'], f['heights'], f['n_peaks']
        if redo and file_exists(hdf5_path):
            remove_file(hdf5_path)
        times, heights = [], []
        simplefilter('ignore', RankWarning)
        wave_forms, trigger_cells = self.WF.get_all(), self.WF.get_trigger_cells()
        self.Ana.PBar.start(trigger_cells.size)
        for i in range(trigger_cells.size):
            j, t, p = self.find(wave_forms[i], trigger_cells[i], thresh=thresh)
            times.append(self.fit_landau(i, j) if fit else t)
            heights.append(p)
            self.Ana.PBar.update()
        find_n_peaks = vectorize(size, otypes=['u2'])
        f = h5py.File(hdf5_path, 'w')
        f.create_dataset('times', data=concatenate(times).astype('f2'))
        f.create_dataset('heights', data=concatenate(heights).astype('f2'))
        f.create_dataset('n_peaks', data=find_n_peaks(heights))
        return f['times'], f['heights'], f['n_peaks']

    def find(self, values, trigger_cell, thresh=None):
        peaks = find_peaks(values, height=self.Threshold if thresh is None else thresh, distance=self.Distance, prominence=self.Prominence)
        return array([peaks[0], array([self.Ana.Waveform.get_calibrated_time(trigger_cell, value) for value in peaks[0]]), peaks[1]['peak_heights']])

    def find_from_event(self, event=1819, fit=False):
        self.Tree.GetBranch('wf{}'.format(self.Channel)).GetEntry(event)
        self.Tree.GetBranch('trigger_cell').GetEntry(event)
        values = self.Polarity * array(getattr(self.Tree, 'wf{}'.format(self.Channel)))
        return self.find(values, self.Tree.trigger_cell, fit)
    # endregion FIND
    # ----------------------------------------

    # ----------------------------------------
    # region cft
    def find_cft(self, values=None, times=None, peaks=None, peak_times=None, ind=None, thresh=.5, show=False):
        x, y, p, t = (times, values, peaks, peak_times) if ind is None else self.get_event(ind)
        cfts = []
        for k, (ip, it) in enumerate(zip(p, t)):
            i = where(x == it)[0][0]
            v = y[max(0, i - 20): i]
            j = argmax(v > ip * thresh)  # find next index greater than threshold
            cft = interpolate_x(x[i + j - 21], x[i + j - 20], v[j - 1], v[j], ip * thresh)
            cfts.append(cft)
            if show:
                # self.WF.draw_single(ind=ind, x_range=self.Ana.get_signal_range(), draw_opt='alp') if not k else do_nothing()
                self.WF.draw_single(ind=ind, draw_opt='alp') if not k else do_nothing()
                Draw.horizontal_line(ip * thresh, 0, 2000, color=4)
                Draw.vertical_line(cft, -1000, 1000, color=2)
        return cfts

    def find_cft0(self, values=None, times=None, fac=.5, delay=None, show=False, i=None):
        delay = choose(delay, int(self.WF.get_average_rise_time() / self.BinWidth))  # [number of bins]
        s, e = self.Ana.SignalRegion
        v = (values if i is None else self.WF.get_all()[i])[s:e]
        t = (times if i is None else self.WF.get_all_times()[i])[s:e]
        v -= roll(v, -delay) * fac
        ineg = argmax(v < 0)  # find index with the  first negative value
        i = argmax(v[ineg:] > 0) + ineg  # find the next index greater than 0 and ignore the first positive ones
        x0 = interpolate_x(t[i - 1], t[i], v[i - 1], v[i], 0)
        if show:
            g = Draw.make_tgrapherrors('g', 'g', x=t, y=v)
            self.Draw(g, draw_opt='apl')
            Draw.horizontal_line(0, 0, 1000)
            Draw.vertical_line(x0, -500, 500)
            update_canvas()
        return x0

    def find_all_cft0(self, fac=.5, delay=None, redo=False):
        def f():
            values, times = self.WF.get_all(), self.WF.get_all_times()
            self.PBar.start(values.shape[0])
            cfts = []
            for i in range(values.shape[0]):
                cfts.append(self.find_cft0(values[i], times[i], fac, delay))
                self.PBar.update()
            return cfts
        return do_hdf5(self.make_simple_hdf5_path('cft0', '{:.0f}'.format(fac * 100)), f, redo=redo)

    def find_all_cft(self, thresh=.5, redo=False):
        def f():
            self.info('calculating constant fraction discrimination times ...')
            values, times, peaks, peak_times = self.WF.get_all(), self.WF.get_all_times(), self.get_heights(), self.get()
            self.PBar.start(values.shape[0])
            cfts = []
            for i in range(values.shape[0]):
                cfts.append(self.find_cft(values[i], times[i], peaks[i], peak_times[i], thresh=thresh))
                self.PBar.update()
            return concatenate(cfts).astype('f2')
        return do_hdf5(self.make_simple_hdf5_path('cft', '{:.0f}{}'.format(thresh * 100, self.Cut.get_name())), f, redo=redo)

    def draw_cft(self, thresh=.5, bin_size=.5, show=True, draw_ph=False, x=None, y=None, x_range=None, y_range=None, smear=None):
        times = choose(x, self.get_all_cft, thresh=thresh)
        self.smear_times(times, smear)
        title = '{:.0f}% Constrant Fraction Times'.format(thresh * 100)
        h = self.Draw.distribution(times, self.Ana.get_t_bins(bin_size), title, lm=.13, rm=.12 if draw_ph else None, show=show, x_tit='Constant Fraction Time [ns]', y_off=1.8)
        c = self.draw_ph(get_last_canvas(), bin_size, times, y, x_range, y_range, show=draw_ph)
        format_statbox(h, all_stat=True, c=c)
        return h

    def draw_height_vs_cft(self, bin_size=None, show=True):
        p = TProfile('phcft', 'Peak Height vs. Constant Fraction Time', *self.Ana.get_t_bins(bin_size))
        x, y = array(self.find_all_cft()), array(self.get_heights(flat=True))
        fill_hist(p, x=x, y=y)
        format_histo(p, x_tit='Constrant Fraction Time [ns]', y_tit='Peak Height [mV]', y_off=1.4)
        self.Draw(p, show, lm=.12)

    def draw_cft_vs_time(self, bin_size=.2, signal=False, show=True):
        h = TH2F('hcftt', 'Constant Fraction vs. Peak Time', *(self.Ana.get_t_bins(bin_size) + self.Ana.get_t_bins(bin_size)))
        x = array(self.get_from_tree()) if signal else self.get(flat=True)
        y = self.get_all_cft() if signal else array(self.find_all_cft()).astype('d')
        h.FillN(x.size, x.astype('d'), y.astype('d'), ones(x.size))
        format_histo(h, y_tit='Constrant Fraction Time [ns]', x_tit='Peak Time [ns]', y_off=1.3)
        self.Draw(h, show=show, lm=.11, draw_opt='colz', rm=.12, stats=set_statbox(entries=True, w=.2))
    # endregion cft
    # ----------------------------------------

    # ----------------------------------------
    # region TOT
    def calc_all_tot(self, thresh=None, fixed=True, redo=False):
        def f():
            self.info('calculating time over threshold ...')
            values, times, peaks, peak_times = self.WF.get_all(), self.WF.get_all_times(), self.get_heights(), self.get()
            self.PBar.start(values.shape[0])
            tots = []
            for i in range(values.shape[0]):
                tots.append(self.calc_tot(values[i], times[i], peaks[i], peak_times[i], thresh=thresh, fixed=fixed))
                self.PBar.update()
            return concatenate(tots).astype('f2')
        suffix = '' if thresh is None else '{:.0f}'.format(thresh if fixed else thresh * 100)
        return do_hdf5(self.make_simple_hdf5_path('TOT', suffix + self.Cut.get_name()), f, redo=redo)

    def calc_tot(self, values=None, times=None, peaks=None, peak_times=None, ind=None, thresh=None, fixed=True, show=False):
        x, y, p, t = (times, values, peaks, peak_times) if ind is None else self.get_event(ind)
        tot = []
        for j, (ip, it) in enumerate(zip(p, t)):
            thresh = ip * thresh if not fixed else self.Threshold * .75 if thresh is None else thresh
            i = where(x == it)[0][0]
            vl, vr = y[max(0, i - 20):i], y[i:i + 40]  # get left and right side of the peak
            l, r = argmax(vl > thresh), argmax(vr < thresh)  # find indices crossing the threshold
            tl = interpolate_x(x[i + l - 21], x[i + l - 20], vl[l - 1], vl[l], thresh)
            tr = interpolate_x(x[i + r - 1], x[i + r], vr[r - 1], vr[r], thresh)
            v = tr - tl
            tot.append(v if v < 1000 else -1)
            if show:
                self.WF.draw_single(ind=ind) if not j else do_nothing()
                Draw.horizontal_line(thresh, 0, 2000, color=4)
                Draw.vertical_line(tl, -1000, 1000, color=2)
                Draw.vertical_line(tr, -1000, 1000, color=2)
        return tot

    def draw_tot(self, thresh=None, fixed=True, show=True):
        values = array(self.get_all_tot(thresh, fixed))
        m, s = mean_sigma(sorted(values[(values > 0) & (values < 1e5)])[100:-100])
        thresh = '{:.0f}{}'.format(self.Threshold * .75 if thresh is None else thresh if not fixed else 100 * thresh, '%' if not fixed else 'mV')
        h = TH1F('htot', 'Time over {} threshold'.format(thresh), 200, *ax_range(m - 3 * s, m + 3 * s, .3, .3))
        h.FillN(values.size, values.astype('d'), ones(values.size))
        format_histo(h, x_tit='ToT [ns]', y_tit='Number of Entries', y_off=1.5, fill_color=Draw.FillColor)
        self.Draw(h, show=show, lm=.13)
    # endregion TOT
    # ----------------------------------------

    # ----------------------------------------
    # region MODEL
    def find_scale(self, n=20, redo=False):
        def f():
            times, heights = self.get_signal_times()[:n], self.get_signal_heights()[:n]
            c = []
            for i in range(n):
                g = self.WF.draw_single(ind=i, show=False)
                fit = g.Fit('landau', 'qs0', '', times[i] - 4, times[i] + 5)
                c.append(fit2u(FitRes(fit), par=0))
            g = Draw.make_tgrapherrors('gsm', 'Model Scale', x=heights, y=c)
            fit = g.Fit('pol1', 'qs0')
            return fit.Parameter(1)
        return do_pickle(self.make_simple_pickle_path('ModelScale'), f, redo=redo)

    @staticmethod
    def _signal(x, height, peak_time, rise_time, rise_fac):
        x0, y0 = peak_time, height
        x1, x2 = x0 - rise_time, x0 + rise_fac * rise_time
        p1 = get_p1(x0, x1 if x < x0 else x2, y0, 0)
        return p1 * x + get_p0(x0, y0, p1) if x1 <= x <= x2 else 0

    def signal1(self, height, peak_time, scale, landau_width=3, noise=4):
        x0, x1 = ax_range(peak_time - 2 * landau_width, peak_time + 4 * landau_width, .5, .5)
        x = arange(x0, x1, self.BinWidth, dtype='d')
        y = array([height * scale * TMath.Landau(ix, peak_time, landau_width) for ix in x])
        return x, y + normal(scale=noise, size=x.size)

    def signal0(self, height, peak_time, rise_time=None, rise_fac=3, noise=None):
        noise = self.Ana.Pedestal.get_raw_noise().n if noise is None else noise
        rise_time = self.WF.get_average_rise_time() if rise_time is None else rise_time
        x0, x1 = ax_range(peak_time - rise_time, peak_time + rise_time * rise_fac, .5, .5)
        x = arange(x0, x1, self.BinWidth, dtype='d')
        y = array([self._signal(ix, height, peak_time, rise_time, rise_fac) for ix in x])
        return x, y + normal(scale=noise, size=x.size)

    def get_signal(self, n, *args, **kwargs):
        return self.signal0(*args, **kwargs) if not n else self.signal1(*args, **kwargs)

    def draw_model_signal(self, model=0, *args, **kwargs):
        x, y = self.get_signal(model, *args, **kwargs)
        g = Draw.make_tgrapherrors('gms', 'Model Signal', x=x, y=y)
        format_histo(g, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, markersize=.5)
        self.Draw(g, lm=.073, rm=.045, bm=.18, w=1.5, h=.5, grid=1)

    def draw_model_signal1(self, height=None, peak_time=None, landau_width=3):
        scale = self.find_scale()
        height = self.get_mpv_sigma_heights()[0] if height is None else height
        peak_time = self.get_mean_sigma()[0] if peak_time is None else peak_time
        self.draw_model_signal(1, height, peak_time, scale, landau_width, self.Ana.Pedestal.get_raw_noise().n)

    def draw_raw_model(self, n=1e6, model=1, show=True, *args, **kwargs):
        times = []
        self.PBar.start(n)
        for _ in range(int(n)):
            x, y = self.get_signal(model, *args, **kwargs)
            times.append(x[argmax(y)])
            self.PBar.update()
        times = array(times)
        h = TH1F('hrm', 'Raw Model Peak Times', *self.Ana.get_t_bins())
        h.FillN(times.size, times, ones(times.size))
        format_histo(h, x_tit='Signal Peak Time [ns]', y_tit='Number of Entries', y_off=1.8, fill_color=Draw.FillColor)
        self.Draw(h, lm=.13, show=show, stats=set_statbox(entries=True, w=.2))

    def model1(self, n=1e6, redo=False, scale=None, landau_width=3, cft=False):
        scale = self.find_scale() if scale is None else scale
        return self.model(n, model=1, cft=cft, redo=redo, scale=scale, landau_width=landau_width)

    def model0(self, n=1e6, redo=False, rise_time=None, rise_fac=3, sigma=None):
        return self.model(n, 0, redo, )

    def model(self, n=1e6, model=1, noise=None, cft=False, redo=False, *args, **kwargs):
        n = int(n)
        hdf5_path = self.make_simple_hdf5_path('M', '{}_{}_{}'.format(n, model, int(cft)))
        if file_exists(hdf5_path) and not redo:
            f = h5py.File(hdf5_path, 'r')
            return f['times'], f['heights']
        if redo and file_exists(hdf5_path):
            remove_file(hdf5_path)
        Landau(5, 2)  # first value is always crap
        noise = self.Ana.Pedestal.get_raw_noise().n if noise is None else noise
        m, s = self.get_mean_sigma()
        mpv, sl = self.get_mpv_sigma_heights()
        t, v = [], []
        self.PBar.start(n)
        peak_times = normal(scale=s, size=n) + m
        for i in range(int(n)):
            x, y = self.get_signal(model, min(500, Landau(mpv, sl)), peak_times[i], noise=noise, *args, **kwargs)
            peak_time = x[argmax(y)]
            peak_height = max(y)
            t.append(peak_time if not cft else self.find_cft(y, x, [peak_height], [peak_time]))
            v.append(peak_height)
            self.PBar.update()
        t, v = array(t), array(v)
        f = h5py.File(hdf5_path, 'w')
        f.create_dataset('times', data=array(t).astype('f2'))
        f.create_dataset('heights', data=array(v).astype('f4'))
        return f['times'], f['heights']

    def draw_model(self, n=1e6, model=1, cft=False, draw_ph=False, show=True):
        x, y = self.model1(n, cft=cft) if model == 1 else self.model0(n)
        self.draw_signal(x=array(x), y=array(y), draw_ph=draw_ph, show=show)

    # endregion MODEL
    # ----------------------------------------

    def fit_landau(self, i, peak_indices):
        values, times = self.WF.get_all()[i], self.WF.get_all_times()[i]
        t = []
        f = TF1('f', 'landau', 0, 512)
        for ip in peak_indices.astype('i2'):
            g = Draw.make_tgrapherrors('g', 'g', x=times[max(0, ip - 6):ip + 8], y=values[max(0, ip - 6):ip + 8])
            g.Fit(f, 'q0')
            t.append(f.GetMaximumX())
            g.Delete()
        f.Delete()
        return t

    def fit(self, values, peaks, trigger_cell):
        peaks = peaks[(peaks > 10) & (peaks < 1014)]
        t = self.Ana.Waveform.get_calibrated_times(trigger_cell)
        time_peaks = []
        for peak in peaks:
            y = savgol_filter(values[peak - 10:peak + 12], 15, 3)
            p_max = where(y == y.max())[0][0]
            if p_max < 5 or p_max > 15:
                continue
            t0 = t[peak - 10:peak + 12]
            p2, p1, p0 = polyfit(t0[p_max - 1:p_max + 2], y[p_max - 1:p_max + 2], deg=2)
            time_peaks.append(-p1 / 2 / p2)
        return time_peaks

    def compare_signal_distributions(self, bins, bin_size=None, x_range=None):
        histos = [self.Ana.draw_signal_distribution(show=False, cut=self.get_tbin_cut(ibin, bin_size) + self.Ana.Cut()) for ibin in bins]
        stack = THStack('ssdt', 'Time Comparison;Time [ns];Number of Entries')
        leg = Draw.make_legend(nentries=2, w=.25)
        l_names = ['low pulse height', 'high pulse height']
        histos.reverse()
        for i, h in enumerate(histos):
            color = self.Draw.get_color(2)
            h.Scale(1 / h.GetMaximum())
            format_histo(h, stats=0, color=color, fill_color=color, opacity=.6)
            leg.AddEntry(h, l_names[i], 'l')
            stack.Add(h)
        format_histo(stack, draw_first=True, x_range=x_range, y_off=1.4)
        self.Draw(stack, 'TimingComparison', draw_opt='nostack', leg=leg, lm=.12)

    def compare_hit_maps(self, bins, res=2, show=True):
        h0, h1 = [self.Ana.draw_hitmap(cut=self.get_tbin_cut(ibin) + self.Cut(), res=res, show=False) for ibin in bins]
        h2 = TH2F()
        h0.Copy(h2)
        h2.Divide(h1)
        self.Draw(h2, show=show, draw_opt='colz', rm=.15)
        self.Ana.draw_fid()
    
    @staticmethod
    def smear_times(times, width=2.5, n=5, gaus=False):
        if width is None:
            return
        if gaus:  # only Gaussian smear
            times += normal(0, width, times.size) if width else 0  # gaussian width
        else:  # flat plus Gaussian smear
            times += rand(times.size) * width - width / 2 if width else 0
            times[::n] += normal(0, width / 2, times.size // 5 + 1)[:times[::n].size] if width else 0
