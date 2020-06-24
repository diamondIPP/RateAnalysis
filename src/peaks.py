#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import *
from ROOT import TH1F, TCut, TProfile, THStack, TH2F, TMath
from ROOT.gRandom import Landau
from scipy.signal import find_peaks, savgol_filter
from numpy import polyfit, pi, RankWarning, vectorize, size, split, ones, ceil, repeat, linspace, argmax, insert, roll
from numpy.random import normal, rand
from warnings import simplefilter
from InfoLegend import InfoLegend


class PeakAnalysis(Analysis):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.Channel
        self.DUT = self.Ana.DUT
        self.Tree = self.Ana.Tree
        self.WF = self.Ana.Waveform
        self.NoiseThreshold = self.calc_noise_threshold()
        self.Threshold = max(self.NoiseThreshold, self.Ana.get_min_signal(self.Ana.get_signal_name(peak_integral=1)))
        self.Cut = self.Ana.Cut()
        self.InfoLegend = InfoLegend(pad_analysis)
        self.StartAdditional = self.get_start_additional()
        self.NBunches = self.calc_n_bunches()
        self.BinWidth = self.Ana.DigitiserBinWidth
        self.BunchSpacing = self.Ana.BunchSpacing
        self.set_pickle_sub_dir('Peaks')

    # ----------------------------------------
    # region INIT
    def calc_noise_threshold(self):
        """ return peak threshold, 5 times the raw noise + mean of the noise. """
        return abs(self.Ana.Pedestal.get_raw_mean() + 5 * self.Ana.Pedestal.get_raw_noise()).n

    def get_start_additional(self, bunch=2):
        """Set the start of the additional peaks 2.5 bunches after the signal peak to avoid the biased bunches after the signal."""
        return int(self.Run.IntegralRegions[self.DUT.Number - 1]['signal_a'][0] + self.Ana.BunchSpacing * (bunch + 0.5) / self.Ana.DigitiserBinWidth)

    def calc_n_bunches(self):
        return int((self.Run.NSamples - self.StartAdditional) * self.Ana.DigitiserBinWidth / self.Ana.BunchSpacing)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_binning(self, bin_size=.5):
        return self.Ana.Waveform.get_binning(bin_size)

    def get_from_tree(self):
        return do_hdf5(self.make_simple_hdf5_path('Peaks', self.get_cut_name()), self.Run.get_root_vec, var=self.Ana.PeakName, cut=self.Cut, dtype='f2')

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

    def get(self, flat=False, fit=False):
        times, heights, n_peaks = self.find_all(fit=fit)
        return array(times) if flat else array(split(times, cumsum(n_peaks)[:-1]))

    def get_heights(self, flat=False):
        times, heights, n_peaks = self.find_all()
        return array(heights) if flat else array(split(heights, cumsum(n_peaks)[:-1]))

    def get_tbin_cut(self, ibin, bin_size=None, corr=True, fine_corr=False):
        cut_name = self.Ana.Timing.get_peak_name(corr, fine_corr)
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

    def get_corrected_times(self, times, n_peaks=None, events=None):
        signal_peak_times = repeat(self.get_from_tree(), n_peaks) if events is None else array(self.get_from_tree())[events]
        return times - (signal_peak_times - self.get_from_tree()[0])

    def get_n_times(self, n=2, ret_indices=False):
        """ :returns all times with exactly [n] additional peaks. """
        times = self.get()
        for i in range(times.size):
            times[i] = times[i][(times[i] > self.StartAdditional * self.BinWidth)]
        return where(array([lst.size for lst in times]) == n)[0] if ret_indices else times[array([lst.size for lst in times]) == n]  # select all events with two additional peaks

    def get_n(self, n=1):
        return self.get_n_times(n, ret_indices=True)

    def get_flux(self, n_peaks=None, lam=None, redo=False, prnt=True):
        def f():
            n = self.find_n_additional() if n_peaks is None and lam is None else n_peaks
            lambda_ = ufloat(mean(n), sqrt(mean(n) / n.size)) if lam is None else lam
            flux = lambda_ / (self.Ana.BunchSpacing * self.NBunches * self.get_area()) * 1e6
            return flux
        value = do_pickle(self.make_simple_pickle_path('Flux'), f, redo=redo)
        self.info('Estimated Flux by number of peaks: {}'.format(make_flux_string(value)), prnt=prnt)
        return value

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

    def get_cut_name(self):
        return self.Cut.GetName() if not self.Cut.GetName().startswith('All') else ''
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw(self, corr=True, scale=False, split_=1, thresh=None, y_range=None, show=True, redo=False):
        def f():
            times, heights, n_peaks = self.find_all(redo=redo, thresh=thresh)
            times = self.get_corrected_times(times, n_peaks) if corr else times
            times = split(times, split_)
            hs = [TH1F('hp{}{}'.format(self.Run.Number, i), 'Peak Times', 512 * 2, 0, 512) for i in range(split_)]
            for i in range(split_):
                v = times[i]
                hs[i].FillN(v.size, array(v, 'd'), ones(v.size))
            return hs[0] if split_ == 1 else hs
        suffix = '{}{}{}'.format(int(corr), '' if split_ == 1 else '_{}'.format(split_), '' if thresh is None else '_{:1.0f}'.format(thresh))
        h = do_pickle(self.make_simple_pickle_path('Histo', suffix), f, redo=redo)
        if scale:
            h.Sumw2()
            h.Scale(1e5 / self.Ana.Waveform.get_from_tree().shape[0])
        if show:
            self.format_statbox(entries=True)
            format_histo(h, x_tit='Time [ns]', y_tit='Number of Peaks', y_off=1.3, fill_color=self.FillColor, y_range=y_range)
            self.draw_histo(h, lm=.12, show=show, x=1.5, y=0.75, logy=True)
        return h

    def draw_signal(self, bin_size=.5, ind=None, fit=False, y=None, x=None, x_range=None, y_range=None, show=True, draw_ph=False, smear=None):
        self.format_statbox(all_stat=True, x=.86 if draw_ph else .95)
        times = choose(x, self.get_signal_times, fit=fit, ind=ind)
        self.smear_times(times, smear)
        h = self.draw_disto(times, 'Signal Peak Times', self.Ana.get_t_bins(bin_size), lm=.13, rm=.12 if draw_ph else None, show=show, x_tit='Signal Peak Time [ns]', y_off=1.8, x_range=x_range)
        self.draw_ph(get_last_canvas(), bin_size, times, y, x_range, y_range, draw_ph)
        return h

    def draw_ph(self, c, bin_size, x, y, x_range, y_range, show):
        if show:
            p = self.Ana.draw_signal_vs_peaktime(show=False, bin_size=bin_size, x=x, y=y)
            values = get_hist_vec(p, err=False)
            format_histo(p, title=' ', stats=0, x_tit='', l_off_x=1, y_range=choose(y_range, increased_range([min(values[values > 0]), max(values)], .3, .3)), x_range=x_range)
            c.cd()
            self.draw_tpad('psph', transparent=True, lm=.13, rm=.12)
            p.Draw('y+')
            update_canvas(c)

    def draw_heights_vs_time(self, bin_size=.5, corr=True, cft=False, show=True):
        times, heights, n_peaks = self.find_all()
        times = self.get_corrected_times(times, n_peaks) if corr else self.find_all_cft() if cft else times
        p = TProfile('pph', 'Peak Heights', *self.get_binning(bin_size))
        p.FillN(times.size, array(times).astype('d'), array(heights).astype('d'), ones(times.size))
        format_histo(p, x_tit='Time [ns]', y_tit='Peak Height [mV]', y_off=1.3, stats=0, fill_color=self.FillColor)
        self.draw_histo(p, lm=.12, show=show, x=1.5, y=0.75)

    def draw_signal_height_vs_time(self, bin_size=None, show=True):
        self.format_statbox(entries=True)
        p = TProfile('pspt', 'Peak Height vs. Constant Fraction Time', *self.Ana.get_t_bins(bin_size))
        fill_hist(p, x=self.get(flat=True), y=self.get_heights(flat=True))
        format_histo(p, x_tit='Constrant Fraction Time [ns]', y_tit='Peak Height [mV]', y_off=1.4)
        self.draw_histo(p, show, lm=.12)
        return p

    def draw_combined_heights(self, hist=False, show=True):
        s1_indices = self.get_n()
        times, heights = concatenate(self.get()[s1_indices]), concatenate(self.get_heights()[s1_indices])
        x = self.get_corrected_times(times[times > self.StartAdditional * self.BinWidth] % self.BunchSpacing, events=s1_indices)
        y = heights[times > self.StartAdditional * self.BinWidth]
        p = TProfile('pcph', 'Combined Bunch Pulse Heights', int(ceil(self.BunchSpacing)) * 2, 0, ceil(self.BunchSpacing))
        p = TH1F('hcph', 'Combined Number of Peaks for all Bunches', int(ceil(self.BunchSpacing)) * 2, 0, ceil(self.BunchSpacing)) if hist else p
        p.FillN(x.size, x.astype('d'), ones(x.size)) if hist else p.FillN(x.size, x.astype('d'), y.astype('d'), ones(x.size))
        format_histo(p, x_tit='Time [ns]', y_tit='Number of Peaks' if hist else 'Peak Height [mV]', y_off=1.6, stats=0, fill_color=self.FillColor)
        self.draw_histo(p, lm=.12, show=show)

    def draw_n(self, do_fit=False, show=True):
        n_peaks = self.find_n_additional()
        h = TH1F('h_pn', 'Number of Peaks', 10, 0, 10)
        h.FillN(n_peaks.size, n_peaks.astype('d'), ones(n_peaks.size))
        self.format_statbox(only_fit=True, w=.3) if do_fit else self.format_statbox(entries=True)
        if do_fit:
            fit_poissoni(h, show=show)
        format_histo(h, x_tit='Number of Peaks', y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor, lw=2)
        self.save_histo(h, 'PeakNumbers{}'.format('Fit' if do_fit else ''), show, logy=True, lm=.11)
        self.get_flux(n_peaks)
        return h

    def draw_n_bunch(self, b=0, y_range=None, show=True):
        bunches = self.find_bunches()
        times = self.get_n_times(n=2)
        times = times[[any((bunches[b] < lst) & (lst < bunches[b + 1])) for lst in times]]  # select events with a peak in bunch b
        times = concatenate(times)
        h = TH1F('h2a', '2 Additional Peaks for Bunch {}'.format(b), *self.get_binning())
        h.FillN(times.size, times.astype('d'), ones(times.size))
        self.format_statbox(entries=True)
        format_histo(h, x_tit='Time [ns]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor, y_range=y_range)
        self.draw_histo(h, lm=.12, show=show, x=1.5, y=0.75, logy=True)
        return times

    def draw_bunch_systematics(self, n=None, show=True):
        n = self.NBunches - 1 if n is None else n
        bunches = self.find_bunches()
        times = self.get_n_times(n=2)
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
        peaks_per_bunch = [ufloat(v, sqrt(v)) for v in peaks_per_bunch] / arange(self.NBunches, self.NBunches - n - 1, -1)  # add errors and normalise
        n_peaks = peaks_per_bunch[1:]  # exclude the signal peak
        g = self.make_tgrapherrors('g', 'Bunch Systematics', x=arange(1, n_peaks.size + 1), y=n_peaks)
        format_histo(g, x_tit='Bunch after a Signal', y_tit='Average Number of Peaks', y_off=1.3)
        self.draw_histo(g, lm=.12, show=show)

    def draw_height_disto(self, show=True):
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.get_signal_heights(), 'Signal Heights', self.Ana.Bins.get_pad_ph(), lm=.12, show=show, x_tit='Peak Height [mV]', y_off=1.5)

    def draw_flux_vs_threshold(self, steps=20):
        x = linspace(self.NoiseThreshold, self.Threshold, steps)
        y = array([self.get_flux(lam=self.get_n_additional(thresh=ix), redo=1) for ix in x]) / 1000.
        g = self.make_tgrapherrors('gft', 'Flux vs. Peak Threshold', x=x, y=y)
        format_histo(g, x_tit='Peak Finding Threshold [mV]', y_tit='Flux [MHz/cm^{2}]', y_off=1.3)
        self.draw_histo(g, draw_opt='ap', lm=.12)

    def draw_peak_spacing(self, overlay=False, bin_size=.5, show=True):
        bf = self.BunchSpacing
        values = concatenate(self.get() - self.get_signal_times())
        values = values[values != 0]
        values = ((values + bf / 2) % bf + bf / 2) if overlay else values
        m, w = mean(values), self.Ana.get_t_bins()[1][-1] - self.Ana.get_t_bins()[1][0]
        bins = arange(m - w / 2, m + w / 2, bin_size) if overlay else arange(0, self.Run.NSamples * self.BinWidth, bin_size)
        h = TH1F('hbs', 'Peak Spacing', bins.size - 1, bins)
        h.FillN(values.size, values.astype('d'), ones(values.size))
        self.format_statbox(entries=True)
        format_histo(h, x_tit='Peak Distance [ns]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor)
        x, y = [None, None] if overlay else [1.5, .75]
        self.draw_histo(h, show=show, lm=.11, x=x, y=y)
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region FIND
    def find_additional(self, h=None, scale=False, show=True):
        values = get_hist_vec(self.draw(scale=scale, show=False) if h is None else h)[self.StartAdditional:]
        peaks = find_peaks([v.n for v in values], height=max(values).n / 2., distance=self.Ana.BunchSpacing)
        g = self.make_tgrapherrors('ga', 'Additional Peak Heights', x=(peaks[0] + self.StartAdditional) / 2., y=values[peaks[0]])
        self.format_statbox(fit=True)
        g.Fit('pol0', 'qs')
        self.draw_histo(g, show=show, x=1.5, y=.75, gridy=1)
        return mean(values[peaks[0]])

    def find_bunches(self, center=False):
        values = get_hist_vec(self.draw(show=False))[self.StartAdditional:]
        bunches = (find_peaks([v.n for v in values], height=max(values).n / 2., distance=self.Ana.BunchSpacing)[0] + self.StartAdditional) * self.BinWidth
        bunches -= self.BunchSpacing / 2. if not center else 0
        return concatenate([bunches, [bunches[-1] + self.BunchSpacing]])

    def find_n_additional(self, start_bunch=None, end_bunch=None, thresh=None):
        times, heights, n_peaks = self.find_all(thresh=thresh)
        times = array(split(times, cumsum(n_peaks)[:-1]))
        start = (self.StartAdditional if start_bunch is None else self.get_start_additional(start_bunch)) * self.Ana.DigitiserBinWidth
        end = (self.Run.NSamples if end_bunch is None else self.get_start_additional(end_bunch)) * self.Ana.DigitiserBinWidth
        for i in range(times.size):
            times[i] = times[i][(times[i] > start) & (times[i] < end)]
        return array([lst.size for lst in times], dtype='u2')

    def draw_additional_disto(self, show=True):
        hs = self.draw(split_=4)
        h0 = TH1F('ht', 'Peak Heights', 20, 170, 260)
        for h in hs:
            for v in self.find_additional(h):
                h0.Fill(v.n)
        format_histo(h0, x_tit='Peak Height', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.draw_histo(h0, show=show)

    def find_all(self, redo=False, thresh=None, fit=False):
        suf = '' if thresh is None and not fit else '{:1.0f}_{}'.format(thresh, int(fit)) if thresh is not None else int(fit)
        hdf5_path = self.make_simple_hdf5_path(suf=suf + self.get_cut_name(), dut=self.Channel)
        if file_exists(hdf5_path) and not redo:
            f = h5py.File(hdf5_path, 'r')
            return f['times'], f['heights'], f['n_peaks']
        if redo and file_exists(hdf5_path):
            remove_file(hdf5_path)
        times, heights = [], []
        simplefilter('ignore', RankWarning)
        wave_forms, trigger_cells = self.WF.get_all(), self.WF.get_trigger_cells()
        self.Ana.PBar.start(trigger_cells.size)
        for i in xrange(trigger_cells.size):
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
        peaks = find_peaks(values, height=self.Threshold if thresh is None else thresh, distance=10, prominence=20)
        return array([peaks[0], array([self.Ana.Waveform.get_calibrated_time(trigger_cell, value) for value in peaks[0]]), peaks[1]['peak_heights']])

    def find_from_event(self, event=1819, fit=False):
        self.Tree.GetBranch('wf{}'.format(self.Channel)).GetEntry(event)
        self.Tree.GetBranch('trigger_cell').GetEntry(event)
        values = self.Ana.Polarity * array(getattr(self.Tree, 'wf{}'.format(self.Channel)))
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
                self.draw_horizontal_line(ip * thresh, 0, 2000, name='thresh{}'.format(k), color=4)
                self.draw_vertical_line(cft, -1000, 1000, name='cft{}'.format(k), color=2)
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
            g = self.make_tgrapherrors('g', 'g', x=t, y=v)
            self.draw_histo(g, draw_opt='apl')
            self.draw_horizontal_line(0, 0, 1000, name='h')
            self.draw_vertical_line(x0, -500, 500, name='v')
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
        return do_hdf5(self.make_simple_hdf5_path('cft', '{:.0f}{}'.format(thresh * 100, self.get_cut_name())), f, redo=redo)

    def draw_cft(self, thresh=.5, bin_size=.5, show=True, draw_ph=False, x=None, y=None, x_range=None, y_range=None, smear=None):
        self.format_statbox(all_stat=True, x=.86 if draw_ph else .95)
        times = choose(x, self.get_all_cft, thresh=thresh)
        self.smear_times(times, smear)
        title = '{:.0f}% Constrant Fraction Times'.format(thresh * 100)
        h = self.draw_disto(times, title, self.Ana.get_t_bins(bin_size), lm=.13, rm=.12 if draw_ph else None, show=show, x_tit='Constant Fraction Time [ns]', y_off=1.8)
        self.draw_ph(get_last_canvas(), bin_size, times, y, x_range, y_range, show=draw_ph)
        return h

    def draw_height_vs_cft(self, bin_size=None, show=True):
        p = TProfile('phcft', 'Peak Height vs. Constant Fraction Time', *self.Ana.get_t_bins(bin_size))
        x, y = array(self.find_all_cft()), array(self.get_heights(flat=True))
        fill_hist(p, x=x, y=y)
        format_histo(p, x_tit='Constrant Fraction Time [ns]', y_tit='Peak Height [mV]', y_off=1.4)
        self.draw_histo(p, show, lm=.12)

    def draw_cft_vs_time(self, bin_size=.2, signal=False, show=True):
        h = TH2F('hcftt', 'Constant Fraction vs. Peak Time', *(self.Ana.get_t_bins(bin_size) + self.Ana.get_t_bins(bin_size)))
        x = array(self.get_from_tree()) if signal else self.get(flat=True)
        y = self.get_all_cft() if signal else array(self.find_all_cft()).astype('d')
        h.FillN(x.size, x.astype('d'), y.astype('d'), ones(x.size))
        format_histo(h, y_tit='Constrant Fraction Time [ns]', x_tit='Peak Time [ns]', y_off=1.3)
        self.format_statbox(entries=True, x=.86)
        self.draw_histo(h, show=show, lm=.11, draw_opt='colz', rm=.12)
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
        return do_hdf5(self.make_simple_hdf5_path('TOT', suffix + self.get_cut_name()), f, redo=redo)

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
                self.draw_horizontal_line(thresh, 0, 2000, name='thresh', color=4)
                self.draw_vertical_line(tl, -1000, 1000, name='l{}'.format(j), color=2)
                self.draw_vertical_line(tr, -1000, 1000, name='r{}'.format(j), color=2)
        return tot

    def draw_tot(self, thresh=None, fixed=True, show=True):
        values = array(self.get_all_tot(thresh, fixed))
        m, s = mean_sigma(sorted(values[(values > 0) & (values < 1e5)])[100:-100])
        thresh = '{:.0f}{}'.format(self.Threshold * .75 if thresh is None else thresh if not fixed else 100 * thresh, '%' if not fixed else 'mV')
        h = TH1F('htot', 'Time over {} threshold'.format(thresh), 200, *increased_range([m - 3 * s, m + 3 * s], .3, .3))
        h.FillN(values.size, values.astype('d'), ones(values.size))
        format_histo(h, x_tit='ToT [ns]', y_tit='Number of Entries', y_off=1.5, fill_color=self.FillColor)
        self.draw_histo(h, show=show, lm=.13)
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
                c.append(make_ufloat(FitRes(fit), par=0))
            g = self.make_tgrapherrors('gsm', 'Model Scale', x=heights, y=c)
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
        x0, x1 = increased_range([peak_time - 2 * landau_width, peak_time + 4 * landau_width], .5, .5)
        x = arange(x0, x1, self.BinWidth, dtype='d')
        y = array([height * scale * TMath.Landau(ix, peak_time, landau_width) for ix in x])
        return x, y + normal(scale=noise, size=x.size)

    def signal0(self, height, peak_time, rise_time=None, rise_fac=3, noise=None):
        noise = self.Ana.Pedestal.get_raw_noise().n if noise is None else noise
        rise_time = self.WF.get_average_rise_time() if rise_time is None else rise_time
        x0, x1 = increased_range([peak_time - rise_time, peak_time + rise_time * rise_fac], .5, .5)
        x = arange(x0, x1, self.BinWidth, dtype='d')
        y = array([self._signal(ix, height, peak_time, rise_time, rise_fac) for ix in x])
        return x, y + normal(scale=noise, size=x.size)

    def get_signal(self, n, *args, **kwargs):
        return self.signal0(*args, **kwargs) if not n else self.signal1(*args, **kwargs)

    def draw_model_signal(self, model=0, *args, **kwargs):
        x, y = self.get_signal(model, *args, **kwargs)
        g = self.make_tgrapherrors('gms', 'Model Signal', x=x, y=y)
        format_histo(g, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, markersize=.5)
        self.draw_histo(g, lm=.073, rm=.045, bm=.18, x=1.5, y=.5, grid=1)

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
        format_histo(h, x_tit='Signal Peak Time [ns]', y_tit='Number of Entries', y_off=1.8, fill_color=self.FillColor)
        self.format_statbox(entries=1)
        self.draw_histo(h, lm=.13, show=show)

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
            g = self.make_tgrapherrors('g', 'g', x=times[max(0, ip - 6):ip + 8], y=values[max(0, ip - 6):ip + 8])
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
        leg = self.make_legend(nentries=2, w=.25)
        l_names = ['low pulse height', 'high pulse height']
        histos.reverse()
        for i, h in enumerate(histos):
            color = get_color_gradient(2)[i]
            h.Scale(1 / h.GetMaximum())
            format_histo(h, stats=0, color=color, fill_color=color, opacity=.6)
            leg.AddEntry(h, l_names[i], 'l')
            stack.Add(h)
        format_histo(stack, draw_first=True, x_range=x_range, y_off=1.4)
        self.draw_histo(stack, 'TimingComparison', draw_opt='nostack', leg=leg, lm=.12)

    def compare_hit_maps(self, bins, res=2, show=True):
        h0, h1 = [self.Ana.draw_hitmap(cut=self.get_tbin_cut(ibin) + self.Cut, res=res, show=False) for ibin in bins]
        h2 = TH2F()
        h0.Copy(h2)
        h2.Divide(h1)
        self.draw_histo(h2, show=show, draw_opt='colz', rm=.15)
        self.Ana.draw_fid_cut()
    
    @staticmethod
    def smear_times(times, width=2.5, n=5, gaus=False):
        if width is None:
            return
        if gaus:  # only Gaussian smear
            times += normal(0, width, times.size) if width else 0  # gaussian width
        else:  # flat plus Gaussian smear
            times += rand(times.size) * width - width / 2 if width else 0
            times[::n] += normal(0, width / 2, times.size // 5 + 1)[:times[::n].size] if width else 0
