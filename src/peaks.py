#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import *
from ROOT import TH1F, TCut
from scipy.signal import find_peaks, savgol_filter
from numpy import polyfit, pi, RankWarning, vectorize, size, split, ones
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
        self.NoiseThreshold = self.calc_threshold()
        self.Cut = self.Ana.Cut()
        self.InfoLegend = InfoLegend(pad_analysis)
        self.StartAdditional = self.get_start_additional()
        self.NBunches = self.calc_n_bunches()

    def calc_threshold(self):
        """ return peak threshold, 6 times the raw noise or the minimum singnal, whatever is higher. """
        return max(abs(self.Ana.Pedestal.get_raw_mean() + 6 * self.Ana.Pedestal.get_raw_noise()), self.Ana.get_min_signal())

    def get_start_additional(self):
        """Set the start of the additional peaks 2.5 bunches after the signal peak to avoid the biased bunches after the signal."""
        return int(self.Run.IntegralRegions[self.DUT.Number - 1]['signal_a'][0] + self.Ana.BunchSpacing * 2.5 / self.Ana.DigitiserBinWidth)

    def calc_n_bunches(self):
        return int((self.Run.NSamples - self.StartAdditional) * self.Ana.DigitiserBinWidth / self.Ana.BunchSpacing)

    def get_all(self):
        return do_hdf5(self.make_hdf5_path('Peaks', 'V1', self.Ana.RunNumber, self.Channel), self.Run.get_root_vec, var=self.Ana.PeakName, cut=self.Cut, dtype='f2')

    def draw(self, corr=True, scale=False, fit=False, split_=1, y_range=None, show=True, redo=False):
        def f():
            values, n_peaks = self.find_all(redo=redo, fit=fit)
            if corr:
                values = array(split(values, cumsum(n_peaks)[:-1]))
                peaks = self.get_all()
                for i in xrange(values.size):
                    values[i] -= peaks[i] - peaks[0]
                values = concatenate(values)
            values = split(values, split_)
            hs = [TH1F('hp{}{}'.format(self.Ana.RunNumber, i), 'Peak Times', 512 * 2, 0, 512) for i in range(split_)]
            for i in range(split_):
                v = values[i]
                hs[i].FillN(v.size, array(v, 'd'), ones(v.size))
            return hs[0] if split_ == 1 else hs
        h = do_pickle(self.make_pickle_path('Peaks', 'Histo', self.Ana.RunNumber, self.Channel, suf=None if split_ == 1 else split_), f, redo=redo)
        if scale:
            h.Scale(1e5 / self.Ana.Waveform.get_all().shape[0])
        if show:
            self.format_statbox(entries=True)
            format_histo(h, x_tit='Time [ns]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor, y_range=y_range)
            self.draw_histo(h, lm=.12, show=show, x=1.5, y=0.75, logy=True)
        return h

    def find_additional(self, h=None, scale=False, show=True):
        values = get_hist_vec(self.draw(scale=scale, show=False) if h is None else h)[self.StartAdditional:]
        peaks = find_peaks([v.n for v in values], height=max(values).n / 2., distance=self.Ana.BunchSpacing)
        g = self.make_tgrapherrors('ga', 'Additional Peak Heights', x=(peaks[0] + self.StartAdditional) / 2., y=values[peaks[0]])
        self.format_statbox(fit=True)
        g.Fit('pol0', 'qs')
        self.draw_histo(g, show=show, x=1.5, y=.75, gridy=1)
        return mean(values[peaks[0]])

    def find_n_additional(self):
        values, n_peaks = self.find_all()
        values = array(split(values, cumsum(n_peaks)[:-1]))
        for i in range(values.size):
            values[i] = values[i][values[i] > self.StartAdditional * self.Ana.DigitiserBinWidth]
        return array([lst.size for lst in values], dtype='u2')

    def draw_additional_disto(self, show=True):
        hs = self.draw(split_=4)
        h0 = TH1F('ht', 'Peak Heights', 20, 170, 260)
        for h in hs:
            for v in self.find_additional(h):
                h0.Fill(v.n)
        format_histo(h0, x_tit='Peak Height', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.draw_histo(h0, show=show)

    def find_all(self, redo=False, fit=False):
        hdf5_path = self.make_hdf5_path('Peaks', run=self.Ana.RunNumber, ch=self.Channel, suf=int(fit))
        if file_exists(hdf5_path) and not redo:
            f = h5py.File(hdf5_path, 'r')
            return f['values'], f['n_peaks']
        peak_times = []
        simplefilter('ignore', RankWarning)
        wfs = self.Ana.Waveform.get_all()
        self.Ana.PBar.start(wfs.shape[0])
        for wf, tc in zip(wfs, self.Ana.Waveform.get_trigger_cells()):
            peak_times.append(self.find(wf, tc, fit=fit))
            self.Ana.PBar.update()
        find_n_peaks = vectorize(size, otypes=['u2'])
        values, n_peaks = concatenate(peak_times), find_n_peaks(peak_times)
        f = h5py.File(hdf5_path, 'w')
        f.create_dataset('values', data=values)
        f.create_dataset('n_peaks', data=n_peaks)
        return f['values'], f['n_peaks']

    def find(self, values, trigger_cell, fit=False):
        peak_values = find_peaks(values, height=self.NoiseThreshold, distance=2, prominence=20)[0]
        if fit:
            return self.fit(values, peak_values, trigger_cell)
        return array([self.Ana.Waveform.get_calibrated_time(trigger_cell, value) for value in peak_values])

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

    def find_from_event(self, event=1819, fit=False):
        self.Tree.GetBranch('wf0').GetEntry(event)
        self.Tree.GetBranch('trigger_cell').GetEntry(event)
        values = self.Ana.Polarity * array(getattr(self.Tree, 'wf{}'.format(self.Channel)))
        return self.find(values, self.Tree.trigger_cell, fit)

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

    def get_flux(self, n_peaks=None):
        n_peaks = self.find_n_additional() if n_peaks is None else n_peaks
        lambda_ = ufloat(mean(n_peaks), sqrt(mean(n_peaks) / n_peaks.size))
        flux = lambda_ / (self.Ana.BunchSpacing * self.NBunches * self.get_area()) * 1e6
        info('Estimated Flux by number of peaks: {}'.format(make_flux_string(flux)))
        return flux

    def get_area(self, bcm=False):
        return self.get_bcm_area() if bcm else .35 ** 2

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

    def draw_positions(self, cut='', corr=False, show=True):
        if not self.Ana.has_branch('peak_positions'):
            warning('The peak_positions branch does not exist!')
            return
        h = TH1F('h_pt', 'Peak {m}'.format(m='Timings' if corr else 'Positions'), 1024, 0, 512 if corr else 1024)
        self.Tree.Draw('peak_{p}[{c}]>>h_pt'.format(c=self.Channel, p='positions' if not corr else 'times'), TCut(cut) + TCut('!pulser'), 'goff')
        format_histo(h, x_tit='Time [ns]' if corr else 'Digitiser Bin', y_tit='Number of Entries', y_off=.4, fill_color=self.FillColor, lw=1, tit_size=.05, stats=0)
        self.save_histo(h, 'PeakTimings', show, logy=True, lm=.045, rm=.045, x=4, y=.5)
        return h

    def draw_timings(self, cut='', show=True):
        self.draw_positions(cut, corr=True, show=show)

    def draw_max_position(self, cut='', corr=False, show=True):
        h = TH1F('h_pt', 'Max Peak {m}'.format(m='Timings' if corr else 'Positions'), 1024, 0, 512 if corr else 1024)
        cut = TCut(cut) + TCut('!pulser') if 'pulser' not in cut else TCut(cut)
        self.Tree.Draw('max_peak_{p}[{c}]>>h_pt'.format(c=self.Channel, p='position' if not corr else 'time'), cut, 'goff')
        format_histo(h, x_tit='Time [ns]' if corr else 'Digitiser Bin', y_tit='Number of Entries', y_off=.4, fill_color=self.FillColor, lw=1, tit_size=.07, stats=0, lab_size=.06)
        self.save_histo(h, 'MaxPeak{m}'.format(m='Timings' if corr else 'Positions'), show, logy=True, lm=.073, rm=.045, x=2, y=.5)

    def draw_max_timing(self, cut='', show=True):
        self.draw_max_position(cut, corr=True, show=show)

    def check_peaks(self):
        h, n = self.Ana.draw_waveforms(t_corr=False)
        self.Tree.GetEntry(self.Ana.StartEvent + self.Count - 1)
        print n, self.Ana.StartEvent, self.Ana.count
        a = self.Tree.Draw('peak_positions[{ch}]:n_peaks[{ch}]'.format(ch=self.Channel), '', 'goff', 1, self.Ana.StartEvent + self.Ana.count - 1)
        print [self.Tree.GetV1()[j] / 2 for j in xrange(a)], [self.Tree.GetV2()[j] for j in xrange(a)][0]
        a = self.Tree.Draw('peaks3_x', '', 'goff', 1, self.Ana.StartEvent + self.Ana.count - 1)
        print [i / 2. for i in [self.Tree.GetV1()[j] for j in xrange(a)]]
