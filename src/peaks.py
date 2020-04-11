#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import *
from ROOT import TH1F, TCut, gROOT
from scipy.signal import find_peaks, savgol_filter
from numpy import polyfit, pi, RankWarning, vectorize, size, split
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

    def calc_threshold(self):
        """ return peak threshold, 6 times the raw noise or the minimum singnal, whatever is higher. """
        return max(abs(self.Ana.Pedestal.get_raw_mean() + 6 * self.Ana.Pedestal.get_raw_noise()), self.Ana.get_min_signal())

    def get_all(self):
        return do_hdf5(self.make_hdf5_path('Peaks', 'V1', self.Ana.RunNumber, self.Channel), self.Run.get_root_vec, var=self.Ana.PeakName, cut=self.Ana.Cut(), dtype='f2')

    def draw(self, corr=True, show=True, redo=False, fit=True):
        h = TH1F('hp', 'Peak Times', 512 * 2, 0, 512)
        values, n_peaks = self.find_all(redo=redo, fit=fit)
        if corr:
            values = array(split(values, cumsum(n_peaks)[:-1]))
            peaks = self.get_all()
            for i in xrange(values.size):
                values[i] -= peaks[i] - peaks[0]
            values = concatenate(values)
        h.FillN(values.size, array(values, 'd'), full(values.size, 1, 'd'))
        self.format_statbox(entries=True)
        format_histo(h, x_tit='Time [ns]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor, stats=0)
        self.draw_histo(h, lm=.12, show=show, x=1.5, y=0.75, logy=True)

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
        find_n_peaks = vectorize(size)
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

    def draw_n_peaks(self, spec=False, show=True, do_fit=False):
        h = TH1F('h_pn', 'Number of Peaks', 10, 0, 10)
        draw_var = '@peaks{ch}_x.size()>>h_pn' if spec and self.Run.has_branch('peaks{ch}_x'.format(ch=self.Channel)) else 'n_peaks[{ch}] - 1>>h_pn'
        self.Tree.Draw(draw_var.format(ch=self.Channel), self.Ana.Cut(), 'goff')
        self.format_statbox(only_fit=True, w=.3) if do_fit else self.format_statbox(entries=True)
        format_histo(h, x_tit='Number of Peaks', y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor, lw=2)
        self.save_histo(h, 'PeakNumbers', show, logy=True, lm=.11)
        return h

    def draw_n_peaks_fit(self, show=True, spec=False):
        h = self.draw_n_peaks(spec, show=False, do_fit=True)
        format_histo(h, 'Fit Result')
        self.draw_histo(h, '', show, logy=True, lm=.11)
        fit_func = fit_poissoni(h, show=show)
        self.save_histo(h, 'PeakNumbersFit', show, logy=True, lm=.11, draw_opt='e1same', canvas=gROOT.GetListOfCanvases()[-1])
        self.get_flux(fit_func=fit_func)
        return fit_func

    def get_flux(self, fit_func=None):
        fit_func = self.draw_n_peaks_fit(show=False) if fit_func is None else fit_func
        lambda_ = fit_func.GetParameter(1)
        err = fit_func.GetParError(1)
        n_bunches = 25
        proc_frequency = 50.6e6
        flux = lambda_ * proc_frequency / (n_bunches * self.get_area()) / 1000.
        err *= proc_frequency / (n_bunches * self.get_area()) / 1000.
        print 'Estimated Flux by number of peaks: {f:0.4f} kHz/cm2'.format(f=flux)
        return flux, err

    def get_area(self):
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
