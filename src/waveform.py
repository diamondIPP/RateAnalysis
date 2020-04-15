#!/usr/bin/env python
# --------------------------------------------------------
#       waveform analysis
# created on May 13th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import *
from InfoLegend import InfoLegend
from ROOT import TCut, TH2F, TProfile, TH1F, TProfile2D
from numpy import ones


class Waveform(Analysis):
    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.Channel
        self.Tree = self.Ana.Tree
        self.Cut = self.Ana.Cut()
        self.set_save_directory(self.Ana.SubDir)
        self.Polarity = self.Ana.Polarity
        self.DUT = self.Ana.DUT.Number
        self.RunNumber = self.Ana.RunNumber
        self.InfoLegend = InfoLegend(pad_analysis)

        self.Count = 0
        self.StartEvent = self.Ana.StartEvent
        self.BinWidth = self.Ana.DigitiserBinWidth

    def draw(self, n=1, cut=None, start_event=None, t_corr=True, channel=None, show=True, x_range=None, y_range=None, grid=False, x=None):
        """ Draws a stack of n waveforms. """
        title = '{n}{tc} Waveform{p}'.format(n=n, tc=' Time Corrected' if t_corr else '', p='s' if n > 1 else '')
        h = TH2F('h_wf', title, 1024, 0, 512, 2048, -512, 512)
        start_count = deepcopy(self.Count)
        values, times = self.get_values(n, cut, start_event, t_corr, channel)
        values = values if x is None else x
        for v, t in zip(values, times):
            h.Fill(t, v)
        y_range = increased_range([min(values), max(values)], .1, .2) if y_range is None else y_range
        h = self.make_tgrapherrors('g_cw', title, x=times, y=values) if n == 1 else h
        format_histo(h, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, y_range=y_range, markersize=.5, x_range=x_range)
        self.draw_histo(h, 'WaveForms{n}'.format(n=n), show=show, draw_opt='col' if n > 1 else 'apl', lm=.073, rm=.045, bm=.18, x=1.5, y=.5, gridy=grid, gridx=grid)
        return h, self.Count - start_count

    def draw_all(self, corr=True, n=-1, x_range=None, y_range=None, show=True):
        h = TH2F('h_wf', 'All Waveforms', 1024, 0, 512, 2048, -512, 512)
        times = concatenate(self.get_all_times(corr=corr)).astype('d')[:n]
        values = concatenate(self.get_all()).astype('d')[:n]
        h.FillN(values.size, times, values, ones(values.size))
        y_range = increased_range([min(values), max(values)], .1, .2) if y_range is None else y_range
        format_histo(h, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, markersize=.5, x_range=x_range, y_range=y_range)
        self.draw_histo(h, 'WaveForms{n}'.format(n='e'), show=show, draw_opt='col', lm=.073, rm=.045, bm=.18, x=1.5, y=.5, grid=1, logz=True)

    def get_trigger_cells(self):
        return self.Run.get_root_vec(var='trigger_cell', cut=self.Cut, dtype='i2')

    def get_all(self, redo=False):
        """ extracts all dut waveforms after all cuts from the root tree and saves it as an hdf5 file """
        def f():
            waveforms = []
            events = self.Ana.get_events(cut=self.Cut)
            self.Ana.PBar.start(events.size)
            for ev in events:
                n = self.Tree.Draw('wf0', '', 'goff', 1, ev)
                waveforms.append(self.Ana.Polarity * self.Run.get_root_vec(n, dtype='f2'))
                self.Ana.PBar.update()
            return array(waveforms)
        return do_hdf5(self.make_hdf5_path('WF', run=self.RunNumber, ch=self.Channel), f, redo=redo)

    def get_all_times(self, corr=False):
        times = array([self.get_calibrated_times(trigger_cell) for trigger_cell in self.get_trigger_cells()])
        peaks = self.Ana.Peaks.get_all() if corr else []
        return times - (peaks - peaks[0]).reshape(peaks.size, 1) if corr else times

    def draw_single(self, cut='', event=None, show=True, show_noise=False):
        h, n = self.draw(n=1, start_event=event, cut=cut, t_corr=True, show=show, grid=True)
        if show_noise:
            self.__draw_noise()
        return h

    def draw_all_single(self, n=1, cut='', start_event=None):
        activated_wfs = [wf for wf in xrange(4) if self.Run.wf_exists(wf)]
        print 'activated wafeforms:', activated_wfs
        wfs = [self.draw(n=n, start_event=start_event, cut=cut, show=False, channel=wf)[0] for wf in activated_wfs]
        n_wfs = len(activated_wfs)
        c = self.make_canvas('c_wfs', 'Waveforms', 2, n_wfs * .5)
        c.Divide(1, n_wfs)
        for i, wf in enumerate(wfs):
            wf.SetTitle('{nam} WaveForm'.format(nam=self.Run.DigitizerChannels[activated_wfs[i]]))
            c.cd(i + 1)
            wf.Draw('aclp')

    def draw_average(self, n=100, cut=None, align_peaks=True, show=True, show_noise=False):
        p = TProfile('pawf', 'Averaged Waveform', 2000, 0, 500)
        cut = self.Ana.Cut(cut)
        values, times = self.get_values(n, cut)
        if align_peaks:
            self.Tree.Draw(self.Ana.PeakName, cut, 'goff')
            peak_times = [self.Tree.GetV1()[i] for i in xrange(n)]
            for i, t in enumerate(peak_times):
                for j in xrange(1024):
                    times[j + i * 1024] = times[j + i * 1024] + peak_times[0] - peak_times[i]
        for t, v in zip(times, values):
            p.Fill(t, v)
        format_histo(p, x_tit='Time [ns]', y_tit='Pulse Height [mv]', y_off=1.2, stats=0, markersize=.5)
        self.draw_histo(p, show=show)
        if show_noise:
            self.__draw_noise()

    def __draw_noise(self):
        c = get_last_canvas()
        mean_noise = self.Ana.Pedestal.get_mean()
        legend = self.make_legend(.8, .4, nentries=1, scale=2)
        c.cd()
        legend.AddEntry(self.draw_horizontal_line(self.Polarity * mean_noise.n, 0, 700, w=2, style=7, color=2), 'mean pedestal', 'l')
        legend.Draw()
        c.Update()

    def get_start_event(self, start_event):
        return self.Count + self.StartEvent if start_event is None else start_event

    def get_values(self, n=1, cut=None, start_event=None, t_corr=True, channel=None):
        """ return lists of the values and times of the waveform. """
        cut = self.Ana.Cut(cut)
        channel = self.Channel if channel is None else channel
        if not self.Run.wf_exists(channel):
            return
        start_event = self.get_start_event(start_event)
        self.info('Starting at event {}'.format(start_event))
        n_events = self.Run.find_n_events(n, cut, start_event)
        self.Tree.SetEstimate(n * 1024)
        n_entries = self.Tree.Draw('wf{ch}:trigger_cell'.format(ch=channel), cut, 'goff', n_events, start_event)
        values = self.Run.get_root_vec(n_entries)
        times = [self.BinWidth * i for i in xrange(1024)] * n
        if t_corr:
            times = [v for lst in [self.get_calibrated_times(int(self.Tree.GetV2()[1024 * i])) for i in xrange(n)] for v in lst]
        self.Tree.SetEstimate()
        self.Count += n_events
        return values, times

    def get_calibrated_times(self, trigger_cell):
        # TODO: try directly with TCalSum
        if all(self.Run.TCal == self.Run.TCal[0]):
            return [self.BinWidth * i for i in xrange(self.Run.NSamples)]
        return cumsum(concatenate([[0], self.Run.TCal[trigger_cell:], self.Run.TCal[0:trigger_cell - 1]])) if trigger_cell else cumsum(concatenate([[0], self.Run.TCal[:-1]]))

    def get_calibrated_time_old(self, trigger_cell, bin_nr):
        return sum(self.Run.TCal[trigger_cell+1:trigger_cell + bin_nr]) + (sum(self.Run.TCal[:bin_nr - (1024 - trigger_cell) + 1]) if trigger_cell + bin_nr > 1024 else 0)

    def get_calibrated_time(self, t, b):
        return self.Run.TCalSum[t + b] - self.Run.TCalSum[t]

    def reset(self):
        self.Count = 0

    def calc_rise_time(self, cut=None, start_event=0):
        values, times = self.get_values(1, cut, start_event)
        pedestal = self.Ana.Pedestal.draw_disto_fit(cut=cut, show=False)
        noise, sigma = (abs(pedestal.Parameter(i)) for i in [1, 2])
        rise_time = []
        tmin, tmax = [t * self.BinWidth for t in self.Ana.SignalRegion]
        data = OrderedDict([(t, v) for t, v in zip(times, values) if tmin - 5 < t < tmax + 5])
        for t, y in data.iteritems():
            if abs(y) == max(abs(v) for v in data.itervalues()):  # stop at the highest point
                break
            if abs(y) > noise + 4 * sigma:
                rise_time.append(t)
        return rise_time[-1] - rise_time[0]

    def draw_rise_time(self, cut=None, show=True):
        h = TH1F('hrt', 'Signal Rise Time', 100, 0, 10)
        self.Tree.Draw('rise_time[{}]>>hrt'.format(self.Channel), self.Ana.Cut(cut), 'goff')
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Rise Time [ns]', y_tit='Number of Entries', y_off=1.4)
        self.save_histo(h, 'RiseTime', lm=.12, show=show)

    def draw_fall_time(self, cut=None, show=True):
        h = TH1F('hft', 'Signal Fall Time', 200, 0, 20)
        self.Tree.Draw('fall_time[{}]>>hft'.format(self.Channel), self.Ana.Cut(cut), 'goff')
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Fall Time [ns]', y_tit='Number of Entries', y_off=1.4)
        self.save_histo(h, 'FallTime', lm=.12, show=show)

    def draw_rise_time_map(self, res=sqrt(12), cut=None, show=True):
        p = TProfile2D('prtm', 'Rise Time Map', *self.Ana.Bins.get_global(res))
        cut = self.Ana.Cut.generate_custom(exclude='fiducial') if cut is None else TCut(cut)
        self.Tree.Draw('rise_time[{}]:{}:{}>>prtm'.format(self.Channel, *self.Ana.Cut.get_track_vars(self.DUT.Number - 1)), cut, 'goff')
        self.Ana.set_dia_margins(p)
        # self.Ana.set_z_range(p, n_sigma=1)
        self.format_statbox(entries=True, x=.84)
        format_histo(p, x_tit='track x [cm]', y_tit='track y [cm]', y_off=1.4, z_off=1.3, z_tit='Rise Time [ns]', ncont=20, ndivy=510, ndivx=510)
        self.draw_histo(p, show=show, draw_opt='colz', rm=.14)
        self.Ana.draw_fiducial_cut()
        self.save_plots('RiseTimeMap')
