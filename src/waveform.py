#!/usr/bin/env python
# --------------------------------------------------------
#       waveform analysis
# created on May 13th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import *
from InfoLegend import InfoLegend
from ROOT import TCut, TH2F, TProfile, TH1F, TProfile2D
from numpy import fft


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
        self.InfoLegend = InfoLegend(pad_analysis)

        self.Count = 0
        self.StartEvent = self.Ana.StartEvent
        self.BinWidth = self.Ana.DigitiserBinWidth
        self.set_pickle_sub_dir('WF')

    def get_binning(self, bin_size=.5):
        bins = arange(0, self.Run.NSamples * self.Ana.DigitiserBinWidth, bin_size)
        return [bins.size - 1, bins]

    def draw(self, n=1, cut=None, start_event=None, t_corr=True, channel=None, show=True, x_range=None, y_range=None, grid=False, x=None, raw=False):
        """ Draws a stack of n waveforms. """
        title = '{n}{tc} Waveform{p}'.format(n=n, tc=' Time Corrected' if t_corr else '', p='s' if n > 1 else '')
        h = TH2F('h_wf', title, 1024, 0, 512, 2048, -512, 512)
        start_count = deepcopy(self.Count)
        values, times = self.get_tree_values(n, cut, start_event, t_corr, channel, raw)
        values = values if x is None else x
        for v, t in zip(values, times):
            h.Fill(t, v)
        y_range = increased_range([min(values), max(values)], .1, .2) if y_range is None else y_range
        h = self.make_tgrapherrors('g_cw', title, x=times, y=values) if n == 1 else h
        format_histo(h, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, y_range=y_range, markersize=.5, x_range=x_range)
        self.save_histo(h, 'WaveForms{n}'.format(n=n), show, draw_opt='col' if n > 1 else 'apl', lm=.073, rm=.045, bm=.18, x=1.5, y=.5, gridy=grid, gridx=grid)
        return h, self.Count - start_count

    def draw_all(self, corr=True, n=None, x_range=None, y_range=None, ind=None, channel=None, draw_opt=None, show=True):
        n = -1 if n is None else 1024 * n
        values, times = self.get_values(ind, channel)[:n], self.get_times(corr, ind)[:n]
        if values.size > self.Run.NSamples:
            h = TH2F('hwf', 'All Waveforms', 1024, 0, 512, 2048, -512, 512)
            h.FillN(values.size, times.astype('d'), values.astype('d'), ones(values.size))
        else:
            h = self.make_tgrapherrors('gaw', 'Waveform', x=times, y=values)
        y_range = increased_range([min(values), max(values)], .1, .2) if y_range is None else y_range
        format_histo(h, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, markersize=.5, x_range=x_range, y_range=y_range)
        draw_opt = draw_opt if draw_opt is not None else 'col' if values.size > self.Run.NSamples else 'ap'
        self.draw_histo(h, show, draw_opt=draw_opt, lm=.073, rm=.045, bm=.18, x=1.5, y=.5, grid=1, logz=True)
        return h, n

    def draw_single(self, cut='', event=None, ind=None, x_range=None, y_range=None, draw_opt=None, show=True, show_noise=False):
        h, n = self.draw(n=1, start_event=event, cut=cut, t_corr=True, show=show, grid=True) if ind is None else self.draw_all(False, 1, x_range, y_range, ind, draw_opt=draw_opt, show=show)
        if show_noise:
            self.__draw_noise()
        return h

    def draw_fft(self, ind):
        # TODO: finish
        values = fft.fft(self.get_values(ind))

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

    def fit_average(self, fit_range=None, n=3, ind=None, show=True):
        max_x = self.Ana.Timing.draw_peaks(show=0).GetListOfFunctions()[1].GetParameter(1)
        h = self.draw_all_average(show=show, ind=ind)
        fit_range = [max_x - 15, max_x + 4] if fit_range is None else fit_range
        from fit import ErfLand
        c = ErfLand(h, fit_range=fit_range)
        c.fit(n, show)
        format_histo(h, x_range=increased_range(fit_range, .5, .5), stats=0)
        return c

    def get_average_rise_time(self, p=.1, ind=None, x_range=None, y_range=None, show=False, redo=False):
        def f():
            h = self.draw_all_average(show=show, ind=ind, x_range=x_range, y_range=y_range)
            maxval = h.GetBinContent(h.GetMaximumBin()) - self.Ana.get_pedestal().n
            bins = [h.FindFirstBinAbove(ip * maxval) for ip in [1 - p, p]]
            coods = [(h.GetBinCenter(ib), h.GetBinContent(ib), h.GetBinCenter(ib - 1), h.GetBinContent(ib - 1), ib) for ib in bins]
            f1, f2 = [interpolate_two_points(*icood) for icood in coods]
            if show:
                self.add(f1, f2)
                self.draw_vertical_line(f1.GetX((1 - p) * maxval), -100, 1e4, name='1')
                self.draw_vertical_line(f2.GetX(p * maxval), -100, 1e4, name='2')
                f1.Draw('same')
                f2.Draw('same')
            return f1.GetX((1 - p) * maxval) - f2.GetX(p * maxval)
        return do_pickle(self.make_simple_pickle_path('RT'), f, redo=redo)

    def draw_all_average(self, corr=True, n=-1, ind=None, prof=True, x_range=None, y_range=None, show=True, show_noise=False, redo=False):
        def f():
            p1 = TProfile('paawf', 'Averaged Waveform', 1024, 0, 512) if prof else TH2F('haawf', 'Overlayed Waveforms', 1024, 0, 512, *self.Ana.Bins.get_pad_ph(4))
            values = self.get_values(ind)[:n]
            p1.FillN(values.size, self.get_times(corr, ind).astype('d')[:n], values.astype('d'), ones(values.size))
            return p1
        p = do_pickle(self.make_simple_pickle_path('AWF', '{}_{}'.format(len(ind) if ind is not None else '', int(prof))), f, redo=redo)
        x_range = increased_range(self.Ana.SignalRegion * self.BinWidth, 0, .3) if x_range is None else x_range
        format_histo(p, x_tit='Time [ns]', y_tit='Pulse Height [mV]', y_off=1.2, stats=0, markersize=.5, x_range=x_range, y_range=y_range)
        self.draw_histo(p, show=show, draw_opt='' if prof else 'col')
        if show_noise:
            self.__draw_noise(pol=False)
        return p

    def compare_averages(self, ind1, ind2, x_range=None, normalise=False, show=True):
        g0, g1 = [self.make_graph_from_profile(self.draw_all_average(show=False, ind=i)) for i in [ind1, ind2]]
        if normalise:
            scale_graph(g0, 1 / max(get_graph_y(g0)).n)
            scale_graph(g1, 1 / max(get_graph_y(g1)).n)
        format_histo(g0, color=get_color(2, 1), x_tit='Peak Time [ns]', y_tit='Signal [mV]', x_range=x_range, y_off=1.2)
        format_histo(g1, color=get_color(2, 0))
        leg = self.make_legend(nentries=2, w=.25)
        leg.AddEntry(g0, 'high pulse height', 'l')
        leg.AddEntry(g1, 'low pulse height', 'l')
        self.draw_histo(g0, draw_opt='ac', show=show, leg=leg, lm=.11)
        g1.Draw('c')

    def get_trigger_cells(self, redo=False):
        return do_hdf5(self.make_simple_hdf5_path('TC', self.get_cut_name()), self.Run.get_root_vec, redo=redo, var='trigger_cell', cut=self.Cut, dtype='i2')

    def get_all(self, channel=None, redo=False):
        """ extracts all dut waveforms after all cuts from the root tree and saves it as an hdf5 file """
        def f():
            waveforms = []
            events = self.Ana.get_events(cut=self.Cut)
            self.Ana.PBar.start(events.size)
            for ev in events:
                n = self.Tree.Draw('wf{}'.format(self.Channel if channel is None else channel), '', 'goff', 1, ev)
                waveforms.append(self.Ana.Polarity * self.Run.get_root_vec(n, dtype='f2'))
                self.Ana.PBar.update()
            return array(waveforms)
        return do_hdf5(self.make_simple_hdf5_path(suf=self.get_cut_name(), dut=self.Channel if channel is None else channel), f, redo=redo)

    def get_cut_name(self):
        return self.Cut.GetName() if not self.Cut.GetName().startswith('All') else ''

    def get_values(self, ind=None, channel=None):
        return array(self.get_all(channel=channel))[ind].flatten()

    def get_times(self, corr=True, ind=None):
        return array(self.get_all_times(corr))[ind].flatten()

    def get_all_times(self, corr=False, redo=False):
        def f():
            self.PBar.start(self.get_trigger_cells().size)
            times = []
            for i, tc in enumerate(self.get_trigger_cells()):
                times.append(self.get_calibrated_times(tc).astype('f2'))
                if not i % 100:
                    self.PBar.update(i)
            self.PBar.finish()
            return array(times)
        t = do_hdf5(self.make_simple_hdf5_path('Times', self.get_cut_name()), f, redo=redo)
        peaks = self.Ana.Peaks.get_from_tree() if corr else []
        return array(t) - (peaks - peaks[0]).reshape(peaks.size, 1) if corr else t

    def draw_average(self, n=100, cut=None, align_peaks=True, show=True, show_noise=False):
        p = TProfile('pawf', 'Averaged Waveform', 2000, 0, 500)
        cut = self.Ana.Cut(cut)
        values, times = self.get_tree_values(n, cut)
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

    def __draw_noise(self, pol=True):
        c = get_last_canvas()
        mean_noise = self.Ana.Pedestal.get_mean()
        legend = self.make_legend()
        c.cd()
        legend.AddEntry(self.draw_horizontal_line(mean_noise.n * (self.Polarity if pol else 1), 0, 700, w=2, style=7, color=2), 'mean pedestal', 'l')
        legend.Draw()
        c.Update()

    def get_start_event(self, start_event):
        return self.Count + self.StartEvent if start_event is None else start_event

    def get_tree_values(self, n=1, cut=None, start_event=None, t_corr=True, channel=None, raw=False):
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
        times = arange(self.Run.NSamples, dtype='u2') * (1 if raw else self.BinWidth)
        if t_corr:
            times = [v for lst in [self.get_calibrated_times(int(self.Tree.GetV2()[1024 * i])) for i in xrange(n)] for v in lst]
        self.Tree.SetEstimate()
        self.Count += n_events
        return values, times

    def get_calibrated_times(self, trigger_cell):
        return self.Run.TCalSum[trigger_cell:trigger_cell + self.Run.NSamples] - self.Run.TCalSum[trigger_cell]

    def get_calibrated_time_old(self, trigger_cell, bin_nr):
        return sum(self.Run.TCal[trigger_cell+1:trigger_cell + bin_nr]) + (sum(self.Run.TCal[:bin_nr - (1024 - trigger_cell) + 1]) if trigger_cell + bin_nr > 1024 else 0)

    def get_calibrated_time(self, t, b):
        return self.Run.TCalSum[t + b] - self.Run.TCalSum[t]

    def get_bin(self, tc, t):
        lst = where(self.Run.TCalSum == t + self.Run.TCalSum[tc])[0]
        return lst[0] - tc if lst.size else None

    def reset(self):
        self.Count = 0

    def calc_rise_time(self, cut=None, start_event=0):
        values, times = self.get_tree_values(1, cut, start_event)
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
