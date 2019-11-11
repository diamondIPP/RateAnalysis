#!/usr/bin/env python
# --------------------------------------------------------
#       waveform analysis
# created on May 13th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import *
from InfoLegend import InfoLegend
from ROOT import TCut, TH2F, TProfile, TH1F, TProfile2D


class Waveform(Analysis):
    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.Channel
        self.Tree = self.Ana.Tree
        self.Cut = self.Ana.Cut
        self.set_save_directory(self.Ana.SubDir)
        self.Polarity = self.Ana.Polarity
        self.DUTName = self.Ana.DUTName
        self.DUTNumber = self.Ana.DUTNumber
        self.RunNumber = self.Ana.RunNumber
        self.InfoLegend = InfoLegend(pad_analysis)

        self.Count = 0
        self.StartEvent = self.Ana.StartEvent
        self.BinWidth = self.Ana.DigitiserBinWidth

    def draw(self, n=1, cut=None, start_event=None, t_corr=True, channel=None, show=True, x_range=None, y_range=None, grid=False):
        """ Draws a stack of n waveforms. """
        title = '{n}{tc} Waveform{p}'.format(n=n, tc=' Time Corrected' if t_corr else '', p='s' if n > 1 else '')
        h = TH2F('h_wf', title, 1024, 0, 512, 2048, -512, 512)
        start_count = deepcopy(self.Count)
        values, times = self.get_values(n, cut, start_event, t_corr, channel)
        for v, t in zip(values, times):
            h.Fill(t, v)
        y_range = increased_range([min(values), max(values)], .1, .2) if y_range is None else y_range
        h = self.make_tgrapherrors('g_cw', title, x=times, y=values) if n == 1 else h
        format_histo(h, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, y_range=y_range, markersize=.5, x_range=x_range)
        self.draw_histo(h, 'WaveForms{n}'.format(n=n), show=show, draw_opt='col' if n > 1 else 'apl', lm=.073, rm=.045, bm=.18, x=1.5, y=.5, gridy=grid, gridx=grid)
        return h, self.Count - start_count

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
        cut = self.Cut.AllCut if cut is None else TCut(cut)
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
        cut = self.Cut() if cut is None else TCut(cut)
        channel = self.Channel if channel is None else channel
        if not self.Run.wf_exists(channel):
            return
        start_event = self.get_start_event(start_event)
        self.info('Starting at event {}'.format(start_event))
        n_events = self.Run.find_n_events(n, cut, start_event)
        self.Tree.SetEstimate(n * 1024)
        n_entries = self.Tree.Draw('wf{ch}:trigger_cell'.format(ch=channel), cut, 'goff', n_events, start_event)
        values = [self.Tree.GetV1()[i] for i in xrange(n_entries)]
        times = [self.BinWidth * i for i in xrange(1024)] * n
        if t_corr:
            times = [v for lst in [self.get_calibrated_times(self.Tree.GetV2()[1024 * i]) for i in xrange(n)] for v in lst]
        self.Tree.SetEstimate()
        self.Count += n_events
        return values, times

    def get_calibrated_times(self, trigger_cell):
        if all(self.Run.TCal[0] == t for t in self.Run.TCal):
            return [self.BinWidth * i for i in xrange(self.Run.NSamples)]
        t = [self.Run.TCal[int(trigger_cell)]]
        for i in xrange(1, self.Run.NSamples):
            t.append(self.Run.TCal[(int(trigger_cell) + i) % self.Run.NSamples] + t[-1])
        return t

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
        self.Tree.Draw('rise_time[{}]>>hrt'.format(self.Channel), self.Cut.AllCut if cut is None else TCut(cut), 'goff')
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Rise Time [ns]', y_tit='Number of Entries', y_off=1.4)
        self.save_histo(h, 'RiseTime', lm=.12, show=show)

    def draw_fall_time(self, cut=None, show=True):
        h = TH1F('hft', 'Signal Fall Time', 200, 0, 20)
        self.Tree.Draw('fall_time[{}]>>hft'.format(self.Channel), self.Cut.AllCut if cut is None else TCut(cut), 'goff')
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Fall Time [ns]', y_tit='Number of Entries', y_off=1.4)
        self.save_histo(h, 'FallTime', lm=.12, show=show)

    def draw_rise_time_map(self, res=sqrt(12), cut=None, show=True):
        p = TProfile2D('prtm', 'Rise Time Map', *self.Ana.Bins.get_global(res))
        cut = self.Cut.generate_special_cut(excluded='fiducial') if cut is None else TCut(cut)
        self.Tree.Draw('rise_time[{}]:{}:{}>>prtm'.format(self.Channel, *self.Cut.get_track_vars(self.DUTNumber - 1)), cut, 'goff')
        self.Ana.set_dia_margins(p)
        # self.Ana.set_z_range(p, n_sigma=1)
        self.format_statbox(entries=True, x=.84)
        format_histo(p, x_tit='track x [cm]', y_tit='track y [cm]', y_off=1.4, z_off=1.3, z_tit='Rise Time [ns]', ncont=20, ndivy=510, ndivx=510)
        self.draw_histo(p, show=show, draw_opt='colz', rm=.14)
        self.Ana.draw_fiducial_cut()
        self.save_plots('RiseTimeMap')
