#!/usr/bin/env python
# --------------------------------------------------------
#       waveform analysis
# created on May 13th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TCut, TMultiGraph
from numpy import fft, argmax, array_split, sum
from src.sub_analysis import PadSubAnalysis
from helpers.draw import *
from helpers.fit import ErfLand


class Waveform(PadSubAnalysis):
    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='WF')
        self.NChannels = self.Run.NChannels
        self.NSamples = self.Run.NSamples

        self.Count = 0
        self.StartEvent = self.Ana.StartEvent
        self.BinWidth = self.Ana.DigitiserBinWidth
        self.Exists = self.has_branch(f'wf{self.Channel}')

    def __getitem__(self, item):
        return self.get_all()[item]

    def reset(self):
        self.Count = 0

    # ----------------------------------------
    # region GET
    def get_active(self):
        return array([wf for wf in range(self.NChannels) if self.Run.wf_exists(wf)])

    def get_binning(self, bin_size=None):
        return make_bins(0, self.Run.NSamples * self.BinWidth, choose(bin_size, self.BinWidth))

    def get_cut(self, cut=None, n=None):
        cut = self.Ana.get_event_cut(self.Cut.to_string(cut)) if type(cut) in [TCut, str] else cut
        cut = array(cut) if is_iter(cut) else self.Ana.get_event_cut() if cut is None else full(self.Run.NEvents, True)
        if n is not None and n < count_nonzero(cut):
            cut[where(cut)[0][n:]] = False
        return cut

    def get_trigger_cells(self, cut=...):
        return self.get_tree_vec('trigger_cell', dtype='i2')[self.get_cut(cut)]

    def get_trigger_cell(self, i):
        return self.get_tree_vec('trigger_cell', '', 'i2', 1, i)[0]

    def get(self, i):
        return self.get_all()[i]

    def get_all(self, channel=None, redo=False):
        """ extracts all dut waveforms after all cuts from the root tree and saves it as an hdf5 file """
        def f():
            with Pool() as pool:
                self.info('Saving signal waveforms to hdf5 ...')
                result = pool.starmap(self._get_all, [(i, channel) for i in array_split(arange(self.Run.NEvents), cpu_count())])
                return concatenate(result)
        return do_hdf5(self.make_simple_hdf5_path(dut=choose(channel, self.Channel)), f, redo=redo)

    def _get_all(self, ind, ch=None):
        var = f'{"-" if self.Ana.Polarity < 0 else ""}wf{choose(ch, self.Channel)}'
        tree, pbar = self.Run.load_rootfile(False), PBar(ind.size) if not ind[0] else None
        return array([self.get_from_tree(tree, ev, var, pbar) for ev in ind])

    @staticmethod
    def get_from_tree(t, ev, var, pbar=None):
        if pbar is not None:
            if ev % 1000:
                pbar.update(ev)
        return get_tree_vec(t, var, '', 'f2', 1, ev)

    def get_values(self, cut=None, channel=None, n=None):
        cut = self.get_cut(cut, n)
        imax = max(where(cut)[0]) + 1
        return array(self.get_all(channel)[:imax])[cut[:imax]].flatten()

    def get_times(self, signal_corr=True, cut=None, n=None, raw=False):
        cut = self.get_cut(cut, n)
        if raw:
            return tile(arange(self.NSamples) * self.BinWidth, count_nonzero(cut))
        t = self.get_all_calibrated_times(cut)
        return (self.correct_times(t, cut) if signal_corr else t).flatten()

    def get_all_calibrated_times(self, cut=None):
        t = array([self.get_calibrated_times(tc) for tc in range(self.NSamples)])
        return t[self.get_trigger_cells(cut)]

    def get_peak_times(self, cut=None):
        return self.Ana.get_peak_times(self.get_cut(cut))

    def correct_times(self, t, cut=None):
        pt = self.get_peak_times(cut)
        return t - (pt - mean(pt)).reshape(pt.size, 1)

    def get_all_cal_times(self):
        return array([self.get_calibrated_times(tc) for tc in range(self.NSamples)])

    def get_calibrated_times(self, trigger_cell):
        return self.Run.TCalSum[trigger_cell:trigger_cell + self.Run.NSamples] - self.Run.TCalSum[trigger_cell]

    def get_calibrated_time(self, t, b):
        return self.Run.TCalSum[t + b] - self.Run.TCalSum[t]

    def get_tree_values(self, n=1, cut=None, t_corr=True, channel=None, raw=False):
        """ return lists of the values and times of the waveform. """
        channel = choose(channel, self.Channel)
        np = self.Count
        if not self.Run.wf_exists(channel):
            return
        events = self.Ana.get_events(cut)[np:np + n]
        self.info('Starting at event {}'.format(events[0]))
        values = concatenate([self.get_tree_vec(['wf{}'.format(channel)], nentries=1, firstentry=ev) for ev in events])
        times = arange(self.Run.NSamples, dtype='u2') * (1 if raw else self.BinWidth)
        if t_corr:
            times = array([self.get_calibrated_times(tc) for tc in self.get_trigger_cells(cut)[np:np + n]]).flatten()
        self.Count += n
        return values, times

    def get_average_rise_time(self, p=.1, redo=False):
        return do_pickle(self.make_simple_pickle_path('RT', int(p * 100)), self.draw_average_rise_time, redo=redo, p=p, show=False)

    def get_max(self, h, region=None):
        x, y = get_graph_vecs(h, err=False)
        xmin, xmax = self.Ana.get_region(region=region) * self.BinWidth
        cut = (x >= xmin) & (x <= xmax)
        i_max = abs(y[cut]).argmax()
        return x[cut][i_max], y[cut][i_max]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region WAVEFORMS
    def draw(self, n=1, cut=None, t_corr=True, channel=None, show=True, x_range=None, y_range=None, grid=False, x=None, raw=False):
        """ Draws a stack of n waveforms. """
        start_count = deepcopy(self.Count)
        title = '{n}{tc} Waveform{p}'.format(n=n, tc=' Time Corrected' if t_corr else '', p='s' if n > 1 else '')
        values, times = self.get_tree_values(n, cut, t_corr, channel, raw)
        values = choose(x, values)
        h = self.Draw.histo_2d(times, values, [1024, 0, 512, 2048, -512, 512], title, show=False)
        y_range = ax_range(min(values), max(values), .1, .2) if y_range is None else y_range
        h = Draw.make_tgrapherrors(times, values, title=title) if n == 1 else h
        format_histo(h, x_tit='Time [ns]', y_tit='Signal [mV]', y_off=.5, stats=0, tit_size=.07, lab_size=.06, y_range=y_range, markersize=.5, x_range=x_range)
        self.Draw(h, 'WaveForms{n}'.format(n=n), show, draw_opt='col' if n > 1 else 'apl', lm=.073, rm=.045, bm=.18, w=1.5, h=.5, gridy=grid, gridx=grid)
        return h, self.Count - start_count

    def draw_all(self, signal_corr=False, raw=False, n=None, x_range=None, y_range=None, cut=None, channel=None, **kwargs):
        n = choose(n, 100000)
        y, x, bins = self.get_values(cut, channel, n), self.get_times(signal_corr, cut, n, raw), self.get_binning() + make_bins(-512, 512.1, .5)
        rx, ry = choose(x_range, [0, 512]), choose(y_range, ax_range(min(y), max(y), .1, .2)),
        h = self.Draw.histo_2d(x, y, bins, f'{n} Waveforms', **Draw.mode(3), **kwargs, x_range=rx, y_range=ry, stats=False, grid=True, logz=True, draw_opt='col', file_name='WaveForms100k')
        return h, n

    def draw_single(self, cut=None, ind=None, x_range=None, y_range=None, t_corr=True, grid=True, raw=False, show=True, show_noise=False, draw_opt='ap'):
        if ind is None:
            g, n = self.draw(n=1, cut=cut, t_corr=t_corr, show=show, grid=grid, x_range=x_range, y_range=y_range)
        else:
            x, y = arange(self.NSamples) if raw else self.get_calibrated_times(self.get_trigger_cell(ind)), self.get(ind)
            g = self.Draw.graph(x, y, 'Single Waveform', **Draw.mode(3), grid=grid, gridy=True, draw_opt=draw_opt, show=show)
        if show_noise:
            self.__draw_noise()
        return g

    def draw_all_single(self, n=1, cut=''):
        activated_wfs = self.get_active()
        self.info('activated wafeforms: {}'.format(activated_wfs))
        wfs = [self.draw(n=n, cut=cut, show=False, channel=wf)[0] for wf in activated_wfs]
        c = Draw.canvas('Waveforms', w=2, h=len(wfs) * .5, divide=(1, len(wfs)))
        for i, wf in enumerate(wfs):
            format_histo(wf, title='{} WaveForm'.format(self.Run.DigitizerChannels[activated_wfs[i]]))
            Draw.histo(wf, canvas=c.cd(i + 1), draw_opt='aclp')

    def draw_bucket(self, show=True, t_corr=True):
        good = self.draw(1, show=False, t_corr=t_corr)[0]
        bucket = self.draw(1, cut=self.Cut.generate_custom(invert='bucket'), show=False, t_corr=t_corr)[0]
        bad_bucket = self.draw(1, cut=self.Cut.get_bad_bucket(), show=False, t_corr=t_corr)[0]
        mg = TMultiGraph('mg_bw', 'Bucket Waveforms')
        leg = Draw.make_legend(.85, .4, nentries=3)
        names = ['good wf', 'bucket wf', 'both wf']
        for i, gr in enumerate([good, bucket, bad_bucket]):
            format_histo(gr, color=self.Draw.get_color(3), markersize=.5)
            mg.Add(gr, 'lp')
            leg.AddEntry(gr, names[i], 'lp')
        self.Draw(mg, show=show, draw_opt='A', w=1.5, h=0.75, lm=.07, rm=.045, bm=.2, leg=leg)
        format_histo(mg, x_range=ax_range(self.Ana.get_region(region='e') * self.BinWidth, None, 0, .3), y_off=.7, x_tit='Time [ns]', y_tit='Signal [mV]')
    # endregion WAVEFORMS
    # ----------------------------------------

    # ----------------------------------------
    # region AVERAGE
    def draw_all_average(self, corr=True, n=100000, cut=None, x_range=None, y_range=None, show=True, show_noise=False, redo=False, draw_opt='ap'):
        def f():
            t, v = self.get_times(corr, cut, n), self.get_values(cut, n=n)
            return Draw.make_graph_from_profile(self.Draw.profile(t, v, self.get_binning(), 'Averaged Waveform', show=False))
        g = do_pickle(self.make_simple_pickle_path('AWF', f'{count_nonzero(self.get_cut(cut))}'), f, redo=redo)
        x_range = ax_range(self.Ana.get_signal_range(), fl=.5, fh=1) if x_range is None else x_range
        format_histo(g, x_tit='Time [ns]', y_tit='Pulse Height [mV]', y_off=1.2, lw=2, stats=0, markersize=.4, x_range=x_range, y_range=y_range)
        self.Draw(g, show=show, draw_opt=draw_opt, gridy=True)
        if show_noise:
            self.__draw_noise(pol=False)
        return g

    def fit_average(self, fit_range=None, n=3, ind=None, show=True):
        max_x = self.Ana.Timing.draw_peaks(show=0).GetListOfFunctions()[1].GetParameter(1)
        h = self.draw_all_average(show=show, cut=ind)
        fit_range = [max_x - 15, max_x + 4] if fit_range is None else fit_range
        c = ErfLand(h, fit_range=fit_range)
        c.fit(n, show)
        format_histo(h, x_range=ax_range(fit_range, fl=.5, fh=.5), stats=0)
        return c

    def fit_landau(self, i, tc, peak_i=None):
        y, t = self.get_all()[i], self.get_calibrated_times(tc)
        imin, imax = self.Ana.SignalRegion
        ip = choose(peak_i, argmax(y[imin:imax])) + imin
        g = Draw.make_tgrapherrors(t[max(0, ip - 6):ip + 8], y[max(0, ip - 6):ip + 8], x_tit='Time [s]', y_tit='Signal [mV]')
        f = TF1('f', 'landau', 0, 512)
        g.Fit(f, 'q')
        m = f.GetMaximumX()
        return f.Integral(m - 4, m + 6) / 10

    def get_fits(self):
        def f():
            events = self.Ana.get_events()
            self.PBar.start(events.size)
            v = []
            tcs = self.get_trigger_cells()
            for ev, tc in zip(events, tcs):
                v.append(self.fit_landau(ev, tc))
                self.PBar.update()
            return array(v).astype('f4')
        return do_hdf5(self.make_simple_hdf5_path('LFits'), f)

    def draw_average_rise_time(self, p=.1, ind=None, x_range=None, y_range=None, show=True):
        x, y = get_graph_vecs(self.draw_all_average(show=show, cut=ind, x_range=x_range, y_range=y_range))
        ymax = max(y).n - self.Ana.get_pedestal().n
        i0, i1 = [next(i for i, v in enumerate(y) if v > ip * ymax) for ip in [p, 1 - p]]
        x0, x1 = [get_x(x[i - 1], x[i], y[i - 1], y[i], ip * ymax) for i, ip in [(i0, p), (i1, 1 - p)]]
        [Draw.vertical_line(x.n) for x in [x0, x1]]
        return x1 - x0

    def compare_averages(self, ind1=None, ind2=None, cut=None, x_range=None, normalise=False, show=True):
        """Compare the average waveform at two subsets, choose index ranges acordingly."""
        cut = self.get_cut(cut)
        n = cut.size // 2
        ind1, ind2 = choose(ind1, concatenate([cut[:n], zeros(cut.size - n, dtype=bool)])), choose(ind2, concatenate([zeros(cut.size - n, dtype=bool), cut[n:]]))
        leg = Draw.make_legend(nentries=2, w=.25)
        graphs = [self.Draw.make_graph_from_profile(self.draw_all_average(show=False, cut=ind, n=None)) for ind in [ind1, ind2]]
        for i, g in enumerate(graphs):
            scale_graph(g, 1 / max(get_graph_y(g)).n) if normalise else do_nothing()
            x_range = ax_range(self.Ana.get_signal_range(), fl=0, fh=1) if x_range is None else x_range
            format_histo(g, color=self.Draw.get_color(2), x_tit='Peak Time [ns]', y_tit='Signal [mV]', x_range=x_range, y_off=1.2)
            leg.AddEntry(g, 'ind{}'.format(i + 1), 'l')
            g.Draw('c') if i else self.Draw(g, draw_opt='ac', show=show, leg=leg, lm=.11)
        update_canvas()

    def draw_average(self, n=100, cut=None, align_peaks=True, show=True, show_noise=False):
        (values, times), peaks = self.get_tree_values(n, self.Cut(cut)), self.get_tree_vec(var=self.Ana.PeakName, cut=self.Cut(cut))[:n]
        times -= (peaks[0] - peaks).repeat(self.NSamples) if align_peaks else 0
        self.Draw.profile(times, values, self.get_binning(), 'Averaged Waveform', x_tit='Time [ns]', y_tit='Pulse Height [mv]', stats=0, markersize=.5, show=show)
        if show_noise:
            self.__draw_noise()
    # endregion AVERAGE
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW INTEGRATION
    def draw_avrg(self, c, cut=None, s=None, n=50000, t_off=0):
        p = self.draw_all_average(cut=cut, n=n, show=0, redo=1)
        if s is not None:
            p.Scale(s / p.GetMaximum())
        x, y = get_hist_vecs(p, err=0)
        g = self.Draw.make_tgrapherrors(x + t_off, y, color=2)
        c.cd()
        self.Draw.legend([g], ['average'], 'l')
        g.Draw('l')

    def draw_signal_integral(self, ind=None, draw_reg=False):
        w = self.draw_single(ind=ind)
        r = [self.draw_sig_region(opacity=.1)] if draw_reg else []
        self.draw_sig_int(w)
        o = [self.draw_peak_pos(w)] + self.draw_sig_peakint(w) + r
        Draw.legend(o, ['peak time', 'left width', 'right width'] + (['signal region'] if draw_reg else []), 'l', scale=1.4)

    def draw_peak_pos(self, h):
        x, y = self.get_max(h)
        return Draw.vertical_line(x, color=4, w=2)

    def draw_buckets(self, start=0, n=8, tl=.04, ts=.05):
        c = get_last_canvas()
        x0 = self.Ana.Peaks.get_bunch_range(start)[0]
        for i in range(n + 1):
            x = x0 + i * self.BunchSpacing
            Draw.vertical_line(x, style=3)
            Draw.vertical_line(x, c.GetUymax() * (1 - tl) + tl * c.GetUymin())
            Draw.tlatex(x + self.BunchSpacing / 2, c.GetUymax() * (1 + tl) + tl * c.GetUymin(), str(i + start), align=21, font=42, size=ts)
        Draw.tlatex(c.GetUxmax(), c.GetUymax() * (1 + 4 * tl) + 4 * tl * c.GetUymin(), 'Bucket Number', font=42, align=31, size=ts)

    def draw_region(self, region=None, lw=2, show_leg=True, fill=False, opacity=None):
        regions = [self.Ana.load_region_name(region=region) for region in make_list(region)]
        lines = []
        for region in regions:
            x1, x2 = self.Ana.get_region(region=region) * self.BinWidth
            color = self.Draw.get_color(len(regions))
            lines.append(Draw.box(x1, -1e6, x2, 1e6, line_color=color, style=2, width=lw, fillcolor=color if fill else None, opacity=choose(opacity, .2)))
        if show_leg:
            Draw.legend(Draw.Objects[-len(regions):], regions, 'l', scale=1.5, w=.1)
        return lines

    def draw_sig_region(self, opacity=None):
        return self.draw_region(show_leg=False, fill=True, opacity=opacity)[0]

    def draw_sig_ped(self):
        s, p = self.draw_sig_region(), self.draw_ped_time()
        Draw.legend([s, p], ['Signal Region', 'Pedestal Time'], scale=2, x2=.95)

    @staticmethod
    def draw_integral(h, xmin, xmax, color=None):
        x, y = get_graph_vecs(h, err=False)
        cut = (x > xmin) & (x < xmax)
        i0, i1 = where(cut)[0][[0, -1]]  # find first and last index fulfilling the condition
        ymin, ymax = interpolate_y(x[i0 - 1], x[i0], y[i0 - 1], y[i0], xmin), interpolate_y(x[i1], x[i1 + 1], y[i1], y[i1 + 1], xmax)
        x = concatenate([[xmin] * 2, x[cut], [xmax] * 2])
        y = concatenate([[0, ymin], y[cut], [ymax, 0]])
        Draw.polygon(x=x, y=y, line_color=2, fill_color=choose(color, Draw.FillColor), opacity=.5)

    def draw_peakint(self, x, peakint=None, y=None):
        y = choose(y, -20 * self.Polarity)
        imin, imax = self.Ana.get_peak_integral(peakint) * self.BinWidth
        return [Draw.arrow(x + v, x, y, y, col=col, width=2, opt='<', size=.02) for v, col in [(-imin, 618), (imax, 434)]]

    def draw_sig_peakint(self, h, peakint=None):
        return self.draw_peakint(self.get_max(h)[0], peakint)

    def draw_sig_int(self, h, peakint=None, color=None):
        xmin, xmax = self.Ana.get_peak_integral(peakint) * [-1, 1] * self.BinWidth + self.get_max(h)[0]
        format_histo(h, x_range=ax_range(xmin, xmax, 2, 2))
        self.draw_integral(h, xmin, xmax, color)

    def draw_ped_time(self, color=4, lw=2):
        return Draw.vertical_line(self.Ana.Pedestal.Region * self.BinWidth, -1000, 1000, color=color, w=lw)

    def draw_pedestal(self, peakint=None):
        x, y = self.Ana.get_region('pedestal')[0] * self.BinWidth, 20 * self.Polarity
        print(x, y)
        Draw.vertical_line(x, -1000, 1000, color=4, w=3)
        Draw.tlatex(x + 1, -y, 'Pedestal', color=4, align=12)
        self.draw_peakint(x, peakint)

    def draw_ped_int(self, h, peakint=None):
        xmin, xmax = (self.Ana.get_peak_integral(peakint) * [-1, 1] + self.Ana.get_region('pedestal')[0]) * self.BinWidth
        format_histo(h, y_range=[-30, 30], x_range=ax_range(xmin, xmax, .3, 1.3))
        self.draw_integral(h, xmin, xmax, color=4)
    # endregion SHOW INTEGRATION
    # ----------------------------------------

    def integrate(self, i, r1=None, r2=None, t1=4, t2=6):
        t, v = self.get_calibrated_times(self.get_trigger_cell(i)), self[i] * self.Polarity
        r1, r2 = (r1, r2) if r1 is not None else self.Ana.SignalRegion * self.BinWidth
        return self._integrate(t, v, r1, r2, t1, t2)

    @update_pbar
    def _integrate(self, t, v, r1, r2, t1, t2):
        ir = where((t > r1) & (t < r2))[0]
        tmax = t[argmax(v[ir]) + ir[0]]
        i = where((t >= tmax - t1) & (t <= tmax + t2))[0]
        i = arange(i[0] - 1, i[-1] + 2)  # include values left and right
        t, v = t[i], v[i]
        v = concatenate([[get_y(t[0], t[1], v[0], v[1], tmax - t1)], v[1:-1], [get_y(t[-2], t[-1], v[-2], v[-1], tmax + t2)]])
        t = concatenate([[tmax - t1], t[1:-1], [tmax + t2]])
        return sum((t[1:] - t[:-1]) * (v[:-1] + v[1:])) / 2 / (t1 + t2)

    def get_integrals(self, sig_region=None, peak_int=None, redo=False):
        def f():
            t, v, tc = array([self.get_calibrated_times(tc) for tc in range(self.NSamples)]), array(self.get_all()), self.get_trigger_cells()
            (r1, r2), (t1, t2) = choose(sig_region, self.Ana.SignalRegion * self.BinWidth), choose(peak_int, array(self.Ana.PeakIntegral) * self.BinWidth)
            self.info(f'averaging waveforms around peak in region [{r1:.1f}, {r2:.1f}] with a range of [{t1:.0f}, {t2:.0f}]')
            self.PBar.start(tc.size)
            return array([self._integrate(t[tc[i]], v[i], r1, r2, t1, t2) for i in range(tc.size)]).astype('f4')
        return array(do_hdf5(self.make_simple_hdf5_path('Int', f'{sig_region}_{peak_int}'), f, redo=redo))

    def __draw_noise(self, pol=True):
        c = get_last_canvas()
        mean_noise = self.Ana.Pedestal.get_mean()
        legend = Draw.make_legend(w=.1)
        c.cd()
        legend.AddEntry(Draw.horizontal_line(mean_noise.n * (self.Polarity if pol else 1), 0, 700, w=2, style=7, color=2), 'mean pedestal', 'l')
        legend.Draw()
        update_canvas()

    def draw_rise_time(self, cut=None, show=True):
        values = self.get_tree_vec(var='rise_time[{}]'.format(self.Channel), cut=self.Cut(cut))
        return self.Draw.distribution(values, make_bins(0, 10, .1), 'Signal Rise Time', x_tit='Rise Time [ns]', file_name='RiseTime', show=show)

    def draw_fall_time(self, cut=None, show=True):
        values = self.get_tree_vec(var='fall_time[{}]'.format(self.Channel), cut=self.Cut(cut))
        return self.Draw.distribution(values, make_bins(0, 20, .1), 'Signal Fall Time', x_tit='Fall Time [ns]', show=show)

    def draw_rise_time_map(self, res=None, cut=None, show=True):
        cut = self.Ana.Cut.generate_custom(exclude='fiducial') if cut is None else TCut(cut)
        rt, x, y = self.get_tree_vec(var=['rise_time[{}]'.format(self.Channel)] + self.Ana.get_track_vars(), cut=cut)
        p = self.Draw.prof2d(x, y, rt, self.Ana.Bins.get_global(res), 'Rise Time Map', show=show, draw_opt='colz')
        format_histo(p, x_tit='track x [cm]', y_tit='track y [cm]', z_tit='Rise Time [ns]', ncont=20, ndivy=510, ndivx=510)
        self.Ana.draw_fid_cut()
        self.Draw.save_plots('RiseTimeMap')

    def draw_fft(self, ind):
        # TODO: finish
        values = fft.fft(self.get_values(ind))
        print(values)
