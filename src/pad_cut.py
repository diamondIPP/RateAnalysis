# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TCut, TF1
from src.cut import Cut, CutString, Bins
from helpers.draw import *
from numpy import quantile, histogram2d


class PadCut(Cut):
    """ Generates all cut strings which correspond to a single DUT with digitiser. """

    def __init__(self, analysis):
        Cut.__init__(self, analysis)
        self.Channel = self.Ana.get_channel()
        self.generate_dut()
        self.ConsecutiveCuts = self.get_consecutive()

    # ----------------------------------------
    # region GET
    def get_raw_pedestal(self, sigma=False, redo=False):
        def f():
            return mean_sigma(self.get_tree_vec(var=self.Ana.get_pedestal_name(), cut=self()))
        return do_pickle(self.make_simple_pickle_path('RawPed'), f, redo=redo)[1 if sigma else 0]

    def get_raw_noise(self, redo=False):
        return self.get_raw_pedestal(sigma=True, redo=redo)

    def get_raw_snr(self):
        return self.get_raw_pulse_height() / self.get_raw_noise()

    def get_pre_bucket(self):
        return CutString('!pre bucket', self.generate_custom(include=['pulser', 'fiducial', 'event range'], prnt=False) + self.generate_pre_bucket().invert())()

    def get_bucket(self):
        return CutString('!bucket', self.generate_custom(include=['pulser', 'fiducial', 'event range'], prnt=False) + self.get('bucket', invert=True))()

    def get_bad_bucket(self):
        return self.generate_custom(exclude='bucket') + self.generate_pre_bucket().invert() + TCut('{} > {}'.format(self.Ana.get_raw_signal_var(), self.get_bucket_threshold()))

    def get_pulser(self, beam_on=True):
        cut = self.generate_custom(include=['ped sigma', 'event range'], prnt=False) + Cut.invert(self.get('pulser'))
        cut += self.get('beam stops', warn=False) if beam_on else Cut.invert(self.generate_jump_cut())
        return CutString('PulserBeam{}'.format('On' if beam_on else 'Off'), cut)

    def get_timing(self):
        return CutString('timing', self.generate_custom(exclude='timing', prnt=False), 'All cuts expect timing.')
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE
    def generate_dut(self):

        # -- WAVEFORM --
        self.CutStrings.register(self.generate_pulser(), 2)
        self.CutStrings.register(self.generate_saturated(), level=3)

        # -- SIGNAL --
        self.CutStrings.register(self.generate_pedestal_sigma(), 30)
        # self.CutStrings.register(self.generate_threshold(), 31)
        # self.CutStrings.register(self.generate_timing(), 35)
        # self.CutStrings.register(self.generate_cft(), 36)

        # -- FIDUCIAL --
        self.CutStrings.register(self.generate_fiducial(), 90)

        # --BUCKET --
        self.CutStrings.register(self.generate_bucket(), 36)

        # --SIGNAL DROPS--
        self.update('event range', self.generate_event_range(None, self.find_zero_ph_event()).Value)  # update event range if drop is found

    def generate_saturated(self):
        cut_string = '!is_saturated[{ch}]'.format(ch=self.Channel)
        description = 'exclude {:.2f}% saturated events'.format(100. * self.find_n_saturated(Cut.invert(cut_string)) / self.Run.NEvents)
        return CutString('saturated', cut_string, description)

    def generate_pulser(self):
        return CutString('pulser', '!pulser', 'exclude {:.1f}% pulser events'.format(100. * self.find_n_pulser('pulser') / self.Run.NEvents))

    def generate_pedestal_sigma(self, sigma=None):
        sigma = choose(sigma, self.Config.get_value('CUT', 'pedestal sigma', dtype=int))
        ped_range = self.calc_pedestal(sigma)
        description = 'sigma <= {} -> [{:1.1f}mV, {:1.1f}mV]'.format(sigma, *ped_range)
        return CutString('ped sigma', '{ped} > {} && {ped} < {}'.format(ped=self.Ana.get_pedestal_name(), *ped_range), description)

    def generate_threshold(self):
        thresh = self.calc_threshold(show=False)
        return CutString('threshold', '{} > {}'.format(self.Ana.get_raw_signal_var(), thresh), thresh)

    def generate_pre_bucket(self):
        """select only events when the signal in the signal and consecutive bucket are the same."""
        return CutString('pre-bucket', '{} == {}'.format(self.Ana.get_signal_name(region='e'), self.Ana.get_signal_name()))

    def generate_bucket(self, thresh=None):
        """select bucket events below threshold and then negate it"""
        thresh = choose(thresh, self.get_bucket_threshold())
        string = self.generate_pre_bucket()() if thresh is None else (CutString('', self.generate_pre_bucket().invert()) + '{} < {}'.format(self.Ana.get_raw_signal_var(), thresh)).invert()
        return CutString('bucket', string, 'bucket events with{}'.format('out threshold' if thresh is None else ' threshold < {:.1f}mV'.format(thresh)))

    def generate_cft(self, n_sigma=3):
        fit = self.calc_cft()
        m, s = fit.Parameter(1), fit.Parameter(2)
        description = '{:1.1f}ns < constant fraction time < {:.1f}ns'.format(m - n_sigma * s, m + n_sigma * s)
        var = self.Ana.Timing.get_cft_var()
        return CutString('cft', '{} < {v} && {v} < {}'.format(m - n_sigma * s, m + n_sigma * s, v=var), description)

    def generate_timing(self, n_sigma=3):
        t_correction, fit = self.calc_timing_range()
        if fit is None:
            return TCut('')
        corrected_time = '{peak} - {t_corr}'.format(peak=self.Ana.Timing.get_peak_var(corr=True), t_corr=t_correction)
        string = 'TMath::Abs({cor_t} - {mp}) / {sigma} < {n_sigma}'.format(cor_t=corrected_time, mp=fit.GetParameter(1), sigma=fit.GetParameter(2), n_sigma=n_sigma)
        description = '{:1.1f}ns < peak timing < {:.1f}ns'.format(fit.GetParameter(1) - fit.GetParameter(2), fit.GetParameter(1) + fit.GetParameter(2))
        return CutString('timing', string, description)
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region COMPUTE
    def find_n_pulser(self, cut, redo=False):
        return do_pickle(self.make_simple_pickle_path('NPulser', dut=''), self.Tree.GetEntries, None, redo, self(cut).GetTitle())

    def find_n_saturated(self, cut, redo=False):
        return do_pickle(self.make_simple_pickle_path('NSaturated'), self.Tree.GetEntries, None, redo, self(cut).GetTitle())

    def find_fid_cut(self, thresh=.93, show=True):
        h = self.Ana.draw_signal_map(show=show)
        px, py = h.ProjectionX(), h.ProjectionY()

        def find_fid_margins(p, t):
            t0 = p.GetMaximum() * t
            xbin1, xbin2 = p.FindFirstBinAbove(t0), p.FindLastBinAbove(t0)
            f1 = interpolate_two_points(p.GetBinCenter(xbin1), p.GetBinContent(xbin1), p.GetBinCenter(xbin1 - 1), p.GetBinContent(xbin1 - 1))
            f2 = interpolate_two_points(p.GetBinCenter(xbin2), p.GetBinContent(xbin2), p.GetBinCenter(xbin2 + 1), p.GetBinContent(xbin2 + 1))
            return [f1.GetX(t0) / 10, f2.GetX(t0) / 10]
        x, y = find_fid_margins(px, thresh), find_fid_margins(py, thresh)
        self.draw_fid()
        Draw.box(*array([x[0], y[0], x[1], y[1]]) * 10, show=show, width=2)
        return '"{}": [{}]'.format(self.Ana.DUT.Name, ', '.join('{:0.3f}'.format(i) for i in x + y))

    def find_beam_interruptions(self, bin_width=100, max_thresh=.6):
        """ Looking for the beam interruptions by investigating the pulser rate. """
        t = self.info('Searching for beam interruptions of run {r} ...'.format(r=self.Run.Number), endl=False)
        x, y = self.get_tree_vec(var=['Entry$', 'pulser'], dtype='i4')
        rates, x_bins, y_bins = histogram2d(x, y, bins=[arange(0, x.size, bin_width, dtype=int), 2])
        rates = rates[:, 1] / bin_width
        thresh = min(max_thresh, mean(rates) + .2)
        events = x_bins[:-1][rates > thresh] + bin_width / 2
        not_connected = where(concatenate([[False], events[:-1] != events[1:] - bin_width]))[0]  # find the events where the previous event is not related to the event (more than a bin width away)
        events = split(events, not_connected)  # events grouped into connecting events
        interruptions = [(ev[0], ev[0]) if ev.size == 1 else (ev[0], ev[-1]) for ev in events] if events[0].size else []
        self.add_to_info(t)
        return array(interruptions, 'i4')

    def get_bucket_threshold(self):
        """get the bucket threshold from the run with the highest rate"""
        high_rate_run = self.Run.get_high_rate_run()
        picklepath = self.make_simple_pickle_path('BucketThreshold', run=high_rate_run)
        if not file_exists(picklepath) and not self.Run.Number == high_rate_run:
            from src.pad_analysis import PadAnalysis
            return PadAnalysis(high_rate_run, self.DUT.Number, self.TCString).Cut.calc_bucket_threshold()
        return self.calc_bucket_threshold(run=high_rate_run)

    def calc_bucket_threshold(self, run=None, show=False):
        """calculated the threshold of the bucket pedestal so that the mean keeps unchanged."""
        def f():
            t = self.Ana.info('Calculating signal threshold for bucket cut of run {run} and {d} ...'.format(run=self.Run.Number, d=self.DUT.Name), endl=False)
            h = self.draw_bucket_histo(self.Ana.get_raw_signal_var(), show=show)
            if h.GetEntries() < 1000:  # cannot fit with too little events
                self.Ana.add_to_info(t)
                info('Not enough bucket events ({:.0f})'.format(h.GetEntries()))
                return
            if self.get_raw_snr() < 4:  # impossible to separate for the small signal to noise ratios...
                warning('Signal to noise ration is too low... -> taking raw pedestal + 2 * raw noise for bucket cut!')
                return (self.get_raw_pedestal() + self.get_raw_noise() * 2).n
            fit = self.fit_bucket(h)  # extract fit functions
            if fit is None or any([abs(fit.GetParameter(i)) < 20 for i in [0, 3]]) or fit.GetParameter(1) < fit.GetParameter(4) or fit.GetParameter(1) > 500:
                self.Ana.add_to_info(t)
                return warning('bucket cut fit failed')
            ped = Draw.make_f('ped', 'gaus(0) + gaus(3)', -50, 300, [fit.GetParameter(i) for i in range(3, 9)])
            h_sig = deepcopy(h)
            h_sig.Add(ped, -1)
            format_histo(h_sig, x_range=Bins.PadPHRange)
            values = self.get_tree_vec(var=self.Ana.get_raw_signal_var(), cut=self.get_pre_bucket())
            f0 = Draw.make_tf1('Mean', lambda x: mean(values[values > x]), -30, fit.GetParameter(1))
            self.add_to_info(t)
            threshold = f0.GetX(h_sig.GetMean())
            Draw.vertical_line(threshold, -50, 1e6)
            self.Draw(h_sig, show=show)
            format_histo(f0, x_tit='Threshold', y_tit='Mean')
            self.Draw(f0, show=show)
            Draw.horizontal_line(h_sig.GetMean(), -100, 100)
            Draw.vertical_line(threshold, -50, 1e6)
            return threshold
        return do_pickle(self.make_simple_pickle_path('BucketThreshold', run=run), f, redo=show or show)

    def draw_error_bucket_thresh(self, show=True, t_range=None):
        h = self.draw_bucket_histo(self.Ana.get_raw_signal_var(), show=False)
        fit = self.fit_bucket(h)  # extract fit functions
        sig = Draw.make_f('sig', 'gaus', -50, 300, [fit.GetParameter(i) for i in range(3)], line_style=3, color=4)
        ped = Draw.make_f('ped', 'gaus(0) + gaus(3)', -50, 300, [fit.GetParameter(i) for i in range(3, 9)], line_style=2)
        h_sig = deepcopy(h)
        h_sig.Add(ped, -1)  # subtract pedestal fit from real signal distribution

        t_range = choose(t_range, [-30, fit.GetParameter(1)])
        e_s = Draw.make_tf1('#varepsilon_{sig}', lambda x: 1 - sig.Integral(-50, x) / h_sig.Integral() / h_sig.GetBinWidth(1), *t_range, w=2)
        e_b = Draw.make_tf1('#varepsilon_{bg}', lambda x: ped.Integral(-50, x) / ped.Integral(-50, 500), *t_range, color=1, w=2)
        roc = Draw.make_tf1('ROC Curve', lambda x: e_s(e_b.GetX(x)), 0, 1, color=1, w=2)
        s = Draw.make_tf1('Signal Events', lambda x: sig.Integral(x, fit.GetParameter(1)) + h_sig.Integral(h_sig.FindBin(fit.GetParameter(1)), h_sig.GetNbinsX()) * h_sig.GetBinWidth(1), *t_range)
        err = Draw.make_tf1('Signal Error', lambda x: s(x) / sqrt(s(x) + ped.Integral(x, 200)), *t_range, w=2)
        t = err.GetMaximumX()

        c = Draw.canvas('Signal Threshold Overview', divide=(2, 2), w=1.5, h=1.5)
        # Bucket cut plot
        self.Draw(h, show=show, lm=.135, canvas=c.cd(1), leg=[sig, ped])
        Draw.y_axis(t, h.GetYaxis().GetXmin(), h.GetMaximum(), 'threshold  ', off=.3, line=True)

        # Efficiency plot
        format_histo(e_s, title='Efficiencies', x_tit='Threshold', y_tit='Efficiency',)
        l2 = Draw.make_legend(.3, .35, scale=1.5)
        [l2.AddEntry(p, p.GetName(), 'l') for i, p in enumerate([e_s, e_b])]
        self.Draw(e_s, show=show, leg=[l2, e_b], canvas=c.cd(2), lm=.12)

        # ROC Curve
        format_histo(roc, x_tit=e_b.GetName(), y_tit=e_s.GetName(), y_off=1.2)
        self.Draw(roc, show=show, gridx=True, gridy=True, canvas=c.cd(3))
        g = self.Draw.graph([e_b(t)], [e_s(t)], color=2, markersize=2, c=c.cd(3), draw_opt='p')
        Draw.tlatex(g.GetX()[0] - .04, g.GetY()[0], 'Working Point', color=2, size=.04, align=32)

        # Error Function plot
        format_histo(err, x_tit='Threshold', y_tit='1 / error', y_off=1.4)
        self.Draw(err, show=show, gridx=True, canvas=c.cd(4), prnt=show)
        self.Draw.save_plots('ErrorBucket', canvas=c.cd())
        return t

    def calc_pedestal(self, n_sigma):
        def f():
            t = self.Ana.info('generating pedestal cut for {dia} of run {run} ...'.format(run=self.Run.Number, dia=self.DUT.Name), endl=False)
            h0 = self.Draw.distribution(self.get_tree_vec(var=self.Ana.get_pedestal_name(), cut=self()), Bins.make(-100, 100, .5), show=False)
            fit = fit_fwhm(h0, show=False)
            self.Ana.add_to_info(t)
            return fit.Parameter(1), fit.Parameter(2)
        m, s = do_pickle(self.make_simple_pickle_path('Pedestal'), f)
        return [m - n_sigma * s, m + n_sigma * s]

    def calc_threshold(self, redo=False, show=True):
        def f():
            h = self.Ana.draw_signal_distribution(show=show, cut='', off_corr=False, evnt_corr=False)
            peaks = find_maxima(h, 3, sigma=2.5, sort_x=True)[:, 0]
            h.GetXaxis().SetRangeUser(peaks[0], peaks[-1])
            threshold = h.GetBinCenter(h.GetMinimumBin())
            Draw.vertical_line(threshold, 0, 1e7)
            return threshold
        return do_pickle(self.make_simple_pickle_path('Threshold'), f, redo=redo or show)

    def calc_cft(self, show=False, redo=False):
        def f():
            t = self.Ana.info('generating cft cut for {dia} of run {run} ...'.format(run=self.Run.Number, dia=self.Ana.DUT.Name), endl=False)
            cut = self.generate_custom(exclude=['cft'], prnt=False, name='cft_cut')
            h = self.Ana.Timing.draw_cft(cut=cut, show=show)
            fit = h.Fit('gaus', 'qs')
            self.Ana.add_to_info(t)
            return FitRes(fit)
        return do_pickle(self.make_simple_pickle_path('CFTFit'), f, redo=redo)

    def calc_timing_range(self, redo=False):
        def f():
            t = self.Ana.info('generating timing cut for {dia} of run {run} ...'.format(run=self.Run.Number, dia=self.Ana.DUT.Name), endl=False)
            cut = self.generate_custom(exclude=['timing'], prnt=False, name='timing_cut')
            t_correction = self.Ana.Timing.get_fine_correction(redo=redo)
            h = self.Ana.Timing.draw_peaks(show=False, cut=cut, fine_corr=t_correction != '0', prnt=False, redo=redo)
            fit = h.GetListOfFunctions()[1]
            x_min, x_max = self.Ana.SignalRegion * self.Ana.DigitiserBinWidth
            if fit.GetParameter(2) > 15 or x_min + 3 > fit.GetParameter(1) or x_max - 3 < fit.GetParameter(1):  # fit failed
                fit.SetParameter(1, h.GetBinCenter(h.GetMinimumBin()))
                fit.SetParameter(2, 15)
            self.Ana.add_to_info(t)
            self.Ana.info('Peak Timing: Mean: {0}, sigma: {1}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            return t_correction, fit
        return do_pickle(self.make_simple_pickle_path('TimingRange'), f, redo=redo)
    # endregion COMPUTE
    # ----------------------------------------

    # ----------------------------------------
    # region ANA
    def draw_bucket_histo(self, sig_var=None, show=True):
        x = self.get_tree_vec(var=choose(sig_var, self.Ana.get_signal_var()), cut=self.get_pre_bucket())
        xmax = quantile(x, .95)
        return self.Draw.distribution(x, Bins.get_pad_ph(max(1, xmax / 40)), 'Bucket Events', y_tit='Pulse Height [mV]', show=show, x_range=[-50, xmax])

    def print_bucket_events(self, prnt=True, redo=False):
        def f():
            info('getting number of bucket events for run {run} and {dia} ... '.format(run=self.Run.Number, dia=self.DUT.Name))
            return self.Ana.get_n_entries(self.get_pre_bucket()), self.Ana.get_n_entries(self.get_bucket())
        n_pre, n_new = do_pickle(self.make_simple_pickle_path('BucketEvents'), f, redo=redo)
        info('Pre Bucket: {: 5d} / {} = {:4.2f}%'.format(n_pre, self.Run.NEvents, n_pre / self.Run.NEvents * 100), prnt=prnt)
        info('New Bucket: {: 5d} / {} = {:4.2f}%'.format(n_new, self.Run.NEvents, n_new / self.Run.NEvents * 100), prnt=prnt)
        return n_pre, n_new, self.Run.NEvents

    def draw_bucket_map(self, pre=False, res=None, show=True):
        h = self.Ana.draw_hitmap(cut=self.get_pre_bucket() if pre else self.get_bucket(), show=show, res=res)
        format_histo(h, title='Bucket Hit Map', **{n: v for n, v in zip(['x_range', 'y_range'], ax_range(1, 1, .1, .1, h))})

    def draw_bucket_pedestal(self, show=True, corr=True, additional_cut=''):
        xbins = Bins.make(*array(self.Run.IntegralRegions[self.DUT.Number - 1]['signal_e']) * self.Ana.DigitiserBinWidth, bin_width=.5)
        cut_string = self.generate_custom(exclude='bucket') + TCut(additional_cut)
        self.Ana.draw_ph_peaktime(region='e', cut=cut_string, fine_corr=corr, prof=False, show=show, xbins=xbins, logz=True)

    def draw_bucket_distributions(self):
        cuts = [self.generate_custom(exclude='bucket', name='nobucket'), self.generate_custom(exclude='bucket', name='prebucket') + self.generate_pre_bucket()(), None]
        histos = [self.Ana.draw_signal_distribution(cut=cut, show=False) for cut in cuts]
        return self.Draw.stack(histos, 'Bucket Distributions', ['no bucket', 'pre bucket', 'new bucket'])

    def draw_bucket_means(self):
        thresh = self.calc_bucket_threshold()
        args = [('no bucket', 'bucket', None), ('pre bucket', 'bucket', self.generate_pre_bucket()), ('with threshold', None, None), ('thresh + 5', 'bucket', self.generate_bucket(thresh + 5)),
                ('thresh - 5', 'bucket', self.generate_bucket(thresh - 5))]
        cuts = [self.generate_custom(name=n, exclude=e, add=a) for n, e, a in args]
        self.draw_means(cuts=cuts, names=[cut.GetName() for cut in cuts])

    def compare_single_cuts(self, scale=True, redo=False):
        histos = [self.Ana.draw_signal_distribution(cut(), show=False, redo=redo, save=False) for cut in self.get_strings()]
        self.Draw.stack(histos, 'Single Cuts', self.get_names(), scale=scale)

    def compare_consecutive_cuts(self, scale=False, short=False, x_range=None, redo=False):
        cuts = self.get_consecutive(short)
        histos = [self.Ana.draw_signal_distribution(cut=cut, show=False, redo=redo, x_range=x_range, save=False) for cut in cuts.values()]
        self.Draw.stack(histos, 'Signal Distribution with Consecutive Cuts', cuts.keys(), scale)

    def draw_means(self, short=False, redo=False, show=True, cuts=None, names=None):
        y = [self.Ana.get_pulse_height(cut=cut, redo=redo) for cut in choose(cuts, self.get_consecutive(short).values())]
        g = self.Draw.graph(arange(len(y)), y, title='Pulse Height for Consecutive Cuts', y_tit='Pulse Height [mV]', draw_opt='ap', show=show, bm=.26, gridy=True, x_range=[-1, len(y)])
        set_bin_labels(g, choose(names, self.get_consecutive(short).keys()))
        update_canvas()

    def draw_signal_vs_signale(self, show=True):
        x, y = self.get_tree_vec(var=[self.Ana.get_signal_var(), self.Ana.get_signal_var(region='e')], cut=self.generate_custom(exclude='bucket'))
        names = self.Ana.SignalRegionName.title().replace('_', ' '), 'Signal E'
        self.Draw.histo_2d(x, y, self.Bins.get_pad_ph(mean_ph=mean(x)) * 2, '{} vs. {}'.format(*names), x_tit='{} [mV]'.format(names[0]), y_tit='{} [mV]'.format(names[1]), pal=53, stats=0,
                           show=show,
                           z_off=1.4, logz=True)

    def draw_cut_vars(self, normalise=False, consecutive=False):
        values = [self.Ana.get_ph_values(cut=(self.generate_threshold() + cut).Value) for cut in (self.get_consecutive().values() if consecutive else self.get_strings(with_raw=True))]
        v = array([[mean_sigma(lst)[0], quantile(lst, .5), get_fw_center(self.Draw.distribution(lst, self.Bins.get_pad_ph(mean_ph=mean(lst)), show=False))] for lst in values])
        if normalise:
            v = v / mean(v, axis=0) + arange(3) / 100
        names = self.get_names(with_raw=True)
        graphs = [self.Draw.make_tgrapherrors(arange(v.shape[0]), v[:, i]) for i in arange(v.shape[1])]
        return self.Draw.multigraph(graphs, '{}Cut Variables'.format('Consecutive ' if consecutive else ''), ['mean', 'median', 'fwc'], names, x_tit='Cut Name',
                                    y_tit='{}Pulse Height [mV]'.format('Norm ' if normalise else ''), gridy=True, lm=.12)
    # endregion ANA
    # ----------------------------------------

    @staticmethod
    def fit_bucket(h, show=True):
        maxima = concatenate(find_maxima(h, n=3, sigma=2.5, sort_x=True))  # n one larger than the expected otherwise buffer full
        if maxima.size < 3 or maxima[1] < 20 or maxima[1] > 1e10:
            return  # didn't find pedestal peak!
        x1, y1, x2, y2 = maxima
        d = x2 - x1
        fit = TF1('fit', 'gaus(0) + gaus(3) + gaus(6)', h.GetXaxis().GetXmin(), x2 + 5)
        fit.SetParameters(y2, x2, 10, y1, x1, 3, min(y1, y2) / 4, x1 + d / 4, 5)
        fit.SetParLimits(1, x2 - 5, x2 + 5)         # signal mean
        fit.SetParLimits(3, 1, y1 * 2)              # pedestal height
        fit.SetParLimits(4, x1 - 10, x1 + 10)       # pedestal mean
        fit.SetParLimits(6, 1, min(y1, y2) / 2)     # guard ring height
        fit.SetParLimits(7, x1, x1 + d / 2)         # guard ring mean
        h.Fit(fit, 'qs{0}'.format('' if show else '0'), '', -50, x2 + 5)
        return fit
