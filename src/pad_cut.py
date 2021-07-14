# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TCut, TF1
from src.cut import Cut, CutString, Bins, low_rate, high_rate
from helpers.draw import *
from numpy import quantile, histogram2d, argmax


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

    @save_pickle('SNR', print_dur=True, low_rate=True)
    @low_rate
    def get_raw_snr(self, _redo=False):
        return self.get_raw_pulse_height() / self.get_raw_noise()

    def get_full_bucket(self):
        return self.generate_custom(exclude='trigger phase', add=self.generate_bucket() + self.generate_b2(4), name='bful', prnt=False)

    def get_bucket(self, all_cuts=False):
        cut = self.generate_custom(exclude=['bucket', 'bucket2'], prnt=False) if all_cuts else self.generate_custom(include=['pulser', 'ped sigma', 'event range'], prnt=False)
        return Cut.make('!bucket', cut + self.generate_bucket().invert())

    def get_tp(self):
        return self.generate_custom(exclude=['bucket', 'bucket2'], add=self.generate_trigger_phase(), prnt=False, name='tp')

    def get_pulser(self, beam_on=True):
        cut = self.generate_custom(include=['ped sigma', 'event range'], prnt=False) + Cut.invert(self.get('pulser'))
        cut += self.get('beam stops', warn=False) if beam_on else Cut.invert(self.generate_jump_cut())
        return CutString('PulserBeam{}'.format('On' if beam_on else 'Off'), cut)

    def get_timing(self):
        return CutString('timing', self.generate_custom(exclude='timing', prnt=False), 'All cuts expect timing.')

    def get_ped_sigma(self, sigma=None):
        return choose(sigma, self.Config.get_value('CUT', 'pedestal sigma', dtype=int))
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE
    def generate_dut(self):

        # -- WAVEFORM --
        self.CutStrings.register(self.generate_pulser(), 2)
        self.CutStrings.register(self.generate_saturated(), level=3)

        # -- SIGNAL --
        self.CutStrings.register(self.generate_pedestal_bucket(), 29)
        self.CutStrings.register(self.generate_pedestal_sigma(), 30)
        # self.CutStrings.register(self.generate_threshold(), 31)
        # self.CutStrings.register(self.generate_timing(), 35)
        # self.CutStrings.register(self.generate_cft(), 36)

        # -- FIDUCIAL --
        self.CutStrings.register(self.generate_fiducial(), 23)

        # --BUCKET --
        if self.get_raw_snr() > 5:
            self.CutStrings.register(self.generate_bucket(), 91)
            if 'voltage' not in self.Run.Info['runtype']:
                self.CutStrings.register(self.generate_b2(n_sigma=4), 92)
        else:
            self.CutStrings.register(self.generate_trigger_phase(), 91)

        # --SIGNAL DROPS--
        self.update('event range', self.generate_event_range(None, self.find_zero_ph_event()).Value)  # update event range if drop is found

    def generate_saturated(self):
        cut_string = '!is_saturated[{ch}]'.format(ch=self.Channel)
        description = 'exclude {:.2f}% saturated events'.format(100. * self.find_n_saturated(Cut.invert(cut_string)) / self.Run.NEvents)
        return CutString('saturated', cut_string, description)

    def generate_pulser(self):
        return CutString('pulser', '!pulser', 'exclude {:.1f}% pulser events'.format(100. * self.find_n_pulser('pulser') / self.Run.NEvents))

    def generate_pedestal_sigma(self, sigma=None):
        sigma = self.get_ped_sigma(sigma)
        ped_range = self.calc_pedestal(sigma)
        description = 'sigma <= {} -> [{:1.1f}mV, {:1.1f}mV]'.format(sigma, *ped_range)
        return CutString('ped sigma', '{ped} > {} && {ped} < {}'.format(ped=self.Ana.get_pedestal_name(), *ped_range), description)

    def generate_threshold(self):
        thresh = self.calc_threshold(show=False)
        return CutString('threshold', '{} > {}'.format(self.Ana.get_raw_signal_var(), thresh), thresh)

    def generate_pedestal_bucket(self):
        """exclude events with a peak in bucket -1 (the pre-pedstal bucket)"""
        return CutString('ped bucket', f'!ped_bucket[{self.DUT.Number - 1}]', 'pedestal bucket events')

    def generate_bucket(self):
        """exclude events with a peak in bucket 2 and no peak in the signal region"""
        return CutString('bucket', f'!bucket[{self.DUT.Number - 1}]', 'bucket events')

    def generate_b2(self, n_sigma):
        fit, b1, b2, noise = self.get_b2_fit(), self.Ana.get_raw_signal_var(), self.Ana.get_b2_var(), self.calc_pedestal(n_sigma)[1]
        string = f'{b2} < {fit[0].n} + {fit[1].n} * {b1} +  {fit[2].n} * pow({b1}, 2) + {noise}'
        descr = f'signals above 3 sigma in bucket 2 with pedestal shape {fit[0].n:.2f} + {fit[1].n:.2f} * x + {fit[2].n:.4f} * x²'
        return CutString('bucket2', string, descr)

    def generate_trigger_phase(self):
        tps = self.get_trigger_phases()
        return CutString('trigger phase', ' && '.join(f'trigger_phase != {i}' for i in tps), f'trigger phase != {tps}')

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

    def find_fiducial(self, thresh=.93, show=True):
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

    def get_b2_fit(self, redo=False):
        """ fit the b2/b1 profile of the lowest rate run with pol2 """
        def f():
            if not self.Run.get_low_rate_run() == self.Run.Number:
                from src.pad_analysis import PadAnalysis
                return PadAnalysis(self.Run.get_low_rate_run(), self.DUT.Number, self.TCString, prnt=False).Cut.get_b2_fit(redo)
            t = self.info(f'generating bucket cut for {self.Run} and {self.DUT}', endl=False)
            h = self.Ana.draw_bucket_profile(self(), show=False)
            xmax = h.GetBinCenter(int(h.GetNbinsX() - argmax(get_h_entries(h)[::-1] > 10)))  # find first bin with more than 10 entries from the back
            fit, ped = FitRes(h.Fit('pol2', 'qs', '', 10, xmax)), self.calc_pedestal(1)
            fit.Pars[0] -= self.Ana.get_polarity() * sum(ped) / 2  # subtract baseline
            self.add_to_info(t)
            return fit
        return do_pickle(self.make_simple_pickle_path('B2Fit', run=self.Run.get_low_rate_run()), f, redo=redo)

    @save_pickle('TP', print_dur=True, high_rate=True)
    @high_rate
    def get_trigger_phase(self, show=False, _redo=False):
        h = self.Draw.distribution(self.get_tree_vec('trigger_phase', cut=self.get_bucket() + self.generate_fiducial()()), make_bins(11), show=show)
        return argmax(get_hist_vec(h))

    def get_trigger_phases(self, n=1):
        return (arange(2 * n + 1) - 1 + self.get_trigger_phase()) % 10  # select also the [n] trigger phases around the MPV

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
    def compare_single_cuts(self, scale=True, redo=False):
        histos = [self.Ana.draw_signal_distribution(cut(), show=False, redo=redo, save=False) for cut in self.get_strings()]
        self.Draw.stack(histos, 'Single Cuts', self.get_names(), scale=scale)

    def compare_consecutive_cuts(self, scale=False, short=False, x_range=None, redo=False):
        cuts = self.get_consecutive(short)
        histos = [self.Ana.draw_signal_distribution(cut=cut, show=False, redo=redo, x_range=x_range, save=False) for cut in cuts.values()]
        self.Draw.stack(histos, 'Signal Distribution with Consecutive Cuts', cuts.keys(), scale)

    @quiet
    def draw_means(self, short=False, cuts=None, names=None, normalise=True, redo=False, **kwargs):
        cuts, labels = choose(cuts, list(self.get_consecutive(short).values())), choose(names, self.get_consecutive(short).keys())
        self.Ana.PBar.start(len(cuts), counter=True) if redo or not file_exists(self.make_simple_pickle_path('Fit', f'{cuts[-1].GetName()}_1', 'PH')) else do_nothing()
        x, y = arange(len(cuts)), array([self.Ana.get_pulse_height(cut=cut, redo=redo) for cut in cuts])
        y /= y[-1] if normalise else 1
        return self.Draw.graph(x, y, title='Pulse Height for Consecutive Cuts', y_tit='Pulse Height [mV]', **prep_kw(kwargs, draw_opt='ap', gridy=True, x_range=[-1, len(y)], bin_labels=labels))

    def draw_cut_vars(self, normalise=False, consecutive=False):
        values = [self.Ana.get_ph_values(cut=(self.generate_threshold() + cut).Value) for cut in (self.get_consecutive().values() if consecutive else self.get_strings(with_raw=True))]
        v = array([[mean_sigma(lst)[0], quantile(lst, .5), get_fw_center(self.Draw.distribution(lst, Bins.get_pad_ph(), show=False))] for lst in values])
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
