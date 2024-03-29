# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TCut, TF1
from src.cut import Cut, CutString, Bins, low_rate, high_rate
from src.telescope import Telescope
from plotting.draw import *
from helpers.utils import *
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
        return self.generate_custom(exclude='trigger phase', add=self.generate_bucket() + self.generate_b2(5), name='bful', prnt=False)

    @property
    def bucket(self):
        return self.generate_custom(include=['bucket', 'bucket2'])

    def inv_bucket(self, all_cuts=False):
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

    def no_bucket(self, cut=None):
        return self.exclude('bucket', 'bucket2', name='nobuc') if cut is None else self(cut)

    @property
    def has_tp(self):
        return self.has('trigger phase')
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE
    def generate_dut(self):
        try:
            # -- WAVEFORM --
            self.CutStrings.register(self.generate_pulser(), 2)
            self.CutStrings.register(self.generate_saturated(), level=3)

            # -- FIDUCIAL --
            self.CutStrings.register(self.generate_fiducial(), 23)

            # -- SIGNAL --
            self.CutStrings.register(self.generate_pedestal_bucket(), 29)
            self.CutStrings.register(self.generate_pedestal_sigma(), 30)
            # self.CutStrings.register(self.generate_threshold(), 31)
            # self.CutStrings.register(self.generate_timing(), 35)
            # self.CutStrings.register(self.generate_cft(), 36)

            # --BUCKET --
            if self.Run.is_volt_scan() or self.get_raw_snr() > 5 or abs(self.Ana.DUT.Bias) < 10:
                self.CutStrings.register(self.generate_bucket(), 91)
                if not self.Run.is_volt_scan() or self.Ana.DUT.Bias < 10:
                    self.CutStrings.register(self.generate_b2(n_sigma=5), 92)
            else:
                self.CutStrings.register(self.generate_trigger_phase(), 91)

            # --SIGNAL DROPS--
            self.update('event range', self.generate_event_range(None, self.find_signal_drops()).Value)  # update event range if drop is found
        except Exception as err:
            self.show_contributions(redo=True)
            critical(err)

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
        wide_range = diff(self.Ana.get_region())[0] * self.Run.DigitiserBinWidth > 1.5 * self.BunchSpacing
        return CutString() if wide_range else CutString('ped bucket', f'!ped_bucket[{self.DUT.Number - 1}]', 'pedestal bucket events')

    def generate_bucket(self):
        """exclude events with a peak in bucket 2 and no peak in the signal region"""
        return CutString('bucket', f'!bucket[{self.DUT.Number - 1}]', 'bucket events')

    def generate_b2(self, n_sigma):
        fit, b1, b2, noise = self.get_b2_fit(), self.Ana.get_raw_signal_var(), self.Ana.get_b2_var(), self.calc_pedestal_()[1] * n_sigma
        string = f'{b2} < {fit[0].n} + {fit[1].n} * {b1} + {fit[2].n} * pow({b1}, 2) + {noise}'
        descr = f'signals above 3 sigma in bucket 2 with pedestal shape {fit[0].n:.2f} + {fit[1].n:.2f} * x + {fit[2].n:.4f} * x²'
        return CutString('bucket2', string, descr)

    def generate_trigger_phase(self):
        tps = self.get_trigger_phases()
        return CutString('trigger phase', ' && '.join(f'{Telescope.tp_var(tree=self.Tree)} != {i}' for i in tps), f'trigger phase != {tps}')

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
    def find_signal_drops(self, thresh=.6, pol=None, _redo=False):
        return super(PadCut, self).find_signal_drops(thresh, choose(pol, self.Ana.get_polarity()), _redo=_redo)

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

    @save_pickle('BeamStops', suf_args='all')
    def find_beam_interruptions(self, bw=100, thresh=.2):
        """ Looking for the beam interruptions by investigating the pulser rate. """
        x, y = self.get_tree_vec(var=['Entry$', 'pulser'], dtype='i4')
        rates, x_bins, y_bins = histogram2d(x, y, bins=[arange(0, x.size, bw, dtype=int), 2])
        rates = rates[:, 1] / bw
        thresh = min(1, mean(rates) + thresh)
        events = x_bins[:-1][rates > thresh] + bw / 2
        not_connected = where(concatenate([[False], events[:-1] != events[1:] - bw]))[0]  # find the events where the previous event is not related to the event (more than a bin width away)
        events = split(events, not_connected)  # events grouped into connecting events
        interruptions = [(ev[0], ev[0]) if ev.size == 1 else (ev[0], ev[-1]) for ev in events] if events[0].size else []
        return array(interruptions, 'i4')

    @save_pickle('B2Fit', print_dur=True, low_rate=True)
    @low_rate
    def get_b2_fit(self, _redo=False):
        """ fit the b2/b1 profile of the lowest rate run with pol2 """
        h = self.Ana.draw_b2_profile(self.generate_custom(include=['pulser', 'ped sigma', 'event range', 'fiducial'], prnt=False), show=False)
        xmax = h.GetBinCenter(int(h.GetNbinsX() - argmax(bins.entries(h)[::-1] > 10)))  # find first bin with more than 10 entries from the back
        fit = FitRes(h.Fit('pol2', 'qs', '', 10, xmax))
        # fit.Pars[0] -= self.Ana.get_polarity() * self.calc_pedestal_[0]  # subtract baseline -> THIS IS WRONG!
        return fit

    def get_fb2(self, redo=False, **kw):
        fit = self.get_b2_fit(_redo=redo)
        return Draw.make_f(None, f'{fit[0].n} + {fit[1].n} * x + {fit[2].n} * pow(x, 2)', -50, 500, **kw)

    @save_pickle('TP', print_dur=True, high_rate=True)
    @high_rate
    def get_trigger_phase(self, show=False, _redo=False):
        h = self.Draw.distribution(self.get_tree_vec(Telescope.tp_var(tree=self.Tree), cut=self.inv_bucket(all_cuts=True)), w=1, x0=0, show=show)
        return argmax(h) - 1  # subtract zero bin

    def get_trigger_phases(self, n=1):
        return (arange(2 * n + 1) - 1 + self.get_trigger_phase()) % 10  # select also the [n] trigger phases around the MPV

    @save_pickle('Pedestal', print_dur=True)
    def calc_pedestal_(self, _redo=False):
        x = self.get_tree_vec(var=self.Ana.get_pedestal_name(), cut=self())
        h0 = self.Draw.distribution(x, Bins.make(-100, 100, max(.1, Bins.find_width(x))), show=False)
        return fit_fwhm(h0, show=False).get_pars(err=False)[1:]

    def calc_pedestal(self, n_sigma, redo=False):
        m, s = self.calc_pedestal_(_redo=redo)
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
    def draw_cut_vars(self, normalise=False, consecutive=False):
        values = [self.Ana.get_ph_values(cut=(self.generate_threshold() + cut).Value) for cut in (self.get_consecutive().values() if consecutive else self.get_strings(with_raw=True))]
        v = array([[mean_sigma(lst)[0], quantile(lst, .5), get_fw_center(self.Draw.distribution(lst, Bins.get_pad_ph(), show=False))] for lst in values])
        if normalise:
            v = v / mean(v, axis=0) + arange(3) / 100
        names = self.get_names(with_raw=True)
        graphs = [self.Draw.make_tgraph(arange(v.shape[0]), v[:, i]) for i in arange(v.shape[1])]
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
