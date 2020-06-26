# --------------------------------------------------------
#       cut sub class to handle all the cut strings for the DUTs with digitiser
# created in 2015 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from utils import *
from ROOT import TCut, TH1F, TF1, TSpectrum
from cut import Cut, CutString, loads, invert
from draw import format_histo, fit_bucket


class PadCut(Cut):
    """ Generates all cut strings which correspond to a single DUT with digitiser. """

    def __init__(self, analysis):
        Cut.__init__(self, analysis)
        self.Channel = self.Analysis.Channel
        self.DUT = analysis.DUT

        self.update_config()
        self.generate_dut()
        self.ConsecutiveCuts = self.generate_consecutive()

    @classmethod
    def from_parent(cls, parent):
        return cls(parent.Analysis)

    def update_config(self):
        self.CutConfig['pedestal_sigma'] = self.Config.getint('CUT', 'pedestal sigma')
        self.CutConfig['fiducial'] = self.load_fiducial()
        self.CutConfig['threshold'] = self.load_dut_config('threshold', store_true=True)
        self.CutConfig['trigger_cell'] = loads(self.Config.get('CUT', 'trigger cell')) if self.Config.has_option('CUT', 'trigger cell') else None

    # ----------------------------------------
    # region GET
    def get_raw_pedestal(self):
        n = self.Analysis.Tree.Draw(self.Analysis.get_signal_name(sig_type='pedestal'), self(), 'goff')
        return make_ufloat(mean_sigma(self.Run.get_root_vec(n)))

    def get_raw_snr(self):
        ped = self.get_raw_pedestal()
        return self.get_raw_pulse_height() / ped.s, ped

    def get_bucket(self):
        return CutString('bucket', self.generate_custom(include=['pulser', 'fiducial', 'event_range'], prnt=False) + invert(self.generate_pre_bucket()))()

    def get_pulser(self, beam_on=True):
        cut = self.generate_custom(include=['ped_sigma', 'event_range'], prnt=False) + invert(self.get('pulser'))
        cut += self.get('beam_interruptions') if beam_on else invert(self.generate_jump_cut())
        return CutString('PulserBeam{}'.format('On' if beam_on else 'Off'), cut)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE
    def generate_dut(self):

        # -- WAVEFORM --
        self.CutStrings.register(self.generate_saturated(), level=2)
        self.CutStrings.register(self.generate_pulser(), 3)

        # -- SIGNAL --
        self.CutStrings.register(self.generate_pedestal_sigma(), 30)
        self.CutStrings.register(self.generate_threshold(), 31)
        # self.CutStrings.register(self.generate_timing(), 35)
        # self.CutStrings.register(self.generate_cft(), 36)

        # -- FIDUCIAL --
        self.CutStrings.register(self.generate_fiducial(), 90)

        # --BUCKET --
        self.CutStrings.register(self.generate_bucket(), 36)

        # --SIGNAL DROPS--
        self.update('event_range', self.generate_event_range(None, self.find_zero_ph_event()).Value)  # update event range if drop is found

    def generate_saturated(self):
        cut_string = '!is_saturated[{ch}]'.format(ch=self.Channel)
        description = 'exclude {:.1f}% saturated events'.format(100. * self.find_n_saturated(invert(cut_string)) / self.Run.NEntries)
        return CutString('saturated', cut_string, description)

    def generate_pulser(self):
        return CutString('pulser', '!pulser', 'exclude {:.1f}% pulser events'.format(100. * self.find_n_pulser('pulser') / self.Run.NEntries))

    def generate_pedestal_sigma(self, sigma=None):
        sigma = self.CutConfig['pedestal_sigma'] if sigma is None else sigma
        ped_range = self.__calc_pedestal_range(sigma)
        description = 'sigma <= {} -> [{:1.1f}mV, {:1.1f}mV]'.format(sigma, *ped_range)
        return CutString('ped_sigma', '{ped}>{}&&{ped}<{}'.format(ped=self.Analysis.PedestalName, *ped_range), description)

    def generate_threshold(self):
        description = self.calc_threshold(show=False)
        return CutString('threshold', '{sig}>{thresh}'.format(sig=self.Analysis.SignalName, thresh=self.calc_threshold(show=False)) if self.CutConfig['threshold'] else '', description)

    def generate_pre_bucket(self):
        try:
            sig2 = self.Analysis.get_signal_name('e', 2)
            string = '{sig2}=={sig1}'.format(sig2=sig2, sig1=self.Analysis.SignalName)
            return string
        except ValueError as err:
            print err
            return ''

    def generate_bucket(self, threshold=None):
        sig = self.Analysis.SignalName
        threshold = self.calc_signal_threshold(show=False) if threshold is None else threshold
        string = '!(!({old_buck}) && ({sig} < {thres}))'.format(sig=sig, thres=threshold, old_buck=self.generate_pre_bucket())
        description = 'bucket events with threshold < {:.1f}mV'.format(threshold)
        return CutString('bucket', string if self.generate_pre_bucket() else '', description)

    def generate_cft(self, n_sigma=3):
        fit = self.fit_cft()
        m, s = fit.Parameter(1), fit.Parameter(2)
        description = '{:1.1f}ns < constant fraction time < {:.1f}ns'.format(m - n_sigma * s, m + n_sigma * s)
        var = self.Analysis.Timing.get_cft_name()
        return CutString('cft', '{} < {v} && {v} < {}'.format(m - n_sigma * s, m + n_sigma * s, v=var), description)

    def generate_timing(self, n_sigma=3):
        t_correction, fit = self.calc_timing_range()
        if fit is None:
            return TCut('')
        corrected_time = '{peak} - {t_corr}'.format(peak=self.Analysis.Timing.get_peak_name(corr=True), t_corr=t_correction)
        string = 'TMath::Abs({cor_t} - {mp}) / {sigma} < {n_sigma}'.format(cor_t=corrected_time, mp=fit.GetParameter(1), sigma=fit.GetParameter(2), n_sigma=n_sigma)
        description = '{:1.1f}ns < peak timing < {:.1f}ns'.format(fit.GetParameter(1) - fit.GetParameter(2), fit.GetParameter(1) + fit.GetParameter(2))
        return CutString('timing', string, description)

    def generate_trigger_cell(self, low=None, high=None):
        tc_range = [low, high] if low is not None else self.CutConfig['trigger_cell']
        return CutString('trigger_cell', 'trigger_cell < {} && trigger_cell >= {}'.format(*tc_range) if tc_range else '', '{} < trigger celll < {}'.format(*tc_range) if tc_range else '')
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region COMPUTE
    def find_n_pulser(self, cut, redo=False):
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'NPulser', self.Run.Number)
        return int(do_pickle(pickle_path, self.Analysis.Tree.GetEntries, None, redo, str(cut)))

    def find_n_saturated(self, cut, redo=False):
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'NSaturated', self.Run.Number)
        return int(do_pickle(pickle_path, self.Analysis.Tree.GetEntries, None, redo, str(cut)))

    def find_fid_cut(self, thresh=.93, show=True):
        h = self.Analysis.draw_signal_map(show=False)
        px = h.ProjectionX()
        format_histo(px, title='Projection X of the Signal Map', y_tit='Number of Entries', y_off=1.5)
        self.Analysis.draw_histo(px, lm=.12, show=show)
        py = h.ProjectionY()
        return '"{}": [{}]'.format(self.Analysis.DUT.Name, ', '.join('{:0.3f}'.format(i) for i in self.find_fid_margins(px, thresh) + self.find_fid_margins(py, thresh)))

    @staticmethod
    def find_fid_margins(proj, thresh):
        thresh = proj.GetMaximum() * thresh
        xbin1, xbin2 = proj.FindFirstBinAbove(thresh), proj.FindLastBinAbove(thresh)
        f1 = interpolate_two_points(proj.GetBinCenter(xbin1), proj.GetBinContent(xbin1), proj.GetBinCenter(xbin1 - 1), proj.GetBinContent(xbin1 - 1))
        f2 = interpolate_two_points(proj.GetBinCenter(xbin2), proj.GetBinContent(xbin2), proj.GetBinCenter(xbin2 + 1), proj.GetBinContent(xbin2 + 1))
        return [f1.GetX(thresh) / 10, f2.GetX(thresh) / 10]

    def calc_signal_threshold(self, use_bg=False, show=True, show_all=False):
        run = self.HighRateRun if self.HighRateRun is not None else self.Run.Number
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'SignalThreshold', run, self.DUT.Number)
        show = False if show_all else show

        def f():
            t = self.Analysis.info('Calculating signal threshold for bucket cut of run {run} and {d} ...'.format(run=self.Run.Number, d=self.DUT.Name), next_line=False)
            h = TH1F('h', 'Bucket Cut', 200, -50, 150)
            self.Analysis.Tree.Draw('{name}>>h'.format(name=self.Analysis.SignalName), self.get_bucket(), 'goff')
            format_histo(h, x_tit='Pulse Height [mV]', y_tit='Entries', y_off=1.8, stats=0, fill_color=self.Analysis.FillColor)
            if h.GetEntries() / self.Run.NEntries < .01:
                log_info('Not enough bucket events ({:.1f}%)'.format(h.GetEntries() / self.Run.NEntries))
                self.Analysis.add_to_info(t)
                return -30
            snr, ped = self.get_raw_snr()
            if snr < 4:  # impossible to separate for the small signal to noise ratios...
                warning('Signal to noise ration is too low... -> taking noise based value for bucket cut!')
                return ped.n + 2 * ped.s
            # extract fit functions
            fit = fit_bucket(h)
            if fit is None or any([abs(fit.GetParameter(i)) < 20 for i in [0, 3]]) or fit.GetParameter(1) < fit.GetParameter(4) or fit.GetParameter(1) > 500:
                warning('bucket cut fit failed')
                self.Analysis.draw_histo(h, show=show)
                self.Analysis.add_to_info(t)
                return -30
            sig_fit = TF1('f1', 'gaus', -50, 300)
            sig_fit.SetParameters(fit.GetParameters())
            ped_fit = TF1('f2', 'gaus(0) + gaus(3)', -50, 300)
            ped_fit.SetParameters(*[fit.GetParameter(i) for i in xrange(3, 9)])
            set_root_output(True)

            # real data distribution without pedestal fit
            signal = deepcopy(h)
            signal.Add(ped_fit, -1)

            gr1 = self.Analysis.make_tgrapherrors('gr1', '#varepsilon_{bg}', marker_size=0.2)
            gr2 = self.Analysis.make_tgrapherrors('gr2', '#varepsilon_{sig}', marker_size=0.2, color=2)
            gr3 = self.Analysis.make_tgrapherrors('gr3', 'ROC Curve', marker_size=0.2)
            gr4 = self.Analysis.make_tgrapherrors('gr4', 'Signal Error', marker_size=0.2)
            xs = arange(-30, sig_fit.GetParameter(1), .1)
            errors = {}
            for i, x in enumerate(xs):
                ped = ped_fit.Integral(-50, x) / ped_fit.Integral(-50, 500)
                sig = 1 - sig_fit.Integral(-50, x) / signal.Integral()
                err = ped_fit.Integral(-50, x) / (sqrt(sig_fit.Integral(-50, x) + ped_fit.Integral(-50, x)))
                s, bg = signal.Integral(h.FindBin(x), signal.GetNbinsX() - 1), ped_fit.Integral(x, 200)
                err1 = s / sqrt(s + bg)
                errors[err1 if not use_bg else err] = x
                gr1.SetPoint(i, x, ped)
                gr2.SetPoint(i, x, sig)
                gr3.SetPoint(i, sig, ped)
                gr4.SetPoint(i, x, err1 if not use_bg else err)
            if len(errors) == 0:
                print ValueError('errors has a length of 0')
                self.Analysis.add_to_info(t)
                return -30
            max_err = max(errors.items())[1]
            c = None
            if show_all:
                set_root_output(True)
                c = self.Analysis.make_canvas('c_all', 'Signal Threshold Overview', divide=(2, 2))
            # Bucket cut plot
            self.Analysis.draw_histo(h, '', show or show_all, lm=.135, canvas=c.cd(1) if show_all else None)
            self.Analysis.draw_y_axis(max_err, h.GetYaxis().GetXmin(), h.GetMaximum(), 'threshold  ', off=.3, line=True)
            ped_fit.SetLineStyle(2)
            ped_fit.Draw('same')
            sig_fit.SetLineColor(4)
            sig_fit.SetLineStyle(3)
            sig_fit.Draw('same')
            self.Analysis.save_plots('BucketCut', canvas=c.cd(1) if show_all else get_last_canvas(), prnt=show)

            # Efficiency plot
            format_histo(gr1, title='Efficiencies', x_tit='Threshold', y_tit='Efficiency', markersize=.2)
            l2 = self.Analysis.make_legend(.78, .3)
            tits = ['#varepsilon_{bg}', gr2.GetTitle()]
            [l2.AddEntry(p, tits[i], 'l') for i, p in enumerate([gr1, gr2])]
            self.Analysis.draw_histo(gr1, '', show_all, draw_opt='apl', leg=l2, canvas=c.cd(2) if show_all else None)
            self.Analysis.draw_histo(gr2, show=show_all, draw_opt='same', canvas=c.cd(2) if show_all else get_last_canvas())
            self.Analysis.save_plots('Efficiencies', canvas=c.cd(2) if show_all else get_last_canvas(), prnt=show)

            # ROC Curve
            format_histo(gr3, y_tit='background fraction', x_tit='excluded signal fraction', markersize=0.2, y_off=1.2)
            self.Analysis.draw_histo(gr3, '', show_all, gridx=True, gridy=True, draw_opt='apl', canvas=c.cd(3) if show_all else None)
            p = self.Analysis.make_tgrapherrors('gr', 'working point', color=2)
            p.SetPoint(0, 1 - sig_fit.Integral(-50, max_err) / signal.Integral(), ped_fit.Integral(-50, max_err) / ped_fit.Integral(-50, 200))
            sleep(.1)
            latex = self.Analysis.draw_tlatex(p.GetX()[0], p.GetY()[0] + .01, 'Working Point', color=2, size=.04)
            p.GetListOfFunctions().Add(latex)
            self.Analysis.draw_histo(p, show=show_all, canvas=c.cd(3) if show_all else get_last_canvas(), draw_opt='p')
            self.Analysis.save_plots('ROC_Curve', canvas=c.cd(3) if show_all else get_last_canvas(), prnt=show)

            # Error Function plot
            format_histo(gr4, x_tit='Threshold', y_tit='1 / error', y_off=1.4)
            self.Analysis.save_histo(gr4, 'ErrorFunction', show_all, gridx=True, draw_opt='al', canvas=c.cd(4) if show_all else None, prnt=show)

            self.Analysis.Objects.append([sig_fit, ped_fit, gr2, c])

            self.Analysis.add_to_info(t)
            return max_err

        return do_pickle(pickle_path, f, redo=show or show_all)

    def __calc_pedestal_range(self, sigma_range):
        picklepath = self.Analysis.make_pickle_path('Pedestal', 'Cut', self.Run.Number, self.Channel)

        def func():
            t = self.Analysis.info('generating pedestal cut for {dia} of run {run} ...'.format(run=self.Run.Number, dia=self.Analysis.DUT.Name), next_line=False)
            h1 = TH1F('h_pdc', 'Pedestal Distribution', 600, -150, 150)
            self.Analysis.Tree.Draw('{name}>>h_pdc'.format(name=self.Analysis.PedestalName), '', 'goff')
            fit_pars = fit_fwhm(h1, do_fwhm=True, draw=False)
            self.Analysis.add_to_info(t)
            return fit_pars

        fit = do_pickle(picklepath, func)
        sigma = fit.Parameter(2)
        mean_ = fit.Parameter(1)
        self.PedestalFit = fit
        return [mean_ - sigma_range * sigma, mean_ + sigma_range * sigma]

    def calc_threshold(self, show=True):
        pickle_path = self.Analysis.make_pickle_path('Cuts', 'Threshold', self.Run.Number, self.Channel)

        def func():
            self.Analysis.Tree.Draw(self.Analysis.SignalName, '', 'goff', 5000)
            xvals = sorted([self.Analysis.Tree.GetV1()[i] for i in xrange(5000)])
            x_range = [xvals[0] - 5, xvals[-5]]
            h = self.Analysis.draw_signal_distribution(show=show, cut=self.generate_fiducial()(), x_range=x_range, bin_width=1)
            s = TSpectrum(3)
            s.Search(h)
            peaks = [s.GetPositionX()[i] for i in xrange(s.GetNPeaks())]
            h.GetXaxis().SetRangeUser(peaks[0], peaks[-1])
            x_start = h.GetBinCenter(h.GetMinimumBin())
            h.GetXaxis().UnZoom()
            fit = TF1('fit', 'landau', x_start, h.GetXaxis().GetXmax())
            h.Fit(fit, 'q{0}'.format(0 if not show else ''), '', x_start, h.GetXaxis().GetXmax())
            return fit.GetX(.1, 0, peaks[-1])

        threshold = func() if show else None
        return do_pickle(pickle_path, func, threshold)

    def fit_cft(self, show=False, redo=False):
        def f():
            t = self.Analysis.info('generating cft cut for {dia} of run {run} ...'.format(run=self.Run.Number, dia=self.Analysis.DUT.Name), next_line=False)
            cut = self.generate_custom(exclude=['cft'], prnt=False, name='cft_cut')
            h = self.Analysis.Timing.draw_cft(cut=cut, show=show)
            fit = h.Fit('gaus', 'qs')
            self.Analysis.add_to_info(t)
            return FitRes(fit)
        return do_pickle(self.Analysis.make_simple_pickle_path('CFTFit', sub_dir='Cuts'), f, redo=redo)

    def calc_timing_range(self, redo=False):
        def f():
            t = self.Analysis.info('generating timing cut for {dia} of run {run} ...'.format(run=self.Run.Number, dia=self.Analysis.DUT.Name), next_line=False)
            cut = self.generate_custom(exclude=['timing'], prnt=False, name='timing_cut')
            t_correction = self.Analysis.Timing.calc_fine_correction(redo=redo)
            h = self.Analysis.Timing.draw_peaks(show=False, cut=cut, fine_corr=t_correction != '0', prnt=False, redo=redo)
            fit = h.GetListOfFunctions()[1]
            x_min, x_max = self.Analysis.SignalRegion * self.Analysis.DigitiserBinWidth
            if fit.GetParameter(2) > 15 or x_min + 3 > fit.GetParameter(1) or x_max - 3 < fit.GetParameter(1):  # fit failed
                fit.SetParameter(1, h.GetBinCenter(h.GetMinimumBin()))
                fit.SetParameter(2, 15)
            self.Analysis.add_to_info(t)
            self.Analysis.info('Peak Timing: Mean: {0}, sigma: {1}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            return t_correction, fit

        return do_pickle(self.Analysis.make_simple_pickle_path('TimingRange', sub_dir='Cuts'), f, redo=redo)
    # endregion COMPUTE
    # ----------------------------------------
