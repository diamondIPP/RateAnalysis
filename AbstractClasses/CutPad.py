from Cut import Cut
from Extrema import Extrema2D
from functools import partial
from ROOT import TCut, TH1F, TF1, TCanvas, TCutG, TSpectrum
from Utils import *
from json import loads
from numpy import array
from ConfigParser import NoOptionError

__author__ = 'micha'


# ==============================================
# MAIN CLASS
# ==============================================
class CutPad(Cut):
    """The ChannelCut contains all cut settings which corresponds to a single diamond in a single run. """

    def __init__(self, analysis, channel=0):
        Cut.__init__(self, analysis, skip=True)
        self.__dict__.update(analysis.Cut.__dict__)
        self.channel = channel
        self.RunNumber = self.analysis.RunNumber
        self.DiamondName = analysis.DiamondName
        self.DiamondNumber = analysis.DiamondNumber

        self.save_dir = analysis.save_dir
        self.TelSaveDir = analysis.TelSaveDir

        self.load_channel_config()

        self.PedestalFit = None
        self.ped_range = None

        self.generate_channel_cutstrings()
        self.all_cut = self.generate_all_cut()
        self.CutStrings['AllCuts'] = self.all_cut

        self.ConsecutiveCuts = self.generate_consecutive_cuts()

    # ==============================================
    # region GET CONFIG
    def load_channel_config(self):
        self.CutConfig['absMedian_high'] = self.load_config_data('absMedian_high')
        self.CutConfig['pedestalsigma'] = self.load_config_data('pedestalsigma')
        self.CutConfig['fiducial'] = self.load_fiducial()
        self.CutConfig['threshold'] = self.load_dia_config('threshold', store_true=True)

    def load_fiducial(self):
        if self.ana_config_parser.has_option('SPLIT', 'fiducial'):
            split_runs = loads(self.ana_config_parser.get('SPLIT', 'fiducial')) + [int(1e10)]
            return next(self.load_dia_config('fid_cuts{n}'.format(n=i if i else '')) for i in xrange(len(split_runs)) if self.RunNumber <= split_runs[i])
        return self.load_dia_config('fid_cuts')

    def load_config_data(self, name):
        value = self.ana_config_parser.getint('CUT', name)
        return value if value > 0 else None

    def load_dia_config(self, name, store_true=False):
        try:
            conf = loads(self.ana_config_parser.get('CUT', name))
            dia = self.analysis.DiamondName
            dia = '{}*{}'.format(dia, self.analysis.DiamondNumber) if '{}*1'.format(dia) in conf else dia
            if store_true:
                return dia in conf
            return conf[dia]
        except (KeyError, NoOptionError):
            log_warning('No option {0} in the analysis config for {1}!'.format(name, make_tc_str(self.TESTCAMPAIGN)))

    # endregion

    # ==============================================
    # region SET CUTS
    def set_cut(self, name, value=None):
        if name not in self.CutStrings:
            log_warning('There is no cut with the name "{name}"!'.format(name=name))
            return
        self.reset_cut(name)
        self.CutStrings[name] += self.generate_cut(name, value)
        self.update_all_cut()

    def set_abs_median_high(self, high=None):
        self.set_cut('median', high)

    def set_pedestal_sigma(self, sigma=None):
        self.set_cut('pedestalsigma', sigma)

    def set_signal_peak_pos(self, x_min, x_max):
        self.set_cut('signal_peak_pos', [x_min, x_max])

    def set_signal_peak_time(self, x_min, x_max):
        self.set_cut('signal_peak_time', [x_min, x_max])

    def set_trigger_cell(self, x_min, x_max):
        self.set_cut('trigger_cell', [x_min, x_max])

    def set_bucket(self, value):
        self.set_cut('bucket', value)
    # endregion

    # ==============================================
    # region GENERATE CUT STRINGS
    def generate_cut(self, name, value):
        dic = {'median': self.generate_median,
               'pedestalsigma': self.generate_pedestalsigma,
               'signal_peak_pos': self.generate_signal_peak_pos,
               'signal_peak_time': self.generate_signal_peak_time,
               'trigger_cell': self.generate_trigger_cell,
               'bucket': self.generate_bucket,
               'chi2X': partial(self.generate_chi2, 'x'),
               'chi2Y': partial(self.generate_chi2, 'y')}
        return dic[name](value)

    def generate_median(self, high=None):
        value = self.CutConfig['absMedian_high'] if high is None else high
        string = ''
        if value is not None:
            assert value > 0, 'The median has to be a positive number!'
            string = 'abs(median[{ch}])<{high}'.format(ch=self.channel, high=float(high))
            self.EasyCutStrings['absMedian_high'] = '|median|<{high}'.format(high=value)
        return TCut(string)

    def generate_pedestalsigma(self, sigma=None):
        sigma = self.CutConfig['pedestalsigma'] if sigma is None else sigma
        string = ''
        if sigma is not None:
            assert sigma > 0, 'The sigma has to be a positive number!'
            ped_range = self.__calc_pedestal_range(sigma)
            self.ped_range = ped_range
            string = '{ped}>{min}&&{ped}<{max}'.format(ped=self.analysis.PedestalName, min=ped_range[0], max=ped_range[1])
            self.EasyCutStrings['pedestalsigma'] = 'PedSigma<{0}'.format(sigma)
        return TCut(string)

    def generate_region(self, signal_histo, mean_histo):
        extrema = Extrema2D(signal_histo, mean_histo)
        extrema.region_scan()
        extrema.show_voting_histos()
        all_string = ''
        nr = self.DiamondNumber - 1
        for col in xrange(extrema.cols):
            all_val = [bool(extrema.VotingHistos['max'].GetBinContent(col, row)) for row in xrange(extrema.rows)]
            # print col, all_val
            if True not in all_val:
                continue
            all_string += '||' if all_string else ''
            xmin = extrema.VotingHistos['max'].GetXaxis().GetBinLowEdge(col)
            xmax = extrema.VotingHistos['max'].GetXaxis().GetBinUpEdge(col)
            all_string += '(dia_track_x[{nr}]>{xmin}&&dia_track_x[{nr}]<{xmax})&&'.format(nr=nr, xmin=xmin, xmax=xmax)
            y_string = ''
            cont = True
            for row in xrange(extrema.rows + 1):
                val = extrema.VotingHistos['max'].GetBinContent(col, row) if not row == extrema.rows else 0
                last_val = extrema.VotingHistos['max'].GetBinContent(col, row - 1) if row else 0
                if val and not last_val:
                    y = extrema.VotingHistos['max'].GetYaxis().GetBinLowEdge(row)
                    if y < abs(1e-10):
                        cont = False
                        continue
                    cont = True
                    y_string += '||' if y_string else '('
                    y_string += 'dia_track_y[{nr}]>{y}&&'.format(nr=nr, y=y)
                elif not val and last_val and cont:
                    y_string += 'dia_track_y[{nr}]<{y}'.format(nr=nr, y=extrema.VotingHistos['max'].GetYaxis().GetBinUpEdge(row))
            y_string += ')'
            all_string += y_string
        self.region_cut += all_string
        return extrema

    def generate_signal_peak_pos(self, min_max):
        assert 0 <= min_max[0] <= 1024, 'min signal peak has to be in [0, 1024], not "{min}"'.format(min=min_max[0])
        assert 0 <= min_max[1] <= 1024, 'max signal peak has to be in [0, 1024], not "{max}"'.format(max=min_max[1])
        self.EasyCutStrings['SignalPeakPos'] = 'Signal Peak in {0}'.format(min_max)
        return TCut('IntegralPeaks[{num}] < {max} && IntegralPeaks[{num}] >= {min}'.format(num=self.analysis.SignalNumber, min=min_max[0], max=min_max[1]))

    def generate_signal_peak_time(self, min_max):
        assert 0 <= min_max[0] <= 1024, 'min signal peak time has to be in [0, 1024], not "{min}"'.format(min=min_max[0])
        assert 0 <= min_max[1] <= 1024, 'max signal peak time has to be in [0, 1024], not "{max}"'.format(max=min_max[1])
        self.EasyCutStrings['SignalPeakTime'] = 'Signal Peak Time in {0}'.format(min_max)
        return TCut('IntegralPeakTime[{num}] < {max} && IntegralPeakTime[{num}] >= {min}'.format(num=self.analysis.SignalNumber, min=min_max[0], max=min_max[1]))

    def generate_trigger_cell(self, min_max):
        assert 0 <= min_max[0] <= 1024, 'min trigger cell has to be in [0, 1024], not "{min}"'.format(min=min_max[0])
        assert 0 <= min_max[1] <= 1024, 'max trigger cell has to be in [0, 1024], not "{max}"'.format(max=min_max[1])
        self.EasyCutStrings['TriggerCell'] = 'Trigger Cell in {0}'.format(min_max)
        return TCut('trigger_cell < {max} && trigger_cell >= {min}'.format(min=min_max[0], max=min_max[1]))

    def generate_old_bucket(self):
        # only generate the cut if the region e2 exists! todo: find a smarter solution for that!
        try:
            sig2 = self.analysis.get_signal_name('e', 2)
            string = '{sig2}=={sig1}'.format(sig2=sig2, sig1=self.analysis.SignalName)
            return TCut(string)
        except ValueError as err:
            print err
            return TCut('')

    def generate_bucket(self, threshold=None):
        # TODO: bucket cut for high irradiation (low signals)
        sig = self.analysis.SignalName
        threshold = self.calc_signal_threshold(show=False) if threshold is None else threshold
        string = '!(!({old_buck}) && ({sig} < {thres}))'.format(sig=sig, thres=threshold, old_buck=self.CutStrings['old_bucket'])
        # string = self.CutStrings['old_bucket'] if threshold == -30 else string
        cut = TCut(string) if self.CutStrings['old_bucket'].GetTitle() else TCut('')
        return cut

    def generate_timing(self, n_sigma=3):
        t_correction, fit = self.calc_timing_range()
        if fit is None:
            return TCut('')
        # corrected_time = '{peak} - {t_corr}'.format(peak=self.analysis.Timing.get_peak_name(corr=True, region='e'), t_corr=t_correction)  # correction for bucket cut
        corrected_time = '{peak} - {t_corr}'.format(peak=self.analysis.Timing.get_peak_name(corr=True), t_corr=t_correction)
        string = 'TMath::Abs({cor_t} - {mp}) / {sigma} < {n_sigma}'.format(cor_t=corrected_time, mp=fit.GetParameter(1), sigma=fit.GetParameter(2), n_sigma=n_sigma)
        return TCut(string)

    def generate_threshold(self):
        return TCut('{sig}>{thresh}'.format(sig=self.analysis.SignalName, thresh=self.calc_threshold(show=False))) if self.CutConfig['threshold'] else TCut('')

    def generate_fiducial(self):
        xy = self.CutConfig['fiducial']
        cut = None
        if xy is not None:
            cut = TCutG('fid{}'.format(self.RunNumber), 5, array([xy[0], xy[0], xy[1], xy[1], xy[0]], 'd'), array([xy[2], xy[3], xy[3], xy[2], xy[2]], 'd'))
            nr = self.analysis.DiamondNumber - 1
            cut.SetVarX(self.get_track_var(nr, 'x'))
            cut.SetVarY(self.get_track_var(nr, 'y'))
            self.ROOTObjects.append(cut)
            cut.SetLineWidth(3)
        return TCut(cut.GetName() if cut is not None else '')

    def find_fid_cut(self, thresh=.93, show=True):
        h = self.analysis.draw_signal_map(show=False)
        px = h.ProjectionX()
        self.format_histo(px, title='Projection X of the Signal Map', y_tit='Number of Entries', y_off=1.5)
        self.draw_histo(px, lm=.12, show=show)
        py = h.ProjectionY()
        return '"{}": [{}]'.format(self.analysis.DiamondName, ', '.join('{:0.3f}'.format(i) for i in self.find_fid_margins(px, thresh) + self.find_fid_margins(py, thresh)))

    @staticmethod
    def find_fid_margins(proj, thresh):
        thresh = proj.GetMaximum() * thresh
        xbin1, xbin2 = proj.FindFirstBinAbove(thresh), proj.FindLastBinAbove(thresh)
        f1 = interpolate_two_points(proj.GetBinCenter(xbin1), proj.GetBinContent(xbin1), proj.GetBinCenter(xbin1 - 1), proj.GetBinContent(xbin1 - 1))
        f2 = interpolate_two_points(proj.GetBinCenter(xbin2), proj.GetBinContent(xbin2), proj.GetBinCenter(xbin2 + 1), proj.GetBinContent(xbin2 + 1))
        return [f1.GetX(thresh) / 10, f2.GetX(thresh) / 10]

    def get_fid_area(self):
        conf = self.CutConfig['fiducial']
        return (conf[1] - conf[0]) * (conf[3] - conf[2])

    def draw_fid_cut(self, scale=1):
        cut = get_object('fid{}'.format(self.RunNumber))
        if cut:
            cut = deepcopy(cut)
            cut.SetName('fid{}'.format(scale))
            for i in xrange(cut.GetN()):
                cut.SetPoint(i, scale * cut.GetX()[i], scale * cut.GetY()[i])
            cut.Draw()
            self.ROOTObjects.append(cut)

    # special cut for analysis
    def generate_pulser_cut(self, beam_on=True):
        cut = self.CutStrings['ped_sigma'] + self.CutStrings['event_range']
        cut.SetName('Pulser{0}'.format('BeamOn' if beam_on else 'BeamOff'))
        cut += self.CutStrings['beam_interruptions'] if beam_on else '!({0})'.format(self.JumpCut)
        cut += '!({0})'.format(self.CutStrings['pulser'])
        return cut

    def get_bucket_cut(self):
        cut = self.CutStrings['fiducial'] + self.CutStrings['pulser'] + TCut('!({})'.format(self.CutStrings['old_bucket']))
        cut.SetName('Bucket')
        return cut

    def generate_channel_cutstrings(self):

        # --THRESHOLD --
        self.CutStrings['threshold'] += self.generate_threshold()

        # -- PULSER CUT --
        self.CutStrings['pulser'] += '!pulser'

        # -- SATURATED CUT --
        self.CutStrings['saturated'] += '!is_saturated[{ch}]'.format(ch=self.channel)

        # -- MEDIAN CUT --
        self.CutStrings['median'] += self.generate_median()

        # -- PEDESTAL SIGMA CUT --
        self.CutStrings['ped_sigma'] += self.generate_pedestalsigma()

        # -- FIDUCIAL --
        self.CutStrings['fiducial'] += self.generate_fiducial()

        # --PEAK POSITION TIMING--
        self.CutStrings['timing'] += self.generate_timing()

        # --BUCKET --
        self.CutStrings['old_bucket'] += self.generate_old_bucket()
        self.CutStrings['bucket'] += self.generate_bucket()

    # endregion

    # ==============================================
    # HELPER FUNCTIONS

    def calc_signal_threshold(self, use_bg=False, show=True, show_all=False):
        pickle_path = self.make_pickle_path('Cuts', 'SignalThreshold', self.analysis.highest_rate_run, self.DiamondNumber)
        show = False if show_all else show

        def func():
            t = self.log_info('Calculating signal threshold for bucket cut of run {run} and {d} ...'.format(run=self.analysis.RunNumber, d=self.DiamondName), next_line=False)
            h = TH1F('h', 'Bucket Cut', 200, -50, 150)
            draw_string = '{name}>>h'.format(name=self.analysis.SignalName)
            fid = self.CutStrings['fiducial']
            cut_string = '!({buc})&&{pul}{fid}'.format(buc=self.CutStrings['old_bucket'], pul=self.CutStrings['pulser'], fid='&&{}'.format(fid.GetTitle()) if fid.GetTitle() else '')
            self.analysis.tree.Draw(draw_string, cut_string, 'goff')
            entries = h.GetEntries()
            if entries < 2000:
                self.add_info(t)
                return -30
            # extract fit functions
            set_root_output(False)
            fit = self.fit_bucket(h)
            if fit is None:
                self.add_info(t)
                return -30
            if fit is None or any(abs(fit.GetParameter(i)) < 20 for i in [0, 3]) or fit.GetParameter(1) < fit.GetParameter(4) or fit.GetParameter(1) > 500:
                warning('bucket cut fit failed')
                self.draw_histo(h, show=show)
                self.add_info(t)
                return -30
            sig_fit = TF1('f1', 'gaus', -50, 300)
            sig_fit.SetParameters(fit.GetParameters())
            ped_fit = TF1('f2', 'gaus(0) + gaus(3)', -50, 300)
            ped_fit.SetParameters(*[fit.GetParameter(i) for i in xrange(3, 9)])
            set_root_output(True)

            # real data distribution without pedestal fit
            signal = deepcopy(h)
            signal.Add(ped_fit, -1)

            gr1 = self.make_tgrapherrors('gr1', '#varepsilon_{bg}', marker_size=0.2)
            gr2 = self.make_tgrapherrors('gr2', '#varepsilon_{sig}', marker_size=0.2, color=2)
            gr3 = self.make_tgrapherrors('gr3', 'ROC Curve', marker_size=0.2)
            gr4 = self.make_tgrapherrors('gr4', 'Signal Error', marker_size=0.2)
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
                self.add_info(t)
                return -30
            max_err = max(errors.items())[1]
            c = None
            if show_all:
                set_root_output(True)
                c = TCanvas('c_all', 'Signal Threshold Overview', self.Res, self.Res)
                c.Divide(2, 2)
            # Bucket cut plot
            self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.8, stats=0)
            self.draw_histo(h, '', show or show_all, lm=.135, canvas=c.cd(1) if show_all else None)
            self.draw_y_axis(max_err, h.GetYaxis().GetXmin(), h.GetMaximum(), 'threshold  ', off=.3, line=True)
            ped_fit.SetLineStyle(2)
            ped_fit.Draw('same')
            sig_fit.SetLineColor(4)
            sig_fit.SetLineStyle(3)
            sig_fit.Draw('same')
            self.save_plots('BucketCut', sub_dir=self.analysis.save_dir, canvas=c.cd(1) if show_all else get_last_canvas(), prnt=show)

            # Efficiency plot
            self.format_histo(gr1, title='Efficiencies', x_tit='Threshold', y_tit='Efficiency', markersize=.2)
            l2 = self.make_legend(.78, .3)
            tits = ['#varepsilon_{bg}', gr2.GetTitle()]
            [l2.AddEntry(p, tits[i], 'l') for i, p in enumerate([gr1, gr2])]
            self.draw_histo(gr1, '', show_all, draw_opt='apl', l=l2, canvas=c.cd(2) if show_all else None)
            self.draw_histo(gr2, show=show_all, draw_opt='same', canvas=c.cd(2) if show_all else get_last_canvas())
            self.save_plots('Efficiencies', canvas=c.cd(2) if show_all else get_last_canvas(), prnt=show)

            # ROC Curve
            self.format_histo(gr3, y_tit='background fraction', x_tit='excluded signal fraction', markersize=0.2, y_off=1.2)
            self.draw_histo(gr3, '', show_all, gridx=True, gridy=True, draw_opt='apl', canvas=c.cd(3) if show_all else None)
            p = self.make_tgrapherrors('gr', 'working point', color=2)
            p.SetPoint(0, 1 - sig_fit.Integral(-50, max_err) / signal.Integral(), ped_fit.Integral(-50, max_err) / ped_fit.Integral(-50, 200))
            sleep(.1)
            l = self.draw_tlatex(p.GetX()[0], p.GetY()[0] + .01, 'Working Point', color=2, size=.04)
            p.GetListOfFunctions().Add(l)
            self.draw_histo(p, show=show_all, canvas=c.cd(3) if show_all else get_last_canvas(), draw_opt='p')
            self.save_plots('ROC_Curve', sub_dir=self.analysis.save_dir, canvas=c.cd(3) if show_all else get_last_canvas(), prnt=show)

            # Error Function plot
            self.format_histo(gr4, x_tit='Threshold', y_tit='1 / error', y_off=1.4)
            self.save_histo(gr4, 'ErrorFunction', show_all, gridx=True, draw_opt='al', canvas=c.cd(4) if show_all else None, prnt=show)

            self.Objects.append([sig_fit, ped_fit, gr2, c])

            self.add_info(t)
            return max_err

        threshold = func() if show or show_all else None
        threshold = do_pickle(pickle_path, func, threshold)
        return threshold

    @staticmethod
    def fit_bucket(histo, show=True):
        set_root_warnings(0)
        h = histo
        fit = TF1('fit', 'gaus(0) + gaus(3) + gaus(6)', h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        s = TSpectrum(3)
        n = s.Search(h, 2.5)
        points = [(s.GetPositionX()[i], s.GetPositionY()[i]) for i in [0, 1 if n == 2 else 2]]
        x1, x2 = (p[0] for p in sorted(points))
        y1, y2 = (p[1] for p in sorted(points))
        if y1 < 20 or y1 > 1e10:
            return  # didn't find pedestal peak!
        diff = x2 - x1
        fit.SetParameters(*[y2, x2, 10, y1, x1, 3, min(y1, y2) / 4, x1 + diff / 4, 5])
        # signal
        fit.SetParLimits(1, x2 - 5, x2 + 5)
        # pedestal
        fit.SetParLimits(3, 1, y1 * 2)
        fit.SetParLimits(4, x1 - 10, x1 + 10)
        # middle ped
        fit.SetParLimits(6, 1, min(y1, y2) / 2)
        fit.SetParLimits(7, x1, x1 + diff / 2)
        for i in xrange(1):
            h.Fit(fit, 'qs{0}'.format('' if show else '0'), '', -50, x2 + 5)
        set_root_warnings(1)
        return fit

    def __calc_pedestal_range(self, sigma_range):
        picklepath = self.make_pickle_path('Pedestal', 'Cut', self.RunNumber, self.channel)

        def func():
            t = self.log_info('generating pedestal cut for {dia} of run {run} ...'.format(run=self.analysis.RunNumber, dia=self.analysis.DiamondName), next_line=False)
            h1 = TH1F('h_pdc', 'Pedestal Distribution', 600, -150, 150)
            self.analysis.tree.Draw('{name}>>h_pdc'.format(name=self.analysis.PedestalName), '', 'goff')
            fit_pars = self.fit_fwhm(h1, do_fwhm=True, draw=False)
            self.add_info(t)
            return fit_pars

        fit = do_pickle(picklepath, func)
        sigma = fit.Parameter(2)
        mean_ = fit.Parameter(1)
        self.PedestalFit = fit
        return [mean_ - sigma_range * sigma, mean_ + sigma_range * sigma]

    def calc_threshold(self, show=True):
        pickle_path = self.make_pickle_path('Cuts', 'Threshold', self.analysis.RunNumber, self.channel)

        def func():
            self.analysis.tree.Draw(self.analysis.SignalName, '', 'goff', 5000)
            xvals = sorted([self.analysis.tree.GetV1()[i] for i in xrange(5000)])
            x_range = [xvals[0] - 5, xvals[-5]]
            h = self.analysis.draw_signal_distribution(show=show, cut=self.generate_fiducial(), x_range=x_range, bin_width=1)
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

    def calc_timing_range(self, redo=False):
        pickle_path = self.make_pickle_path('Cuts', 'TimingRange', self.RunNumber, self.DiamondNumber)

        def func():
            t = self.log_info('Generating timing cut for {dia} of run {run} ...'.format(run=self.analysis.RunNumber, dia=self.analysis.DiamondName), next_line=False)
            cut = self.generate_special_cut(excluded=['timing'], prnt=False, name='timing_cut')
            t_correction = self.analysis.Timing.calc_fine_correction(redo=redo)
            h = self.analysis.Timing.draw_peaks(show=False, cut=cut, fine_corr=t_correction != '0', prnt=False, redo=redo)
            fit = h.GetListOfFunctions()[2]
            if fit.GetParameter(2) > 15:  # fit failed
                fit.SetParameter(1, h.GetBinCenter(h.GetMinimumBin()))
                fit.SetParameter(2, 15)
            self.add_info(t)
            self.log_info('Peak Timing: Mean: {0}, sigma: {1}'.format(fit.GetParameter(1), fit.GetParameter(2)))
            return t_correction, fit

        return do_pickle(pickle_path, func, redo=redo)

    def generate_consecutive_cuts(self):
        cuts = OrderedDict([('raw', TCut('0', ''))])
        for i, (key, value) in enumerate([(key, value) for key, value in self.CutStrings.iteritems() if str(value) and key != 'AllCuts' and not key.startswith('old')], 1):
            new_cut = cuts.values()[i - 1] + value
            key = 'beam_stops' if 'beam' in key else key
            cuts[key] = TCut('{n}'.format(n=i), str(new_cut))
        return cuts
