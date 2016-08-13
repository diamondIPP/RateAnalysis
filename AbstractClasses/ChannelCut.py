from Cut import Cut
from Extrema import Extrema2D
from copy import deepcopy
from ROOT import TCut, TH1F, TH2F, TF1, TCanvas, TLegend, gROOT, TProfile, THStack
from collections import OrderedDict
from Utils import *
from ConfigParser import NoOptionError

__author__ = 'micha'


# ==============================================
# MAIN CLASS
# ==============================================
class ChannelCut(Cut):
    """The ChannelCut contains all cut settings which corresponds to a single diamond in a single run. """

    def __init__(self, analysis, channel=0):
        Cut.__init__(self, analysis, skip=True)
        self.__dict__.update(analysis.Cut.__dict__)
        self.channel = channel

        self.load_channel_config()

        self.PedestalFit = None
        self.ped_range = None

        self.generate_channel_cutstrings()
        self.all_cut = self.generate_all_cut()
        self.CutStrings['all_cuts'] = self.all_cut
        self.run_number = self.analysis.run_number

        self.ConsecutiveCuts = self.load_consecutive_cuts()

    # ==============================================
    # region GET CONFIG
    def load_channel_config(self):
        self.CutConfig['absMedian_high'] = self.load_config_data('absMedian_high')
        self.CutConfig['pedestalsigma'] = self.load_config_data('pedestalsigma')

    def load_config_data(self, name):
        value = self.ana_config_parser.getint('CUT', name)
        return value if value > 0 else None

    def load_dia_config(self, name, store_true=False):
        try:
            conf = loads(self.ana_config_parser.get('CUT', name))
            if not store_true:
                return conf[self.analysis.diamond_name] if self.analysis.diamond_name in conf else None
            else:
                return True if self.analysis.diamond_name in conf else False
        except NoOptionError:
            log_warning('No option {0} in the analysis config for {1}!'.format(name, make_tc_str(self.TESTCAMPAIGN)))

    # endregion

    # ==============================================
    # region SET CUTS
    def set_cut(self, name, value):
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

    def load_consecutive_cuts(self):
        dic = OrderedDict()
        for key, value in self.CutStrings.iteritems():
            if (str(value) or key == 'raw') and key != 'all_cuts' and not key.startswith('old'):
                dic[key] = value
        return dic
    # endregion

    # ==============================================
    # region GENERATE CUT STRINGS
    def generate_cut(self, name, value):
        if name == 'median':
            return self.generate_median(value)
        if name == 'pedestalsigma':
            return self.generate_pedestalsigma(value)
        if name == 'signal_peak_pos':
            return self.generate_signal_peak_pos(value)
        if name == 'signal_peak_time':
            return self.generate_signal_peak_time(value)
        if name == 'trigger_cell':
            return self.generate_trigger_cell(value)
        if name == 'bucket':
            return self.generate_bucket(value)

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
        nr = 2 if self.channel else 1
        for col in xrange(extrema.cols):
            all_val = [bool(extrema.VotingHistos['max'].GetBinContent(col, row)) for row in xrange(extrema.rows)]
            # print col, all_val
            if True not in all_val:
                continue
            all_string += '||' if all_string else ''
            xmin = extrema.VotingHistos['max'].GetXaxis().GetBinLowEdge(col)
            xmax = extrema.VotingHistos['max'].GetXaxis().GetBinUpEdge(col)
            all_string += '(diam{nr}_track_x>{xmin}&&diam{nr}_track_x<{xmax})&&'.format(nr=nr, xmin=xmin, xmax=xmax)
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
                    y_string += 'diam{nr}_track_y>{y}&&'.format(nr=nr, y=y)
                elif not val and last_val and cont:
                    y_string += 'diam{nr}_track_y<{y}'.format(nr=nr, y=extrema.VotingHistos['max'].GetYaxis().GetBinUpEdge(row))
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
        except AssertionError as err:
            print err
            return TCut('')

    def generate_bucket(self, threshold=None):
        sig = self.analysis.SignalName
        threshold = self.calc_signal_threshold(show=False) if threshold is None else threshold
        string = '!(!({old_buck}) && ({sig} < {thres}))'.format(sig=sig, thres=threshold, old_buck=self.CutStrings['old_bucket'])
        cut = TCut(string) if self.CutStrings['old_bucket'].GetTitle() else TCut('')
        return cut

    def generate_timing(self, n_sigma=3):
        dic = self.calc_timing_range(show=False)
        num = self.analysis.SignalNumber
        t_correction = '({p1}* trigger_cell + {p2} * trigger_cell*trigger_cell)'.format(p1=dic['t_corr'].GetParameter(1), p2=dic['t_corr'].GetParameter(2))
        corrected_time = 'IntegralPeakTime[{num}] - {t_corr}'.format(num=num, t_corr=t_correction)
        try:
            string = 'TMath::Abs({cor_t} - {mp}) / {sigma} < {n_sigma}'.format(cor_t=corrected_time, mp=dic['timing_corr'].GetParameter(1), sigma=dic['timing_corr'].GetParameter(2), n_sigma=n_sigma)
        except: 
            print dic['timing_corr']
            raise Exception()
        return TCut(string), corrected_time, t_correction

    # special cut for analysis
    def generate_pulser_cut(self, beam_on=True):
        cut = self.CutStrings['ped_sigma'] + self.CutStrings['event_range'] + self.CutStrings['saturated']
        cut.SetName('Pulser{0}'.format('BeamOn' if beam_on else 'BeamOff'))
        cut += self.CutStrings['beam_interruptions'] if beam_on else '!({0})'.format(self.JumpCut)
        cut += '!({0})'.format(self.CutStrings['pulser'])
        return cut

    def generate_channel_cutstrings(self):

        # -- PULSER CUT --
        self.CutStrings['pulser'] += '!pulser'

        # -- SATURATED CUT --
        self.CutStrings['saturated'] += '!is_saturated[{ch}]'.format(ch=self.channel)

        # -- MEDIAN CUT --
        self.CutStrings['median'] += self.generate_median()

        # -- PEDESTAL SIGMA CUT --
        self.CutStrings['ped_sigma'] += self.generate_pedestalsigma()

        # --PEAK POSITION TIMING--
        cut, corrected_peak_time, time_correction = self.generate_timing()
        self.analysis.CorrectedTime = corrected_peak_time
        self.analysis.TimeCorrection = time_correction
        self.CutStrings['timing'] += cut

        # --BUCKET --
        self.CutStrings['old_bucket'] += self.generate_old_bucket()
        self.CutStrings['bucket'] += self.generate_bucket()

    # endregion

    # ==============================================
    # HELPER FUNCTIONS

    def calc_signal_threshold(self, bg=False, show=True, show_all=False):
        pickle_path = self.analysis.PickleDir + 'Cuts/SignalThreshold_{tc}_{run}_{ch}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.highest_rate_run, ch=self.channel)
        show = False if show_all else show

        def func():
            print 'calculating signal threshold for bucket cut of run {run} and ch{ch}...'.format(run=self.analysis.run_number, ch=self.channel)
            h = TH1F('h', 'Bucket Cut', 100, -50, 150)
            draw_string = '{name}>>h'.format(name=self.analysis.SignalName)
            cut_string = '!({buc})&&{pul}'.format(buc=self.CutStrings['old_bucket'], pul=self.CutStrings['pulser'])
            self.analysis.tree.Draw(draw_string, cut_string, 'goff')
            entries = h.GetEntries()
            if entries < 2000:
                return 30
            h.Rebin(2) if entries < 5000 else self.do_nothing()
            # extract fit functions
            self.set_root_output(False)
            fit = self.triple_gauss_fit(h)
            sig_fit = TF1('f1', 'gaus', -50, 300)
            sig_fit.SetParameters(fit.GetParameters())
            ped_fit = TF1('f2', 'gaus(0) + gaus(3)', -50, 300)
            pars = [fit.GetParameter(i) for i in xrange(3, 9)]
            ped_fit.SetParameters(*pars)
            self.set_root_output(True)

            # real data distribution without pedestal fit
            signal = deepcopy(h)
            signal.Add(ped_fit, -1)

            gr1 = self.make_tgrapherrors('gr1', '#varepsilon_{bg}', marker_size=0.2)
            gr2 = self.make_tgrapherrors('gr2', '#varepsilon_{sig}', marker_size=0.2, color=2)
            gr3 = self.make_tgrapherrors('gr3', 'ROC Curve', marker_size=0.2)
            gr4 = self.make_tgrapherrors('gr4', 'Signal Error', marker_size=0.2)
            xs = [i / 10. for i in xrange(-300, int(sig_fit.GetMaximumX()) * 10)]
            errors = {}
            for i, x in enumerate(xs):
                ped = ped_fit.Integral(-50, x) / ped_fit.Integral(-50, 300)
                sig = 1 - sig_fit.Integral(-50, x) / signal.Integral()
                err = ped_fit.Integral(-50, x) / (sqrt(sig_fit.Integral(-50, x) + ped_fit.Integral(-50, x)))
                err1 = signal.Integral(h.FindBin(x), h.FindBin(300)) / sqrt(signal.Integral(h.FindBin(x), h.FindBin(300)) + ped_fit.Integral(x, 300))
                errors[err1 if not bg else err] = x
                gr1.SetPoint(i, x, ped)
                gr2.SetPoint(i, x, sig)
                gr3.SetPoint(i, sig, ped)
                gr4.SetPoint(i, x, err1 if not bg else err)
            if len(errors) == 0:
                print ValueError('errors has a length of 0')
                return -30
            max_err = errors[max(errors.keys())]

            c = None
            if show_all:
                self.set_root_output(True)
                c = TCanvas('c_all', 'Signal Threshold Overview', self.Res, self.Res)
                c.Divide(2, 2)

            # Bucket cut plot
            self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.8, stats=0)
            self.draw_histo(h, '', show, lm=.135, canvas=c.cd(1) if show_all else None)
            self.draw_y_axis(max_err, h.GetYaxis().GetXmin(), h.GetMaximum(), 'threshold  ', off=.3, line=True)
            ped_fit.SetLineStyle(2)
            ped_fit.Draw('same')
            sig_fit.SetLineColor(4)
            sig_fit.SetLineStyle(3)
            sig_fit.Draw('same')
            self.save_plots('BucketCut', sub_dir=self.analysis.save_dir)

            # Efficiency plot
            self.format_histo(gr1, title='Efficiencies', x_tit='Threshold', y_tit='Efficiency', markersize=.2)
            l2 = self.make_legend(.78, .3)
            tits = ['#varepsilon_{bg}', gr2.GetTitle()]
            [l2.AddEntry(p, tits[i], 'l') for i, p in enumerate([gr1, gr2])]
            self.draw_histo(gr1, '', False, draw_opt='apl', l=l2, canvas=c.cd(2) if show_all else None)
            gr2.Draw('pl')
            self.save_plots('Efficiencies')

            # ROC Curve
            self.format_histo(gr3, y_tit='background fraction', x_tit='excluded signal fraction', markersize=0.2, y_off=1.2)
            self.draw_histo(gr3, '', False, gridx=True, gridy=True, draw_opt='apl', canvas=c.cd(3) if show_all else None)
            p = self.make_tgrapherrors('gr', 'working point', color=2)
            p.SetPoint(0, 1 - sig_fit.Integral(-50, max_err) / signal.Integral(), ped_fit.Integral(-50, max_err) / ped_fit.Integral(-50, 300))
            l = self.draw_tlatex(p.GetX()[0], p.GetY()[0] + .01, 'Working Point', color=2, size=.04)
            p.GetListOfFunctions().Add(l)
            p.Draw('p')
            self.save_plots('ROC_Curve', sub_dir=self.analysis.save_dir)

            # Error Function plot
            self.format_histo(gr4, x_tit='Threshold', y_tit='1 / error', y_off=1.4)
            self.save_histo(gr4, 'ErrorFunction', False, gridx=True, draw_opt='al', canvas=c.cd(4) if show_all else None)

            self.RootObjects.append([sig_fit, ped_fit, gr2, c])

            return max_err

        threshold = func() if show or show_all else None
        threshold = self.do_pickle(pickle_path, func, threshold)
        return threshold

    def find_ped_range(self):
        self.analysis.tree.Draw(self.analysis.PedestalName, '', 'goff', 1000)
        return calc_mean([self.analysis.tree.GetV1()[i] for i in xrange(1000)])

    def __calc_pedestal_range(self, sigma_range):
        ped_range = self.find_ped_range()
        x_range = [ped_range[0] - 5 * ped_range[1], ped_range[0] + 10 * ped_range[1]]
        fit = self.analysis.show_pedestal_histo(region=self.analysis.PedestalRegion, peak_int=self.analysis.PeakIntegral, save=False, cut='', show=False, x_range=x_range)
        sigma = fit.Parameter(2)
        mean = fit.Parameter(1)
        self.PedestalFit = fit
        return [mean - sigma_range * sigma, mean + sigma_range * sigma]

    def calc_timing_range(self, show=True, n_sigma=4):
        pickle_path = self.analysis.PickleDir + 'Cuts/TimingRange_{tc}_{run}_{ch}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run_number, ch=self.channel)

        def func():
            print 'generating timing cut for {dia} of run {run}...'.format(run=self.analysis.run_number, dia=self.analysis.diamond_name)

            gROOT.SetBatch(1) if not show else self.do_nothing()
            num = self.analysis.SignalNumber
            cut = self.generate_special_cut(excluded_cuts=['bucket', 'timing'])

            # estimate timing
            draw_string = 'IntegralPeakTime[{num}]'.format(num=num)
            self.analysis.tree.Draw(draw_string, cut, 'goff')
            h1 = gROOT.FindObject('htemp')
            fit1 = TF1('fit1', 'gaus', -50, 1024)
            max_bin = h1.GetBinCenter(h1.GetMaximumBin())
            fit1.SetParLimits(1, max_bin - 5, max_bin + 5)
            h1.Fit(fit1, 'q0')
            h1.GetListOfFunctions().Add(fit1)
            original_mpv = fit1.GetParameter(1)
            print 'mean: {0}, sigma: {1}'.format(original_mpv, fit1.GetParameter(2))

            # extract timing correction
            h2 = TProfile('tcorr', 'Original Peak Position vs Trigger Cell', 1024, 0, 1024)
            self.analysis.tree.Draw('IntegralPeakTime[{num}]:trigger_cell>>tcorr'.format(num=num), cut, 'goff')
            fit2 = TF1('fit2', 'pol2', -50, 1024)
            h2.Fit(fit2, 'q0',)
            h2.GetListOfFunctions().Add(fit2)
            self.format_histo(h2, x_tit='trigger cell', y_tit='signal peak time', y_off=1.5)
            self.RootObjects.append(self.save_histo(h2, 'OriPeakPosVsTriggerCell', False, self.analysis.save_dir, lm=.12))
            t_correction = '({p1}* trigger_cell + {p2} * trigger_cell*trigger_cell)'.format(p1=fit2.GetParameter(1), p2=fit2.GetParameter(2))

            # get time corrected sigma
            h3 = TH1F('h3', 'Corrected Timing', 80, int(original_mpv - 10), int(original_mpv + 10))
            self.analysis.tree.Draw('(IntegralPeakTime[{num}] - {t_corr}) >> h3'.format(num=num, t_corr=t_correction), cut, 'goff')
            fit3 = TF1('fit3', 'gaus', -50, 1024)
            h3.Fit(fit3, 'q0')
            h3.GetListOfFunctions().Add(fit3)
            self.format_histo(h3, x_tit='time [ns]', y_tit='entries', y_off=2.1)
            self.RootObjects.append(self.save_histo(h3, 'TimingCorrection', False, self.analysis.save_dir, lm=.15))
            gROOT.SetBatch(0)

            if show:
                corrected_time = 'IntegralPeakTime[{num}] - {t_corr}'.format(num=num, t_corr=t_correction)
                t_cut = TCut('TMath::Abs({cor_t} - {mp}) / {sigma} < {n_sigma}'.format(cor_t=corrected_time, mp=fit3.GetParameter(1), sigma=fit3.GetParameter(2), n_sigma=n_sigma))
                # print results
                c = TCanvas('c_timing', 'Timing Cut Results', 1000, 1000)
                c.Divide(2, 2)
                # fit for correction
                c.cd(1)
                h2.Draw()
                # corrected timing
                c.cd(2)
                h4 = TProfile('h4', 'Corrected Peak Position vs Trigger Cell', 512, 0, 1024)
                h5 = TProfile('h5', 'Corrected Peak Position vs Trigger Cell with Cut', 512, 0, 1024)
                self.analysis.tree.Draw('{cor}:trigger_cell>>h4'.format(num=num, cor=corrected_time), cut, 'goff')
                self.analysis.tree.Draw('{cor}:trigger_cell>>h5'.format(num=num, cor=corrected_time), cut + t_cut, 'goff')
                self.format_histo(h4, x_tit='trigger cell', y_tit='signal peak times [ns]', y_off=1.6, color=self.get_color(), markersize=.5)
                self.format_histo(h5, color=self.get_color(), markersize=.5)
                h4.SetLineColor(1)
                h5.SetLineColor(1)
                self.reset_colors()
                h4.SetStats(0)
                h4.Draw()
                h5.Draw('same')
                # compare distributions
                c.cd(3)
                h6 = TH1F('h6', 'Corrected Timing with Cut', 80, int(original_mpv - 10), int(original_mpv + 10))
                self.analysis.tree.Draw('(IntegralPeakTime[{num}] - {t_corr}) >> h6'.format(num=num, t_corr=t_correction), cut + t_cut, 'goff')
                stack = THStack('stack', 'Time Comparison;time [ns];entries')
                mu = 0
                l = TLegend(.6, .78, .88, .88)
                l_names = ['before', 'after', 'after']
                for i, h in enumerate([h1, h3, h6]):
                    h.SetStats(0)
                    h.SetLineColor(self.get_color())
                    if len(h.GetListOfFunctions()):
                        fit = deepcopy(h.GetListOfFunctions()[-1])
                        fit.SetLineColor(h.GetLineColor())
                        mu = fit.GetParameter(1)
                        sig = fit.GetParameter(2)
                        l.AddEntry(fit, 'sigma {nam} corr.: {sig:1.2} ns'.format(nam=l_names[i], sig=sig), 'l')
                        fit.SetParameter(1, 0)
                    xax = h.GetXaxis()
                    xax.SetLimits(xax.GetXmin() - mu, xax.GetXmax() - mu)
                    stack.Add(h)
                stack.Draw('nostack')
                stack.GetXaxis().SetRangeUser(-4, 4)
                l.Draw()

                self.RootObjects.append([c, h4, h5, h6, h1, stack, l])

            return {'t_corr': fit2, 'timing_corr': fit3}
        fits = func() if show else None
        fits = self.do_pickle(pickle_path, func, fits)
        return fits

    def generate_timing_cut(self, sigma=4, show=False):
        # picklepath = 'Configuration/Individual_Configs/TimingCut/{tc}_{run}_{mod}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run.run_number, mod=mode.title())
        # if not bPlot:
        #     return TCut('')
        print 'generate_timing_cut with %s sigma' % sigma
        cut = self.generate_special_cut(excluded_cuts=['timing','bucket'])
        # Estimate Timing
        print ' * Estimate Timing',
        # hTiming = TH1F('hTiming','hTiming',4096,0,512)
        draw_string = 'IntegralPeakTime[{num}]>>hTiming'.format(num=self.analysis.SignalNumber)
        self.analysis.tree.Draw(draw_string,cut,'goff')
        hTiming = gROOT.FindObject('hTiming')
        mp = hTiming.GetBinCenter(hTiming.GetMaximumBin())
        fGaus = TF1('fGaus','gaus',mp-10,mp+10)
        hTiming.Fit(fGaus,'Q','goff',mp-10,mp+10)
        orignal_mp = fGaus.GetParameter(1)
        print '  - mean',orignal_mp,'sigma',fGaus.GetParameter(2)

        # Extract Timing Correction
        print ' * Extract Timing Correction'
        hTimingVsTriggerCell_pfy = TProfile('hTimingVsTriggerCell_pfy','Orig. Peak Position vs Trigger Cell',1024,0,1024)
        draw_string = 'IntegralPeakTime[{num}]:trigger_cell>>hTimingVsTriggerCell_pfy'.format(num=self.analysis.SignalNumber)
        self.analysis.tree.Draw(draw_string,cut,'pfy goff')
        fPol2 = TF1('fPol2','pol2',0,1024)
        gROOT.GetListOfGlobalFunctions().Add(fPol2)
        hTimingVsTriggerCell_pfy.Fit(fPol2,'Q','goff')


        # Define Timing Correction
        timing_correction = '({p0} + {p1}* trigger_cell + {p2} * trigger_cell*trigger_cell)'.format(
            p0=fPol2.GetParameter(0),
            p1=fPol2.GetParameter(1),
            p2=fPol2.GetParameter(2)
        )

        print ' * Time Correction: ',
        print fPol2.GetParameter(0), fPol2.GetParameter(1), fPol2.GetParameter(2)
        print '                  ',timing_correction

        # Get Time Corrected Sigma
        hNew = TH1F('hNewTiming','Corrected Timing; corrected Peak Time /ns;number of entries',100,-10,10)
        hNew.SetLineColor(2)
        draw_string = 'IntegralPeakTime[{num}]-{cor}>>hNewTiming'.format(num=self.analysis.SignalNumber,cor= timing_correction)
        self.analysis.tree.Draw(draw_string,cut,'goff')
        hNew.Fit(fGaus,'Q','goff')
        hNew.SetTitle('Corrected Timing: #sigma=%.3f'%fGaus.GetParameter(2))
        print 'New Sigma:',fGaus.GetParameter(2),'@',fGaus.GetParameter(1)


        # Define Cut
        lowEdge = fGaus.GetParameter(1) - sigma * fGaus.GetParameter(2)
        upEdge  = fGaus.GetParameter(1) + sigma * fGaus.GetParameter(2)
        peakTime = 'IntegralPeakTime[{num}]'.format(num=self.analysis.SignalNumber)
        timingCut = TCut('{lowEdge} < {peakTime} - {correctedTime} &&{peakTime} - {correctedTime} < {upEdge}'.format(
            lowEdge=lowEdge,upEdge=upEdge,correctedTime=timing_correction,peakTime=peakTime))
        timingCut.SetName('timing')

        gROOT.SetBatch(0)
        if show:
            self.analysis.histos.append(hTimingVsTriggerCell_pfy)
            self.analysis.histos.append(hNew)
            # Get Original Timing Distribution
            hOrigTiming = TH1F('hOrigTiming','Original Timing; Peak Time - MP / ns; number of entries',100,-10,10)
            c2 = TCanvas()
            print 'orignal_mp',orignal_mp
            draw_string = 'IntegralPeakTime[{num}]-{mp}>>hOrigTiming'.format(num=self.analysis.SignalNumber,mp=orignal_mp)
            print draw_string
            self.analysis.tree.Draw(draw_string,cut,'')

            self.analysis.canvases['cTiming2'] = c2
            hOrigTiming.Draw()
            fGaus2= TF1('fGaus2','gaus',-5,5)
            fGaus2.SetLineColor(3)
            hOrigTiming.Fit(fGaus2,'Q','')
            print 'Original Sigma:',fGaus2.GetParameter(2),'@',fGaus2.GetParameter(1), 'with ',hOrigTiming.GetEntries()
            hOrigTiming.SetTitle('Original Timing: #sigma=%.3f'%fGaus2.GetParameter(2))
            self.analysis.histos.append(hOrigTiming)
            c2.Update()

            # Get New Timing with applied Cut
            hNewTiming2 = hNew.Clone('hNewTiming2')
            hNewTiming2.SetTitle('Corrected Timing w Timing Cut')
            hNewTiming2.SetLineColor(4)
            fGaus3 = fGaus.Clone('fGaus3')
            draw_string = 'IntegralPeakTime[{num}]-{cor}>>hNewTiming2'.format(num=self.analysis.SignalNumber,cor= timing_correction)
            self.analysis.tree.Draw(draw_string,cut+timingCut,'goff')
            fGaus3.SetLineColor(5)
            fGaus3.SetLineStyle(2)
            hNewTiming2.Fit(fGaus3,'Q')
            self.analysis.histos.append(hNewTiming2)


            hTimingVsTriggerCellCorrected_pfy = TProfile('hTimingVsTriggerCellCorrected_pfy','Corrected Peak Position vs Trigger Cell',1024,0,1024)
            draw_string = 'IntegralPeakTime[{num}]-{pol}:trigger_cell>>hTimingVsTriggerCellCorrected_pfy'.format(num=self.analysis.SignalNumber,pol=timing_correction)
            self.analysis.tree.Draw(draw_string,cut,'colz goff')
            self.analysis.histos.append(hTimingVsTriggerCellCorrected_pfy)

            hTimingVsTriggerCellCorrected2_pfy = TProfile('hTimingVsTriggerCellCorrected2_pfy','Corrected Peak Position vs Trigger Cell',1024,0,1024)
            draw_string = 'IntegralPeakTime[{num}]-{pol}:trigger_cell>>hTimingVsTriggerCellCorrected2_pfy'.format(num=self.analysis.SignalNumber,pol=timing_correction)
            self.analysis.tree.Draw(draw_string,cut+timingCut,'colz goff')
            hTimingVsTriggerCellCorrected2_pfy.SetLineColor(3)
            self.analysis.histos.append(hTimingVsTriggerCellCorrected2_pfy)


            draw_string = 'IntegralPeakTime[{num}]-{pol}:trigger_cell>>hTimingVsTriggerCellCorrected'.format(num=self.analysis.SignalNumber,pol=timing_correction)
            hTimingVsTriggerCellCorrected = TH2F('hTimingVsTriggerCellCorrected','Peak Position vs Trigger Cell',1024,0,1024,100,-10,10)
            self.analysis.tree.Draw(draw_string,cut,'colz goff')
            self.analysis.histos.append(hTimingVsTriggerCellCorrected)
            c1 = TCanvas('cTiming','cTiming',1000,1000)
            c1.Divide(2,2)

            # Draw Orig Time vs cell with fit
            c1.cd(1)
            hTimingVsTriggerCell_pfy.Draw('')
            c1.Update()

            # Draw Stack of histos: hNewTiming, hOriginal, hNewTimingWithCut
            c1.cd(3)
            stack = THStack('stack','Time Comparison;time /ns;number of entries')
            stack.Add(hNew)
            stack.Add(hNewTiming2)
            stack.Add(hOrigTiming)
            stack.Draw('nostack')
            c1.cd(3).BuildLegend()
            c1.Update()
            self.analysis.histos.append(stack)

            # Draw Corrected Time vs cell with and without time cut
            c1.cd(2)
            hTimingVsTriggerCellCorrected_pfy.Draw()
            hTimingVsTriggerCellCorrected2_pfy.Draw('same')
            c1.Update()

            #Draw Corrected Time vs cell  Profile
            c1.cd(4)
            hTimingVsTriggerCellCorrected.Draw('colz')
            c1.Update()
            self.analysis.canvases['cTiming'] = c1
        return timingCut
