from Cut import Cut
from Extrema import Extrema2D
from copy import deepcopy
from ROOT import TCut, TH1F, TH2F, TF1, TCanvas, TLegend, gROOT, TProfile, THStack
from time import sleep
from numpy import sqrt

__author__ = 'micha'


# ==============================================
# MAIN CLASS
# ==============================================
class ChannelCut(Cut):
    """ The ChannelCut contains all cut settings which corresponds to a single diamond in a single run. """

    def __init__(self, analysis, channel=0):
        Cut.__init__(self, analysis, skip=True)
        self.__dict__.update(analysis.Cut.__dict__)
        self.channel = channel

        self.load_channel_config()

        self.PedestalFit = None

        self.generate_channel_cutstrings()
        self.all_cut = self.generate_all_cut()
        self.CutStrings['all_cuts'] = self.all_cut

    # ==============================================
    # region GET CONFIG
    def load_channel_config(self):
        self.CutConfig['absMedian_high'] = self.load_config_data('absMedian_high')
        self.CutConfig['pedestalsigma'] = self.load_config_data('pedestalsigma')

    def load_config_data(self, name):
        value = self.ana_config_parser.getint('CUT', name)
        return value if value > 0 else None

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

    def set_trigger_cell(self, x_min, x_max):
        self.set_cut('trigger_cell', [x_min, x_max])

    def set_bucket(self, value):
        self.set_cut('bucket', value)
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
        assert 0 <= min_max[1] <= 1024, 'max trigger cell has to be in [0, 1024], not "{max}"'.format(min=min_max[1])
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

    # special cut for analysis
    def generate_pulser_cut(self, beam_on=True):
        cut = self.CutStrings['ped_sigma'] + self.CutStrings['event_range'] + self.CutStrings['saturated']
        cut.SetName('Pulser{0}'.format('BeamOn' if beam_on else 'BeamOff'))
        cut += self.CutStrings['beam_interruptions'] if beam_on else '!({0})'.format(self.JumpCut)
        cut += '!({0})'.format(self.CutStrings['pulser'])
        return cut

    def generate_channel_cutstrings(self):
        # -- SATURATED CUT --
        self.CutStrings['saturated'] += '!is_saturated[{ch}]'.format(ch=self.channel)

        # -- MEDIAN CUT --
        self.CutStrings['median'] += self.generate_median()

        # -- PEDESTAL SIGMA CUT --
        self.CutStrings['ped_sigma'] += self.generate_pedestalsigma()

        # --PEAK POSITION TIMING--
        # todo: add a method that fits the real time disto and sets the cut to 4 sigma!
        # self.CutStrings['signal_peak_time'] += self.generate_signal_peak_time()
        self.CutStrings['timing'] = TCut('timing','')  # self.generate_timing_cut()

        # --BUCKET --
        self.CutStrings['old_bucket'] += self.generate_old_bucket()
        self.CutStrings['bucket'] += self.generate_bucket()

    # endregion

    # ==============================================
    # HELPER FUNCTIONS

    def calc_signal_threshold(self, bg=False, show=True):
        pickle_path = self.analysis.PickleDir + 'Cuts/SignalThreshold_{tc}_{run}_{ch}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.highest_rate_run, ch=self.channel)

        def func():
            print 'calculating signal threshold for bucket cut of run {run} and ch{ch}...'.format(run=self.analysis.run_number, ch=self.channel)
            h = TH1F('h', 'Bucket Cut', 250, -50, 300)
            draw_string = '{name}>>h'.format(name=self.analysis.SignalName)
            cut_string = '!({buc})&&{pul}'.format(buc=self.CutStrings['old_bucket'], pul=self.CutStrings['pulser'])
            self.analysis.tree.Draw(draw_string, cut_string, 'goff')
            entries = h.GetEntries()
            if entries < 1000:
                return 30
            h.Rebin(2) if entries < 5000 else self.do_nothing()
            # extract fit functions
            fit = self.triple_gauss_fit(h)
            sig_fit = TF1('f1', 'gaus', -50, 300)
            sig_fit.SetParameters(fit.GetParameters())
            ped_fit = TF1('f2', 'gaus(0) + gaus(3)', -50, 300)
            pars = [fit.GetParameter(i) for i in xrange(3, 9)]
            ped_fit.SetParameters(*pars)

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
            max_err = errors[max(errors.keys())]
            if show:
                c1 = TCanvas('c1', 'c', 1000, 1000)
                c1.SetLeftMargin(.127)
                self.format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.8)
                h.Draw()
                sleep(.1)
                a = self.make_tgaxis(max_err, c1.GetUymin(), c1.GetUymax(), 'threshold', offset=.3)
                a.Draw()
                # add subfunction to the plot
                ped_fit.SetLineStyle(2)
                ped_fit.Draw('same')
                sig_fit.SetLineColor(4)
                sig_fit.SetLineStyle(3)
                sig_fit.Draw('same')
                c1.Update()
                self.save_plots('BucketCut', sub_dir=self.analysis.save_dir)
                c2 = TCanvas('c2', 'c', 1000, 1000)
                self.format_histo(gr1, title='Efficiencies', x_tit='Threshold', y_tit='Efficiency', markersize=.2)
                gr1.Draw('apl')
                gr2.Draw('pl')
                leg = TLegend(.75, .8, .9, .9)
                gr5 = deepcopy(gr1)
                gr5.SetTitle('#varepsilon_{bg}')
                [leg.AddEntry(gr, gr.GetTitle(), 'l') for gr in [gr5, gr2]]
                leg.Draw()
                self.save_plots('Efficiencies', sub_dir=self.analysis.save_dir)

                c3 = TCanvas('c3', 'c', 1000, 1000)
                c3.SetGrid()
                self.format_histo(gr3, y_tit='background fraction', x_tit='excluded signal fraction', markersize=0.2, y_off=1.2)
                gr3.Draw('apl')
                gr = self.make_tgrapherrors('gr', 'working point', color=2)
                gr.SetPoint(0, 1 - sig_fit.Integral(-50, max_err) / signal.Integral(), ped_fit.Integral(-50, max_err) / ped_fit.Integral(-50, 300))
                l = self.make_tlatex(gr.GetX()[0], gr.GetY()[0] + .01, 'Working Point', color=2, size=.04)
                gr.GetListOfFunctions().Add(l)
                gr.Draw('p')
                self.save_plots('ROC_Curve', sub_dir=self.analysis.save_dir)
                self.save_plots('ROC_Curve', file_type='root', sub_dir=self.analysis.save_dir)

                c4 = TCanvas('c4', 'c', 1000, 1000)
                c4.SetGridx()
                self.format_histo(gr4, x_tit='Threshold', y_tit='1 / error', y_off=1.4)
                gr4.Draw('al')
                self.save_plots('ErrorFunction', sub_dir=self.analysis.save_dir)

                self.histos[0] = [h, gr1, gr2, c1, gr3, c2, c3, c4, gr4, a, leg, gr]
            return max_err

        threshold = func() if show else 0
        threshold = self.do_pickle(pickle_path, func, threshold)
        return threshold

    def __calc_pedestal_range(self, sigma_range):
        fit = self.analysis.show_pedestal_histo(region=self.analysis.PedestalRegion, peak_int=self.analysis.PeakIntegral, draw=False, cut='')
        sigma = fit.Parameter(2)
        mean = fit.Parameter(1)
        self.PedestalFit = fit
        return [mean - sigma_range * sigma, mean + sigma_range * sigma]

    def calc_timing_range(self, sigma):
        num = self.analysis.SignalNumber
        # estimate timing
        cut = self.generate_special_cut(excluded_cuts=['bucket'])
        draw_string = 'IntegralPeakTime[{num}]>>h1'.format(num=num)
        self.analysis.tree.Draw(draw_string, cut, 'goff')
        h1 = gROOT.FindObject('h1')
        fit1 = self.fit_fwhm(h1)
        original_mpv = fit1.Parameter(1)
        print 'mean: {0}, sigma: {1}'.format(original_mpv, fit1.Parameter(2))

        # extract timing correction
        h2 = TProfile('tcorr', 'Original Peak Position vs Trigger Cell', 1024, 0, 1024)
        self.analysis.tree.Draw('IntegralPeakTime[{num}]:trigger_cell>>tcorr'.format(num=num), cut, 'goff')
        fit2 = h2.Fit('pol2', 'qs')
        t_correction = '({p1}* trigger_cell + {p2} * trigger_cell*trigger_cell)'.format(p1=fit2.Parameter(1), p2=fit2.Parameter(2))

        # get time corrected sigma
        h3 = TH1F('h3','Corrected Timing', 100, original_mpv - 10, original_mpv + 10)
        self.analysis.tree.Draw('(IntegralPeakTime[{num}] + {t_corr}) >> h3'.format(num=num, t_corr=t_correction), cut, 'goff')
        fit3 = self.fit_fwhm(h3, draw=1)
        self.data.append(self.draw_histo(h3, 'bla', 1, self.save_directory))


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
        print ' * Time Correction: ',timing_correction

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
