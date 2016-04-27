from Cut import Cut
from Extrema import Extrema2D
from copy import deepcopy
from ROOT import TCut, TH1F, TF1, TCanvas, TLegend
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

    def load_abs_median_high(self, value):
        if value > 0:
            self.EasyCutStrings['absMedian_high'] = '|median|<{high}'.format(high=value)
            return value
        else:
            return -1

    def set_abs_median_high(self, high):
        self.CutConfig['absMedian_high'] = self.load_abs_median_high(high)

    def load_pedestal_sigma(self, value):
        if value > 0:
            self.EasyCutStrings['pedestalsigma'] = 'PedSigma' + str(value)
            return value
        else:
            self.EasyCutStrings['pedestalsigma'] = ''
            return -1

    def set_pedestal_sigma(self, sigma=-1):
        self.CutConfig['pedestalsigma'] = self.load_pedestal_sigma(sigma)

    def set_peak_value_pos(self, x_min, x_max):
        assert 0 <= x_min <= 1024, 'min signal peak has to be in [0, 1024], not "{min}"'.format(min=x_min)
        assert 0 <= x_max <= 1024, 'max trigger cell has to be in [0, 1024], not "{max}"'.format(max=x_max)
        if x_max < x_min:
            x_min, x_max = x_max, x_min
        self.CutConfig['signal_peak_pos'] = [x_min, x_max]
        self.EasyCutStrings['SignalPeakPos'] = 'Signal Peak in {0}'.format([x_min, x_max])

    def set_trigger_cell(self, x_min, x_max):
        assert 0 <= x_min <= 1024, 'min trigger cell has to be in [0, 256], not "{min}"'.format(min=x_min)
        assert 0 <= x_max <= 1024, 'max trigger cell has to be in [0, 256], not "{max}"'.format(max=x_max)
        if x_max < x_min:
            x_min, x_max = x_max, x_min
        self.CutConfig['trigger_cell'] = [x_min, x_max]
        self.EasyCutStrings['TriggerCell'] = 'Trigger Cell in {0}'.format([x_min, x_max])
    # endregion

    # ==============================================
    # region GENERATE CUT STRINGS
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

    def generate_signal_peak_pos(self):
        lst = self.CutConfig['signal_peak_pos']
        if lst:
            self.CutStrings['signal_peak_pos'] += 'IntegralPeaks[{num}] < {max} && IntegralPeaks[{num}] >= {min}'.format(num=self.analysis.SignalNumber, min=lst[0], max=lst[1])

    def generate_trigger_cell(self):
        lst = self.CutConfig['trigger_cell']
        if lst:
            self.CutStrings['trigger_cell'] += 'trigger_cell < {max} && trigger_cell >= {min}'.format(min=lst[0], max=lst[1])

    def add_signal_peak_pos_cut(self, xmin, xmax):
        self.set_peak_value_pos(xmin, xmax)
        self.CutStrings['signal_peak_pos'].SetTitle('')
        self.generate_signal_peak_pos()
        self.all_cut = self.generate_all_cut()

    def add_trigger_cell_cut(self, xmin, xmax):
        self.set_trigger_cell(xmin, xmax)
        self.CutStrings['trigger_cell'].SetTitle('')
        self.generate_trigger_cell()
        self.all_cut = self.generate_all_cut()

    def update_all_cut(self):
        self.all_cut = self.generate_all_cut()

    def update_bucket(self, threshold=None):
        threshold = self.calc_signal_threshold(show=False) if threshold is None else threshold
        self.reset_cut('bucket')
        self.CutStrings['bucket'] += self.generate_bucket(threshold)
        self.update_all_cut()

    def generate_old_bucket(self):
        # only generate the cut if the region e2 exists! todo: find a smarter solution for that!
        try:
            sig2 = self.analysis.get_signal_name('e', 2)
            string = '{sig2}=={sig1}'.format(sig2=sig2, sig1=self.analysis.SignalName)
            return TCut(string)
        except:
            return TCut('')

    def generate_bucket(self, threshold=None):
        sig = self.analysis.SignalName
        threshold = self.calc_signal_threshold(show=False) if threshold is None else threshold
        string = '!(!({old_buck})&&({sig}<{thres}))'.format(sig=sig, thres=threshold, old_buck=self.CutStrings['old_bucket'])
        cut =  TCut(string) if self.CutStrings['old_bucket'].GetTitle() else TCut('')
        return cut

    def add_signal_peak_time_cut(self, low, up):
        self.reset_cut('signal_peak_time')
        self.CutStrings['signal_peak_time'] += self.generate_signal_peak_time(low, up)
        self.update_all_cut()

    def generate_signal_peak_time(self, low = 67.068, up = 72.750):
        string = '{low}<IntegralPeakTime[{num}]&&IntegralPeakTime[{num}]<{up}'.format(low=low, up=up, num=self.analysis.SignalNumber)
        return TCut(string)

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

    def __calc_pedestal_range(self):
        fit = self.analysis.show_pedestal_histo(region=self.analysis.PedestalRegion, peak_int=self.analysis.PeakIntegral, draw=False, cut='')
        sigma = fit.Parameter(2)
        mean = fit.Parameter(1)
        self.PedestalFit = fit
        sigma_range = self.CutConfig['pedestalsigma']
        return [mean - sigma_range * sigma, mean + sigma_range * sigma]

    def generate_channel_cutstrings(self):
        # -- SATURATED CUT --
        self.CutStrings['saturated'] += '!is_saturated[{ch}]'.format(ch=self.channel)

        # -- SPREAD LOW CUT --
        if self.CutConfig['spread_low'] > 0:
            self.CutStrings['spread_low'] += 'sig_spread[{ch}]>{low}'.format(ch=self.channel, low=int(self.CutConfig['spread_low']))
            self.EasyCutStrings['spread_low'] = 'spread>{low}'.format(low=int(self.CutConfig['spread_low']))

        # -- MEDIAN CUT --
        if self.CutConfig['absMedian_high'] > 0:
            print 'CONFIG:', self.CutConfig['absMedian_high'], type(self.CutConfig['absMedian_high'])
            self.CutStrings['median'] += 'abs(median[{ch}])<{high}'.format(ch=self.channel, high=int(self.CutConfig['absMedian_high']))
            self.EasyCutStrings['absMedian_high'] = '|median|<{high}'.format(high=int(self.CutConfig['absMedian_high']))

        # -- PEDESTAL SIGMA CUT --
        if self.CutConfig['pedestalsigma'] > 0:
            ped_range = self.__calc_pedestal_range()
            self.CutStrings['ped_sigma'] += '{ped}>{min}&&{ped}<{max}'.format(ped=self.analysis.PedestalName, min=ped_range[0], max=ped_range[1])
            self.EasyCutStrings["pedestalsigma"] = "PedSigma<" + str(self.CutConfig['pedestalsigma'])

        # --PEAK POSITION TIMING--
        # todo: add a method that fits the real time disto and sets the cut to 4 sigma!
        # self.CutStrings['signal_peak_time'] += self.generate_signal_peak_time()

        # --BUCKET --
        self.CutStrings['old_bucket'] += self.generate_old_bucket()
        self.CutStrings['bucket'] += self.generate_bucket()

    # endregion

    def generate_pulser_cut(self, beam_on=True):
        cut = self.CutStrings['ped_sigma'] + self.CutStrings['event_range'] + self.CutStrings['saturated']
        cut.SetName('Pulser{0}'.format('BeamOn' if beam_on else 'BeamOff'))
        cut += self.CutStrings['beam_interruptions'] if beam_on else '!({0})'.format(self.JumpCut)
        cut += '!({0})'.format(self.CutStrings['pulser'])
        return cut
