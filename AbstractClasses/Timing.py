#!/usr/bin/env python
# --------------------------------------------------------
#       Pedestal analysis of the pad waveforms
# created on March 1st 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from Elementary import Elementary
from InfoLegend import InfoLegend
from Utils import *
from ROOT import TH1F, TF1, TCut, TH2F, TProfile, THStack
from numpy import pi


class TimingAnalysis(Elementary):
    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Elementary.__init__(self, verbose=self.Ana.verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.channel
        self.Tree = self.Ana.tree
        self.Cut = self.Ana.Cut
        self.TimingCut = self.Cut.generate_special_cut(excluded=['timing'], prnt=False, name='Timing')
        self.save_dir = self.Ana.save_dir
        self.Polarity = self.Ana.Polarity
        self.DiamondName = self.Ana.DiamondName
        self.DiamondNumber = self.Ana.DiamondNumber
        self.RunNumber = self.Ana.RunNumber
        self.InfoLegend = InfoLegend(pad_analysis)

    def draw_all(self):
        self.draw_peaks(show=False, prnt=False, show_cut=True)
        self.draw_peaks_tc(cut=self.get_raw_cut(), show=False, prnt=False)
        self.draw_comparison(show=False, prnt=False)
        self.draw_fine_correction(show=False, prnt=False)

    def reload_cut(self):
        self.TimingCut = self.Cut.generate_special_cut(excluded=['timing'], prnt=False, name='Timing')

    # --------------------------
    # region RUN CONFIG

    def draw_raw_peaks(self, xmin=100, xmax=400, bin_width=1, ch=None, corr=False, cut='', show=True):
        h = TH1F('h_pt', 'Max Peak Timings', int((xmax - xmin) / bin_width), xmin, xmax)
        channel = self.Ana.channel if ch is None else ch
        self.Tree.Draw('max_peak_{p}[{c}]>>h_pt'.format(c=channel, p='position' if not corr else 'time'), cut, 'goff')
        self.format_histo(h, x_tit='Time [ns]' if corr else 'Digitiser Bin', y_tit='Number of Entries', fill_color=self.FillColor)
        self.draw_histo(h, show=show)
        return h

    def create_run_config(self, ch=0, off=0):
        h0 = self.draw_raw_peaks(max(5, 100 - off), 450 - off, ch=ch, show=False)
        h0.SetName('h0')
        h1 = self.draw_raw_peaks(500 - off, 900 - off, ch=ch, show=False)
        sig_max_x = h0.GetBinCenter(h0.GetMaximumBin())
        self.format_histo(h0, x_range=[sig_max_x - 20, sig_max_x + 20])
        m = h0.GetMean()
        print 'signal_a_region = [{}, {}]'.format(int(round(m)), int(round(m)))
        bunch_width = int(20 / self.Ana.DigitiserBinWidth)
        b_start = int(round(m - bunch_width / 2, -1))
        print 'signal_b_region = [{}, {}]'.format(b_start, b_start + bunch_width)
        print 'signal_e_region = [{}, {}]'.format(b_start - bunch_width, b_start + bunch_width * 2)
        self.format_histo(h0, x_range=[b_start - bunch_width, b_start + bunch_width * 2])
        print 'pedestal_ab_region = [{}, {}]'.format(int(round(m)) - bunch_width, int(round(m)) - bunch_width)
        pul_max = h1.GetMaximum()
        x1, x2 = (h1.GetBinCenter(ibin) for ibin in [h1.FindFirstBinAbove(pul_max / 2) - 2, h1.FindLastBinAbove(pul_max / 2) + 2])
        self.format_histo(h1, x_range=[round_down_to(x1 - bunch_width / 2, 10), round_up_to(x2, 10)])
        print 'pulser_region = [{}, {}]'.format(round_down_to(x1, 10), round_up_to(x2, 10))
        print 'pedestal_ac_region = [{}, {}]'.format(round_down_to(x1 - bunch_width / 2, 10), round_down_to(x1 - bunch_width / 2, 10))
        c = self.make_canvas(x=2, y=.7)
        c.Divide(2, 1)
        for i, h in enumerate((h0, h1), 1):
            c.cd(i)
            h.Draw()

    # endregion
    # --------------------------

    # --------------------------
    # region PEAK TIMING

    def get_peak_name(self, corr, fine_corr=False, cut=None, region=None, redo=False):
        fine_corr = ' - {}'.format(self.make_fine_correction_str(cut, redo=redo)) if fine_corr else ''
        return '{}{}{}'.format(self.Ana.get_peak_name(t_corr=corr, region=region), '*{}'.format(self.Ana.DigitiserBinWidth) if not corr else '', fine_corr)

    def draw_peaks(self, fit=True, cut=None, corr=True, fine_corr=True, show=True, prnt=True, redo=False, save=True, show_cut=False):

        cut = self.TimingCut if cut is None else TCut(cut)
        pickle_path = self.make_pickle_path('Timing', 'Peak', self.RunNumber, self.DiamondNumber, suf='{}_{}{}'.format(cut.GetName(), int(corr), int(fine_corr)))
        xmin, xmax = [value * self.Ana.DigitiserBinWidth for value in self.Ana.SignalRegion]
        name = 'htrc{}{}'.format(int(corr), int(fine_corr))

        def f():
            set_root_warnings(False)
            h1 = TH1F(name, '{}Peak Positions'.format('Time Corrected ' if corr else ''), int((xmax - xmin + 20) * (8 if corr else 1 / self.Ana.DigitiserBinWidth)), xmin - 10, xmax + 10)
            self.Ana.tree.Draw('{}>>{}'.format(self.get_peak_name(corr, fine_corr, cut, redo=redo), name), cut, 'goff')
            if not h1.GetEntries():
                return
            return h1

        h = do_pickle(pickle_path, f, redo=redo)
        if h is None:
            return
        self.set_statbox(fit=fit, n_entries=5, w=.2, all_stat=not fit)
        if fit:
            self.fit_peaks(h)
        self.format_histo(h, x_tit='Time [ns]', y_tit='Number of Entries', y_off=2.0, fill_color=self.FillColor)
        set_drawing_range(h, thresh=.05 * h.GetMaximum(), lfac=.1, rfac=.7)
        prefix = 'Raw' if not corr else 'Fine' if fine_corr else ''
        self.draw_histo(h, lm=.13, show=show, prnt=prnt)
        if show_cut:
            b = self.__draw_cut(h)
            h.Draw('same')
            b.Draw('l')
        self.save_plots('{}PeakPos{}'.format(prefix, 'Fit' if fit else ''), show=show, prnt=prnt, save=save)
        return h

    def __draw_cut(self, h):
        fit = h.GetListOfFunctions()[2]
        xmin, xmax = fit.GetParameter(1) - 3 * fit.GetParameter(2), fit.GetParameter(1) + 3 * fit.GetParameter(2)
        b = self.draw_box(xmin, -10, xmax, 1e7, color=2, width=2, fillstyle=3001, name='timing', style=7)
        legend = self.make_legend(.59, y2=.43, nentries=1, margin=.45, name='la', scale=1.25)
        legend.AddEntry(b, 'cut (3 sigma)', 'lf')
        legend.Draw()
        return b

    def draw_comparison(self, show=True, prnt=True):
        stack = THStack('stack', 'Time Comparison;Time [ns];Number of Entries')
        l = self.make_legend(nentries=3, scale=.7, x1=.6)
        l_names = ['raw', 'time', 'fine']
        histos = [self.draw_peaks(show=False, corr=c, fine_corr=fc, fit=f, prnt=False) for c, fc, f in [(0, 0, 0), (1, 0, 1), (1, 1, 1)]]
        h0 = histos[0]
        h0.Scale(.25)
        l.AddEntry(h0, l_names[0], 'l')
        for i, h in enumerate(histos):
            self.format_histo(h, stats=0, color=self.get_color(), fill_color=0)
            if list(h.GetListOfFunctions()):
                fit = h.GetListOfFunctions()[-1]
                fit.SetLineColor(h.GetLineColor())
                fit.SetRange(-50, 50)
                mu = fit.GetParameter(1)
                sig = fit.GetParameter(2)
                l.AddEntry(fit, 'sigma {nam} corr.: {sig:1.2} ns'.format(nam=l_names[i], sig=sig), 'l')
                fit.SetParameter(1, 0)
            else:
                mu = h.GetMean()
            xax = h.GetXaxis()
            xax.SetLimits(xax.GetXmin() - mu, xax.GetXmax() - mu)
            stack.Add(h)
        self.format_histo(stack, x_range=[h0.GetBinCenter(h0.FindFirstBinAbove(0)), h0.GetBinCenter(h0.FindLastBinAbove(0))], draw_first=True, y_off=2.1)
        self.save_histo(stack, 'TimingComparison',  draw_opt='nostack', l=l, show=show, prnt=prnt, lm=.14)
        self.reset_colors()

    @staticmethod
    def fit_peaks(h):
        max_val = h.GetBinCenter(h.GetMaximumBin())
        fit1 = h.Fit('gaus', 'qs0', '', max_val - 2, max_val + 2)
        mean_, sigma = fit1.Parameter(1), fit1.Parameter(2)
        fit = TF1('f', 'gaus', mean_ - sigma, mean_ + sigma)
        fit.SetParameters(*fit1.Parameters())
        h.Fit('f', 'q0', '', mean_ - sigma, mean_ + sigma)
        fit2 = TF1('f1', 'gaus', mean_ - 5 * sigma, mean_ + 5 * sigma)
        fit2.SetParameters(fit.GetParameters())
        fit2.SetLineStyle(2)
        h.GetListOfFunctions().Add(fit)
        h.GetListOfFunctions().Add(fit2)
        return fit

    def draw_peaks_tc(self, corr=True, fit=True, cut=None, show=True, prnt=True, save=True, redo=False):

        cut = self.Cut.generate_special_cut(excluded=['timing'], prnt=prnt) if cut is None else TCut(cut)
        pickle_path = self.make_pickle_path('Timing', 'PeakTC', self.RunNumber, self.DiamondNumber, suf=cut.GetName())

        def f():
            set_root_warnings(False)
            h = TProfile('tcorr', '{}Peak Position vs Trigger Cell'.format('Corrected ' if corr else ''), 256, 0, 1024)
            self.Ana.tree.Draw('{}:trigger_cell>>tcorr'.format(self.get_peak_name(corr)), cut, 'goff')
            if fit:
                fit_func = TF1('fcor', '[0]*TMath::Sin([1]*(x - [2])) + [3]', -50, 1024)
                fit_func.SetParNames('Scale', 'Period', 'Phase', 'Offset')
                fit_func.SetParameters(1, 1 / 1024. * 2 * pi, 100, mean([h.GetBinContent(i) for i in xrange(h.GetNbinsX())]))
                fit_func.FixParameter(1, 1 / 1024. * 2 * pi)
                h.Fit(fit_func, 'q0')
                h.GetListOfFunctions().Add(fit_func)
            return h

        histo = do_pickle(pickle_path, f, redo=redo)
        self.format_histo(histo, x_tit='Trigger Cell', y_tit='Signal Peak Time [ns]', y_off=1.8, stats=fit)
        set_statbox(only_fit=fit, x=.7)
        self.save_histo(histo, 'OriPeakPosVsTriggerCell', show, lm=.13, prnt=prnt, save=save)
        return histo

    def get_raw_cut(self, cut=None):
        h1 = self.draw_peaks(show=False, fine_corr=False, prnt=False, cut=self.TimingCut, save=False)
        fit1 = h1.GetListOfFunctions()[2]
        return TCut('RawTiming', '({} - {}) / {} < 3'.format(self.Ana.PeakName, fit1.GetParameter(1), fit1.GetParameter(2))) + (self.TimingCut if cut is None else cut)

    def calc_fine_correction(self, cut=None, redo=False):
        h = self.draw_peaks_tc(show=False, prnt=False, cut=self.get_raw_cut(cut), redo=redo)
        fit = h.GetListOfFunctions()[0]
        fine_corr = fit.GetChisquare() / fit.GetNDF() < 100 and abs(fit.GetParameter(0)) < 10  # require decent chi2 and a meaningful scaling of the sin(x)
        return self.make_fine_correction_str(cut=self.TimingCut if cut is None else cut, fit=fit) if fine_corr else '0'

    def make_fine_correction_str(self, cut=None, fit=None, redo=False):
        fit = self.draw_peaks_tc(show=False, prnt=False, cut=cut, redo=redo).GetListOfFunctions()[0] if fit is None else fit
        return '({0} * TMath::Sin({1} * (trigger_cell - {2})))'.format(*[fit.GetParameter(i) for i in xrange(3)])

    def draw_fine_correction(self, canvas=None, show=True, prnt=True):
        set_root_warnings(False)
        p1 = TProfile('p1', 'Corrected Peak Position vs Trigger Cell', 64, 0, 1024)
        p2 = TProfile('p2', 'Corrected Peak Position vs Trigger Cell with Cut', 64, 0, 1024)
        self.Ana.tree.Draw('{}:trigger_cell>>p1'.format(self.Ana.PeakName), self.TimingCut, 'goff')
        self.Ana.tree.Draw('{} - {}:trigger_cell>>p2'.format(self.Ana.PeakName, self.calc_fine_correction()), self.get_raw_cut(), 'goff')
        y_range = increased_range([min(p1.GetMinimum(), p2.GetMinimum()), max(p1.GetMaximum(), p2.GetMaximum())], .1, .1)
        self.format_histo(p1, x_tit='Trigger Cell', y_tit='Signal Peak Times [ns]', y_off=1.8, color=self.get_color(), markersize=.5, y_range=y_range, line_color=1, stats=0, )
        self.format_histo(p2, color=self.get_color(), markersize=.5, line_color=1, stats=0)
        l = self.make_legend(nentries=2)
        l.AddEntry(p1, 'Uncorrected')
        l.AddEntry(p2, 'Fine Correction')
        self.draw_histo(p1, '', show, l=l, canvas=canvas)
        self.draw_histo(p2, '', show, draw_opt='same', canvas=get_last_canvas(), lm=.14)
        self.save_plots('FineCorrection', canvas=get_last_canvas(), prnt=prnt)

    # endregion
    # --------------------------

    # --------------------------
    # region TRIGGER CELL

    def draw_trigger_cell(self, show=True, cut=None):
        h = TH1F('tc', 'Trigger Cell', 1024, 0, 1024)
        self.Ana.tree.Draw('trigger_cell>>tc', self.Cut.all_cut if cut is None else cut, 'goff')
        self.set_statbox(fit=True, y=.4)
        self.format_histo(h, x_tit='trigger cell', y_tit='Entries', y_off=1.7, fill_color=self.FillColor, y_range=[0, h.GetMaximum() * 1.05])
        self.draw_histo(h, show=show, lm=.11)
        h.Fit('pol0', 'qs')
        self.save_plots('TriggerCell')

    def draw_intlength_vs_triggercell(self, show=True, bin_width=4):
        h = TProfile('hltc', 'Integral Length vs. Triggercell', 1024 / bin_width, 0, 1024)
        self.Ana.tree.Draw('IntegralLength[{num}]:trigger_cell>>hltc'.format(num=self.Ana.SignalNumber), self.Cut.all_cut, 'goff')
        self.set_statbox(only_fit=True, n_entries=5, w=.25)
        self.format_histo(h, x_tit='Triggercell', y_tit='Integral Length [ns]', y_off=1.6, z_tit='Number of Entries')
        self.draw_histo(h, 'IntLengthVsTriggerCell', show, lm=.12)
        fit = TF1('f', '[0] * sin([1] * (x - [2])) + [3]')
        fit.SetParNames('Scale', 'Period', 'Phase', 'Offset')
        fit.SetParLimits(0, .1, 3)
        fit.FixParameter(1, 1 / 1024. * 2 * pi)
        h.Fit(fit, 'q')
        self.save_plots('IntLengthVsTriggerCell')

    def draw_intdiff_vs_triggercell(self, show=True, prof=False):
        h = TH2F('hdtc', 'Difference of the Integral Definitions vs Triggercell', 1024 / 2, 0, 1024, 200, 0, 25)
        hprof = TProfile('hdtc_p', 'Difference of the Integral Definitions vs Triggercell', 1024 / 8, 0, 1024)
        self.Ana.tree.Draw('(TimeIntegralValues[{num}]-IntegralValues[{num}]):trigger_cell>>hdtc'.format(num=self.Ana.SignalNumber), self.Cut.all_cut, 'goff')
        self.Ana.tree.Draw('(TimeIntegralValues[{num}]-IntegralValues[{num}]):trigger_cell>>hdtc_p'.format(num=self.Ana.SignalNumber), self.Cut.all_cut, 'goff')
        gStyle.SetPalette(53)
        self.format_histo(h, x_tit='Triggercell', y_tit='Integral2 - Integral1 [au]', z_tit='Number of Entries', stats=0, y_off=1.4, z_off=1.1)
        self.ROOTObjects.append(self.draw_histo(h, '', show, draw_opt='colz', lm=.12, rm=.15))
        self.format_histo(hprof, lw=3, color=600)
        hprof.Draw('hist same') if prof else do_nothing()
        p = h.ProjectionY()
        h.GetYaxis().SetRangeUser(0, p.GetBinCenter(p.FindLastBinAbove(1)))
        self.ROOTObjects.append(hprof)
        self.save_plots('IntDiffVsTriggerCell', self.save_dir)
        gStyle.SetPalette(1)

    def draw_trigger_cell_vs_forc(self, show=True, cut=None, full_range=False, corr=False):
        if not full_range:
            self.Ana.tree.Draw('forc_pos', 'forc_pos[0]>20', 'goff')
            htemp = gROOT.FindObject('htemp')
            x = [int(htemp.GetBinCenter(htemp.FindFirstBinAbove(5000))) - 10, int(htemp.GetBinCenter(htemp.FindLastBinAbove(5000))) + 10]
        else:
            x = [0, 1024]
        h = TH2D('tcf', 'Trigger Cell vs. FORC Timing', 1024, 0, 1024, x[1] - x[0], x[0] / 2., x[1] / 2.)
        gStyle.SetPalette(55)
        forc = 'forc_pos/2.' if not corr else 'forc_time'
        self.Ana.tree.Draw('{forc}:trigger_cell>>tcf'.format(forc=forc), self.Cut.all_cut if cut is None else cut, 'goff')
        self.format_histo(h, x_tit='trigger cell', y_tit='forc timing [ns]', y_off=1.4, stats=0)
        self.save_histo(h, 'TriggerCellVsFORC{0}'.format('FullRange' if full_range else ''), show, lm=.11, draw_opt='colz', rm=.15)

    def draw_tcal(self, show=True):
        gr = self.make_tgrapherrors('gr_tcal', 'DRS4 Bin Sizes', marker_size=.5, x=range(len(self.Run.TCal)), y=self.Run.TCal)
        gr.Fit('pol0', 'qs') if show else do_nothing()
        self.format_histo(gr, x_tit='Bin number', y_tit='Length [ns]', y_off=1.5, y_range=[0, 1])
        set_statbox(only_fit=True)
        self.save_histo(gr, 'DRSBinSizes', x=1.5, y=.75, show=show)

    def draw_tcal_disto(self, show=True):
        h = TH1F('h_tcal', 'Bin Size Distribution', *self.Ana.Plots.get_tcal_bins())
        for value in self.Run.TCal:
            h.Fill(value)
        self.set_statbox(all_stat=1)
        self.format_histo(h, x_tit='Time [ns]', y_tit='Number of Entries', y_off=1.5, fill_color=self.FillColor)
        self.save_histo(h, 'DRSBinSizeDisto', show)

    # endregion
    # --------------------------

    # --------------------------
    # region BUCKET
    def _add_buckets(self, y1=None, y2=None, x2=512, avr_pos=-2, full_line=False, size=.03, ch=0):
        ymin, ymax = (get_last_canvas().GetUymin(), get_last_canvas().GetUymax()) if y1 is None else y1, y2
        start = self.Run.IntegralRegions[ch]['signal_b'][0] % 40 / 2
        bucket0 = self.Run.IntegralRegions[0]['signal_b'][0] / 40
        y_range = ymax - ymin
        self.draw_tlatex(round(x2, -1), ymax - .05 * y_range, 'Buckets', align=10, color=418, size=size)
        peak_fit = self.Run.IntegralRegions[0]['signal_a'][0] / 2.
        for i, x in enumerate(xrange(start, x2, 20), -bucket0):
            self.draw_vertical_line(x, ymin, ymax, 418, style=3 if full_line else 1, tline=True)
            if x <= x2 - 20:
                self.draw_tlatex(x + 10, ymax - .05 * y_range, str(i), align=20, color=418, size=size)
                if peak_fit:
                    pos = peak_fit % 20
                    p_pos = round_down_to(x, 20) + pos
                    p_pos += 20 if p_pos < x else 0
                    if i == avr_pos:
                        self.draw_tlatex(p_pos, ymin + 0.05 * y_range, 'Average Peak Position', color=807, size=size)
                    self.draw_arrow(p_pos, p_pos, ymin + 1, ymin + 0.04 * y_range, col=807, width=2)

    def draw_bucket_disto(self, fit=False, show=True):
        h = self.Ana.draw_signal_distribution(cut=self.Cut.get_bucket_cut(), show=show)
        # entries = h.GetEntries()
        if fit:
            self.fit_bucket(h, show)

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
    # endregion
    # --------------------------

    def draw_integration(self, event=None, region=None, integral=None, max_peak=True, pedestal=False, ped_int=False, show=True, wide=True):
        h = self.Ana.draw_single_waveform(event=event, show=show if wide else False)
        r = array(self.Ana.SignalRegion if region is None else region) * self.Ana.DigitiserBinWidth
        i = array(self.Ana.PeakIntegral if integral is None or not integral else integral) * self.Ana.DigitiserBinWidth
        if not wide:
            self.format_histo(h, lab_size=.04, tit_size=.04, y_off=1.3, x_range=[r[0] - 20, r[1] + 20])
            self.draw_histo(h, draw_opt='apl', lm=.11, show=show, gridx=True, gridy=True)
        if region or region is None and integral is not None:
            self.draw_box(r[0], -1000, r[1], 1000, color=2, style=7, fillstyle=3002, width=2)
        x_bins = xrange(*[int(v / self.Ana.DigitiserBinWidth) for v in [r[0] - 20, r[1]]])
        values = OrderedDict((j, h.GetY()[j]) for j in x_bins)
        x_vals = [h.GetX()[j] for j in x_bins]
        y_vals = values.values()
        max_x = h.GetX()[max(values, key=lambda val: abs(values[val]))]
        ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
        if max_peak and not pedestal:
            self.draw_vertical_line(max_x, -1000, 1000, color=4, w=3)
            self.draw_tlatex(max_x + 2, ymax - .05 * (ymax - ymin), 'Peak Position', color=4, align=12)
        if integral or integral is None:
            x1, x2 = max_x - i[0], max_x + i[1]
            i1, i2 = next(j for j, val in enumerate(x_vals) if val > x1), next(j for j, val in enumerate(x_vals) if val > x2) - 1
            y1, y2 = interpolate_two_points(x_vals[i1 - 1], y_vals[i1 - 1], x_vals[i1], y_vals[i1])(x1), interpolate_two_points(x_vals[i2], y_vals[i2], x_vals[i2 + 1], y_vals[i2 + 1])(x2)
            x = [x1, x1] + [x_vals[j] for j in xrange(len(x_vals) - 1) if x_vals[j] > x1 and x_vals[j - 1] < x2] + [x2, x2]
            y = [0, y1] + [y_vals[j] for j in xrange(len(x_vals) - 1) if x_vals[j] > x1 and x_vals[j - 1] < x2] + [y2, 0]
            self.draw_n_box(x=x, y=y, color=4, name='a', fillstyle=3344)
            y_pos = ymin + .1 * (ymax - ymin)
            self.draw_arrow(max_x - i[0], max_x, y_pos, y_pos, col=618, width=3, opt='<', size=.02)
            self.draw_arrow(max_x + i[1], max_x, y_pos, y_pos, col=434, width=3, opt='<', size=.02)
        if pedestal:
            x_pos, y_pos = r[0] - 10, -20
            self.draw_vertical_line(x_pos, -1000, 1000, color=4, w=3)
            self.draw_tlatex(r[0] - 8, ymin + .1 * (ymax - ymin), 'Pedestal', color=4, align=12)
            self.draw_arrow(x_pos - i[0], x_pos, y_pos, y_pos, col=618, width=3, opt='<', size=.02)
            self.draw_arrow(x_pos + i[1], x_pos, y_pos, y_pos, col=434, width=3, opt='<', size=.02)
            if ped_int:
                self.format_histo(h, y_range=[-30, 30])
                x1, x2 = x_pos - i[0], x_pos + i[1]
                i1, i2 = next(j for j, val in enumerate(x_vals) if val > x1), next(j for j, val in enumerate(x_vals) if val > x2) - 1
                y1, y2 = interpolate_two_points(x_vals[i1 - 1], y_vals[i1 - 1], x_vals[i1], y_vals[i1])(x1), interpolate_two_points(x_vals[i2], y_vals[i2], x_vals[i2 + 1], y_vals[i2 + 1])(x2)
                x = [x1, x1] + [x_vals[j] for j in xrange(len(x_vals) - 1) if x_vals[j] > x1 and x_vals[j - 1] < x2] + [x2, x2]
                y = [0, y1] + [y_vals[j] for j in xrange(len(x_vals) - 1) if x_vals[j] > x1 and x_vals[j - 1] < x2] + [y2, 0]
                self.draw_n_box(x=x, y=y, color=4, name='b', fillstyle=3344)
        self.save_histo(h, 'Integtation', draw_opt='same', canvas=get_last_canvas())

    def draw_max_peaks(self, cut=None):
        self.set_statbox(entries=True)
        h = self.draw_raw_peaks(0, 512, bin_width=.5, corr=1, show=False, cut=self.TimingCut if cut is None else TCut(cut))
        self.format_histo(h, y_off=.5, tit_size=.07, lab_size=.06)
        self.save_histo(h, 'MaxPeakTimings', logy=True, lm=.073, rm=.045, bm=.18, x=1.5, y=.5)

    def draw_fit_peak_timing(self, show=True):
        xmin, xmax = [t * self.Ana.DigitiserBinWidth for t in self.Ana.SignalRegion]
        h = TH1F('hfpt', 'Fitted Peak Positions', int((xmax - xmin) * 4), xmin, xmax)
        self.Ana.tree.Draw('fit_peak_time[{}]>>hfpt'.format(self.Ana.channel), self.Cut.all_cut, 'goff')
        self.set_statbox(all_stat=True)
        self.format_histo(h, x_tit='Time [ns]', y_tit='Number of Entries', y_off=1.3)
        self.draw_histo(h, show=show, lm=.12)

    def draw_peaking_time(self, show=True):
        h = TH1F('hpt', 'Peaking Time', 400, 0, 20)
        self.Ana.tree.Draw('peaking_time[{}]>>hpt'.format(self.Ana.channel), self.Cut.all_cut, 'goff')
        self.set_statbox(all_stat=True)
        self.format_histo(h, x_tit='Time [ns]', y_tit='Number of Entries', y_off=1.3)
        self.draw_histo(h, show=show, lm=.12)

    def draw_forc_times(self, show=True, corr=False):
        self.Ana.tree.Draw('forc_pos', 'forc_pos[0]>20', 'goff')
        htemp = gROOT.FindObject('htemp')
        x = [int(htemp.GetBinCenter(htemp.FindFirstBinAbove(5000))) - 10, int(htemp.GetBinCenter(htemp.FindLastBinAbove(5000))) + 10]
        h = TH1F('ft', 'FORC Timing', x[1] - x[0], x[0] / 2., x[1] / 2.)
        forc = 'forc_pos/2.' if not corr else 'forc_time'
        self.Ana.tree.Draw('{forc}>>ft'.format(forc=forc), self.Cut.all_cut, 'goff')
        self.format_histo(h, x_tit='Time [ns]', y_tit='Entries', y_off=2, fill_color=self.FillColor)
        self.ROOTObjects.append(self.save_histo(h, 'FORCTiming', show, sub_dir=self.save_dir, lm=.14))

    def fit_rf(self, evnt=0, t_corr=True, show=True):
        gr, n = self.Ana.draw_waveforms(start_event=evnt, t_corr=t_corr, show=show, channel=self.Run.get_rf_channel(), cut='')
        set_statbox(only_fit=True, entries=10)
        fit = TF1('f', '[0] * TMath::Sin((x+[1])*2*pi/[2])+[3]', 0, 500)
        fit.SetParNames('Scale', 'Period', 'Phase', 'Offset')
        fit.SetNpx(1000)
        fit.SetParameters(100, 20, 3, -40)
        fit.SetParLimits(2, -20, 20)
        fit_res = gr.Fit('f', 'qs{}'.format('' if show else 0))
        print evnt, fit.GetChisquare() / fit.GetNDF()
        return FitRes(fit_res)

    def draw_rf_phases(self, n=100, start_event=0, show=True):
        g = self.make_tgrapherrors('grf', 'RF Phases')
        for i in xrange(n):
            fit = self.fit_rf(start_event + i, show=False)
            get_last_canvas().Update()
            g.SetPoint(i, i, fit.Parameter(2))
        self.format_histo(g, x_tit='event', y_tit='Number of Entries')
        self.save_histo(g, 'RFPhases', show=show)

    def draw_rf_period(self, cut=None, show=True):
        cut = self.Cut.all_cut if cut is None else TCut(cut)
        h = TH1F('hrfp', 'Beam Period', 500, 19.7, 19.8)
        self.Tree.Draw('rf_period>>hrfp', cut, 'goff')
        set_statbox(entries=4)
        self.format_histo(h, x_tit='RF Period [ns]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor, ndivx=507)
        self.save_histo(h, 'RFPeriod', lm=.12, show=show)

    def draw_rf_phase(self, cut=None, show=True):
        cut = self.Cut.all_cut + TCut('rf_chi2<100') if cut is None else TCut(cut)
        h = TH1F('hrfph', 'RF Phase', 500, -30, 30)
        self.Tree.Draw('rf_phase>>hrfph', cut, 'goff')
        set_statbox(entries=4)
        self.format_histo(h, x_tit='RF Period [ns]', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.save_histo(h, 'RFPeriod', lm=.12, show=show)

    def draw_signal_peak(self, cut=None, corr=False, show=True):
        h = TH1F('hspt', 'Signal Peak Timings', 2000 * (2 if corr else 1), 0, 500)
        cut = self.Cut.generate_special_cut(excluded=['timing']) if cut is None else TCut(cut)
        corr = '+rf_phase' if corr else ''
        self.Tree.Draw('signal_peak_time[{ch}]{c}>>hspt'.format(ch=self.Ana.channel, c=corr), cut, 'goff')
        set_statbox(w=.3, entries=6)
        x_range = increased_range([h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(5), h.FindLastBinAbove(5)]], .1, .3)
        self.format_histo(h, x_tit='Signal Peak Timing [ns]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor, x_range=x_range)
        self.draw_histo(h, show=show, lm=.12)

    def draw_rf_vs_peak(self, show=True):
        h = TH2F('hsprf', 'RF Phase vs. Peak Timings', 4000, 0, 500, 240, -12, 12)
        cut = self.Cut.all_cut + TCut('rf_chi2<100')
        self.Tree.Draw('rf_phase:signal_peak_time[{ch}]>>hsprf'.format(ch=self.Ana.channel), cut, 'goff')
        x_range = increased_range([h.GetXaxis().GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(1), h.FindLastBinAbove(5)]], .2, .2)
        y_range = increased_range([h.GetYaxis().GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(5, 2), h.FindLastBinAbove(5, 2)]], .2, .2)
        self.format_histo(h, x_tit='Signal Peak Timing [ns]', y_tit='RF Phase [ns]', z_tit='Number of Entries', z_off=1.2, y_off=1.2, stats=0, x_range=x_range, y_range=y_range)
        self.save_histo(h, 'RFSignalPeak', rm=.18, draw_opt='colz', show=show)

    def draw_inflexion_time(self, corr=False, channel=None, show=True):
        channel = self.Ana.channel if channel is None else channel
        h = TH1F('hrt', 'Inflexion Time', 2000 * (2 if corr else 1), 0, 500)
        corr = '+rf_phase' if corr else ''
        self.Tree.Draw('rise_time[{ch}]{c}>>hrt'.format(ch=channel, c=corr), self.Cut.all_cut + TCut('rise_time[{ch}]'.format(ch=channel)), 'goff')
        set_statbox(w=.3, entries=6)
        x_range = increased_range([h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(5), h.FindLastBinAbove(5)]], .1, .3)
        self.format_histo(h, x_tit='Rise Time [ns]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor, x_range=x_range)
        self.draw_histo(h, show=show, lm=.12)

    def draw_scint_inflexion(self, corr=False, show=True):
        return self.draw_inflexion_time(corr=corr, channel=self.get_scint_channel(), show=show)

    def draw_peak_width(self, show=True):
        h = TH1F('hpw', 'Peak Width', 500, 0, 1)
        self.Tree.Draw('abs(rise_width[{ch}])>>hpw'.format(ch=self.Ana.channel), self.Cut.all_cut, 'goff')
        set_statbox(w=.3, entries=6)
        x_range = increased_range([h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(5), h.FindLastBinAbove(5)]], .1, .3)
        self.format_histo(h, x_tit='Rise Width [au]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor, x_range=x_range)
        self.draw_histo(h, show=show, lm=.12)

    def draw_threshold(self, corr=False, channel=None, show=True):
        channel = self.Ana.channel if channel is None else channel
        h = TH1F('hrt', 'Signal over Threshold', 2000 * (2 if corr else 1), 0, 500)
        corr = '+rf_phase' if corr else ''
        self.Tree.Draw('t_thresh[{ch}]{c}>>hrt'.format(ch=channel, c=corr), self.Cut.all_cut + TCut('t_thresh[{ch}] > 10'.format(ch=channel)), 'goff')
        set_statbox(w=.3, entries=6)
        x_range = increased_range([h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(5), h.FindLastBinAbove(5)]], .1, .3)
        self.format_histo(h, x_tit='Threshold Time [ns]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor, x_range=x_range)
        self.draw_histo(h, show=show, lm=.12)

    def draw_scint_threshold(self, corr=False, show=True):
        return self.draw_threshold(corr=corr, channel=self.get_scint_channel(), show=show)

    def draw_inter_dia_corr(self, show=True):
        h = TH2F('hidc', 'Inflextion Times of Diamond Signals', 4000, 0, 500, 4000, 0, 500)
        cut = self.Cut.all_cut + TCut('rise_time[{c1}]'.format(c1=self.Run.Channels[0]))
        self.Tree.Draw('rise_time[{c1}]:rise_time[{c2}]>>hidc'.format(c1=self.Run.Channels[0], c2=self.Run.Channels[1]), cut, 'goff')
        x_range = increased_range([h.GetXaxis().GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(1), h.FindLastBinAbove(1)]], .2, .2)
        y_range = increased_range([h.GetYaxis().GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(1, 2), h.FindLastBinAbove(1, 2)]], .2, .2)
        self.format_histo(h, x_tit='Inflexion Time1 [ns]', y_tit='Inflexion Time2 [ns]', z_tit='Number of Entries', z_off=1.2, y_off=1.2, stats=0, x_range=x_range, y_range=y_range)
        self.save_histo(h, 'InterDiaCorrelation', rm=.18, draw_opt='colz', show=show)

    def draw_inter_dia(self, show=True):
        h = TH1F('hid', 'Inter Diamond Timing Correction', 200, -10, 10)
        cut = self.Cut.all_cut + TCut('rise_time[{c1}]'.format(c1=self.Run.Channels[0]))
        self.Tree.Draw('rise_time[{c1}] - rise_time[{c2}]>>hid'.format(c1=self.Run.Channels[0], c2=self.Run.Channels[1]), cut, 'goff')
        x_range = increased_range([h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(5), h.FindLastBinAbove(5)]], .1, .3)
        set_statbox(entries=4)
        self.format_histo(h, x_tit='Timing Correction [ns]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor, x_range=x_range)
        self.draw_histo(h, show=show, lm=.12)

    def get_channel_nr(self, name):
        try:
            return self.Run.DigitizerChannels.index(next(ch for ch in self.Run.DigitizerChannels if name in ch.lower()))
        except StopIteration:
            log_warning('There is no Digitiser Channel with {n} in it'.format(n=name))

    def get_rf_channel(self):
        return self.get_channel_nr('rf')

    def get_scint_channel(self):
        return self.get_channel_nr('scint')
