#! /usr/bin/env python
from ROOT import TExec

from Extrema import Extrema2D
from PulserCollection import PulserCollection
from analysis_collection import *
from pad_analysis import PadAnalysis


class PadCollection(AnalysisCollection):
    """ Analysis of the various runs of a single runplan. """

    def __init__(self, run_plan, dut_nr, test_campaign=None, load_tree=True, verbose=False):

        self.DUTType = 'pad'
        AnalysisCollection.__init__(self, run_plan, dut_nr, test_campaign, load_tree, verbose)

        self.Channel = self.FirstAnalysis.Channel

        # Sub Classes
        self.Pulser = PulserCollection(self)

    def load_analysis(self, run_number):
        return PadAnalysis(run_number, self.DUT.Number, self.TCString, self.Threads[run_number].Tuple, self.Threads[run_number].Time, self.Verbose, prnt=False)

    @staticmethod
    def load_dummy():
        return PadAnalysis

    def generate_threshold_pickle(self):
        if not file_exists(self.make_pickle_path('Cuts', 'SignalThreshold', run=self.MaxFluxRun, ch=self.DUT.Number)):
            PadAnalysis(self.MaxFluxRun, self.DUT.Number, self.TCString, self.Verbose, prnt=False)

    # ----------------------------------------
    # region RESULTS
    def print_all_off_results(self):
        string = 'Signal\tPedest.\tPulser\n'
        for ana in self.Analysis.itervalues():
            string += '{0}\n'.format(ana.print_off_results(prnt=False))
        print string

    def draw_all(self, redo=False):
        t0 = self.info('Generate all plots ... ')
        old_verbose = self.FirstAnalysis.Verbose
        self.set_verbose(False)
        if self.Type == 'voltage scan':
            self.VoltageScan.draw_all(redo)
        else:
            self.draw_pulse_heights(show=False, redo=redo, prnt=False)
            self.draw_scaled_pulse_heights(show=False)
            self.draw_scaled_pulse_heights(show=False, vs_time=True)
            self.Pulser.draw_pulse_heights(show=False, redo=redo)
            self.draw_pedestals(show=False, redo=redo)
            self.draw_pedestals(show=False, sigma=True, redo=redo)
            self.draw_pulser_pedestals(show=False, redo=redo)
            self.draw_pulser_pedestals(show=False, redo=redo, sigma=True)
        self.draw_currents(show=False, draw_opt='al')
        self.draw_flux(show=False)
        self.draw_ph_with_currents(show=False)
        self.draw_signal_distributions(show=False, redo=redo)
        self.draw_signal_maps(redo=redo)
        self.draw_hitmaps(redo=redo)
        self.draw_run_currents()
        self.draw_chi2s()
        self.draw_angles()
        self.draw_occupancies()
        self.draw_timing()
        self.set_verbose(old_verbose)
        print_elapsed_time(t0)
    # endregion RESULTS
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_pedestals(self, pulser=False, redo=False, runs=None):
        means = self.get_run_values('{}pedestals'.format('pulser ' if pulser else ''), self.Analysis.get_pedestal, runs=runs, pulser=pulser, par=1, redo=redo)
        sigmas = self.get_run_values('', self.Analysis.get_pedestal, runs=runs, pbar=False, pulser=pulser, par=2)
        return array(means), array(sigmas)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL/PEDESTAL
    @staticmethod
    def scale_current_gr(gr):
        vals = [gr.GetY()[i] for i in xrange(gr.GetN())]
        mi, ma, mn = min(vals), max(vals), mean(vals)
        y = [mi, mn + 3 * (mn - mi) if ma - mn > 50 else ma]
        format_histo(gr, name='cur', color=899)
        return increased_range(y, .3, .2)

    def draw_ph_with_currents(self, show=True, scale=1):
        ph = self.get_pulse_height_graph(vs_time=True, first_last=False, bin_width=10000, legend=False)
        self.Currents.set_graphs()
        cur = self.Currents.CurrentGraph.Clone()
        cur_range = self.scale_current_gr(cur)
        pul = self.Pulser.get_pulse_height_graph(vs_time=True, legend=False, show_flux=False)
        pul_fac = pul.GetListOfGraphs()[0].GetY()[0] if scale else ph.GetListOfGraphs()[0].GetY()[0] / pul.GetListOfGraphs()[0].GetY()[0]
        ph_fac = ph.GetListOfGraphs()[0].GetY()[0] if scale else sign(self.DUT.Bias)
        scale_multigraph(pul, ph.GetListOfGraphs()[0].GetY()[0] if scale is None else scale)
        scale_multigraph(ph, scale)

        # legends
        entries = [2, 3, 1]
        lm = .09
        positions = [[lm + .02, .95], [lm + .02, .95], [lm + .02, .95]]
        legends = [self.make_legend(*positions[i], nentries=entries[i], scale=2.5 * (27 / 36. if i == 2 else 1), x2=lm + .22) for i in xrange(3)]
        dummy_gr = self.make_tgrapherrors('g', 'g', width=2)
        legends[0].AddEntry(pul.GetListOfGraphs()[0], 'scaled pulser data', 'p')
        legends[0].AddEntry(0, 'factor {sgn}{fac:4.0f}'.format(fac=pul_fac, sgn='-' if self.FirstAnalysis.PulserPolarity < 0 else ''), '')
        legends[1].AddEntry(ph.GetListOfGraphs()[1], '{sgn:2.0f} * (signal - ped.) data'.format(sgn=ph_fac), 'p')
        legends[1].AddEntry(0, 'flux in kHz/cm^{2}', 'p')
        legends[1].AddEntry(dummy_gr, 'duration', 'l')
        legends[2].AddEntry(cur, 'current', 'l')

        # Drawing
        lab_size = .12
        for f in ph.GetListOfGraphs()[0].GetListOfFunctions():
            f.SetTextSize(.09)
            f.SetY(f.GetY() - .01)
        format_histo(ph, draw_first=True, lab_size=lab_size, tit_size=lab_size, y_off=.33)
        format_histo(cur, x_tit='Time [hh:mm]', lab_size=lab_size * 27 / 36., tit_size=lab_size * 27 / 36., y_off=.33 * 36 / 27.)
        format_histo(pul, color=859, draw_first=True, lab_size=lab_size, tit_size=lab_size, y_off=.33)
        format_histo(pul.GetListOfGraphs()[0], color=859)
        draw_opts = ['a', 'a', '']
        y_tits = ['Pulser ', 'Signal', 'Current [nA] ']
        y_ranges = [increased_range(scale_margins(pul, ph), .3), increased_range(scale_margins(ph, pul), .3), cur_range]
        y_ranges = [[.94, 1.06], [.94, 1.06], cur_range] if scale else y_ranges
        divs = [504, 504, None]
        ypos = [[1, .73], [.73, .46], [.46, .1]]
        margins = [[lm, .05, 0, 0], [lm, .05, 0, 0], [lm, .05, 9 / 36., 0]]
        x_range = increased_range([cur.GetX()[0], cur.GetX()[cur.GetN() - 1]], .1, .1)
        close_last_canvas()

        c = self.make_canvas('c_phc', 'c', 1.5, 1, transp=True, show=show)

        for i, gr in enumerate([pul, ph, cur]):
            format_histo(gr, title='', y_range=y_ranges[i], y_tit=y_tits[i], center_y=True, ndivy=divs[i], t_ax_off=self.Currents.Time[0])
            self.draw_tpad('p{}'.format(i), 'p{}'.format(i), pos=[0, ypos[i][0], 1, ypos[i][1]], gridx=True, gridy=True, margins=margins[i], transparent=True)
            gr.GetXaxis().SetLimits(*x_range)
            gr.Draw(draw_opts[i])
            legends[i].Draw()
            c.cd()

        self.save_plots('PhPulserCurrent', all_pads=False)
        self.Objects.append([ph, cur, pul, c, legends])

    def draw_pedestals(self, vs_time=False, sigma=False, pulser=False, redo=False, show=True):
        y_tit = 'Noise' if sigma else 'Pedestal'
        g = self.make_tgrapherrors('gps', '{} {}'.format('Pulser ' if pulser else '', y_tit), x=self.get_x_var(vs_time), y=self.get_pedestals(pulser, redo=redo)[1 if sigma else 0])
        format_histo(g, color=810, x_tit=self.get_x_tit(vs_time), y_tit='{} [mV]'.format(y_tit), y_off=1.45, x_off=1.3, t_ax_off=self.get_tax_off(vs_time), x_range=self.Bins.FluxRange)
        save_name = '{p}Pedestal{s}{mod}'.format(mod=self.get_mode(vs_time), s='Noise' if sigma else '', p='Pulser' if pulser else '')
        self.save_histo(g, save_name=save_name, show=show, logx=False if vs_time else True, lm=.12, draw_opt='ap')
        return g

    def draw_ped_spread(self, pulser=False, redo=False, show=True):
        values = self.get_pedestals(pulser, redo=redo)[0][self.get_fluxes().argsort()]  # sort pedestal by ascending fluxes
        h = TH1F('hps', 'Relative Pedestal Spread', 20, -.5, .5)
        for lst in split(values, self.get_flux_splits(show=False)):
            for value in lst:
                h.Fill((value - mean(lst)).n)
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Relative Pedestal', y_tit='Number of Entries', y_off=1.2)
        self.save_histo(h, 'PedestalSpread', lm=.11, show=show)

    def draw_noise_spread(self, pulser=False, redo=False, show=True):
        h = TH1F('hns', 'Relative Noise Spread', 20, -.3, .3)
        values = self.get_pedestals(pulser, redo=redo)[1]
        h.FillN(self.NRuns, array([(v - mean(values)).n for v in values], 'd'), full(self.NRuns, 1, 'd'))
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Relative Noise', y_tit='Number of Entries', y_off=1.2)
        self.save_histo(h, 'NoiseSpread', lm=.11, show=show)

    def draw_noise(self, vs_time=False, show=True):
        return self.draw_pedestals(vs_time=vs_time, show=show, sigma=True)

    def draw_pulser_pedestals(self, show=True, redo=False, sigma=False):
        self.draw_pedestals(pulser=True, show=show, redo=redo, sigma=sigma)

    def draw_snrs(self, flux=True, draw=True):
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        gr = self.make_tgrapherrors('gr', 'SNR vs {mode}'.format(mode=mode))
        i = 0
        for key, ana in self.Analyses.iteritems():
            snr = ana.calc_snr()
            x = ana.Run.Flux if flux else key
            gr.SetPoint(i, x, snr[0])
            gr.SetPointError(i, 0, snr[1])
            i += 1
        if draw:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        if flux:
            c.SetLogx()
        gr.Draw('ap')
        self.save_plots('AllSNRs', canvas=c)
        self.Objects.append([gr, c])
        gROOT.SetBatch(0)

    def calc_pedestal_spread(self, flux=150):
        values = self.get_pedestals(runs=self.get_runs_below_flux(flux))[0]
        return max(values) - min(values)

    def calc_full_pedestal_spread(self):
        values = self.get_pedestals()[0]
        return max(values) - min(values)

    def draw_signal_legend(self):
        sig = 'positive' if self.FirstAnalysis.Polarity > 0 else 'negative'
        l1 = self.make_legend(.17, .88, nentries=2, margin=.05, clean=True, x2=.57 if self.Legend else .52, cols=2)
        l1.AddEntry(0, 'Signal Polarity:', '')
        l1.AddEntry(0, sig, '').SetTextAlign(12)
        l1.AddEntry(0, 'Pedestal Substraction:', '')
        l1.AddEntry(0, 'yes', '').SetTextAlign(12)
        l1.Draw()
        self.Objects.append(l1)
        return l1

    def compare_pedestals(self, show=True):
        gr1 = self.draw_pedestals(show=False)
        gr2 = self.draw_pedestals(pulser=True, show=False)
        graphs = [gr1, gr2]
        margins = find_graph_margins(graphs)
        gROOT.SetBatch(0) if show else gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Pedestal Comparison', 1000, 1000)
        c.SetLogx()
        legend = TLegend(.7, .78, .88, .88)
        names = ['Signal', 'Pulser', 'BeamOff']
        for i, gr in enumerate(graphs):
            gr.GetYaxis().SetRangeUser(*margins)
            format_histo(gr, color=self.get_color())
            gr.Draw('lp') if i else gr.Draw('alp')
            legend.AddEntry(gr, names[i], 'pl')
        legend.Draw()
        gROOT.SetBatch(0)
        self.save_plots('PulserPedestalComparison')
        self.Objects.append([c, graphs, legend])

    def compare_signal_vs_peak_height(self, i0=0, i1=-1, ym=.05, cft=False, show=True, redo=False):
        def f():
            func = self.Analysis.draw_signal_vs_cft if cft else self.Analysis.draw_signal_vs_peaktime
            x0, y0, x1, y1 = [j for ind in [i0, i1] for j in get_hist_vecs(func(self.get_ana(ind), show=False))]
            y0 = where(y0 == 0, 1e10, y0)
            return x0, y1 / y0
        x, y = do_pickle(self.make_simple_pickle_path('SigPeakRatio', int(cft), sub_dir='Peaks', dut='{}{}'.format(i0, i1)), f, redo=redo)
        flux0, flux1 = [make_flux_string(self.get_ana(i).get_flux()) for i in [i0, i1]]
        g = self.make_tgrapherrors('gcspt', 'Signal Ratio Vs Peak Time at {} and {}'.format(flux0, flux1), x=x, y=y)
        format_histo(g, x_tit='Signal Peak Time [ns]', y_tit='Signal Ratio', y_off=1.7, y_range=array([-ym, ym]) + 1)
        self.draw_histo(g, show=show, lm=.13, gridy=True)
        y = y[y > .1]
        m, s = mean_sigma(y)
        return ufloat(m, s / sqrt(y.size))

    def compare_signal_vs_cft(self, i0=0, i1=-1, ym=.05, show=True, redo=False):
        return self.compare_signal_vs_peak_height(i0, i1, ym, True, show, redo)

    def compare_all_sig_vs_peakheight(self, ym=.05, show=True):
        values = []
        self.PBar.start(self.NRuns - 1)
        for i in range(1, self.NRuns):
            values.append(self.compare_signal_vs_peak_height(0, i, show=False))
            self.PBar.update()
        g = self.make_tgrapherrors('gasphr', 'Signal Ratio Vs Flux', x=self.get_fluxes()[1:], y=values)
        format_histo(g, y_tit='Signal Ratio', y_off=1.7, y_range=array([-ym, ym]) + 1, **self.get_x_args(False))
        self.draw_histo(g, show=show, lm=.13, logx=True, gridy=True)

    # endregion SIGNAL/PEDESTAL
    # ----------------------------------------

    # ----------------------------------------
    # region CUTS
    def draw_bucket_info(self, flux=True, show=True, mean_=True):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        gr1 = self.make_tgrapherrors('gr1', '', color=self.get_color())
        prefix = 'Number of Bucket Cut Events' if not mean_ else 'Mean Pulse Height with Different Bucket Cuts'
        gr2 = self.make_tgrapherrors('gr2', '{pref} vs {mod}'.format(pref=prefix, mod=mode), color=self.get_color())
        gr3 = self.make_tgrapherrors('gr3', '', color=self.get_color())
        i = 0
        for key, ana in self.Analyses.iteritems():
            x = ana.Run.Flux if flux else key
            if not mean_:
                n = ana.show_bucket_numbers(show=False)
                gr1.SetPoint(i, x, n['new'] / n['all'] * 100)
                gr2.SetPoint(i, x, n['old'] / n['all'] * 100)
            else:
                data = ana.show_bucket_means(show=False, plot_histos=False)
                gr1.SetPoint(i, x, data['new'][0])
                gr2.SetPoint(i, x, data['old'][0])
                gr3.SetPoint(i, x, data['no'][0])
                gr1.SetPointError(i, 0, data['new'][1])
                gr2.SetPointError(i, 0, data['old'][1])
                gr3.SetPointError(i, 0, data['no'][1])
            i += 1
        c = TCanvas('c', 'Bucket Numbers', 1000, 1000)
        c.SetLeftMargin(.13)
        if flux:
            c.SetLogx()
        format_histo(gr2, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='Events [%]' if not mean_ else 'Mean [au]', y_off=1.7, color=None)
        gr2.Draw('apl')
        gr1.Draw('pl')
        if mean_:
            gr3.Draw('pl')
        leg = TLegend(.2, .8, .35, .9)
        leg.AddEntry(gr2, 'old cut', 'pl')
        leg.AddEntry(gr1, 'new cut', 'pl')
        if mean_:
            leg.AddEntry(gr3, 'no bucket', 'pl')
        leg.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('{mode}_' + mode)
        self.Objects.append([c, gr1, gr2, gr3, leg])
    # endregion CUTS
    # ----------------------------------------

    # ----------------------------------------
    # region TIMING
    def draw_peak_timings(self, vs_time=False, redo=False, show=True):
        # TODO add pulser
        """ Shows the means of the signal peak distributions. """
        g = self.make_tgrapherrors('gr', 'Peak Times', x=self.get_x_var(vs_time), y=self.get_values('peak timings', self.Analysis.get_peak_timing, redo=redo))
        format_histo(g, y_tit='Peak Time [ns]', y_off=1.9, **self.get_x_args(vs_time))
        self.save_histo(g, 'PeakTimes', lm=.14, logx=not vs_time, show=show)

    def show_peak_distribution(self, show=True):
        """ Shows the positions of the peaks of the 2D map. """
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        # create an overall VotingHistogram
        ana = self.FirstAnalysis
        ana.draw_mean_signal_distribution(show=False)
        extrema = Extrema2D(ana.SignalMapHisto, ana.MeanSignalHisto)
        h = extrema.create_voting_histo()
        for run, ana in self.Analyses.iteritems():
            if ana.IsAligned:
                ana.find_2d_extrema(histo=h, show=False)
            else:
                print 'Run {run} is not aligned...'.format(run=run)
        c = TCanvas('c', 'Voting Histos', 1600, 800)
        c.Divide(2, 1)
        # new_pal = ar.array('i', [kYellow, kYellow, kOrange, kOrange - 3, kOrange + 7, kRed])
        ex = [TExec('ex1', 'gStyle->SetPalette(56);'), TExec('ex2', 'gStyle->SetPalette(51)')]
        if show:
            gROOT.SetBatch(0)
        for i, histo in enumerate(h.itervalues(), 1):
            c.cd(i)
            histo.Draw('col')
            ex[i - 1].Draw()
            histo.Draw('colz same')
        self.save_plots('PeakDistribution')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.Objects.append([c, ex])
    # endregion TIMING
    # ----------------------------------------

    # ----------------------------------------
    # region PEAKS
    def save_additional(self):
        for i, ana in enumerate(self.get_analyses()):
            ana.Peaks.find_additional(show=False)
            self.save_canvas(get_last_canvas(), '', 'r{}'.format(i), res_dir='', ftype='pdf')

    def draw_peaks_vs_rate(self, normalise=False, log_=True, show=True):
        fluxes = self.get_fluxes()
        peak_heights = array([ana.Peaks.find_additional(scale=True, show=False) for ana in self.get_analyses()]) / (fluxes if normalise else 1)
        peak_heights /= mean(peak_heights) if normalise else 1
        g = self.make_tgrapherrors('gpr', 'Number of Additional Peaks vs Flux', x=fluxes, y=peak_heights)
        if not normalise:
            self.format_statbox(fit=True, x=.52)
            g.Fit('pol1', 'qs')
        format_histo(g, y_tit='Number of Additional Peaks {}'.format('/ Flux' if normalise else ''), y_off=1.5, **self.get_x_args(False))
        self.draw_histo(g, show=show, lm=.12, logy=log_ and not normalise, logx=log_)

    def draw_ph_vs_peaks(self, show=True):
        peak_heights = array([ana.Peaks.find_additional(scale=True, show=False) for ana in self.get_analyses()]) / self.get_fluxes()
        peak_heights /= mean(peak_heights)
        pulse_heights = self.get_pulse_heights()
        g = self.make_tgrapherrors('gpr', 'Pulse Height vs Normalised Peak Height', y=pulse_heights, x=peak_heights)
        format_histo(g, y_tit='Pulse Height', x_tit='Normalised Peak Height')
        self.draw_histo(g, show=show, lm=.12)

    def compare_fluxes(self, normalise=False, fit=True, log_=True, corr=False, avrg=False, redo=False, y_range=None, show=True):
        f0 = self.get_fluxes(corr=corr, avrg=avrg, rel_error=-.07)
        f1 = self.get_values('Peak fluxes', self.Analysis.get_peak_flux, avrg=avrg, prnt=False, redo=redo) / (f0 if normalise else 1)
        f1 /= mean(f1) if normalise else 1.
        g = self.make_tgrapherrors('gff', 'FAST-OR Flux vs Peak Flux', x=f0, y=f1)
        x_range = self.Bins.FluxRange
        y_range = self.Bins.FluxRange if y_range is None else y_range
        format_histo(g, x_tit='FAST-OR Flux [kHz/cm^{2}]', y_tit='Peak Flux {}'.format('/ FAST-OR Flux' if normalise else '[kHz/cm^{2}]'), y_off=1.3, x_off=1.2, y_range=y_range, x_range=x_range)
        if fit:
            self.format_statbox(only_fit=True, w=.2, x=.5)
            g.Fit('pol1', 'qs')
        self.draw_histo(g, show=show, lm=.12, logx=log_, logy=log_ and not normalise)

    def draw_bunch_systematics(self, bunch=0, show=True):
        all_bunches = self.get_values('All bunches', self.Analysis.get_n_peaks) / self.FirstAnalysis.Peaks.NBunches
        single_bunch = self.get_values('Single bunch', self.Analysis.get_n_peaks, start_bunch=bunch, end_bunch=bunch + 1)
        g = self.make_tgrapherrors('gps', 'Systematics of the Number of Peaks of Bunch {}'.format(bunch), x=self.get_fluxes(), y=single_bunch / all_bunches)
        format_histo(g, y_tit='N Peaks in Bunch {} / Average N Peaks per Bunch'.format(bunch), y_off=1.3, **self.get_x_args(False))
        self.draw_histo(g, lm=.12, show=show, logx=True)
    # endregion PEAKS
    # ----------------------------------------

    # ----------------------------------------
    # region 2D SIGNAL MAP
    def draw_flux_comparison(self):
        g1 = self.make_tgrapherrors('g_fp', 'Number of Peaks', color=self.get_color())
        pixel_fluxes = [ana.Run.Flux / 1000. for ana in self.Analyses.values()]
        g2 = self.make_tgrapherrors('g_ff', 'Pixel Fast-OR', x=self.Runs, y=pixel_fluxes, color=self.get_color())
        for i, (run, ana) in enumerate(self.Analyses.iteritems()):
            flux, err = ana.Peaks.get_flux()
            g1.SetPoint(i, run, flux / 1000.)
            g1.SetPointError(i, 0, err / 1000.)
        mg = TMultiGraph('mg_ff', 'Flux Comparison')
        mg.Add(g1, 'pl')
        mg.Add(g2, 'pl')
        leg = self.make_legend(nentries=3, x2=.4)
        leg.AddEntry(g1, g1.GetTitle(), 'pl')
        leg.AddEntry(g2, g2.GetTitle(), 'pl')
        names = ['1x1', '2x1', '2x2', '4x2', '4x4']
        format_histo(mg, x_tit='Pattern', y_tit='Flux [MHz/cm^{2}]', y_off=1.4, draw_first=True)
        for i, run in enumerate(self.Runs):
            bin_x = mg.GetXaxis().FindBin(run)
            mg.GetXaxis().SetBinLabel(bin_x, names[i])
        self.save_histo(mg, 'FluxComparison', draw_opt='a', lm=.12, bm=.2, leg=leg)

    def draw_timing(self):
        self.generate_plots('timing', PadAnalysis.draw_timing)

    def draw_pulser_rates(self):
        self.generate_plots('pulser rates', PadAnalysis.draw_pulser_rate, {'show': False, 'prnt': False})


if __name__ == '__main__':

    p = init_argparser(run=12, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)
    p.add_argument('-r', '--runs', action='store_true')
    p.add_argument('-d', '--draw', action='store_true')
    p.add_argument('-rd', '--redo', action='store_true')
    pargs = p.parse_args()

    z = PadCollection(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
    z.print_loaded()
    if pargs.runs:
        z.Currents.draw_indep_graphs()
        raw_input('Press any button to exit')
    if pargs.draw:
        z.draw_all(pargs.redo)
