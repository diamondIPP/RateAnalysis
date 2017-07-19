# ==============================================
# IMPORTS
# ==============================================
from argparse import ArgumentParser
from collections import OrderedDict
from numpy import log, array, zeros, std, mean
from dispy import JobCluster
from functools import partial

from ROOT import gROOT, TCanvas, TLegend, TExec, gStyle, TMultiGraph, THStack, TF1, TH1F, TH2F

from CurrentInfo import Currents
from Elementary import Elementary
from Extrema import Extrema2D
from PadAnalysis import PadAnalysis
from RunSelection import RunSelection
from TelescopeAnalysis import Analysis
from Run import Run
from Utils import *


# ==============================================
# MAIN CLASS
# ==============================================
class AnalysisCollection(Elementary):
    """
    An object of this class contains several analysis of runs.
    It gives the ability to compare the data from different runs.
    """
    current_run_number = -1

    def __init__(self, run_selection, load_tree=True, verbose=False):
        Elementary.__init__(self, verbose=verbose)
        # dict where all analysis objects are saved
        self.collection = OrderedDict()
        self.selection = run_selection

        self.runs = run_selection.get_selected_runs()
        self.channel = self.load_channel()
        # self.diamonds = self.load_diamonds(diamonds, list_of_runs)
        self.min_max_rate_runs = self.get_high_low_rate_runs()

        if load_tree:
            self.generate_slope_pickle()
            self.generate_threshold_pickle()

        self.add_analyses(load_tree)
        self.FirstAnalysis = self.get_first_analysis()

        self.signalValues = None

        # important information
        self.NRuns = len(self.collection)
        self.DiamondName = self.FirstAnalysis.DiamondName
        self.DiamondNumber = self.FirstAnalysis.DiamondNumber
        self.Bias = self.FirstAnalysis.Bias

        # root stuff
        self.RunPlan = run_selection.SelectedRunplan
        self.Type = run_selection.SelectedType
        self.save_dir = '{dia}/runplan{plan}'.format(tc=self.TESTCAMPAIGN[2:], plan=self.RunPlan, dia=self.DiamondName)
        self.RootObjects = []
        # important plots
        self.FWHM = None
        self.PulseHeight = None
        self.Pedestal = None
        self.PeakDistribution = None

        # current information
        self.StartTime = float(self.FirstAnalysis.run.log_start.strftime('%s'))
        self.Currents = Currents(self)

    def __del__(self):
        print '\ndeleting AnalysisCollection...'
        for nr, ana in self.collection.iteritems():
            ana.__del__()

    def close_files(self):
        for ana in self.collection.itervalues():
            ana.run.tree.Delete()
            ana.run.RootFile.Close()

    # ============================================
    # region INIT
    def load_run_config(self):
        return self.load_run_configs(0)

    def create_analysis(self, run, load_tree):
        self.collection[run] = PadAnalysis(run, self.channel, self.min_max_rate_runs, load_tree=load_tree, verbose=self.verbose)

    def add_analyses_fast(self, load_tree):
        cluster = JobCluster(self.create_analysis)
        jobs = []
        for i, run in enumerate(self.runs):
            job = cluster.submit((run, load_tree,))
            job.id = i
            jobs.append(job)
        for job in jobs:
            job()
            print job.id, job.start_time, job.end_time

    def add_analyses(self, load_tree):
        """ Creates and adds Analysis objects with run numbers in runs. """
        for run in self.runs:
            analysis = PadAnalysis(run, self.selection.SelectedDiamondNr, self.min_max_rate_runs, load_tree=load_tree, verbose=self.verbose)
            self.collection[analysis.run.RunNumber] = analysis
            self.current_run_number = analysis.run.RunNumber

    def load_channel(self):
        binary = self.run_config_parser.getint('ROOTFILE_GENERATION', 'active_regions')
        dia_nr = self.selection.SelectedDiamondNr
        return [i for i in xrange(16) if self.has_bit(binary, i)][dia_nr - 1]

    def get_high_low_rate_runs(self):
        runs = [Run(run_number=run, load_tree=False, test_campaign=self.TESTCAMPAIGN) for run in self.runs]
        fluxes = OrderedDict()
        if self.verbose:
            print 'RUN FLUX [kHz/cm2]'
            for run in runs:
                fluxes[run.Flux] = run.RunNumber
                print '{run:3d} {flux:14.2f}'.format(run=run.RunNumber, flux=run.Flux)
        print '\n'
        return {'min': fluxes[min(fluxes)], 'max': fluxes[max(fluxes)]}

    def generate_slope_pickle(self):
        picklepath = 'Configuration/Individual_Configs/Slope/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.min_max_rate_runs['min'])
        if file_exists(picklepath):
            return
        Analysis(self.min_max_rate_runs['min'])

    def generate_threshold_pickle(self):
        picklepath = 'Configuration/Individual_Configs/Cuts/SignalThreshold_{tc}_{run}_{ch}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.min_max_rate_runs['max'], ch=0)
        if file_exists(picklepath):
            return
        Analysis(self.min_max_rate_runs['max'])

    # endregion

    def create_all_single_run_plots(self):
        old_verbose = self.FirstAnalysis.verbose
        self.set_verbose(False)
        for key, ana in self.collection.iteritems():
            print 'Create Plots for Run ', key
            ana.compare_consecutive_cuts(scale=False)
            ana.compare_consecutive_cuts(scale=True)
            ana.show_cut_contributions(show=False)
            ana.draw_bucket_pedestal(show=False)
        self.set_verbose(old_verbose)

    # ============================================
    # region SIGNAL/PEDESTAL
    def print_all_off_results(self):
        string = 'Signal\tPedest.\tPulser\n'
        for ana in self.collection.itervalues():
            string += '{0}\n'.format(ana.print_off_results(prnt=False))
        print string

    def draw_all(self, save=True):
        old_verbose = self.FirstAnalysis.verbose
        self.set_verbose(False)
        self.draw_ph_with_currents(show=False)
        self.draw_pulse_heights(binning=10000, show=False)
        self.draw_pulser_info(do_fit=False, show=False)
        self.draw_pedestals(show=False, save=save)
        self.draw_noise(show=False)
        self.draw_pulser_pedestals(show=False, save=save)
        self.draw_signal_distributions(show=False)
        self.set_verbose(old_verbose)
        self.print_all_off_results()
    
    def scale_current_gr(self, gr):
        vals = [gr.GetY()[i] for i in xrange(gr.GetN())]
        mi, ma, mn = min(vals), max(vals), mean(vals)
        y = [mi, mn + 3 * (mn - mi) if ma - mn > 50 else ma]
        self.format_histo(gr, name='cur', color=899)
        print y
        return increased_range(y, .3, .2)

    def draw_ph_with_currents(self, show=True):
        ph = self.draw_pulse_heights(show=False, vs_time=True, fl=False, save_comb=False, binning=10000)
        self.Currents.set_graphs()
        cur = self.Currents.CurrentGraph.Clone()
        cur_range = self.scale_current_gr(cur)
        scale = ph.GetListOfGraphs()[0].GetY()[0]
        pul_orginal = self.draw_pulser_info(show=False, do_fit=False, vs_time=True, save_comb=False)
        fac = scale / pul_orginal.GetY()[0]
        pul = self.draw_pulser_info(show=False, do_fit=False, vs_time=True, scale=scale, save_comb=False)
        pul.SetLineColor(859)
        pul.SetMarkerColor(859)

        entries = [2, 3, 1]
        positions = [[.8, .25], [.8, .3], [.8, .2]]
        legends = [self.make_legend(*positions[i], nentries=entries[i], scale=1.7, x2=0.94) for i in xrange(3)]
        legends[1].AddEntry(ph.GetListOfGraphs()[2], '{sgn}(signal - ped.) data'.format(sgn='-1 * ' if self.Bias < 0 else ''), 'p')
        dummy_gr = self.make_tgrapherrors('g', 'g', width=2)
        legends[1].AddEntry(0, 'flux in kHz/cm^{2}', 'p')
        legends[1].AddEntry(dummy_gr, 'duration', 'l')
        legends[0].AddEntry(pul, 'scaled pulser data', 'p')
        legends[0].AddEntry(0, 'factor {sgn}{fac:4.2f}'.format(fac=fac, sgn='-' if self.FirstAnalysis.PulserPolarity < 0 else ''), '')
        legends[2].AddEntry(cur, 'current', 'l')

        gROOT.SetBatch(1) if not show else do_nothing()
        c = TCanvas('c', 'c', 1500, 1000)
        margins = [[.075, .05, 0, .1], [.075, .05, 0, 0], [.075, .05, 0.05, 0]]
        pads = [self.Currents.draw_tpad('p{0}'.format(i + 1), 'p{0}'.format(i + 1), pos=[0, (3 * i + 1) / 10., 1, ((i + 1) * 3 + 1) / 10.], gridx=True, margins=margins[2 - i]) for i in xrange(3)]
        for pad in pads:
            pad.Draw()
        draw_opts = ['pl', '', 'l']
        y_tits = ['Pulser Pulse Height [au] ', 'Signal Pulse Height [au] ', 'Current [nA] ']
        y_ranges = [increased_range(scale_margins(pul, ph), .3), increased_range(scale_margins(ph, pul), .3), cur_range]
        divs = [204, 204, None]

        pad = None
        for i, gr in enumerate([pul, ph, cur], 1):
            pad = pads[3 - i]
            pad.cd()
            pad.SetMargin(*margins[i - 1])
            if i != 3:
                pad.SetGridy()
            ymin, ymax = y_ranges[i - 1]
            self.Currents.draw_frame(pad, ymin, ymax, y_tits[i - 1], div=divs[i - 1])
            if i == 2:
                self.Currents.draw_time_axis(ymin, opt='t')
            if i == 3:
                self.Currents.draw_time_axis(ymax, opt='t')
            gr.Draw(draw_opts[i - 1])
            legends[i - 1].Draw()

        run_info = self.FirstAnalysis.run.get_runinfo(self.FirstAnalysis.DiamondNumber, pad=pad)
        width = len(run_info[0].GetListOfPrimitives()[1].GetLabel()) * 0.0064
        self.FirstAnalysis.run.scale_runinfo_legend(w=width)
        c.cd()
        run_info[0].Draw()
        run_info[1].Draw()
        c.Update()
        gROOT.SetBatch(0)
        self.save_canvas(c, self.save_dir, 'PhPulserCurrent', show=show)
        self.RootObjects.append([ph, cur, pul, c, legends, pads])
        self.FirstAnalysis.run.reset_info_legend()

    def draw_ph_vs_voltage(self, binning=10000, pulser=False, redo=False):
        gr1 = self.make_tgrapherrors('gStatError', 'stat. error', self.get_color())
        gStyle.SetEndErrorSize(4)
        gr_first = self.make_tgrapherrors('gFirst', 'first run', marker=22, color=2, marker_size=2)
        gr_last = self.make_tgrapherrors('gLast', 'last run', marker=23, color=2, marker_size=2)
        gr_errors = self.make_tgrapherrors('gFullError', 'stat. + repr. error', marker=0, color=602, marker_size=0)

        # flux_errors = self.draw_ph_distributions_below_flux(flux=80, show=False, save_plot=False)
        # rel_sys_error = flux_errors[1] / flux_errors[0]
        rel_sys_error = 0
        i, j = 0, 0
        for key, ana in self.collection.iteritems():
            fit1 = ana.draw_pulse_height(binning, evnt_corr=True, save=redo, show=False) if not pulser else ana.Pulser.draw_distribution_fit(show=False, save=False)
            x = ana.run.RunInfo['dia{nr}hv'.format(nr=self.DiamondNumber)]
            print x, '\t',
            s, e = (fit1.Parameter(0), fit1.ParError(0)) if not pulser else (fit1.Parameter(1), fit1.ParError(1))
            gr1.SetPoint(i, x, s)
            print s, e
            gr1.SetPointError(i, 0, e)
            gr_errors.SetPoint(i, x, s)
            gr_errors.SetPointError(i, 0, e + rel_sys_error * s)
            # set special markers for the first and last run
            if i == 0:
                gr_first.SetPoint(0, x, s)
            if j == len(self.collection) - 1:
                gr_last.SetPoint(0, x, s)
            i += 1
            j += 1
        graphs = [gr_errors, gr1]
        gr_line = gr1.Clone()
        self.format_histo(gr_line, name='gLine', color=920)
        graphs += [gr_first, gr_last]
        legend = self.make_legend(.65, .35, nentries=len(graphs))
        # gr1.SetName('data') if len(graphs) < 5 else self.do_nothing()

        mg = TMultiGraph('mg_ph', '' + self.DiamondName)
        mg.Add(gr_line, 'l')
        for gr in graphs:
            if gr.GetName().startswith('gFull'):
                legend.AddEntry(gr, gr.GetTitle(), 'l')
            else:
                legend.AddEntry(gr, gr.GetTitle(), 'p')
            mg.Add(gr, 'p')
        self.format_histo(mg, x_tit='Voltage [V]', y_tit='Pulse Height [au]', y_off=1.3, draw_first=True)
        self.save_histo(mg, '{s}VoltageScan'.format(s='Signal' if not pulser else 'Pulser'), draw_opt='a', lm=.12)

    def draw_slope_vs_voltage(self, show=True, gr=False):
        h = TH1F('hSV', 'PH Slope Distribution', 10, -1, 1) if not gr else self.make_tgrapherrors('gSV', 'PH Slope vs. Voltage')

        self.start_pbar(self.NRuns)
        for i, (key, ana) in enumerate(self.collection.iteritems()):
            ana.draw_ph(corr=True, show=False)
            fit = ana.PulseHeight.Fit('pol1', 'qs')
            x = ana.run.RunInfo['dia{nr}hv'.format(nr=self.FirstAnalysis.run.channels.index(self.channel) + 1)]
            s, e = fit.Parameter(1), fit.ParError(1)
            if gr:
                h.SetPoint(i, x, s)
                h.SetPointError(i, 0, e)
            else:
                set_statbox(entries=4, only_fit=True)
                h.Fill(s)
                h.Fit('gaus', 'q')
            self.ProgressBar.update(i + 1)
        self.format_histo(h, x_tit='Slope [mV/min]', y_tit='Number of Entries', y_off=1.3)
        self.draw_histo(h, show=show, draw_opt='alp' if gr else '')

    def draw_pulse_heights(self, binning=10000, flux=True, raw=False, all_corr=False, show=True, save_plots=True, vs_time=False, fl=True, save_comb=True, y_range=None, redo=False):

        pickle_path = self.make_pickle_path('Ph_fit', 'PulseHeights', self.RunPlan, ch=self.DiamondName, suf=binning)
        flux = False if vs_time else flux

        def func(y_ran):

            mode = self.get_mode(flux, vs_time)
            prefix = 'Pulse Height vs {mod} - '.format(mod=mode)
            marker_size = 2

            gr1 = self.make_tgrapherrors('gStatError', 'stat. error', self.get_color(), marker_size=marker_size)
            gr2 = self.make_tgrapherrors('gBinWise', prefix + 'binwise correction', self.get_color(), marker_size=marker_size)
            gr3 = self.make_tgrapherrors('gMeanPed', prefix + 'mean correction', self.get_color(), marker_size=marker_size)
            gr4 = self.make_tgrapherrors('gRaw', prefix + 'raw', self.get_color(), marker_size=marker_size)
            gr5 = self.make_tgrapherrors('gFlux', 'bla', 1, width=1, marker_size=0)
            gStyle.SetEndErrorSize(4)
            gr_first = self.make_tgrapherrors('gFirst', 'first run', marker=22, color=2, marker_size=marker_size * 2)
            gr_last = self.make_tgrapherrors('gLast', 'last run', marker=23, color=2, marker_size=marker_size * 2)
            gr_errors = self.make_tgrapherrors('gFullError', 'stat. + repr. error', marker=0, color=602, marker_size=0)
            not_found_for = False in [coll.run.FoundForRate for coll in self.collection.itervalues()]

            flux_errors = self.get_repr_errors(80, False)
            log_message('Getting pulse heights{0}'.format(' vs time' if vs_time else ''))
            rel_sys_error = flux_errors[1] / flux_errors[0]
            i, j = 0, 0
            self.start_pbar(self.NRuns)
            for key, ana in self.collection.iteritems():
                fit1 = ana.draw_pulse_height(binning, evnt_corr=True, save=False)
                if all_corr:
                    fit2 = ana.draw_pulse_height(binning, bin_corr=True, save=False)
                    fit3 = ana.draw_pulse_height(binning, off_corr=True, save=False, evnt_corr=False)
                    fit4 = ana.draw_pulse_height(binning, evnt_corr=False, save=False)
                x = key
                if flux:
                    x = ana.run.RunInfo['measured flux'] if not_found_for else ana.run.Flux
                if vs_time:
                    self.set_root_output(False)
                    x_err = ana.run.duration.seconds / 2.
                    x = int(ana.run.log_start.strftime('%s')) + x_err - self.StartTime
                    gr5.SetPoint(i, x, fit1.Parameter(0))
                    gr5.SetPointError(i, x_err, 0)
                    l1 = self.draw_tlatex(gr5.GetX()[i] - x_err, gr5.GetY()[i] + .03, '{0:5.0f}'.format(ana.run.Flux), color=1, align=10, size=.04)
                    gr1.GetListOfFunctions().Add(l1)
                if fit1.Parameter(0) > 10:
                    gr1.SetPoint(i, x, fit1.Parameter(0))
                    gr1.SetPointError(i, .1 * x if flux else 0, fit1.ParError(0))
                    gr_errors.SetPoint(i, x, fit1.Parameter(0))
                    gr_errors.SetPointError(i, .1 * x if flux else 0, fit1.ParError(0) + rel_sys_error * fit1.Parameter(0))
                    if all_corr:
                        gr2.SetPoint(i, x, fit2.Parameter(0))
                        gr3.SetPoint(i, x, fit3.Parameter(0))
                        gr4.SetPoint(i, x, fit4.Parameter(0))
                        gr2.SetPointError(i, 0, fit2.ParError(0))
                        gr3.SetPointError(i, 0, fit3.ParError(0))
                        gr4.SetPointError(i, 0, fit4.ParError(0))
                    # set special markers for the first and last run
                    if i == 0:
                        gr_first.SetPoint(0, x, fit1.Parameter(0))
                    if j == len(self.collection) - 1:
                        gr_last.SetPoint(0, x, fit1.Parameter(0))
                    i += 1
                self.ProgressBar.update(j + 1)
                j += 1
            self.ProgressBar.finish()
            graphs = [gr_errors, gr1]
            gr_line = gr1.Clone()
            self.format_histo(gr_line, name='gLine', color=920)
            if fl:
                graphs += [gr_first, gr_last]
            if all_corr:
                graphs += [gr2, gr3]
            if raw:
                graphs.append(gr4)
            legend = self.make_legend(.65, .35, nentries=len(graphs))
            # gr1.SetName('data') if len(graphs) < 5 else self.do_nothing()

            mg = TMultiGraph('mg_ph', prefix + self.DiamondName)
            mg.Add(gr_line, 'l') if not self.Type == 'random scan' else do_nothing()
            for gr in graphs:
                if gr.GetName().startswith('gFull'):
                    legend.AddEntry(gr, gr.GetTitle(), 'l')
                else:
                    legend.AddEntry(gr, gr.GetTitle(), 'p')
                mg.Add(gr, 'p')

            # small range
            self.format_histo(mg, color=None, x_tit=mode + ' [kHz/cm^{2}]' if flux else '', y_tit='Signal Pulse Height [au]', y_off=1.75, x_off=1.3, draw_first=True)
            ymin, ymax = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
            yrange = increased_range([ymin, ymax], .3, .15) if y_ran is None else y_ran
            mg.GetYaxis().SetRangeUser(*yrange)
            if vs_time:
                mg.Add(gr5, '[]')
                mg.Add(gr5, 'p')
                mg.GetXaxis().SetTimeDisplay(1)
                mg.GetXaxis().SetTimeFormat('%H:%M%F2000-02-28 23:00:00')
                mg.GetXaxis().SetLabelSize(.03)
            x_vals = sorted([gr1.GetX()[i] for i in xrange(gr1.GetN())])
            mg.GetXaxis().SetLimits(x_vals[0] * 0.8, x_vals[-1] * 1.2) if flux else self.do_nothing()
            self.save_histo(mg, 'PulseHeight{mod}'.format(mod=mode.title()), False, self.save_dir, lm=.14, draw_opt='A', l=legend, logx=True if flux else 0, gridy=1 if vs_time else 0,
                            gridx=True if vs_time else 0)

            # no zero suppression
            mg1 = mg.Clone()
            mg1.SetName('mg1_ph')
            mg1.GetListOfGraphs()[0].SetLineColor(self.colors[0])
            mg1.GetYaxis().SetRangeUser(0, ymax * 1.1)
            self.save_histo(mg1, 'PulseHeightZero{mod}'.format(mod=mode.title()), False, self.save_dir, lm=.14, draw_opt='A', l=legend, logx=True if flux else 0)

            self.reset_colors()

            self.PulseHeight = gr1
            if save_comb:
                run_info = self.FirstAnalysis.run.get_runinfo(self)
                y_min = increased_range([ymin, ymax], .3)[0] if y_ran is None else y_ran[0]
                self.save_combined_pulse_heights(mg, mg1, legend, y_min, show=show, run_info=run_info, pulser_leg=self.__draw_signal_legend)
            return mg

        f = partial(func, y_range)
        mg2 = func(y_range) if save_plots or redo else None
        return self.do_pickle(pickle_path, f, mg2)

    def draw_pedestals(self, region='ab', peak_int='2', flux=True, all_regions=False, sigma=False, show=True, cut=None, beam_on=True, save=False):

        pickle_path = self.make_pickle_path('Pedestal', 'AllPedestals', self.RunPlan, self.DiamondName, suf='Sigma' if sigma else 'Mean')
        mode = 'Flux' if flux else 'Run'
        log_message('Getting pedestals')
        self.start_pbar(self.NRuns)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        legend.SetName('l1')
        cut_string = self.get_first_analysis().Cut.generate_pulser_cut(beam_on=beam_on) if cut == 'pulser' else cut

        def func():
            y_val = 'Sigma' if sigma else 'Mean'
            prefix = 'Pulser ' if cut is not None and cut.GetName().startswith('Pulser') else ''
            gr1 = self.make_tgrapherrors('pedestal', '{pre}Pedestal {y} in {reg}'.format(y=y_val, reg=region + peak_int, pre=prefix))
            regions = self.get_first_analysis().run.pedestal_regions
            graphs = [self.make_tgrapherrors('pedestal', '{pre}Pedestal {y} in {reg}'.format(y=y_val, reg=reg + peak_int, pre=prefix), color=self.get_color()) for reg in regions]
            i = 0
            par = 2 if sigma else 1
            for key, ana in self.collection.iteritems():
                fit_par = ana.show_pedestal_histo(region, peak_int, cut=cut_string, save=save, show=False)
                x = ana.run.Flux if flux else key
                gr1.SetPoint(i, x, fit_par.Parameter(par))
                gr1.SetPointError(i, 0, fit_par.ParError(par))
                if all_regions:
                    for reg, gr in zip(regions, graphs):
                        fit_par = ana.show_pedestal_histo(reg, peak_int, save=False)
                        gr.SetPoint(i, x, fit_par.Parameter(par))
                        gr.SetPointError(i, 0, fit_par.ParError(par))
                self.ProgressBar.update(i + 1)
                i += 1
            self.format_histo(gr1, color=None, x_tit=self.make_x_tit(mode, flux), y_tit='Mean Pedestal [au]', y_off=1.45)
            if all_regions:
                for i, gr in enumerate(graphs):
                    legend.AddEntry(gr, str(regions.values()[i]), 'p')
                    gr.Draw('alp') if not i else gr.Draw('lp')
            self.Pedestal = gr1
            self.reset_colors()
            return gr1

        self.ProgressBar.finish()
        graph = func() if show or save else None
        graph = self.do_pickle(pickle_path, func, graph)
        save_name = 'Pedestal{s}{mod}{cut}'.format(mod=mode, cut='' if cut is None else cut_string.GetName(), s='Sigma' if sigma else 'Mean')
        print save_name
        self.save_histo(graph, save_name=save_name, show=show, logx=True if flux else False, l=legend if all_regions else None, lm=.12)
        return

    def draw_noise(self, flux=True, show=True, save=False):
        return self.draw_pedestals(flux=flux, show=show, save=save, sigma=True)

    def draw_pulser_pedestals(self, show=True, save=False):
        self.draw_pedestals(cut=self.FirstAnalysis.Pulser.PulserCut, show=show, save=save)

    def draw_signal_distributions(self, show=True, off=3, redo=False):

        stack = THStack('hsd', 'Pulse Height Distributions')
        legend = self.make_legend(.67, .88, nentries=self.get_number_of_analyses())
        log_message('Generating signal distributions!')
        histos = []
        self.start_pbar(self.NRuns)
        for i, ana in enumerate(self.collection.itervalues(), 1):
            histos.append(ana.show_signal_histo(show=False, redo=redo))
            self.ProgressBar.update(i)
        for i, h in enumerate(histos):
            self.format_histo(h, lw=2, color=self.get_color())
            h.Scale(1 / h.GetMaximum())
            stack.Add(h)
            legend.AddEntry(h, '{0:06.1f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux()), 'l')
        self.format_histo(stack, y_off=1.55, draw_first=True, x_tit='Pulse Height [au]', y_tit='Number of Entries')
        log_stack = stack.Clone()
        log_stack.SetMaximum(off)
        log_stack.SetNameTitle('hsdl', 'Signal Distribution LogY')
        self.save_histo(stack, 'SignalDistributions', False, self.save_dir, lm=.13, draw_opt='nostack', l=legend)
        self.save_histo(log_stack, 'SignalDistributionsLogY', False, self.save_dir, lm=.13, draw_opt='nostack', logy=True, l=legend)
        if show:
            c = TCanvas('c_sd1', 'Signal Distributions', 1500, 750)
            c.Divide(2)
            legends = [legend, legend.Clone()]
            for i, s in enumerate([stack, log_stack], 1):
                pad = c.cd(i)
                pad.SetLeftMargin(.14)
                pad.SetLogy() if i == 2 else self.do_nothing()
                s.Draw('nostack')
                legends[i - 1].Draw()
            self.RootObjects.append([c, legends])
        self.reset_colors()

    def draw_snrs(self, flux=True, draw=True):
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        gr = self.make_tgrapherrors('gr', 'SNR vs {mode}'.format(mode=mode))
        i = 0
        for key, ana in self.collection.iteritems():
            snr = ana.calc_snr()
            x = ana.run.Flux if flux else key
            gr.SetPoint(i, x, snr[0])
            gr.SetPointError(i, 0, snr[1])
            i += 1
        if draw:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        if flux:
            c.SetLogx()
        gr.Draw('ap')
        self.save_plots('AllSNRs', canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        gROOT.SetBatch(0)

    def draw_all_ph_distributions(self, binning=5000, show=False):

        pickle_path = self.FirstAnalysis.PickleDir + 'Ph_fit/Ph_distos_fits_{tc}_{rp}_{dia}_{bin}.pickle'.format(tc=self.TESTCAMPAIGN, rp=self.RunPlan, dia=self.DiamondName, bin=binning)

        def func():
            collimator_settings = [(ana.run.RunInfo['fs11'], ana.run.RunInfo['fsh13']) for key, ana in self.collection.iteritems()]
            collimator_settings = set(collimator_settings)
            fits = {}
            for s in collimator_settings:
                retval = self.draw_ph_distributions(binning, show=show, fs11=s[0], fsh13=s[1])
                fits[s] = retval
            return fits

        res = func() if show else None
        return self.do_pickle(pickle_path, func, res)

    def get_repr_errors(self, flux, show=True):
        runs = self.get_runs_below_flux(flux)
        vals = [self.collection[run].draw_pulse_height(save=False).Parameter(0) for run in runs]
        gr = self.make_tgrapherrors('gr_re', 'Pulse Heights Below {f} kHz/cm^{{2}}'.format(f=flux))
        for i, run in enumerate(runs):
            ana = self.collection[run]
            fit = ana.draw_pulse_height(save=False)
            gr.SetPoint(i, ana.run.Flux, fit.Parameter(0))
            gr.SetPointError(i, 0, fit.ParError(0))
        set_statbox(entries=2, only_fit=True)
        gr.Fit('pol0', 'qs{s}'.format(s='' if show else '0'))
        self.format_histo(gr, x_tit='Flux [kHz/cm^{2}]', y_tit='Mean Pulse Height [au]', y_off=1.7)
        self.save_histo(gr, 'ReprErrors', show, draw_opt='ap', lm=.115, prnt=show)
        return mean(vals), std(vals)

    def draw_ph_distributions(self, binning=5000, fsh13=.5, fs11=65, show=True):
        runs = self.get_runs_by_collimator(fsh13=fsh13, fs11=fs11)
        return self.draw_combined_ph_distributions(runs, binning, show)

    def draw_ph_distributions_below_flux(self, binning=None, flux=150, show=True, save_plot=True):
        binning = int(self.FirstAnalysis.run.n_entries * .3 / 60 if binning is None else binning)
        pickle_path = self.make_pickle_path('Ph_fit', 'PhDistoBel', self.RunPlan, self.DiamondName, suf='{bin}_{flux}'.format(bin=binning, flux=flux))

        def func():
            log_message('Getting representative errors')
            runs = self.get_runs_below_flux(flux)
            return self.draw_combined_ph_distributions(runs, binning, show)

        err = func() if save_plot else None
        return self.do_pickle(pickle_path, func, err)

    def draw_combined_ph_distributions(self, runs, binning=200, show=True):
        stack = THStack('s_phd', 'Pulse Height Distributions')
        self.reset_colors()
        self.start_pbar(len(runs))
        for i, run in enumerate(runs, 1):
            ana = self.collection[run]
            self.set_root_output(False)
            h = ana.draw_ph_distribution(show=False, binning=binning, fit=False, save=False)
            self.set_root_output(True)
            self.format_histo(h, fill_color=4000)
            h.SetStats(False)
            h.SetLineColor(self.get_color())
            stack.Add(h)
            self.ProgressBar.update(i)
        self.ProgressBar.finish()
        self.draw_histo(stack, '', show, draw_opt='')
        summed_stack = stack.GetStack().Last().Clone()
        summed_stack.SetFillStyle(0)

        fit = TF1('fit', 'gaus', 0, 160)
        fit.SetNpx()
        fitptr = summed_stack.Fit(fit, 'SQ0', 'same')
        summed_stack.SetStats(1)
        set_statbox(only_fit=True)
        summed_stack.SetName('Fit Results')
        summed_stack.Draw('sames')
        stack.GetXaxis().SetRange(summed_stack.FindFirstBinAbove(0) - 2, summed_stack.FindLastBinAbove(0) + 2)
        fit.Draw('same')
        self.ROOTObjects.append(summed_stack)
        self.ROOTObjects.append(fit)
        self.save_plots('PulseHeightDistributions')
        return fitptr.Parameter(1), fitptr.Parameter(2), fitptr.Chi2() / fitptr.Ndf()

    def calc_pedestal_spread(self, fsh13=.5, fs11=65):
        runs = self.get_runs_by_collimator(fs11, fsh13)
        values = []
        for run in runs:
            ana = self.collection[run]
            fit = ana.show_pedestal_histo(save=False)
            values.append(fit.Parameter(1))
        return max(values) - min(values)

    def calc_all_pedestal_spreads(self, recalc=False):

        pickle_path = self.FirstAnalysis.PickleDir + 'Pedestal/PedSpread_{tc}_{rp}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, rp=self.RunPlan, dia=self.DiamondName)

        def func():
            collimator_settings = [(ana.run.RunInfo['fs11'], ana.run.RunInfo['fsh13']) for key, ana in self.collection.iteritems()]
            collimator_settings = set(collimator_settings)
            mins = {}
            for s in collimator_settings:
                retval = self.calc_pedestal_spread(s[1], s[0])
                mins[s] = retval
            return mins

        res = func() if recalc else None
        return self.do_pickle(pickle_path, func, res)

    def calc_full_pedestal_spread(self):
        values = []
        for ana in self.collection.values():
            fit = ana.show_pedestal_histo(save=False)
            values.append(fit.Parameter(1))
        return max(values) - min(values)

    # endregion

    # ============================================
    # region PULSER
    def draw_pulser_info(self, flux=True, show=True, mean_=True, corr=True, beam_on=True, vs_time=False, do_fit=True, scale=1, save_comb=True, save=True, ret_mg=False):

        pickle_path = self.make_pickle_path('Pulser', 'PulseHeights', self.RunPlan, self.DiamondName)
        flux = False if vs_time else flux
        mode = self.get_mode(flux, vs_time)
        log_message('Getting pulser info{0}'.format(' vs time' if vs_time else ''))

        def func():
            self.start_pbar(len(self.collection))
            gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
            gr = self.make_tgrapherrors('data', 'pulser data', color=602)
            gr1 = self.make_tgrapherrors('gFirst', 'first run', marker=22, color=2, marker_size=2)
            gr2 = self.make_tgrapherrors('gLast', 'last run', marker=23, color=2, marker_size=2)
            y0 = None
            for i, (key, ana) in enumerate(self.collection.iteritems()):
                x = ana.run.Flux if flux else key
                fit = ana.Pulser.draw_distribution_fit(save=False, corr=corr, beam_on=beam_on)
                par = 1 if mean_ else 2
                cut = ana.Cut.generate_pulser_cut(beam_on)
                ped_fit = ana.show_pedestal_histo(cut=cut, save=False)
                ped_err = ped_fit.ParError(par)
                if vs_time:
                    xerr = ana.run.duration.seconds / 2.
                    x = int(ana.run.log_start.strftime('%s')) + xerr - self.StartTime
                y = fit.Parameter(par)
                y0 = y if y0 is None else y0
                yerr = sqrt(pow(fit.ParError(par), 2) + pow(ped_err, 2))
                if scale != 1:
                    y *= scale / y0
                    yerr *= scale / y0
                gr.SetPoint(i, x, y)
                gr.SetPointError(i, .1 * x if flux else 0, yerr)
                if i == 0:
                    gr1.SetPoint(0, x, y)
                if i == len(self.collection) - 1:
                    gr2.SetPoint(0, x, y)
                self.ProgressBar.update(i + 1)
            self.ProgressBar.finish()
            if vs_time:
                gr.GetXaxis().SetTimeDisplay(1)
                gr.GetXaxis().SetTimeFormat('%H:%M%F2000-02-28 23:00:00')
                gr.GetXaxis().SetLabelSize(.03)
            if do_fit:
                gStyle.SetOptFit(1)
                gr.Fit('pol0', 'q')
            mg = TMultiGraph('mg_pph', 'Pulser Signal vs {mod}'.format(mod=mode))
            gr_line = gr.Clone()
            self.format_histo(gr_line, name='gLine', color=920)
            graphs = [gr]
            if not vs_time:
                graphs += [gr1, gr2]
            l = self.make_legend(.17, .35, nentries=3, x2=.4)
            mg.Add(gr_line, 'l')
            for graph in graphs:
                l.AddEntry(graph, graph.GetTitle(), 'p')
                mg.Add(graph, 'p')

            gROOT.SetBatch(1)
            self.format_histo(mg, x_tit=self.make_x_tit(mode, flux), y_tit='{mean} [au]'.format(mean='Pulser Pulse Height' if mean_ else 'Sigma'), draw_first=True)
            y = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
            mg_y = y[0] * 1.3 - y[1] * .3
            self.format_histo(mg, y_range=[mg_y, y[1] + (y[1] - y[0]) * .3], y_off=1.75, x_off=1.3)
            x_vals = sorted([gr.GetX()[i] for i in xrange(gr.GetN())])
            mg.GetXaxis().SetLimits(x_vals[0] * 0.8, x_vals[-1] * 1.2) if flux else self.do_nothing()
            self.save_histo(mg, 'Pulser{mean}{a}{b}'.format(mean='Mean' if mean_ else 'Sigma', a=corr, b=beam_on), lm=.14, draw_opt='A', logx=True if flux else 0, l=l, show=False)
            mg1 = mg.Clone()
            mg1.GetListOfGraphs()[0].SetLineColor(602)
            self.__draw_pulser_legend()
            if save_comb:
                self.save_combined_pulse_heights(mg, mg1, l, mg_y, show, name='CombinedPulserPulseHeights', pulser_leg=self.__draw_pulser_legend)
                self.ROOTObjects.append(mg1)
            return mg

        gra = func() if save else None
        gra = self.do_pickle(pickle_path, func, gra)
        return gra if ret_mg else gra.GetListOfGraphs()[1]

    def __draw_pulser_legend(self):
        try:
            typ = self.FirstAnalysis.RunInfo['pulser']
            pol = 'positive' if self.FirstAnalysis.PulserPolarity > 0 else 'negative'
            sig = 'positive' if self.FirstAnalysis.Polarity > 0 else 'negative'
            l1 = self.make_legend(.17, .88, nentries=3, margin=.05, felix=True, x2=.5)
            l1.AddEntry(0, 'Pulser Type:', '')
            l1.AddEntry(0, typ, '').SetTextAlign(12)
            l1.AddEntry(0, 'Pulser Polarity:', '')
            l1.AddEntry(0, pol, '').SetTextAlign(12)
            l1.AddEntry(0, 'Signal Polarity:', '')
            l1.AddEntry(0, sig, '').SetTextAlign(12)
            l1.AddEntry(0, 'Pulser Ped. Substr.:', '')
            l1.AddEntry(0, 'yes', '').SetTextAlign(12)
            l1.SetNColumns(2)
            l1.Draw()
            self.ROOTObjects.append(l1)
        except KeyError:
            pass

    def __draw_signal_legend(self):
        sig = 'positive' if self.FirstAnalysis.Polarity > 0 else 'negative'
        l1 = self.make_legend(.17, .88, nentries=2, margin=.05, felix=True, x2=.5)
        l1.AddEntry(0, 'Signal Polarity:', '')
        l1.AddEntry(0, sig, '').SetTextAlign(12)
        l1.AddEntry(0, 'Pedestal Substraction:', '')
        l1.AddEntry(0, 'yes', '').SetTextAlign(12)
        l1.SetNColumns(2)
        l1.Draw()
        self.ROOTObjects.append(l1)

    def draw_pulser_histos(self, show=True, corr=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = {i: ana.show_pulser_histo(show=False, corr=corr) for i, ana in enumerate(self.collection.itervalues()) if ana.IsAligned}
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Histos', 1000, 1000)
        legend = TLegend(.13, .88 - self.get_number_of_analyses() * 0.03, .33, .88)
        histos[0].SetTitle('Pulser Distributions {0}Corrected'.format('Pedestal' if corr else 'Un'))
        for i, h in histos.iteritems():
            h.SetStats(0)
            h.GetXaxis().SetRangeUser(h.GetBinCenter(h.FindFirstBinAbove(2) * 10 / 10 - 20), h.GetBinCenter(h.FindLastBinAbove(2) * 10 / 10 + 10))
            h.Scale(1 / h.GetMaximum())
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux()) + '^{2}', 'l')
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.save_plots('AllPulserHistos{0}'.format('Uncorrected' if not corr else ''), sub_dir=self.save_dir)
        self.RootObjects.append([c, legend] + histos.values())
        z.reset_colors()

    def draw_all_pulser_info(self, mean_=True):
        graphs = [self.draw_pulser_info(show=False, mean_=mean_, corr=x, beam_on=y) for x, y in zip([1, 1, 0, 0], [1, 0, 1, 0])]
        margins = self.find_graph_margins(graphs)
        c = TCanvas('c', 'Pulser Info', 1500, 1500)
        c.Divide(2, 2)
        for i, gr in enumerate(graphs, 1):
            gr.GetYaxis().SetRangeUser(*margins)
            self.format_histo(gr, y_off=1.3)
            pad = c.cd(i)
            pad.SetLogx()
            pad.SetBottomMargin(.15)
            gr.Draw('alp')
        gROOT.SetBatch(0)
        self.RootObjects.append([graphs, c])
        self.save_plots('AllPulserOverview{0}'.format('Mean' if mean_ else 'Sigma'), sub_dir=self.save_dir)

    def compare_pedestals(self, show=True):
        gr1 = self.draw_pedestals(show=False)
        gr2 = self.draw_pedestals(cut='pulser', show=False)
        gr3 = self.draw_pedestals(cut='pulser', show=False, beam_on=False)
        graphs = [gr1, gr2, gr3]
        margins = self.find_graph_margins(graphs)
        gROOT.SetBatch(0) if show else gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulser Pedestal Comparison', 1000, 1000)
        c.SetLogx()
        legend = TLegend(.7, .78, .88, .88)
        names = ['Signal', 'Pulser', 'BeamOff']
        for i, gr in enumerate(graphs):
            gr.GetYaxis().SetRangeUser(*margins)
            self.format_histo(gr, color=self.get_color())
            gr.Draw('lp') if i else gr.Draw('alp')
            legend.AddEntry(gr, names[i], 'pl')
        legend.Draw()
        gROOT.SetBatch(0)
        self.save_plots('PulserPedestalComparison', sub_dir=self.save_dir)
        self.RootObjects.append([c, graphs, legend])

    def draw_pulser_rates(self, show=True, flux=True, real=False):
        mode = self.get_mode(flux)
        gr = self.make_tgrapherrors('gr', 'Pulser Rate vs {mod} '.format(mod=mode))
        for i, (key, ana) in enumerate(self.collection.iteritems()):
            x = ana.run.Flux if flux else key
            fit = ana.Pulser.calc_fraction() if not real else ana.Pulser.calc_real_fraction(), 0
            gr.SetPoint(i, x, fit[0])
            gr.SetPointError(i, 0, fit[1])
        self.format_histo(gr, x_tit=self.make_x_tit(mode, flux), y_tit='Pulser Rate [Hz]')
        self.save_histo(gr, 'PulserRate{0}'.format(mode), show, logx=flux, draw_opt='alp')
        return gr

    def calc_pulser_error(self, fs11=65, fsh13=.5):
        runs = self.get_runs_by_collimator(fs11=fs11, fsh13=fsh13)
        means = []
        errors = []
        for run in runs:
            ana = self.collection[run]
            fit = ana.calc_pulser_fit(show=False)
            means.append(fit.Parameter(1))
            errors.append(fit.ParError(1))
        means_ = calc_mean(means)
        w_means = calc_weighted_mean(means, errors)
        return means_ if means_[1] > w_means[1] else w_means

    def calc_all_pulser_errors(self, recalc=False):
        pickle_path = self.FirstAnalysis.PickleDir + 'Ph_fit/PulserErrros_{tc}_{rp}_{dia}.pickle'.format(tc=self.TESTCAMPAIGN, rp=self.RunPlan, dia=self.DiamondName)

        def func():
            collimator_settings = set([(ana.run.RunInfo['fs11'], ana.run.RunInfo['fsh13']) for key, ana in self.collection.iteritems()])
            fits = {}
            for s in collimator_settings:
                retval = self.calc_pulser_error(fs11=s[0], fsh13=s[1])
                fits[s] = retval
            return fits

        errors = func() if recalc else None
        return self.do_pickle(pickle_path, func, errors)

    # endregion

    # ============================================
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
        for key, ana in self.collection.iteritems():
            x = ana.run.Flux if flux else key
            if not mean_:
                n = ana.show_bucket_numbers(show=False)
                gr1.SetPoint(i, x, n['new'] / n['all'] * 100)
                gr2.SetPoint(i, x, n['old'] / n['all'] * 100)
            else:
                info = ana.show_bucket_means(show=False, plot_histos=False)
                gr1.SetPoint(i, x, info['new'][0])
                gr2.SetPoint(i, x, info['old'][0])
                gr3.SetPoint(i, x, info['no'][0])
                gr1.SetPointError(i, 0, info['new'][1])
                gr2.SetPointError(i, 0, info['old'][1])
                gr3.SetPointError(i, 0, info['no'][1])
            i += 1
        c = TCanvas('c', 'Bucket Numbers', 1000, 1000)
        c.SetLeftMargin(.13)
        if flux:
            c.SetLogx()
        self.format_histo(gr2, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='Events [%]' if not mean_ else 'Mean [au]', y_off=1.7, color=None)
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
        self.save_plots('{mode}_' + mode, sub_dir=self.save_dir)
        self.RootObjects.append([c, gr1, gr2, gr3, leg])

    # endregion

    # ============================================
    # region PEAK VALUES
    def draw_signal_peaks(self, flux=True, draw=True, pulser=False):
        """
        Shows the means of the signal peak distribution.
        :param flux:
        :param draw:
        :param pulser:
        """
        mode = 'Flux' if flux else 'Run'
        signal = 'Pulser' if pulser else 'Signal'
        prefix = 'Mean of {pul} Peaks: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.Bias, pul='Pulser' if pulser else 'Signal')
        gr = self.make_tgrapherrors('gr', prefix)
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.fit_peak_values(draw=False, pulser=pulser)
            x = ana.run.Flux if flux else key
            gr.SetPoint(i, x, fit.Parameter(1))
            gr.SetPointError(i, 0, fit.ParError(1))
            i += 1
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''))
        c = TCanvas('c', 'Mean of {0} Peaks'.format(signal), 1000, 1000)
        c.SetLogx()
        gr.Draw('alp')
        if not draw:
            c.Close()
        self.RootObjects.append([gr, c])

    def draw_signal_fwhm(self, flux=True, draw=True):
        """
        Shows the FWHM of the signal peak distribution.
        :param flux:
        :param draw:
        """
        mode = 'Flux' if flux else 'Run'
        prefix = 'FWHM of Signal Peaks: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.Bias)
        gr = self.make_tgrapherrors('gr1', prefix)
        i = 0
        for key, ana in self.collection.iteritems():
            fwhm = ana.calc_peak_value_fwhm()
            x = ana.run.Flux if flux else key
            gr.SetPoint(i, x, fwhm)
            i += 1
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''))
        c = TCanvas('c', 'FWHM of Signal Peaks', 1000, 1000)
        gr.Draw('alp')
        if not draw:
            c.Close()
        self.RootObjects.append([gr, c])

    # endregion

    # ============================================
    # region 2D SIGNAL MAP
    def draw_mean_fwhm(self, saveplots=True, flux=True, draw=True):
        """
        Creates the FWHM Distribution of all selected MeanSignalHistograms
        :param saveplots: if True saves the plot
        :param flux: draw vs flux if True else vs run
        :param draw:
        """
        if not draw:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        prefix = 'FWHM of Mean Signal Histogram: {dia} @ {bias}V vs {mode} '.format(mode=mode, dia=self.collection.values()[0].diamond_name, bias=self.Bias)
        gr = self.make_tgrapherrors('pedestal', prefix)
        conversion_factor = 2 * sqrt(2 * log(2))  # sigma to FWHM
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.fit_mean_signal_distribution()
            x = ana.run.Flux if flux else key
            gr.SetPoint(i, x, fit.Parameter(2) * conversion_factor)
            gr.SetPointError(i, 0, fit.ParError(2) * conversion_factor)
            i += 1
        c = TCanvas('c', 'FWHM', 1000, 1000)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='Flux [kHz/cm2]', y_tit='FWHM [au]', y_off=1.1)
        gr.Draw()
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        if saveplots:
            self.save_plots('Mean_FWHM_' + mode, canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        self.FWHM = gr

    def draw_fluxes(self):
        g1 = self.make_tgrapherrors('g_fp', 'Number of Peaks', color=self.get_color())
        pixel_fluxes = [ana.run.Flux / 1000. for ana in self.collection.values()]
        g2 = self.make_tgrapherrors('g_ff', 'Pixel Fast-OR', x=self.runs, y=pixel_fluxes, color=self.get_color())
        for i, (run, ana) in enumerate(self.collection.iteritems()):
            flux, err = ana.Peaks.get_flux()
            g1.SetPoint(i, run, flux / 1000.)
            g1.SetPointError(i, 0, err / 1000.)
        mg = TMultiGraph('mg_ff', 'Flux Comparison')
        mg.Add(g1, 'pl')
        mg.Add(g2, 'pl')
        l = self.make_legend(nentries=3, x2=.4)
        l.AddEntry(g1, g1.GetTitle(), 'pl')
        l.AddEntry(g2, g2.GetTitle(), 'pl')
        names = ['1x1', '2x1', '2x2', '4x2', '4x4']
        self.format_histo(mg, x_tit='Pattern', y_tit='Flux [MHz/cm^{2}]', y_off=1.4, draw_first=True)
        for i, run in enumerate(self.runs):
            bin_x = mg.GetXaxis().FindBin(run)
            mg.GetXaxis().SetBinLabel(bin_x, names[i])
        self.save_histo(mg, 'FluxComparison', draw_opt='a', lm=.12, bm=.2, l=l)

    def save_signal_maps(self, hitmap=False):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        graphs = []
        for i, ana in enumerate(self.collection.values()):
            h = ana.draw_signal_map(show=False, hitmap=hitmap, cut='' if hitmap else None)
            # self.format_histo(h, z_range=[0, 140], stats=0) if not hitmap else do_nothing()
            h.SetContour(50)
            graphs.append(h)
        # find min/max
        glob_min = int(min([gr.GetMinimum() for gr in graphs])) / 5 * 5
        glob_max = (int(max([gr.GetMaximum() for gr in graphs])) + 5) / 5 * 5
        print glob_min, glob_max

        c = TCanvas('sig_map', 'Signal Maps', 1000, 1000)
        c.SetTheta(55)
        c.SetPhi(20)
        for i, gr in enumerate(graphs):
            gr.GetZaxis().SetRangeUser(glob_min, glob_max)
            gr.Draw('surf2')
            self.save_plots('map{}'.format(i), canvas=c, sub_dir=self.save_dir, ind=i)
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        # return graphs

    def draw_signal_spreads(self, flux=True, draw=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        mode = 'Flux' if flux else 'Run'
        gr = self.make_tgrapherrors('gr', 'Relative Spread vs {mode}'.format(mode=mode))
        i = 0
        for key, ana in self.collection.iteritems():
            rel_spread = ana.calc_signal_spread()
            x = ana.run.Flux if flux else key
            gr.SetPoint(i, x, rel_spread[0])
            gr.SetPointError(i, 0, rel_spread[1])
            i += 1
        if draw:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='Relative Spread [%]')
        gr.Draw('ap')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('RelativeSpread', canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        gROOT.SetBatch(0)

    def show_peak_distribution(self, show=True):
        """
        Shows the positions of the peaks of the 2D map.
        :param show:
        """
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        # create an overall VotingHistogram
        ana = self.get_first_analysis()
        ana.draw_mean_signal_distribution(show=False)
        extrema = Extrema2D(ana.SignalMapHisto, ana.MeanSignalHisto)
        self.PeakDistribution = extrema.create_voting_histo()
        for run, ana in self.collection.iteritems():
            if ana.IsAligned:
                ana.find_2d_extrema(histo=self.PeakDistribution, show=False)
            else:
                print 'Run {run} is not aligned...'.format(run=run)
        c = TCanvas('c', 'Voting Histos', 1600, 800)
        c.Divide(2, 1)
        # new_pal = ar.array('i', [kYellow, kYellow, kOrange, kOrange - 3, kOrange + 7, kRed])
        ex = [TExec('ex1', 'gStyle->SetPalette(56);'), TExec('ex2', 'gStyle->SetPalette(51)')]
        if show:
            gROOT.SetBatch(0)
        for i, histo in enumerate(self.PeakDistribution.itervalues(), 1):
            c.cd(i)
            histo.Draw('col')
            ex[i - 1].Draw()
            histo.Draw('colz same')
        self.save_plots('PeakDistribution', sub_dir=self.save_dir)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.RootObjects.append([c, ex])

    def draw_signal_map(self, show=True, low=1e20, fid=True, factor=1):
        suffix = '{fid}{fac}'.format(fid='Fid' if fid else '', fac=factor)
        pickle_path = self.make_pickle_path('SignalMaps', 'SigMaps', self.RunPlan, self.channel, suffix)

        def func():
            histos = [ana.draw_signal_map(show=False, marg=False, fid=fid, factor=factor) for ana in self.collection.itervalues() if ana.run.Flux < low]
            h = histos[0]
            sig_map = TH2F('h_sms', 'Combined Signal Maps', h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax(), h.GetNbinsY(), h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax())
            for h in histos:
                sig_map.Add(h)
            sig_map.Scale(1. / len(histos))
            self.format_histo(sig_map, x_tit='track_x [cm]', y_tit='track_y [cm]', y_off=1.4, z_off=1.3, stats=0, z_tit='Pulse Height [au]')
            self.__adjust_sig_map(sig_map)
            self.save_histo(sig_map, 'CombinedSignalMaps', show, lm=.12, rm=.16, draw_opt='colz')
            return sig_map

        hist = self.do_pickle(pickle_path, func)
        if not gROOT.FindObject('h_sms'):
            gStyle.SetPalette(53)
            self.__adjust_sig_map(hist)
            self.draw_histo(hist, '',  show, lm=.12, rm=.16, draw_opt='colz')
        return hist

    def __adjust_sig_map(self, h):
        h.GetZaxis().SetRangeUser(0, 500)
        h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinCenter(h.FindFirstBinAbove(0)), h.GetXaxis().GetBinCenter(h.FindLastBinAbove(0)))
        h.GetYaxis().SetRangeUser(h.GetYaxis().GetBinCenter(h.FindFirstBinAbove(0, 2)), h.GetYaxis().GetBinCenter(h.FindLastBinAbove(0, 2)))
        h1 = TH1F('h_av', 'hav', 100, 1, h.GetBinContent(h.GetMaximumBin()))
        for i in xrange(1, h.GetNbinsX() * h.GetNbinsY() + 1):
            h1.Fill(h.GetBinContent(i))
        hmax = h1.GetMaximum()
        ph_min, ph_max = h1.GetBinCenter(h1.FindFirstBinAbove(hmax * .08)), h1.GetBinCenter(h1.FindLastBinAbove(hmax * .02))
        h.GetZaxis().SetRangeUser(ph_min, ph_max)
        self.draw_histo(h1)

    def draw_hit_map(self, show=True):
        pass

    # endregion

    # ====================================================================================
    # region BEAM PROFILE
    def draw_beam_info(self, mean_=True, flux=True, show=True, direction='x', fit_margin=.6):
        if not show:
            gROOT.SetBatch(1)
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        mode = 'Flux' if flux else 'Run'
        title = 'Mean' if mean_ else 'Sigma'
        gr = self.make_tgrapherrors('gr', '{tit} of the Beam Profile in {dir}'.format(tit=title, dir=direction.title()))
        i = 0
        for key, ana in self.collection.iteritems():
            x = ana.run.Flux if flux else key
            fit = ana.fit_beam_profile(direction, show=False, fit_margin=fit_margin)
            par = 1 if mean_ else 2
            gr.SetPoint(i, x, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
            i += 1
        c = TCanvas('c', 'Beam Profile', 1000, 1000)
        c.SetLeftMargin(.125)
        c.SetLogx() if flux else self.do_nothing()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod=mode, unit=' [kHz/cm2]' if flux else ''), y_tit='{tit} [cm]'.format(tit=title), y_off=1.8)
        gr.Draw('alp')
        gROOT.SetBatch(0)
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('BeamProfile{tit}{dir}{mar}'.format(tit=title, dir=direction.title(), mar=fit_margin * 100), sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])
        return gr

    def draw_xy_profiles(self, flux=True, show=True, fitx=.4, fity=.7):
        gr1 = self.draw_beam_info(mean_=True, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr2 = self.draw_beam_info(mean_=True, flux=flux, show=False, direction='y', fit_margin=fity)
        gr3 = self.draw_beam_info(mean_=False, flux=flux, show=False, direction='x', fit_margin=fitx)
        gr4 = self.draw_beam_info(mean_=False, flux=flux, show=False, direction='y', fit_margin=fity)
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Pulse Height Distribution', 1500, 1500)
        c.Divide(2, 2)
        for i, gr in enumerate([gr1, gr2, gr3, gr4], 1):
            self.format_histo(gr, y_off=1.3)
            pad = c.cd(i)
            pad.SetLogx() if flux else self.do_nothing()
            pad.SetBottomMargin(.15)
            gr.Draw('alp')
        gROOT.SetBatch(0)
        self.save_plots('BeamProfileOverview', sub_dir=self.save_dir)
        self.RootObjects.append([c, gr1, gr2, gr3, gr4])

    def draw_beam_profiles(self, show=True, direction='x'):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [ana.draw_beam_profile(mode=direction, show=False, fit=False) for ana in self.collection.itervalues()]
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        legend = TLegend(.4, .6 - self.get_number_of_analyses() * 0.03, .6, .6)
        for i, h in enumerate(histos):
            h.SetStats(0)
            self.normalise_histo(h)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, '{0:6.2f} kHz/cm'.format(self.collection.values()[i].get_flux()) + '^{2}', 'l')
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        gROOT.SetBatch(0)
        self.save_plots('AllBeamProfiles{mod}'.format(mod=direction.title()), sub_dir=self.save_dir)
        self.RootObjects.append([c, legend] + histos)

    # endregion

    # ====================================================================================
    # region TRACKS
    def show_chi2s(self, mode=None, show=True, disto=False):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        self.reset_colors()
        histos = [ana.show_chi2(mode=mode, show=False, prnt=False) for ana in self.collection.itervalues()]
        yq = zeros(1)
        cuts = []
        for h in histos:
            h.GetQuantiles(1, yq, array([.9]))
            cuts.append(yq[0])
        cut = min(cuts)
        xmax = max(cuts) * 1.1
        legend = self.make_legend(y2=.95, nentries=self.get_number_of_analyses() + 1)
        stack = THStack('hx2', '#chi^{{2}}{mode}'.format(mode=' in ' + mode if mode is not None else ''))
        ymax = 0
        for i, h in enumerate(histos):
            self.format_histo(h, stats=0, color=self.get_color(), lw=2)
            self.normalise_histo(h, to100=True)
            stack.Add(h)
            ymax = max(ymax, h.GetBinContent(h.GetMaximumBin()))
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux()), 'l')
        self.format_histo(stack, x_tit='#chi^{2}', y_tit='Fraction of Events [%]', y_off=1.5, draw_first=True)
        stack.GetXaxis().SetRangeUser(0, xmax)
        mode = '' if mode is None else mode
        self.draw_histo(stack, '', show, self.save_dir, lm=.15, draw_opt='nostack', l=legend)
        l = self.draw_vertical_line(cut, -1e9, 1e9, color=2, style=2)
        legend.AddEntry(l, 'cut: 90% quantile', 'l')
        if disto:
            nominal_chi2 = TF1('f', '[1]*ROOT::Math::chisquared_pdf(x, {ndf})'.format(ndf=4 if not mode else 2), 0, xmax)
            histos[0].Fit(nominal_chi2, 'q0')
            nominal_chi2.SetNpx(1000)
            nominal_chi2.Draw('same')
            self.RootObjects.append(nominal_chi2)
        self.save_plots('AllChi2{mod}'.format(mod=mode.title()))

    def draw_all_chi2s(self, show=True):
        self.show_chi2s(show=show)
        self.show_chi2s('x', show)
        self.show_chi2s('y', show)

    def show_angles(self, mode='x', show=True):
        self.reset_colors()
        histos = [ana.draw_angle_distribution(mode=mode, show=False, print_msg=False) for ana in self.collection.itervalues()]

        legend = self.make_legend(nentries=self.get_number_of_analyses())
        stack = THStack('has', 'Track Angles in {mode}'.format(mode=mode.title()))
        for i, h in enumerate(histos):
            self.format_histo(h, stats=0, color=self.get_color())
            self.normalise_histo(h, to100=True)
            stack.Add(h)
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(self.collection.values()[i].get_flux()), 'l')
        self.format_histo(stack, x_tit='Angle [deg]', y_tit='Fraction of Events [%]', y_off=1.5, draw_first=True)
        stack.GetXaxis().SetRangeUser(-3, 4)
        self.RootObjects.append(self.save_histo(stack, 'AllTrackAngles{mod}'.format(mod=mode.title()), show, self.save_dir, lm=.15, draw_opt='nostack', l=legend))

    def draw_both_angles(self, show=True):
        self.show_angles('x', show)
        self.show_angles('y', show)

    def show_angle_peaks(self, mode='x', sigma=False, flux=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        tit = '{mod} of Track Angle Distribution in {dir}'.format(mod='Sigma' if sigma else 'Mean', dir=mode.title())
        gr = self.make_tgrapherrors('gr', tit)
        i = 0
        for key, ana in self.collection.iteritems():
            fit = ana.calc_angle_fit(mode, show=False)
            x = ana.run.Flux if flux else key
            par = 2 if sigma else 1
            gr.SetPoint(i, x, fit.Parameter(par))
            gr.SetPointError(i, 0, fit.ParError(par))
            i += 1
        c = TCanvas('c', 'Angle Peaks', 1000, 1000)
        c.SetLeftMargin(.12)
        if flux:
            c.SetLogx()
        self.format_histo(gr, x_tit='{mod}{unit}'.format(mod='Flux' if flux else 'Run', unit=' [kHz/cm2]' if flux else ''), y_tit='Sigma [deg]' if sigma else 'Mean [deg]', y_off=1.5)
        gr.Draw('ap')
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.save_plots('TrackAngle{0}s{mod}'.format('Sigma' if sigma else 'Mean', mod=mode.title()), canvas=c, sub_dir=self.save_dir)
        self.RootObjects.append([gr, c])

    # endregion

    def draw_currents(self, v_range=None, rel_time=False):
        self.Currents.draw_indep_graphs(rel_time=rel_time, v_range=v_range)

    def make_flux_table(self):
        # for ana in self.collection.itervalues():
        #     print ana.run.RunInfo
        fluxes = OrderedDict(sorted({cs: [] for cs in set([(ana.run.RunInfo['fs11'], ana.run.RunInfo['fs13']) for ana in self.collection.itervalues()])}.iteritems()))
        for ana in self.collection.itervalues():
            fluxes[(ana.run.RunInfo['fs11'], ana.run.RunInfo['fs13'])].append('{0:3.1f}'.format(ana.run.calc_flux()))
        first_col = [[''] * max([len(col) for col in fluxes.itervalues()])]
        header = ['FS11/FSH13'] + ['{0}/{1}'.format(make_col_str(head[0]), make_col_str(head[1])) for head in fluxes.iterkeys()]
        print make_latex_table(header=header, cols=first_col + fluxes.values(), endline=True),
        print make_latex_table_row(['Mean'] + ['{0:3.1f}'.format(calc_mean(flux)[0]) for flux in fluxes.values()]),
        print make_latex_table_row(['Sigma'] + ['{0:3.1f}'.format(calc_mean(flux)[1]) for flux in fluxes.values()]),
        print '\\bottomrule'

    def make_signal_analysis(self, saveplots=True):
        """
        Run all available signal analyises together and plot them in an overview.
        :param saveplots:
        """
        start_time = time()
        self.draw_pulse_heights(show=False)
        self.draw_pedestals(show=False)
        self.draw_mean_fwhm(draw=False)
        c = TCanvas('c', 'overview', 800, 1000)
        c.Divide(1, 3)
        plots = [self.PulseHeight, self.Pedestal, self.FWHM]
        for i, plot in enumerate(plots, 1):
            pad = c.cd(i)
            pad.SetLogx()
            plot.Draw()
        if saveplots:
            self.save_plots('Overview Plot', sub_dir=self.save_dir, canvas=c)
        self.RootObjects.append(c)

        print '\nThe preanalysis for this selection took', self.print_elapsed_time(start_time)

    def get_systematic_error_table(self, latex=False):
        f = open('PlotsFelix/table_{tc}_{rp}_{dia}.txt'.format(tc=self.TESTCAMPAIGN, rp=self.RunPlan, dia=self.DiamondName), 'w')
        l = '&\t' if latex else ''
        print 'Pulser', 'SignalFlux', 'Signal < 150', 'Signal < 80', 'Full Pulser Spread', 'Single Pulser Spreads'
        pulser_err = self.calc_all_pulser_errors().values()
        p_mean = calc_mean([x[0] for x in pulser_err])
        out = '{0:2.2f}\t{l}{1:2.2f}\t{l}{2:2.2f}\t{l}'.format(min(x[1] / x[0] for x in pulser_err) * 100, max(x[1] / x[0] for x in pulser_err) * 100, p_mean[1] / p_mean[0] * 100, l=l)
        scan_errors = self.draw_all_ph_distributions(show=False).values()
        out += '{0:2.2f}\t{l}{1:2.2f}\t{l}'.format(min(x[1] / x[0] for x in scan_errors) * 100, max(x[1] / x[0] for x in scan_errors) * 100, l=l)
        flux_errors = self.draw_ph_distributions_below_flux(flux=150, show=False)
        out += '{0:2.2f}\t{l}'.format(flux_errors[1] / flux_errors[0] * 100, l=l)
        flux_errors = self.draw_ph_distributions_below_flux(flux=80, show=False)
        out += '{0:2.2f}\t{l}'.format(flux_errors[1] / flux_errors[0] * 100, l=l)
        out += '{0:2.2f}\t{l}'.format(self.calc_full_pedestal_spread(), l=l)
        spreads = self.calc_all_pedestal_spreads().values()
        try:
            min_spread = min(x for x in spreads if x)
        except ValueError:
            min_spread = 0
        out += '{0:2.2f}\t{l}{1:2.2f}'.format(min_spread, max(x for x in spreads), l=l)
        print out
        f.write(out)

    def get_runs_by_collimator(self, fs11=65, fsh13=.5):
        return [key for key, ana in self.collection.iteritems() if ana.run.RunInfo['fs11'] == fs11 and ana.run.RunInfo['fsh13'] == fsh13]

    def get_runs_below_flux(self, flux):
        return [key for key, ana in self.collection.iteritems() if ana.run.Flux < flux]

    def select_runs_in_range(self, start, stop):
        new_collection = OrderedDict()
        for key, ana in self.collection.iteritems():
            if start <= key <= stop:
                new_collection[key] = ana
        if not new_collection:
            print 'You did not select any run! No changes were made!'
        else:
            self.collection = new_collection

    def set_channel(self, ch):
        """
        Sets the channels to be analysed by the SignalAnalysisobjects.
        :param ch: int (0 || 3)
        """
        for ana in self.collection.values():
            ana.set_channel(ch)

    def get_fluxes(self):
        flux = OrderedDict()
        for key, ana in self.collection.iteritems():
            flux[key] = ana.run.get_flux()
        return flux

    @staticmethod
    def get_mode(flux, vs_time=False):
        string = 'Run'
        string = 'Flux' if flux else string
        return 'Time' if vs_time else string

    def get_first_analysis(self):
        return self.collection.values()[0]

    def get_last_analysis(self):
        return self.collection.values()[-1]

    def get_run_numbers(self):
        """ :return: sorted list of run numbers in AnalysisCollection instance """
        return sorted(self.collection.keys())

    def get_number_of_analyses(self):
        """ :return: number of analyses that the analysis collection object contains """
        return len(self.collection)

    def show_information(self):
        print "ANALYSIS COLLECTION INFO:"
        self.get_first_analysis().print_info_header()
        for ana in self.collection.itervalues():
            ana.print_information(header=False)

    def print_loaded(self):
        print '\033[1A\rRuns {0}-{1} were successfully loaded!{2}\n'.format(self.runs[0], self.runs[-1], 20 * ' ')

    def set_verbose(self, status):
        self.verbose = status
        for ana in self.collection.itervalues():
            ana.verbose = status
            ana.Pulser.verbose = status

    @staticmethod
    def make_x_tit(mode, flux):
        return '{mod}{unit}'.format(mod=mode, unit=' [kHz/cm^{2}]' if flux else '')


if __name__ == "__main__":
    st = time()
    main_parser = ArgumentParser()
    main_parser.add_argument('runplan', nargs='?', default=3)
    main_parser.add_argument('dia', nargs='?', default=1, type=int)
    main_parser.add_argument('dia2', nargs='?', default=1, type=int)
    main_parser.add_argument('-tc', '--testcampaign', nargs='?', default='')
    main_parser.add_argument('-t', '--tree', default=True, action='store_false')
    main_parser.add_argument('-r', '--runs', action='store_true')
    args = main_parser.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    run_plan = args.runplan
    diamond = args.dia
    a = Elementary(tc, True, get_resolution())
    sel = RunSelection(testcampaign=tc)
    sel.select_runs_from_runplan(run_plan, ch=diamond) if not args.runs else sel.select_runs([int(args.runplan), int(args.dia if args.dia else args.runplan)], args.dia2 if args.dia2 else 1)
    a.print_banner('STARTING PAD-ANALYSIS COLLECTION OF RUNPLAN {0}'.format(run_plan))
    a.print_testcampaign()

    z = AnalysisCollection(sel, load_tree=args.tree, verbose=True)
    z.print_loaded()
    z.print_elapsed_time(st, 'Instantiation')
    if args.runs:
        z.Currents.draw_indep_graphs()
        raw_input('Press any button to exit')
