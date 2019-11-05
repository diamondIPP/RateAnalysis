#! /usr/bin/env python
from ROOT import THStack, TF1, TProfile2D

from currents import Currents
from InfoLegend import InfoLegend
from run_selection import RunSelection
from VoltageScan import VoltageScan
from analysis import *
from dut_analysis import DUTAnalysis
from telescope_analysis import TelecopeAnalysis
from pad_analysis import PadAnalysis


class AnalysisCollection(Analysis):
    """ Analysis of the various runs of a single runplan. """

    def __init__(self, run_plan, dut_nr, test_campaign=None, load_tree=True, verbose=False):

        Analysis.__init__(self, test_campaign, verbose)
        self.print_start(run_plan)

        # Run Selection Info
        self.RunSelection = RunSelection(self.TCString, run_plan, dut_nr, verbose)
        self.Runs = array(self.RunSelection.get_selected_runs())
        self.NRuns = len(self.Runs)
        self.RunPlan = self.RunSelection.SelectedRunplan
        self.DUTName = self.RunSelection.SelectedDUT
        self.DUTNumber = self.RunSelection.SelectedDUTNr
        self.Bias = self.RunSelection.SelectedBias
        self.Type = self.RunSelection.SelectedType
        self.Fluxes = self.RunSelection.get_selected_fluxes()

        self.print_start_info()

        # Loading the Trees and Time Vectors in Parallel
        self.LoadTree = load_tree
        self.Threads = load_root_files(self.RunSelection, load_tree)

        # Make Common Pickles
        self.MinFluxRun, self.MaxFluxRun = self.get_high_low_rate_runs()
        self.generate_common_pickles()  # make sure the pickles for the common cuts exist

        # Loading the Single Analyses
        self.Analyses = self.load_analyses()
        self.Analysis = self.load_dummy()  # dummy to get access to the methods
        self.FirstAnalysis = self.Analyses.values()[0]
        self.LastAnalysis = self.Analyses.values()[-1]
        self.Bins = self.FirstAnalysis.Bins if load_tree else None
        self.StartTime = self.FirstAnalysis.Run.StartTime if self.LoadTree else time_stamp(self.FirstAnalysis.Run.LogStart)

        # Directory for the Plots
        self.set_save_directory(join(self.DUTName, 'runplan{}'.format(self.RunPlan)))

        # Sub Classes
        self.VoltageScan = VoltageScan(self)
        self.InfoLegend = InfoLegend(self)
        self.Currents = Currents(self)

        self.print_finished()

    def __del__(self):
        print '\n good bye... '

    def draw_all(self, redo=False):
        pass

    def show_information(self):
        print_table(rows=concatenate(self.get_values('', self.Analysis.show_information, pbar=False, prnt=False)), header=self.FirstAnalysis.get_info_header())

    def print_loaded(self):
        print '\033[1A\rRuns {0}-{1} were successfully loaded!{2}\n'.format(self.Runs[0], self.Runs[-1], 20 * ' ')

    # ----------------------------------------
    # region INIT
    def print_start_info(self):
        bias = ['{:+8.0f}'.format(bias) for bias in self.RunSelection.get_selected_biases()]
        times = [datetime.fromtimestamp(duration).strftime('%H:%M').rjust(15) for duration in self.RunSelection.get_selected_durations()]
        rows = zip(self.Runs, ['{:14.1f}'.format(flux.n) for flux in self.Fluxes], bias, times)
        print_table(header=['Run', 'Flux [kHz/cm2]', 'Bias [V]', 'Duration [hh:mm]'], rows=rows, prnt=self.Verbose)

    def get_high_low_rate_runs(self):
        return self.Runs[where(self.Fluxes == self.Fluxes.min())[0]][0], self.Runs[where(self.Fluxes == self.Fluxes.max())[0]][0]

    def generate_common_pickles(self):
        if self.LoadTree:
            self.generate_slope_pickle()
            self.generate_threshold_pickle()

    def generate_slope_pickle(self):
        picklepath = self.make_pickle_path('TrackAngle', 'x', run=self.MinFluxRun)
        if not file_exists(picklepath):
            TelecopeAnalysis(self.MinFluxRun, self.TCString, verbose=self.Verbose)

    def generate_threshold_pickle(self):
        picklepath = self.make_pickle_path('Cuts', 'SignalThreshold', run=self.MaxFluxRun, ch=self.DUTNumber)
        if not file_exists(picklepath) and self.RunSelection.get_type(self.MaxFluxRun) == 'pad':
            PadAnalysis(self.MaxFluxRun, self.DUTNumber, self.TCString, self.Verbose)

    def load_analysis(self, run_number):
        return DUTAnalysis(run_number, self.DUTNumber, self.TCString, self.Threads[run_number].Tuple, self.Threads[run_number].Time, self.Verbose, prnt=False)

    @staticmethod
    def load_dummy():
        return DUTAnalysis

    def load_analyses(self):
        """ Creates and adds Analysis objects with run numbers in runs. """
        analyses = OrderedDict()
        for run in self.Runs:
            analysis = self.load_analysis(run)
            if self.LoadTree:
                analysis.Cut.set_high_low_rate_run(low_run=self.MinFluxRun, high_run=self.MaxFluxRun)
                analysis.Cut.update_all_cut()
            analyses[run] = analysis
        return analyses

    def close_files(self):
        for ana in self.Analyses.itervalues():
            ana.Run.tree.Delete()
            ana.Run.RootFile.Close()

    def delete_trees(self):
        for ana in self.Analyses.itervalues():
            ana.Tree.Delete()
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region BINNING
    def get_binning(self, bin_width=None, rel_time=False):
        bins = concatenate([ana.Bins.get_raw_time(bin_width, t_from_event=True)[1] for ana in self.Analyses.itervalues()])
        return [bins.size - 1, bins - (self.FirstAnalysis.Run.StartTime if rel_time else 0)]

    def get_time_bins(self, bin_width=None, only_edges=False):
        bins = self.fix_t_arrays([ana.Bins.get_time(bin_width)[1] for ana in self.Analyses.itervalues()])
        bins = array([[b[0], b[-1]] for b in bins]).flatten() if only_edges else concatenate(bins)
        return [bins.size - 1, bins]

    def get_raw_time_bins(self, bin_width=None, only_edges=False, t_from_event=False):
        bins = self.fix_t_arrays([ana.Bins.get_raw_time(bin_width, t_from_event=t_from_event)[1] for ana in self.Analyses.itervalues()])
        bins = array([[b[0], b[-1]] for b in bins]).flatten() if only_edges else concatenate(bins)
        return [bins.size - 1, bins]

    def fix_t_arrays(self, t_arrays):
        """ Add the logged time between two runs if the start time of a run is less than the stop time of the previous. """
        for i in xrange(len(t_arrays) - 1):
            diff = t_arrays[i + 1][0] - t_arrays[i][-1]  # check the time difference from the first bin of the next array and the last of the current
            if diff < 0:
                t_arrays[i + 1] += self.get_break_time(i) - diff
        return t_arrays

    def get_break_time(self, ind):
        # because the log stop is not so reliable...
        return (self.get_ana(ind + 1).Run.LogStart - self.get_ana(ind).Run.LogStart - self.get_ana(ind).Run.Duration).total_seconds()

    # endregion BINNING
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_ana(self, ind):
        return self.Analyses.values()[ind]

    def get_analyses(self):
        return self.Analyses.itervalues()

    def get_hv_name(self):
        return self.Currents.Name

    def get_fluxes(self, rel_error=0.):
        return OrderedDict((key, ana.Run.get_flux(rel_error)) for key, ana in self.Analyses.iteritems())

    def get_times(self):
        return OrderedDict((key, ana.Run.get_time()) for key, ana in self.Analyses.iteritems())

    def get_x_var(self, vs_time=False, rel_time=False, rel_error=0.):
        return array(self.get_times().values()) - (time_stamp(self.FirstAnalysis.Run.LogStart) + 3600 if rel_time else 0) if vs_time else array(self.get_fluxes(rel_error).values())

    def get_irradiation(self):
        return self.FirstAnalysis.get_irradiation()

    def get_attenuator(self):
        return self.FirstAnalysis.get_attenuator()

    def get_hv_device(self):
        return self.FirstAnalysis.Currents.Name

    def get_currents(self):
        return OrderedDict((key, ana.Currents.get_current()) for key, ana in self.Analyses.iteritems())

    def get_values(self, string, f, pbar=True, *args, **kwargs):
        return self.generate_plots(string, f, pbar, *args, **kwargs)

    def generate_plots(self, string, f, pbar=True, *args, **kwargs):
        self.info('Generating {} ...'.format(string), prnt=pbar)
        self.PBar.start(self.NRuns) if pbar else do_nothing()
        plots = []
        for i, ana in enumerate(self.get_analyses()):
            plots.append(f(ana, *args, **kwargs))
            self.PBar.update(i) if pbar else do_nothing()
        return plots

    @staticmethod
    def get_mode(vs_time):
        return 'Time' if vs_time else 'Flux'

    def get_sm_std_devs(self, redo=False):
        pickle_path = self.make_pickle_path('Uniformity', 'STD', self.RunPlan, self.DUTNumber)

        def f():
            self.info('Getting STD of Signal Map ... ')
            return OrderedDict((key, ana.get_sm_std(redo=redo)) for key, ana in self.Analyses.iteritems())

        return do_pickle(pickle_path, f, redo=redo)

    def get_sm_std(self, redo=False, low=False, high=False):
        pickle_path = self.make_pickle_path('Uniformity', 'SMSTD', self.RunPlan, self.DUTNumber, suf='{}{}'.format(int(low), int(high)))
        runs = self.get_runs_below_flux(110) if low else self.get_runs_above_flux(2000) if high else self.Runs

        def f():
            return mean_sigma([v for run, v in self.get_sm_std_devs(redo).iteritems() if run in runs]) if runs else make_ufloat((0, 0))

        return do_pickle(pickle_path, f, redo=redo)

    def get_pulse_heights(self, bin_width=None, redo=False):
        phs = OrderedDict()
        self.info('Getting pulse heights ... ')
        self.PBar.start(self.NRuns)
        for i, (key, ana) in enumerate(self.Analyses.iteritems()):
            phs[key] = {'flux': ana.get_flux(), 'ph': ana.get_pulse_height(bin_width, redo=redo), 'time': ana.Run.get_time()}
            self.PBar.update(i)
        return phs

    def get_runs_by_collimator(self, fs11=65, fsh13=.5):
        return [key for key, ana in self.Analyses.iteritems() if ana.Run.RunInfo['fs11'] == fs11 and ana.Run.RunInfo['fs13'] == fsh13]

    def get_runs_below_flux(self, flux):
        return [key for key, ana in self.Analyses.iteritems() if ana.Run.Flux <= flux]

    def get_runs_above_flux(self, flux):
        return [key for key, ana in self.Analyses.iteritems() if ana.Run.Flux >= flux]

    def get_repr_error(self, flux, show=True, redo=False):

        pickle_path = self.make_pickle_path('Errors', 'Repr', self.RunPlan, self.DUTName, suf=flux)

        def f():
            runs = self.get_runs_below_flux(flux)
            if not runs:
                return 0
            data = self.get_pulse_heights(redo=redo)
            values = [data[run]['ph'] for run in runs]
            fluxes = [data[run]['flux'] for run in runs]
            gr = self.make_tgrapherrors('gr_re', 'Pulse Heights Below {f} kHz/cm^{{2}}'.format(f=flux), x=fluxes, y=values)
            self.format_statbox(entries=2, only_fit=True)
            gr.Fit('pol0', 'qs{s}'.format(s='' if show else '0'))
            format_histo(gr, x_tit='Flux [kHz/cm^{2}]', y_tit='Mean Pulse Height [au]', y_off=1.7)
            self.save_histo(gr, 'ReprErrors', show, draw_opt='ap', lm=.14, prnt=show)
            if len(values) == 1:
                return .01  # take 1% if there is only one measurement below the given flux
            m, s = mean_sigma(values)
            return s / m

        return do_pickle(pickle_path, f, redo=redo)

    def get_uniformities(self, bins=10, redo=False, low_flux=False, high_flux=False):
        runs = self.get_runs_below_flux(110) if low_flux else self.get_runs_above_flux(2000) if high_flux else self.Runs
        values = [array(self.Analyses[run].get_uniformity(bins, redo=redo)) for run in runs]
        return array([v if v[0] is not None else full(3, ufloat(0, 0)) for v in values])

    def get_mean_uniformity(self, bins=10, redo=False, low_flux=False, high_flux=False):
        values = self.get_uniformities(bins, redo, low_flux, high_flux)
        return [mean_sigma(values[:, i][where(values[:, i] > 0)[0]]) for i in xrange(values[0].size)]

    @staticmethod
    def get_x_tit(vs_time):
        return '{mod}{unit}'.format(mod='Time' if vs_time else 'Flux', unit=' [hh:mm]' if vs_time else ' [kHz/cm^{2}]')

    def get_tax_off(self, vs_time, rel_time=False):
        return None if not vs_time else self.StartTime if rel_time else 0

    def get_xrange(self, vs_time):
        return None if vs_time else self.Bins.FluxRange

    def get_x_args(self, vs_time, rel_time):
        return {'x_tit': self.get_x_tit(vs_time), 't_ax_off': self.get_tax_off(vs_time, rel_time), 'x_range': self.get_xrange(vs_time)}

    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def draw_full_pulse_height(self, bin_width=10000, show=True, rel_t=True, redo=False, with_flux=True):

        pickle_path = self.make_pickle_path('PulseHeight', 'FullPH', self.RunPlan, ch=self.DUTNumber, suf=bin_width)

        def f():
            pulse_heights = [get_hist_vec(ana.draw_pulse_height(bin_width, corr=True, redo=redo, show=False)[0]) for ana in self.Analyses.itervalues()]
            h1 = TH1F('hfph', 'Pulse Height for Run Plan {n}'.format(n=self.RunPlan), *self.get_time_bins(bin_width))
            i_bin = 1  # the first bin is empty
            for ary in pulse_heights:
                for ph in ary:
                    h1.SetBinContent(i_bin, ph.n)
                    h1.SetBinError(i_bin, ph.s)
                    i_bin += 1
                i_bin += 1  # there is an empty bin after each run
            return h1

        h = do_pickle(pickle_path, f, redo=redo)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Mean Pulse Height [au]', y_off=.8, fill_color=self.FillColor, stats=0, y_range=[0, h.GetMaximum() * 1.05])
        set_time_axis(h, off=self.FirstAnalysis.Run.StartTime if rel_t else 0)
        c = self.draw_histo(h, show=show, draw_opt='hist', x=1.5, y=.75, lm=.065, gridy=True, rm=.1 if with_flux else None)
        if with_flux:
            h = self.draw_flux(rel_time=rel_t, show=False)
            c.cd()
            self.draw_tpad('pr', margins=[.065, .1, .15, .1], transparent=True, logy=True)
            x_range = [h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()]
            format_histo(h, title=' ', fill_color=2, fill_style=3002, lw=1, y_range=[1, h.GetMaximum() * 1.2], stats=0, y_off=1.05, x_range=x_range)
            h.Draw('histy+')
        self.save_plots('FullPulseHeight')

    def get_pulse_height_graph(self, binning=None, vs_time=False, first_last=True, redo=False, legend=True):

        marker_size = 1
        gStyle.SetEndErrorSize(4)
        ph = self.get_pulse_heights(binning, redo)
        y_values = [dic['ph'] for dic in ph.itervalues()]
        x_values = [dic['time' if vs_time else 'flux'] for dic in ph.itervalues()]
        g = self.make_tgrapherrors('g', 'stat. error', self.get_color(), marker_size=marker_size, x=x_values, y=y_values)
        rel_sys_error = self.get_repr_error(105, show=False)
        y_values = [make_ufloat((v.n, v.s + rel_sys_error * abs(v.n))) for v in y_values]
        g_errors = self.make_tgrapherrors('gerr', 'full error', marker=0, color=602, marker_size=0, x=x_values, y=y_values)
        g_first, g_last = [self.make_tgrapherrors('g1', 'first run', marker=22, color=2, marker_size=marker_size, x=[x_values[i].n], y=[y_values[i].n]) for i in [0, -1]]
        graphs = [g, g_errors]
        graphs += [g_first, g_last] if first_last else []
        leg = self.make_legend(.75, .37, nentries=len(graphs))
        mg = TMultiGraph('mg_ph', 'Pulse Height vs {mod} - {dia}'.format(mod='Time' if vs_time else 'Flux', dia=self.DUTName))
        for gr in graphs:
            leg.AddEntry(gr, gr.GetTitle(), 'l' if gr.GetName() == 'gerr' else 'p')
            mg.Add(gr, 'pl')
        if legend:
            mg.GetListOfFunctions().Add(leg)
        self.reset_colors()
        if vs_time:
            g = mg.GetListOfGraphs()[0]
            for i, (ana, x) in enumerate(zip(self.Analyses.itervalues(), x_values)):
                y, ey = g.GetY()[i], g.GetErrorY(i)
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(self.draw_tlatex(x.n, y + ey * 1.2, '{:1.0f}'.format(ana.get_flux().n), color=1, align=21, size=.02, show=0))
        return mg

    def draw_scaled_pulse_heights(self, scale=1, binning=None, vs_time=False, show=True, y_range=None, redo=False, scale_to_low=False):

        mg = self.get_pulse_height_graph(binning, vs_time, first_last=not vs_time, redo=redo)
        scale_multigraph(mg, scale, scale_to_low)
        xtit = 'Time [hh:mm]' if vs_time else 'Flux [kHz/cm^{2}]'
        y_range1 = [.95, 1.05] if y_range is None else y_range
        format_histo(mg, x_tit=xtit, y_tit='Scaled Pulse Height', y_off=1.75, x_off=1.3, draw_first=True, t_ax_off=0 if vs_time else None, y_range=y_range1, ndivx=503, center_y=1)
        mg.GetXaxis().SetLimits(1, 40000) if not vs_time else do_nothing()
        move_legend(mg.GetListOfFunctions()[0], .25, .20)
        self.draw_histo(mg, '', show, lm=.14, draw_opt='a', logx=not vs_time, grid=vs_time, gridy=True, bm=.18)
        self.draw_irradiation(make_irr_string(self.RunSelection.get_irradiation()))
        self.save_plots('ScaledPulseHeights{}'.format(xtit[:4]))
        return mg.GetListOfGraphs()[0]

    def draw_pulse_heights(self, binning=None, vs_time=False, show=True, show_first_last=True, save_comb=True, y_range=None, redo=False):

        mode = 'Time' if vs_time else 'Flux'
        mg = self.get_pulse_height_graph(binning, vs_time, show_first_last, redo)

        # small range
        format_histo(mg, color=None, x_tit='Time [hh:mm]' if vs_time else 'Flux [kHz/cm^{2}]', y_tit='Signal Pulse Height [mV]', y_off=1.75, x_off=1.3, draw_first=True,
                     t_ax_off=0 if vs_time else None)
        ymin, ymax = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
        yrange = increased_range([ymin, ymax], .5, .15) if y_range is None else y_range
        mg.GetYaxis().SetRangeUser(*yrange)
        mg.GetXaxis().SetLimits(1, 40000) if not vs_time else do_nothing()
        self.save_histo(mg, 'PulseHeight{mod}'.format(mod=mode), show=False, lm=.14, draw_opt='A', logx=not vs_time, grid=vs_time)

        # no zero suppression
        mg1 = mg.Clone()
        mg1.GetListOfFunctions().Clear()
        mg1.SetName('mg1_ph')
        mg1.GetListOfGraphs()[0].SetLineColor(self.Colors[0])
        mg1.GetYaxis().SetRangeUser(0, ymax * 1.1)
        self.save_histo(mg1, 'PulseHeightZero{mod}'.format(mod=mode), not save_comb, lm=.14, draw_opt='Al', logx=not vs_time)
        self.reset_colors()

        if save_comb:
            y_min = increased_range([ymin, ymax], .5)[0] if y_range is None else y_range[0]
            # TODO fix vs time and comb plot
            self.save_combined_pulse_heights(mg, mg1, y_min, show=show, pulser_leg=self.draw_signal_legend)

        return mg

    def draw_signal_legend(self):
        pass

    def draw_signal_distributions(self, show=True, off=3, redo=False):

        stack = THStack('hsd', 'Pulse Height Distributions')
        leg = self.make_legend(.67, .88, nentries=self.NRuns)
        self.info('Generating signal distributions ...')
        histos = []
        self.PBar.start(self.NRuns)
        for i, ana in enumerate(self.Analyses.itervalues(), 1):
            h = ana.draw_signal_distribution(show=False, redo=redo, prnt=False)
            format_histo(h, fill_color=0, fill_style=4000)
            histos.append(h)
            self.PBar.update(i)
        for i, h in enumerate(histos):
            format_histo(h, lw=2, color=self.get_color())
            h.Scale(1 / h.GetMaximum())
            stack.Add(h)
            leg.AddEntry(h, '{0:06.1f} kHz/cm^{{2}}'.format(self.Analyses.values()[i].get_flux().n), 'l')
            # legend.AddEntry(h, '{0:2d} V'.format(self.has_collection.values()[i].Bias), 'l')
        format_histo(stack, y_off=1.55, draw_first=True, x_tit='Pulse Height [au]', y_tit='Number of Entries')
        log_stack = stack.Clone()
        log_stack.SetMaximum(off)
        log_stack.SetNameTitle('hsdl', 'Signal Distribution LogY')
        self.save_histo(stack, 'SignalDistributions', show=False, lm=.13, draw_opt='nostack', leg=leg)
        self.save_histo(log_stack, 'SignalDistributionsLogY', show=False, lm=.13, draw_opt='nostack', logy=True, leg=leg)
        if show:
            c = self.make_canvas('csd1', 'Signal Distributions', x=1.5, y=.75, divide=2)
            legends = [leg, leg.Clone()]
            for i, s in enumerate([stack, log_stack], 1):
                pad = c.cd(i)
                pad.SetLeftMargin(.14)
                pad.SetLogy() if i == 2 else do_nothing()
                s.Draw('nostack')
                legends[i - 1].Draw()
            self.Objects.append([c, legends])
        self.reset_colors()

    def draw_all_ph_distributions(self, bin_width=None, show=False):

        pickle_path = self.make_pickle_path('ph_fit', 'PhDistoFits', run=self.RunPlan, ch=self.DUTNumber, suf=self.Bins.BinSize if bin_width is None else bin_width)

        def f():
            collimator_settings = set([(ana.Run.RunInfo['fs11'], ana.Run.RunInfo['fsh13']) for key, ana in self.Analyses.iteritems()])
            return {s: self.draw_ph_distributions(bin_width, show=show, fs11=s[0], fsh13=s[1]) for s in collimator_settings}
        return do_pickle(pickle_path, f)

    def draw_ph_distributions(self, binning=None, fsh13=.5, fs11=65, show=True):
        runs = self.get_runs_by_collimator(fsh13=fsh13, fs11=fs11)
        return self.draw_combined_ph_distributions(runs, binning, show)

    def draw_ph_distributions_below_flux(self, bin_width=None, flux=150, show=True):
        pickle_path = self.make_pickle_path('Ph_fit', 'PhDistoBel', self.RunPlan, self.DUTNumber, suf='{bin}_{flux}'.format(bin=self.Bins(bin_width), flux=flux))

        def f():
            runs = self.get_runs_below_flux(flux)
            return self.draw_combined_ph_distributions(runs, bin_width, show=False)

        h = do_pickle(pickle_path, f, redo=show)
        return h

    def draw_combined_ph_distributions(self, runs, bin_width=.5, show=True):
        s = THStack('s_phd', 'Pulse Height Distributions')
        self.PBar.start(len(runs))
        leg = self.make_legend(nentries=len(runs))
        for i, run in enumerate(runs):
            set_root_output(False)
            h = self.Analyses[run].draw_ph_pull(show=False, bin_width=bin_width, fit=False, save=False)
            format_histo(h, fill_color=4000, stats=0, line_color=self.get_color())
            s.Add(h)
            leg.AddEntry(h, '{:2.1f} kHz/cm^{{2}}'.format(self.Analyses[run].get_flux().n), 'l')
            self.PBar.update(i)
        self.draw_histo(s, show=show, draw_opt='', leg=leg)
        ls = s.GetStack().Last()
        self.format_statbox(only_fit=True)
        format_histo(s, 'Fit Results', stats=True, x_range=[ls.GetBinCenter(ibin) for ibin in [ls.FindFirstBinAbove(0) - 2, ls.FindLastBinAbove(0) + 2]])
        format_histo(ls, y_range=[0, get_hist_vec(ls).max().n * 1.1])
        self.save_plots('PulseHeightDistributions')
        return ls.GetMean(), ls.GetStdDev()
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region MISCELLANEOUS DRAW
    def draw_currents(self, v_range=None, rel_time=False, averaging=1, with_flux=False, c_range=None, f_range=None, draw_opt='ap', show=True):
        self.Currents.draw_indep_graphs(rel_time=rel_time, v_range=v_range, averaging=averaging, with_flux=with_flux, c_range=c_range, f_range=f_range, show=show, draw_opt=draw_opt)

    def draw_iv(self, show=True):
        g = self.Currents.draw_iv(show=False)
        self.save_histo(g, 'IV', draw_opt='ap', logy=True, lm=.12, show=show)

    def draw_fwhm(self, arg=1, vs_time=False, show=True, redo=False):
        y_values = self.get_uniformities(redo=redo)[:, arg]
        x_values = self.get_times() if vs_time else self.get_fluxes()
        g = self.make_tgrapherrors('gfm', 'Full Width Half Maximum', x=x_values.values(), y=y_values)
        y_tit = 'FWHM [mV]' if arg == 1 else 'Most Probable Value [mV]' if not arg else 'FWHM/MPV'
        format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit=y_tit, y_off=1.3, t_ax_off=0 if vs_time else None)
        self.save_histo(g, 'FWHM', show, lm=.12, logx=not vs_time)
        return g

    def draw_mpv(self, vs_time=False, show=True, redo=False):
        return self.draw_fwhm(0, vs_time, show, redo)

    def draw_uniformity(self, vs_time=False, show=True, redo=False):
        return self.draw_fwhm(2, vs_time, show, redo)

    def draw_sm_std(self, vs_time=False, redo=False, show=True):
        y_values = self.get_sm_std_devs(redo).values()
        x_values = self.get_times() if vs_time else self.get_fluxes()
        g = self.make_tgrapherrors('gstd', 'STD of the Signal Map', x=x_values.values(), y=y_values)
        format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='rel. STD', y_off=1.3, t_ax_off=0 if vs_time else None)
        self.save_histo(g, 'STDSigMap', show, lm=.12, logx=not vs_time)

    def draw_raw_flux(self, bin_width=None, rel_t=False, show=True):
        h = TH1F('hr1', 'Flux of Run Plan {r}'.format(r=self.RunPlan), *self.get_time_bins(bin_width, only_edges=True))
        for i, ana in enumerate(self.Analyses.itervalues()):
            h.SetBinContent(i * 2 + 1, ana.Run.Flux.n)
            h.SetBinError(i * 2 + 1, ana.Run.Flux.s)
        format_histo(h, y_tit='Flux [kHz/cm^{2}]', x_tit='Time [hh:mm]', y_off=.8, fill_color=self.FillColor, stats=0, t_ax_off=self.StartTime if rel_t else 0)
        self.draw_histo(h, x=1.5, y=.75, lm=.065, gridy=True, logy=True)
        self.save_histo(h, 'FluxTime', canvas=get_last_canvas(), show=show, draw_opt='esame')
        return h

    def draw_beam_current(self, bin_width=60, rel_t=True, show=True):
        h = TH1F('hr1', 'Beam Current of Run Plan {r}'.format(r=self.RunPlan), *self.get_raw_time_bins(bin_width))
        values = [get_hist_vec(ana.draw_beam_current_prof(bin_width, show=False, save=False)) for ana in self.Analyses.itervalues()]
        # add last flux value for time bin between runs
        values = concatenate([concatenate([run_values, [run_values[-1]]]) for run_values in values])
        for i, value in enumerate(values, 1):
            h.SetBinContent(i, value.n)
            h.SetBinError(i, value.s)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', y_off=.85, fill_color=self.FillColor, stats=0, markersize=.3, t_ax_off=self.StartTime if rel_t else 0)
        self.save_histo(h, 'AllBeamRate', show=show, draw_opt='hist', x=1.5, y=.75, lm=.065)

    def draw_flux(self, bin_width=5, rel_time=True, show=True):

        pickle_path = self.make_pickle_path('Flux', 'FullFlux', self.RunPlan, self.DUTNumber, bin_width)

        def f():
            if not self.FirstAnalysis.has_branch('rate'):
                return self.draw_raw_flux(rel_t=rel_time, show=show)
            else:
                h1 = TH1F('hff', 'Flux Profile', *self.get_raw_time_bins(bin_width))
                values = [get_hist_vec(ana.draw_flux(bin_width, show=False, prnt=False)) for ana in self.Analyses.itervalues()]
                # add extra zero for time bin between runs
                values = concatenate([concatenate([run_values, [ufloat(0, 0)]]) for run_values in values])
                for i, value in enumerate(values, 1):
                    h1.SetBinContent(i, value.n)
                    h1.SetBinError(i, value.s)
                return h1

        hist = do_pickle(pickle_path, f)
        format_histo(hist, x_tit='Time [hh:mm]', y_tit='Flux [kHz/cm^{2}]', t_ax_off=self.InitTime if rel_time else 0, fill_color=self.FillColor, y_range=[1, 20000], stats=0)
        self.save_histo(hist, 'FluxEvo', x=1.5, y=.75, show=show, logy=True, draw_opt='bar' if not self.FirstAnalysis.has_branch('rate') else '')
        return hist

    def draw_current_flux(self, c_range=None, fit=True, show=True):
        currents = self.get_currents().values()
        g = self.make_tgrapherrors('gcf', 'Leakage Current vs. Flux', x=self.get_fluxes(rel_error=.01).values(), y=currents)
        c_range = [.1, max([c.n for c in currents if c is not None]) * 2] if c_range is None else c_range
        format_histo(g, x_tit='Flux [kHz/cm^{2}]', y_tit='Current [nA]', y_off=1.3, y_range=c_range, draw_first=True, x_range=[1, 20000])
        self.draw_histo(g, 'FluxCurrent', show, lm=.13, draw_opt='ap', logx=True, logy=True)
        if fit:
            self.format_statbox(only_fit=True, y=.33, entries=6, w=.22)
            f = TF1('fcf', 'pol1', .1, 1e5)
            f.SetParLimits(0, .1, 5)
            f.SetParLimits(1, 1e-5, 5e-3)
            g.Fit('fcf', 'q')
        self.draw_irradiation(make_irr_string(self.get_irradiation()))
        self.save_plots('FluxCurrent', show=show)
        return g
    # endregion MISCELLANEOUS DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL MAP
    def draw_signal_map(self, chi2=None, fid=False, res=.7, hitmap=False, redo=False, cut=None, show=True):
        self.PBar.start(self.NRuns)
        maps = []
        for i, ana in enumerate(self.Analyses.values()):
            ana.Cut.set_chi2(chi2) if chi2 else do_nothing()
            maps.append(ana.draw_signal_map(fid=fid, cut='' if cut is None and hitmap else cut, res=res, redo=redo, show=False, prnt=False, hitmap=hitmap))
            self.PBar.update(i)
        h_all = TH2I('hchm', 'Cumulative Diamond Hit Map', *self.Bins.get_global(res, mm=True)) if hitmap else TProfile2D('hcsm', 'Cumulative Signal Map', *self.Bins.get_global(res, mm=True))
        for h in maps:
            h_all.Add(h)
        self.draw_histo(h_all, show=show, lm=.12, rm=.16, draw_opt='colz')
        self.FirstAnalysis.draw_fid_cut(scale=10)
        self.save_plots('Cumulative{s}Map'.format(s='Hit' if hitmap else 'Signal'))

    def draw_signal_maps(self, hitmap=False, redo=False):
        f = PadAnalysis.draw_hitmap if hitmap else PadAnalysis.draw_signal_map
        histos = self.generate_plots('{} maps'.format('hit' if hitmap else 'signal'), f, show=False, prnt=False, redo=redo)
        glob_max = (int(max([h.GetMaximum() for h in histos])) + 5) / 5 * 5
        glob_min = int(min([h.GetMinimum() for h in histos])) / 5 * 5
        for i, h in enumerate(histos):
            format_histo(h, z_range=[glob_min, glob_max]) if not hitmap else do_nothing()
            self.save_histo(h, '{n}Map{nr}'.format(nr=str(i).zfill(2), n='Hit' if hitmap else 'Signal'), show=False, ind=i, draw_opt='colz', rm=.16, lm=.12, prnt=False)

    def draw_signal_map_ratio(self, run1, run2, m=10, n=10, grid=True, show=True):
        h1, h2 = [self.Analyses[run].split_signal_map(m, n, show=False)[0] for run in [run1, run2]]
        xbins, ybins = self.Analyses[run1].split_signal_map(m, n, show=False)[1:]
        h = h1.Clone('hnew')
        for i in xrange((h.GetNbinsX() + 2) * (h.GetNbinsY() + 2)):
            v1, v2 = h1.GetBinContent(i), h2.GetBinContent(i)
            h.SetBinEntries(i, 1)
            h.SetBinContent(i, v1 / v2 if v2 else -1)
        format_histo(h, z_range=[0, 3], stats=0, z_tit='Pulse Height Ratio', title='Signal Map Ratio of Run {} & {}'.format(run1, run2))
        self.draw_histo(h, lm=.12, rm=.16, draw_opt='colzsame', show=show)
        self.draw_grid(xbins, ybins, width=2) if grid else do_nothing()
        self.save_plots('SigMapRatio{}{}'.format(run1, run2))

    def draw_signal_spreads(self, vs_time=True, rel_time=True, show=True):
        spreads = [ana.get_signal_spread(prnt=False) for ana in self.get_analyses()]
        g = self.make_tgrapherrors('gss', 'Relative Spread', x=self.get_x_var(vs_time, rel_time=rel_time), y=spreads)
        format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='Relative Spread [%]', y_off=1.2, t_ax_off=self.get_tax_off(vs_time))
        self.save_histo(g, 'RelativeSpread', lm=.12, logx=not vs_time, show=show)
    # endregion SIGNAL MAP
    # ----------------------------------------

    # ----------------------------------------
    # region BEAM PROFILE
    def draw_beam_info(self, mode='x', fit_margin=.5, vs_time=True, rel_time=True, show=True):
        tits = ['Mean', 'Sigma']
        values = self.get_values('beam profile {}'.format(mode), self.Analysis.draw_beam_profile, show=False, fit=True, fit_range=fit_margin, mode=mode, prnt=False)
        values = [[make_ufloat(value, par=par) for value in values] for par in [1, 2]]
        graphs = [self.make_tgrapherrors('gbi{}'.format(tit), '{} of the Beam Profile in {}'.format(tit, mode.title()), x=self.get_x_var(vs_time, rel_time), y=vals) for tit, vals in zip(tits, values)]
        c = self.make_canvas('cbi', 'Beam Infos {}'.format(mode.title()), x=1.5, y=.75, divide=2, show=show)
        for i, g in enumerate(graphs, 1):
            format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='{tit} [cm]'.format(tit=tits[i - 1]), y_off=1.8, t_ax_off=self.get_tax_off(vs_time, rel_time))
            self.save_histo(g, 'BeamProfile{}{}{:1.0f}'.format(tits[i - 1], mode.title(), fit_margin * 100), lm=.125, show=False, logx=not vs_time)
            self.draw_histo(g, canvas=c.cd(i), lm=.125, show=False, logx=not vs_time)
        self.save_plots('BeamProfile{}{:1.0f}'.format(mode.title(), fit_margin * 100), show=show, canvas=c)
        return graphs

    def draw_xy_beam_info(self, vs_time=True, show=True, fitx=.5, fity=.5):
        graphs = concatenate([self.draw_beam_info(vs_time=vs_time, show=False, mode=m, fit_margin=f) for m, f in zip(['x', 'y'], [fitx, fity])])
        c = self.make_canvas('cbpf', 'Beam Profiles', x=1.5, y=1.5, divide=(2, 2), show=show)
        for i, g in enumerate(graphs, 1):
            format_histo(g, y_off=1.3)
            self.draw_histo(g, logx=not vs_time, canvas=c.cd(i))
        self.save_plots('BeamProfileOverview', show=show, canvas=c)

    def draw_beam_profiles(self, mode='x', show=True):
        histos = self.generate_plots('beam profiles in {}'.format(mode), self.Analysis.draw_beam_profile, show=False, prnt=False, fit=False, mode=mode)
        leg = self.make_legend(nentries=self.NRuns)
        stack = THStack('sbp', 'AllBeamProfiles{mod}'.format(mod=mode.title()))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes().values())):
            format_histo(h, lw=2, stats=0, normalise=True, line_color=self.get_color(), sumw2=False, fill_style=4000)
            stack.Add(h)
            leg.AddEntry(h, '{0:6.2f} kHz/cm^{{2}}'.format(flux.n), 'l')
        self.save_histo(stack, 'AllBeamProfiles{mod}'.format(mod=mode.title()), draw_opt='nostack', show=show, leg=leg)
        self.reset_colors()

    # endregion BEAM PROFILE
    # ----------------------------------------

    # ----------------------------------------
    # region TRACKS
    def draw_chi2(self, mode=None, fit=False, show=True):
        mod_str = '' if mode is None else mode
        histos = self.generate_plots('chi squares {}'.format(mod_str), self.Analysis.draw_chi2, show=False, prnt=False, mode=mode)
        cut_val = self.FirstAnalysis.Cut.CutConfig['chi2{}'.format(mod_str.title() if mod_str else 'X')]
        cuts = zeros((self.NRuns, 2))
        for i, h in enumerate(histos):
            h.GetQuantiles(2, cuts[i], array([cut_val / 100., .93]))
        legend = self.make_legend(nentries=self.NRuns + 1, scale=.7, w=.2)
        stack = THStack('hx2', '#chi^{{2}}{}'.format(' in {}'.format(mod_str) if mod_str else ''))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes().itervalues())):
            format_histo(h, stats=0, color=self.get_color(), lw=2, normalise=True, sumw2=False)
            stack.Add(h)
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(flux.n), 'l')
        format_histo(stack, x_tit='#chi^{2}', y_tit='Fraction of Events [%]', y_off=1.5, draw_first=True, x_range=increased_range([0, cuts.max()], fac_top=.4))
        self.draw_histo(stack, '', show, lm=.15, draw_opt='nostack', leg=legend)
        if mod_str:
            line = self.draw_vertical_line(cuts.min(), -1e9, 1e9, color=2, style=2, name=mod_str)
            legend.AddEntry(line, 'cut: {}% quantile'.format(cut_val), 'l')
            histos[0].GetListOfFunctions().Add(line)
        if fit:
            nominal_chi2 = TF1('f', '[1]*ROOT::Math::chisquared_pdf(x, {ndf})'.format(ndf=4 if not mode else 2), 0, cuts.max())
            histos[0].Fit(nominal_chi2, 'q0')
            nominal_chi2.SetNpx(1000)
            stack.GetHists().Last().GetListOfFunctions().Add(nominal_chi2)
            legend.AddEntry(nominal_chi2, '   #chi^{2} fit', 'l')
        histos[0].GetListOfFunctions().Add(legend)
        self.save_plots('AllChi2{mod}'.format(mod=mod_str.title()))
        self.reset_colors()
        return stack

    def draw_all_chi2s(self, show=True, fit=False):
        stacks = [self.draw_chi2(mode, fit, show=False) for mode in [None, 'x', 'y']]
        c = self.make_canvas('ccss', 'AllChiSquares', x=2, y=.66, divide=3, show=show)
        for i, s in enumerate(stacks, 1):
            self.draw_histo(s, draw_opt='nostack', canvas=c.cd(i), show=show)
        self.save_plots('AllChiSquares', show=show, canvas=c)

    def draw_angle(self, mode='x', show=True):
        histos = self.generate_plots('angle distribution {}'.format(mode), self.Analysis.draw_angle_distribution, mode=mode, show=False, prnt=False)
        legend = self.make_legend(nentries=self.NRuns, scale=.7, w=.2)
        stack = THStack('has', 'Track Angles in {mode}'.format(mode=mode.title()))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes().itervalues())):
            format_histo(h, stats=0, color=self.get_color(), normalise=True, sumw2=False)
            stack.Add(h)
            legend.AddEntry(h, '{0: 6.0f} kHz/cm^{{2}}'.format(flux.n), 'l')
        histos[0].GetListOfFunctions().Add(legend)
        format_histo(stack, x_tit='Angle [deg]', y_tit='Fraction of Events [%]', y_off=1.9, draw_first=True, x_range=[-3, 4])
        self.save_histo(stack, 'AllTrackAngles{mod}'.format(mod=mode.title()), show, lm=.15, draw_opt='nostack', leg=legend)
        self.reset_colors()
        return stack

    def draw_both_angles(self, show=True):
        stacks = [self.draw_angle(mode, show=False) for mode in ['x', 'y']]
        c = self.make_canvas('ccss', 'AllAngles', x=1.5, y=.75, divide=2, show=show)
        for i, s in enumerate(stacks, 1):
            self.draw_histo(s, draw_opt='nostack', canvas=c.cd(i), show=show, lm=.15)
        self.save_plots('AllChiAngles', show=show, canvas=c)

    def draw_mean_angles(self, mode='x', vs_time=True, rel_time=True, show=True):
        values = [ufloat(value[0].n, value[1].n / 2) for value in self.get_values('mean angles {}'.format(mode), self.Analysis.get_mean_angle, mode=mode)]
        g = self.make_tgrapherrors('gma', 'Mean Track Angles in {}'.format(mode.title()), x=self.get_x_var(vs_time, rel_time), y=values)
        format_histo(g, y_tit='Mean Angle [deg]', y_off=1.3, **self.get_x_args(vs_time, rel_time))
        self.save_histo(g, 'MeanTrackAngle'.format(mode.title()), show, logx=not vs_time)
    # endregion TRACKS
    # ----------------------------------------

    # region GENERATE PLOTS
    # ----------------------------------------
    def draw_run_currents(self):
        self.get_values('currents', self.Analysis.get_current)
        self.generate_plots('currents', self.Analysis.draw_current, show=False)

    def draw_chi2s(self):
        self.generate_plots('chi2s', self.Analysis.draw_all_chi2, show=False, prnt=False)

    def draw_angles(self):
        self.generate_plots('angles', self.Analysis.draw_both_angles, show=False, prnt=False)

    def draw_occupancies(self):
        self.generate_plots('occupancies', self.Analysis.draw_occupancies, show=False, prnt=False, cluster=True)

    def select_runs_in_range(self, start, stop):
        new_collection = OrderedDict((key, ana) for key, ana in self.Analyses.iteritems() if start <= key <= stop)
        if not new_collection:
            warning('You did not select any run! No changes were made!')
        else:
            self.Analyses = new_collection

    def set_verbose(self, status):
        self.Verbose = status
        self.VoltageScan.Verbose = status
        for ana in self.Analyses.itervalues():
            ana.verbose = status
            ana.Cut.verbose = status


if __name__ == '__main__':

    p = init_argparser(run=10, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)
    p.add_argument('-r', '--runs', action='store_true')
    p.add_argument('-d', '--draw', action='store_true')
    p.add_argument('-rd', '--redo', action='store_true')
    pargs = p.parse_args()

    z = AnalysisCollection(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
    z.print_loaded()
    if pargs.runs:
        z.Currents.draw_indep_graphs()
        raw_input('Press any button to exit')
    if pargs.draw:
        z.draw_all(pargs.redo)
