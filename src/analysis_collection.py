#! /usr/bin/env python
from __future__ import print_function
from ROOT import THStack, TF1, TProfile2D, TMultiGraph, TMath

from currents import Currents
from InfoLegend import InfoLegend
from run_selection import RunSelection
from VoltageScan import VoltageScan
from analysis import *
from dut_analysis import DUTAnalysis
from telescope_analysis import TelecopeAnalysis
from numpy import histogram, cumsum, sort, split, exp, log, ones


class AnalysisCollection(Analysis):
    """ Analysis of the various runs of a single runplan. """

    def __init__(self, run_plan, dut_nr, test_campaign=None, load_tree=True, verbose=False):

        Analysis.__init__(self, test_campaign, verbose)

        # Run Selection Info
        self.RunSelection = RunSelection(self.TCString, run_plan, dut_nr, verbose)
        self.Runs = array(self.RunSelection.get_selected_runs())
        self.NRuns = len(self.Runs)
        self.RunPlan = self.RunSelection.SelectedRunplan
        self.DUT = self.RunSelection.SelectedDUT
        self.print_start(run_plan, prnt=load_tree, dut=self.DUT.Name)
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
        self.set_save_directory(join(self.DUT.Name, 'runplan{}'.format(self.RunPlan)))

        # Sub Classes
        self.VoltageScan = VoltageScan(self)
        self.InfoLegend = InfoLegend(self)
        self.Currents = Currents(self)

        self.print_finished() if verbose else print()

    def __del__(self):
        print('\n good bye... ')

    def draw_all(self, redo=False):
        pass

    def show_information(self):
        print_table(rows=concatenate(self.get_values('', self.Analysis.show_information, pbar=False, prnt=False)), header=self.FirstAnalysis.get_info_header())

    def print_loaded(self):
        print('\033[1A\rRuns {0}-{1} were successfully loaded!{2}\n'.format(self.Runs[0], self.Runs[-1], 20 * ' '))

    def remove_pickles(self):
        files = glob(join(self.PickleDir, '*', '*{}*_{}*'.format(self.TCString, self.RunPlan)))
        self.info('Removing {} pickle files for run plan {}'.format(len(files), self.RunPlan))
        for f in files:
            remove_file(f)

    # ----------------------------------------
    # region INIT
    def print_start_info(self):
        bias = ['{:+8.0f}'.format(bias) for bias in self.RunSelection.get_selected_biases()]
        times = [datetime.fromtimestamp(duration - 3600).strftime('%H:%M').rjust(15) for duration in self.RunSelection.get_selected_durations()]
        rows = zip(self.Runs, ['{:14.1f}'.format(flux.n) for flux in self.Fluxes], bias, self.RunSelection.get_selected_start_times(), times)
        print_table(header=['Run', 'Flux [kHz/cm2]', 'Bias [V]', 'Start', 'Duration [hh:mm]'], rows=rows, prnt=self.Verbose)

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
        pass

    def load_analysis(self, run_number):
        return DUTAnalysis(run_number, self.DUT.Number, self.TCString, self.Threads[run_number].Tuple, self.Threads[run_number].Time, self.Verbose, prnt=False)

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
                analysis.Cut.reload()
            analyses[run] = analysis
        self.Threads = None
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

    def get_analyses(self, runs=None):
        return self.Analyses.values() if runs is None else [ana for key, ana in self.Analyses.iteritems() if key in runs]

    def get_hv_name(self):
        return self.Currents.Name

    def get_fluxes(self, rel_error=0., corr=True, runs=None, avrg=False, pbar=True):
        picklepath = self.make_simple_pickle_path(sub_dir='Flux', run='{}', dut='')
        pbar = False if not self.FirstAnalysis.has_branch('rate') else pbar
        return self.get_run_values('fluxes', DUTAnalysis.get_flux, runs, pbar, avrg=avrg, picklepath=picklepath, rel_error=rel_error, corr=corr)

    def get_flux_splits(self, redo=False, show=True):
        def f():
            values = sort([flux.n for flux in self.get_fluxes()])
            h = TH1F('hmf', 'Fluxes', *log_bins(50, 1, 1e5))
            h.FillN(values.size, values, full(values.size, 1, 'd'))
            self.format_statbox(entries=True, h=.2)
            format_histo(h, x_tit='Flux [kHz/cm^{2}]', y_tit='Number of Entries', y_off=.5, x_off=1.2, lab_size=.05, tit_size=.05)
            self.draw_histo(h, lm=.07, logx=True, show=show, y=.75, x=1.5, bm=.14)
            s = TSpectrum(20)
            s.Search(h, 1)
            bins = sorted(s.GetPositionX()[i] for i in xrange(s.GetNPeaks()))
            split_bins = histogram(values, concatenate(([0], [[ibin / 10 ** .1, ibin * 10 ** .1] for ibin in bins], [1e5]), axis=None))[0]
            return cumsum(split_bins[where(split_bins > 0)])[:-1]
        return do_pickle(self.make_simple_pickle_path('Splits', sub_dir='Flux'), f, redo=redo or show)

    def get_flux_average(self, values):
        values = values[self.get_fluxes().argsort()]  # sort by ascending fluxes
        return array([mean(lst) for lst in split(values, self.get_flux_splits(show=False))])  # split into sub-lists of similar flux and take average

    def get_times(self, runs=None):
        return self.get_run_values('times', DUTAnalysis.get_time, runs, pbar=False)

    def get_x_var(self, vs_time=False, rel_time=False, rel_error=0., avrg=False):
        return array(self.get_times()) - (time_stamp(self.FirstAnalysis.Run.LogStart) + 3600 if rel_time else 0) if vs_time else array(self.get_fluxes(rel_error, avrg=avrg))

    def get_irradiation(self):
        return self.FirstAnalysis.get_irradiation()

    def get_attenuator(self):
        return self.FirstAnalysis.get_attenuator()

    def get_hv_device(self):
        return self.FirstAnalysis.Currents.Name

    def get_currents(self):
        return OrderedDict((key, ana.Currents.get_current()) for key, ana in self.Analyses.iteritems())

    def get_run_values(self, string, f, runs=None, pbar=None, avrg=False, picklepath=None, *args, **kwargs):
        return self.generate_run_plots(string, f, runs, pbar, avrg, picklepath, *args, **kwargs)

    def get_values(self, string, f, pbar=True, avrg=False, picklepath=None, *args, **kwargs):
        return self.generate_run_plots(string, f, None, pbar, avrg, picklepath, *args, **kwargs)

    def generate_run_plots(self, string, f, runs=None, pbar=None, avrg=False, picklepath=None, *args, **kwargs):
        pbar = not all(file_exists(picklepath.format(run)) for run in self.Runs) if picklepath is not None and pbar is None else pbar
        pbar = True if 'redo' in kwargs and kwargs['redo'] else pbar
        self.info('Generating {} ...'.format(string), prnt=pbar)
        self.PBar.start(self.NRuns if runs is None else len(runs)) if pbar else do_nothing()
        plots = []
        for ana in self.get_analyses(runs):
            plots.append(f(ana, *args, **kwargs))
            self.PBar.update() if pbar else do_nothing()
        return array(self.get_flux_average(array(plots)) if avrg else plots)

    def generate_plots(self, string, f, pbar=True, *args, **kwargs):
        return self.generate_run_plots(string, f, runs=None, pbar=pbar, *args, **kwargs)

    @staticmethod
    def get_mode(vs_time):
        return 'Time' if vs_time else 'Flux'

    def get_sm_std_devs(self, redo=False):
        pickle_path = self.make_pickle_path('Uniformity', 'STD', self.RunPlan, self.DUT.Number)

        def f():
            self.info('Getting STD of Signal Map ... ')
            return OrderedDict((key, ana.get_sm_std(redo=redo)) for key, ana in self.Analyses.iteritems())

        return do_pickle(pickle_path, f, redo=redo)

    def get_sm_std(self, redo=False, low=False, high=False):
        pickle_path = self.make_pickle_path('Uniformity', 'SMSTD', self.RunPlan, self.DUT.Number, suf='{}{}'.format(int(low), int(high)))
        runs = self.get_runs_below_flux(110) if low else self.get_runs_above_flux(2000) if high else self.Runs

        def f():
            return mean_sigma([v for run, v in self.get_sm_std_devs(redo).iteritems() if run in runs]) if runs else make_ufloat((0, 0))

        return do_pickle(pickle_path, f, redo=redo)

    def get_pulse_heights(self, bin_width=None, redo=False, runs=None, corr=True, err=True, pbar=None, avrg=False, peaks=False):
        error = self.get_repr_error(110, peaks, redo) if err else 0
        picklepath = self.make_simple_pickle_path('Fit', '{}_eventwise_AllCuts'.format(self.Bins.BinSize), run='{}', sub_dir='Ph_fit')
        pbar = False if peaks else pbar
        return self.get_run_values('pulse heights', self.Analysis.get_pulse_height, runs, pbar, avrg, picklepath, bin_size=bin_width, redo=redo, corr=corr, sys_err=error, peaks=peaks)

    def get_pulse_height(self):
        return mean_sigma(self.get_pulse_heights())

    def get_rate_dependence(self, redo=False):
        values = self.get_pulse_heights(redo=redo, pbar=False)
        return mean_sigma(values)[1] / mean(values), (max(values) - min(values)) / mean(values)

    def print_rate_dependence(self):
        s1, s2 = self.get_rate_dependence()
        print('Rel STD:    {:2.1f}'.format(s1.n * 100))
        print('Rel Spread: {:2.1f} \\pm {:0.1f}'.format(s2.n * 100, s2.s * 100))

    def get_runs_by_collimator(self, fs11=65, fsh13=.5):
        return [key for key, ana in self.Analyses.iteritems() if ana.Run.RunInfo['fs11'] == fs11 and ana.Run.RunInfo['fs13'] == fsh13]

    def get_runs_below_flux(self, flux):
        return [key for key, ana in self.Analyses.iteritems() if ana.Run.Flux <= flux]

    def get_runs_above_flux(self, flux):
        return [key for key, ana in self.Analyses.iteritems() if ana.Run.Flux >= flux]

    def get_repr_error(self, flux=None, peaks=False, redo=False):
        values = self.draw_signal_spread(peaks, redo, show=False, save=False)
        return self.get_repr_error_old(flux, show=False) if values is None else mean_sigma(values)[1]

    def get_repr_error_old(self, flux, show=True, redo=False):

        pickle_path = self.make_pickle_path('Errors', 'Repr', self.RunPlan, self.DUT.Number, suf=flux)

        def f():
            runs = self.get_runs_below_flux(flux)
            if not runs:
                return 0
            values = self.get_pulse_heights(runs=runs, redo=redo, err=False)
            gr = self.make_tgrapherrors('gr_re', 'Pulse Heights Below {f} kHz/cm^{{2}}'.format(f=flux), x=self.get_fluxes(runs=runs), y=values)
            self.format_statbox(entries=2, only_fit=True)
            gr.Fit('pol0', 'qs{s}'.format(s='' if show else '0'))
            format_histo(gr, x_tit='Flux [kHz/cm^{2}]', y_tit='Mean Pulse Height [au]', y_off=1.7)
            self.save_histo(gr, 'ReprErrors', show, draw_opt='ap', lm=.14, prnt=show)
            if len(values) == 1:
                return .01  # take 1% if there is only one measurement below the given flux
            return mean_sigma(values)[1]

        return do_pickle(pickle_path, f, redo=redo)

    def get_uniformities(self, bins=10, redo=False, low_flux=False, high_flux=False):
        runs = self.get_runs_below_flux(110) if low_flux else self.get_runs_above_flux(2000) if high_flux else self.Runs
        values = [array(self.Analyses[run].get_uniformity(bins, redo=redo)) for run in runs]
        return array([v if v[0] is not None else full(3, ufloat(0, 0)) for v in values])

    def get_mean_uniformity(self, bins=10, redo=False, low_flux=False, high_flux=False):
        values = self.get_uniformities(bins, redo, low_flux, high_flux)
        return [mean_sigma(values[:, i][where(values[:, i] > 0)[0]]) for i in xrange(values[0].size)]

    def get_additional_peak_heights(self):
        return self.get_values('Peak Heights', f=self.Analysis.get_additional_peak_height)

    def get_peak_flux(self):
        return self.get_values('Peak Flux', f=self.Analysis.get_peak_flux)

    @staticmethod
    def get_x_tit(vs_time):
        return '{mod}{unit}'.format(mod='Time' if vs_time else 'Flux', unit=' [hh:mm]' if vs_time else ' [kHz/cm^{2}]')

    def get_tax_off(self, vs_time, rel_time=False):
        return None if not vs_time else self.StartTime if rel_time else 0

    def get_xrange(self, vs_time, x_range=None):
        return x_range if vs_time else self.Bins.FluxRange

    def get_x_args(self, vs_time, rel_time=False, x_range=None):
        return {'x_tit': self.get_x_tit(vs_time), 't_ax_off': self.get_tax_off(vs_time, rel_time), 'x_range': self.get_xrange(vs_time, x_range), 'x_off': None if vs_time else 1.2}

    def get_cmd_strings(self, cmd, kwargs):
        return '?'.join(['python analyse.py {} {} -tc {} -d -cmd {} -kw {}'.format(run, self.DUT.Number, self.TCString, cmd, kwargs) for run in self.Runs])
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def draw_full_pulse_height(self, bin_width=10000, show=True, rel_t=True, redo=False, with_flux=True):

        pickle_path = self.make_pickle_path('PulseHeight', 'FullPH', self.RunPlan, ch=self.DUT.Number, suf=bin_width)

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

    def get_pulse_height_graph(self, bin_width=None, vs_time=False, first_last=True, redo=False, legend=True, corr=True, err=True, avrg=False, peaks=False):

        marker_size = 1
        gStyle.SetEndErrorSize(4)
        x = self.get_x_var(vs_time, avrg=avrg)
        g = self.make_tgrapherrors('g', 'stat. error', self.Colors[0], marker_size=marker_size, x=x, y=self.get_pulse_heights(bin_width, redo, corr=corr, err=False, avrg=avrg, peaks=peaks))
        values = self.get_pulse_heights(bin_width, redo, corr=corr, err=err, avrg=avrg, peaks=peaks)
        g_errors = self.make_tgrapherrors('gerr', 'full error', marker=0, color=self.Colors[5], marker_size=0, x=x, y=values)
        g1, g_last = [self.make_tgrapherrors('g{}'.format(i), '{} run'.format('last' if i else 'first'), marker=22 - i, color=2, marker_size=1.5, x=[x[i].n], y=[values[i].n]) for i in [0, -1]]
        graphs = [g, g_errors]
        first_last = first_last and not avrg
        graphs += [g1, g_last] if first_last else []
        leg = self.make_legend(x2=.37, y1=.21, nentries=len(graphs), w=.2)
        mg = TMultiGraph('mg_ph', 'Pulse Height vs {mod} - {dia}'.format(mod='Time' if vs_time else 'Flux', dia=self.DUT.Name))
        for i, gr in enumerate(graphs):
            leg.AddEntry(gr, gr.GetTitle(), 'l' if gr.GetName() in ['gerr', 'g'] else 'p')
            mg.Add(gr, 'p')
        if legend:
            mg.GetListOfFunctions().Add(leg)
        if vs_time:
            g = mg.GetListOfGraphs()[0]
            for i, (ana, ix) in enumerate(zip(self.Analyses.itervalues(), x)):
                y, ey = g.GetY()[i], g.GetErrorY(i)
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(self.draw_tlatex(ix.n, y + ey * 1.2, '{:1.0f}'.format(ana.get_flux().n), color=1, align=21, size=.02, show=0))
        return mg

    def draw_scaled_pulse_heights(self, scale=1, binning=None, vs_time=False, show=True, y_range=None, redo=False, scale_to_low=False, avrg=False, peaks=False):

        mg = self.get_pulse_height_graph(binning, vs_time, first_last=not vs_time, redo=redo, avrg=avrg, peaks=peaks)
        scale_multigraph(mg, scale, scale_to_low)
        y_range = [.95, 1.05] if y_range is None else y_range
        format_histo(mg, y_tit='Scaled Pulse Height', y_off=1.75, draw_first=True, y_range=y_range, ndivx=503, center_y=True, **self.get_x_args(vs_time))
        self.draw_histo(mg, show, lm=.14, draw_opt='a', logx=not vs_time, grid=vs_time, gridy=True, bm=.18)
        self.draw_irradiation(make_irr_string(self.RunSelection.get_irradiation()))
        self.save_plots('ScaledPulseHeights{}'.format('Time' if vs_time else 'Flux'))
        return mg.GetListOfGraphs()[0]

    def draw_pulse_heights(self, bin_width=None, vs_time=False, show=True, show_first_last=True, save_comb=False, y_range=None, corr=True, redo=False, prnt=True, err=True, avrg=False,
                           fit=False, peaks=False):

        mg = self.get_pulse_height_graph(bin_width, vs_time, show_first_last, redo, corr=corr, err=err, avrg=avrg, peaks=peaks)

        # small range
        ymin, ymax = [getattr(mg.GetListOfGraphs()[0].GetYaxis(), 'GetX{}'.format(w))() for w in ['min', 'max']]
        y_range = increased_range([ymin, ymax], .5, .15) if y_range is None else y_range
        format_histo(mg, color=None, y_tit='Signal Pulse Height [mV]', y_off=1.75, draw_first=True, y_range=y_range, **self.get_x_args(vs_time, self.Bins.FluxRange))
        self.save_histo(mg, 'PulseHeight{mod}'.format(mod=self.get_mode(vs_time)), show=False, lm=.14, draw_opt='A', logx=not vs_time, grid=vs_time)

        # no zero suppression
        mg1 = mg.Clone()
        leg = deepcopy(mg1.GetListOfFunctions()[0])
        mg1.GetListOfFunctions().Clear()
        format_histo(mg1, 'mg1_ph', draw_first=True, y_range=[0, ymax * 1.1] if y_range is None else y_range)
        if fit:
            self.format_statbox(only_fit=1, form='2.1f')
            mg1.GetListOfGraphs()[1].Fit('pol0', 'qs')
        self.save_histo(mg1, 'PulseHeightZero{mod}'.format(mod=self.get_mode(vs_time)), not save_comb and show, lm=.14, draw_opt='a', logx=not vs_time, leg=[self.draw_signal_legend(), leg])

        if save_comb:
            y_min = increased_range([ymin, ymax], .5)[0] if y_range is None else y_range[0]
            # TODO fix vs time and comb plot
            self.save_combined_pulse_heights(mg, mg1, y_min, show=save_comb, pulser_leg=self.draw_signal_legend, prnt=prnt)

        return mg

    def show_models(self, a=.001, b=50, c=2000, d=400, e=10, f=20):
        def f1(x, pars):
            v = pars[0] * x[0] + pars[1]
            return v if x[0] < pars[2] else pars[0] * pars[2] + pars[1]

        def f2(x, pars):
            return exp(-x[0] / pars[0] * log(2)) * pars[1] + pars[2]

        def f3(x, pars):
            v = pars[0] * x[0] + pars[1]
            return (v if x[0] < pars[2] else pars[0] * pars[2] + pars[1]) + exp(-x[0] / pars[3] * log(2)) * pars[4] + pars[5]

        f1 = TF1('f1', f1, 0.1, 1e5, 3)
        f2 = TF1('f2', f2, 0.1, 1e5, 3)
        f3 = TF1('f3', f3, 0.1, 1e5, 6)
        f1.SetParameters(a, b, c)
        f2.SetParameters(d, e, f)
        f3.SetParameters(a, b, c, d, e, f)
        self.draw_histo(f1, logx=True)
        format_histo(f1, y_range=[0, 2 * b], lw=2, x_range=[.1, 1e5])
        f2.Draw('same')
        f3.SetLineStyle(7)
        f3.Draw('same')
        self.add(f2, f3)

    def fit_pulse_height(self, logx=True):
        values = self.get_pulse_heights(pbar=False)
        g = self.make_tgrapherrors('gfph', 'Pulse Height Fit', x=self.get_fluxes(), y=values)

        def f1(x, pars):
            v = pars[0] * x[0] + pars[1]
            # return (v if v < pars[2] else pars[2]) + exp(-x[0] * log(2) / pars[3]) + pars[4]
            return (v if x[0] < pars[2] else pars[0] * pars[2] + pars[1]) - pars[3] * TMath.Erf(pars[4] * (x[0] - pars[5]))

        f = TF1('f1', f1, .1, 1e6, 6)
        # f = TF1('fph', '[0] - [1] * exp(-x/[2]*log(2)) + [3] * exp(-x/[4]*log(2))', .1, 1e6)
        f.SetParNames('m', 'c_{1}', 'c_{max}', 'scale', 'width', 'offset')
        f.SetParLimits(0, 0.0001, 0.1)  # slope of the linear part
        f.SetParLimits(1, 100, 130)  # offset of the linear part
        f.SetParLimits(2, 10, 3000)  # constant of the linear part
        f.SetParLimits(3, 1, 50)  # half life for the exp
        f.SetParLimits(4, 1e-6, 1e-3)  # asymptote for of the exp
        f.SetParLimits(5, 100, 5e3)  # asymptote for of the exp
        set_root_output(False)
        g.Fit(f, 'q')
        format_histo(g, x_tit='Flux [kHz/cm^{2}]', y_tit='Pulse Height [mV]', y_off=1.4, x_range=self.Bins.FluxRange if logx else [0, 10100])
        self.format_statbox(only_fit=True, x=.5, y=.45 if self.Title else .38)
        self.draw_histo(g, logx=logx, lm=.12)
        l1 = self.draw_horizontal_line(f.GetParameter(0), .1, 1e6, color=2, style=7, w=2)
        l2 = self.draw_horizontal_line(f.GetParameter(0) - f.GetParameter(1) + f.GetParameter(3), .1, 1e6, color=4, style=7, w=2, name='as')
        leg = self.make_legend(w=.2)
        leg.AddEntry(l1, 'c', 'l')
        leg.AddEntry(l2, 'c - c_{1} + c_{2}', 'l')
        leg.Draw()
        f.SetLineStyle(7)
        f.Draw('same')
        g.GetListOfFunctions()[-1].Draw()  # redraw stats
        return f

    def draw_signal_legend(self):
        pass

    def draw_signal_distributions(self, show=True, off=3, redo=False):

        stack = THStack('hsd', 'Pulse Height Distributions')
        leg = self.make_legend(.67, .88, nentries=self.NRuns)
        self.info('Generating signal distributions ...')
        histos = []
        self.PBar.start(self.NRuns)
        for i, ana in enumerate(self.Analyses.itervalues()):
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

        pickle_path = self.make_pickle_path('ph_fit', 'PhDistoFits', run=self.RunPlan, ch=self.DUT.Number, suf=self.Bins.BinSize if bin_width is None else bin_width)

        def f():
            collimator_settings = set([(ana.Run.RunInfo['fs11'], ana.Run.RunInfo['fsh13']) for key, ana in self.Analyses.iteritems()])
            return {s: self.draw_ph_distributions(bin_width, show=show, fs11=s[0], fsh13=s[1]) for s in collimator_settings}

        return do_pickle(pickle_path, f)

    def draw_ph_distributions(self, binning=None, fsh13=.5, fs11=65, show=True):
        runs = self.get_runs_by_collimator(fsh13=fsh13, fs11=fs11)
        return self.draw_combined_ph_distributions(runs, binning, show)

    def draw_ph_distributions_below_flux(self, bin_width=None, flux=150, show=True):
        pickle_path = self.make_pickle_path('Ph_fit', 'PhDistoBel', self.RunPlan, self.DUT.Number, suf='{bin}_{flux}'.format(bin=self.Bins(bin_width), flux=flux))

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

    def draw_signal_spread(self, peaks=False, redo=False, show=True, save=True):
        values = self.get_pulse_heights(redo=redo, err=False, peaks=peaks)[self.get_fluxes().argsort()]  # sort by ascending fluxes
        rel_values = array([value - mean(lst) for lst in split(values, self.get_flux_splits(show=False)) for value in lst if lst.size > 1])
        if rel_values.size < 2:
            warning('Not enough data for signal spread ...')
            return
        if show:
            h = TH1F('hps', 'Relative Signal Spread', 40, -2, 2)
            h.FillN(rel_values.size, array([v.n for v in rel_values], 'd'), ones(rel_values.size))
            self.format_statbox(all_stat=True)
            format_histo(h, x_tit='Relative Signal', y_tit='Number of Entries', y_off=1.2)
            self.save_histo(h, 'SignalSpread', lm=.11, show=show, save=save)
        return rel_values
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region MISCELLANEOUS DRAW
    def draw_currents(self, v_range=None, rel_time=False, averaging=1, with_flux=False, c_range=None, f_range=None, draw_opt='al', show=True):
        self.Currents.draw_indep_graphs(rel_time=rel_time, v_range=v_range, averaging=averaging, with_flux=with_flux, c_range=c_range, f_range=f_range, show=show, draw_opt=draw_opt)

    def draw_iv(self, show=True):
        g = self.Currents.draw_iv(show=False)
        self.save_histo(g, 'IV', draw_opt='ap', logy=True, lm=.12, show=show)

    def draw_fwhm(self, arg=1, vs_time=False, show=True, redo=False):
        y_values = self.get_uniformities(redo=redo)[:, arg]
        x_values = self.get_times() if vs_time else self.get_fluxes()
        g = self.make_tgrapherrors('gfm', 'Full Width Half Maximum', x=x_values, y=y_values)
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
        g = self.make_tgrapherrors('gstd', 'STD of the Signal Map', x=x_values, y=y_values)
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

    def draw_flux(self, bin_width=5, rel_time=True, show=True, redo=False):

        pickle_path = self.make_pickle_path('Flux', 'FullFlux', self.RunPlan, self.DUT.Number, bin_width)

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

        hist = do_pickle(pickle_path, f, redo=redo)
        format_histo(hist, x_tit='Time [hh:mm]', y_tit='Flux [kHz/cm^{2}]', t_ax_off=self.InitTime if rel_time else 0, fill_color=self.FillColor, y_range=[1, 20000], stats=0)
        self.save_histo(hist, 'FluxEvo', x=1.5, y=.75, show=show, logy=True, draw_opt='bar' if not self.FirstAnalysis.has_branch('rate') else '')
        return hist

    def draw_current_flux(self, c_range=None, fit=True, show=True):
        currents = self.get_currents().values()
        g = self.make_tgrapherrors('gcf', 'Leakage Current vs. Flux', x=self.get_fluxes(rel_error=.01), y=currents)
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
        f = DUTAnalysis.draw_hitmap if hitmap else DUTAnalysis.draw_signal_map
        histos = self.generate_plots('{} maps'.format('hit' if hitmap else 'signal'), f, show=False, prnt=False, redo=redo)
        glob_max = (int(max([h.GetMaximum() for h in histos])) + 5) / 5 * 5
        glob_min = int(min([h.GetMinimum() for h in histos])) / 5 * 5
        for i, h in enumerate(histos):
            format_histo(h, z_range=[glob_min, glob_max]) if not hitmap else do_nothing()
            self.save_histo(h, '{n}Map{nr}'.format(nr=str(i).zfill(2), n='Hit' if hitmap else 'Signal'), show=False, ind=i, draw_opt='colz', rm=.16, lm=.12, prnt=False)

    def draw_hitmaps(self, redo=False):
        self.draw_signal_maps(hitmap=True, redo=redo)

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
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
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
        cut_val = self.FirstAnalysis.Cut.CutConfig['chi2_{}'.format(mod_str if mod_str else 'X')]
        cuts = zeros((self.NRuns, 2))
        for i, h in enumerate(histos):
            h.GetQuantiles(2, cuts[i], array([cut_val / 100., .93]))
        legend = self.make_legend(nentries=self.NRuns + 1, scale=.7, w=.2)
        stack = THStack('hx2', '#chi^{{2}}{}'.format(' in {}'.format(mod_str) if mod_str else ''))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
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
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
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
    p.add_argument('-p', '--prnt', action='store_true')
    p.add_argument('-cmd', '--command', nargs='?', help='method to be executed')
    p.add_argument('-kw', '--kwargs', nargs='?', help='key word arguments as dict {"show": 1}', default='{}')
    pargs = p.parse_args()

    z = AnalysisCollection(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree and not pargs.prnt, pargs.verbose)
    z.print_loaded()
    if pargs.runs:
        z.Currents.draw_indep_graphs()
        raw_input('Press any button to exit')
    if pargs.draw:
        z.draw_all(pargs.redo)
    if pargs.prnt:
        print(z.get_cmd_strings(pargs.command, pargs.kwargs))
