#! /usr/bin/env python
from ROOT import THStack, TF1, TMath
from numpy import sort, log, argmin, argmax

from src.analysis import Analysis, glob
from src.currents import Currents
from src.dut_analysis import DUTAnalysis, Bins
from src.run_selection import RunPlan
from helpers.draw import *


class AnalysisCollection(Analysis):
    """ Analysis of the various runs of a single runplan. """

    def __init__(self, run_plan, dut_nr, testcampaign=None, load_tree=True, verbose=False):

        # RUN SELECTION
        self.RunSelection = RunPlan(testcampaign, run_plan, dut_nr, verbose)
        self.Runs = array(self.RunSelection.get_selected_runs())
        self.NRuns = self.Runs.size
        self.RunPlan = self.RunSelection.SelectedRunplan
        self.DUT = self.RunSelection.SelectedDUT
        self.Type = self.RunSelection.SelectedType
        self.Fluxes = self.RunSelection.get_selected_fluxes()

        super(AnalysisCollection, self).__init__(testcampaign, sub_dir=join(self.DUT.Name, 'RP{}'.format(self.RunPlan)), verbose=verbose)
        self.print_start(run_plan, prnt=load_tree, dut=self.DUT.Name)
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
        self.FirstAnalysis = list(self.Analyses.values())[0]
        self.LastAnalysis = list(self.Analyses.values())[-1]
        self.Bins = self.FirstAnalysis.Bins if load_tree else None
        self.StartTime = self.FirstAnalysis.Run.StartTime if self.LoadTree else time_stamp(self.FirstAnalysis.Run.LogStart)

        # Sub Classes
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

    def make_flux_legend(self, histos, runs=None, x1=.76, y2=.88, draw_style='l'):
        runs = choose(runs, self.Runs)
        legend = self.Draw.make_legend(x1, y2, nentries=len(runs))
        for h, run in zip(histos, runs):
            legend.AddEntry(h, make_flux_string(self.get_fluxes()[list(self.Runs).index(run)].n, prec=0), draw_style)
        return legend

    # ----------------------------------------
    # region INIT
    def print_start_info(self):
        bias = ['{:+8.0f}'.format(bias) for bias in self.RunSelection.get_selected_biases()]
        times = [datetime.fromtimestamp(duration - 3600).strftime('%H:%M').rjust(15) for duration in self.RunSelection.get_selected_durations()]
        rows = array([self.Runs, ['{:14.1f}'.format(flux.n) for flux in self.Fluxes], bias, self.RunSelection.get_selected_start_times(), times]).T
        print_table(header=['Run', 'Flux [kHz/cm2]', 'Bias [V]', 'Start', 'Duration [hh:mm]'], rows=rows, prnt=self.Verbose)

    def get_high_low_rate_runs(self):
        return self.Runs[argmin(self.Fluxes)], self.Runs[argmax(self.Fluxes)]

    def generate_common_pickles(self):
        if self.LoadTree:
            self.generate_slope_pickle()

    def generate_slope_pickle(self):
        picklepath = self.make_pickle_path('TrackAngle', 'x', run=self.MinFluxRun)
        if not file_exists(picklepath):
            DUTAnalysis(self.MinFluxRun, self.DUT.Number, self.TCString, prnt=False)

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
                analysis.Cut.set_low_rate_run(low_run=self.MinFluxRun)
                analysis.Cut.reload()
            analyses[run] = analysis
        self.Threads = None
        return analyses

    def close_files(self):
        for ana in list(self.Analyses.values()):
            ana.Run.tree.Delete()
            ana.Run.RootFile.Close()

    def delete_trees(self):
        for ana in list(self.Analyses.values()):
            ana.Tree.Delete()
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region BINNING
    def get_binning(self, bin_width=None, rel_time=False):
        bins = concatenate([ana.Bins.get_raw_time(bin_width, t_from_event=True)[1] for ana in list(self.Analyses.values())])
        return [bins.size - 1, bins - (self.FirstAnalysis.Run.StartTime if rel_time else 0)]

    def get_time_bins(self, bin_width=None, only_edges=False):
        bins = self.fix_t_arrays([ana.Bins.get_time(bin_width)[1] for ana in list(self.Analyses.values())])
        bins = array([[b[0], b[-1]] for b in bins]).flatten() if only_edges else concatenate(bins)
        return [bins.size - 1, bins]

    def get_raw_time_bins(self, bin_width=None, only_edges=False, t_from_event=False):
        bins = self.fix_t_arrays([ana.Bins.get_raw_time(bin_width, t_from_event=t_from_event)[1] for ana in list(self.Analyses.values())])
        bins = array([[b[0], b[-1]] for b in bins]).flatten() if only_edges else concatenate(bins)
        return [bins.size - 1, bins]

    def fix_t_arrays(self, t_arrays):
        """ Add the logged time between two runs if the start time of a run is less than the stop time of the previous. """
        for i in range(len(t_arrays) - 1):
            delta = t_arrays[i + 1][0] - t_arrays[i][-1]  # check the time difference from the first bin of the next array and the last of the current
            if delta < 0:
                t_arrays[i + 1] += self.get_break_time(i) - delta
        return t_arrays

    def get_break_time(self, ind):
        # because the log stop is not so reliable...
        return (self.get_ana(ind + 1).Run.LogStart - self.get_ana(ind).Run.LogStart - self.get_ana(ind).Run.Duration).total_seconds()
    # endregion BINNING
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_runs(self):
        return self.Runs

    def get_ana(self, ind):
        return list(self.Analyses.values())[ind]

    def get_analyses(self, runs=None):
        return list(self.Analyses.values()) if runs is None else [ana for key, ana in list(self.Analyses.items()) if key in runs]

    def get_hv_name(self):
        return self.Currents.Name

    def get_fluxes(self, rel_error=0., corr=True, runs=None, avrg=False, pbar=True):
        picklepath = self.make_simple_pickle_path(sub_dir='Flux', run='{}', dut='')
        pbar = False if not self.FirstAnalysis.has_branch('rate') else pbar
        return self.get_values('fluxes', DUTAnalysis.get_flux, runs, pbar, avrg=avrg, picklepath=picklepath, rel_error=rel_error, corr=corr)

    def get_flux_strings(self, prec=0, runs=None):
        return [make_flux_string(flux.n, prec=prec) for flux in self.get_fluxes(runs=runs)]

    def get_flux_splits(self, redo=False, corr=True, show=True):
        def f():
            values = sort([flux.n for flux in self.get_fluxes(corr=corr)])
            h = TH1F('hmf', 'Fluxes', *log_bins(50, 1, 1e5))
            h.FillN(values.size, values, full(values.size, 1, 'd'))
            format_histo(h, x_tit='Flux [kHz/cm^{2}]', y_tit='Number of Entries', y_off=.5, x_off=1.2, lab_size=.05, tit_size=.05)
            self.Draw(h, lm=.07, logx=True, show=show, h=.75, w=1.5, bm=.14, stats=set_entries())
            s = TSpectrum(20)
            s.Search(h, 1)
            bins = sorted(s.GetPositionX()[i] for i in range(s.GetNPeaks()))
            split_bins = histogram(values, concatenate(([0], [[ibin / 10 ** .1, ibin * 10 ** .1] for ibin in bins], [1e5]), axis=None))[0]
            return cumsum(split_bins[where(split_bins > 0)])[:-1]
        return do_pickle(self.make_simple_pickle_path('Splits', sub_dir='Flux'), f, redo=redo or show)

    def get_flux_average(self, values):
        values = values[self.get_fluxes().argsort()]  # sort by ascending fluxes
        return array([mean(lst) for lst in split(values, self.get_flux_splits(show=False))])  # split into sub-lists of similar flux and take average

    def get_times(self, runs=None):
        return self.get_values('times', DUTAnalysis.get_time, runs, pbar=False)

    def get_events(self, redo=False):
        return self.get_values('events', DUTAnalysis.get_events, picklepath=self.make_simple_hdf5_path('', 'AllCuts', 'Events', '{}'), redo=redo)

    def get_n_events(self, redo=False):
        return array([e.size for e in self.get_events(redo)])

    def get_x_var(self, vs_time=False, rel_error=0., avrg=False):
        return self.get_times() if vs_time else array(self.get_fluxes(rel_error, avrg=avrg))

    def get_irradiation(self):
        return self.FirstAnalysis.get_irradiation()

    def get_attenuator(self):
        return self.FirstAnalysis.get_attenuator()

    def get_hv_device(self):
        return self.FirstAnalysis.Currents.Name

    def get_currents(self):
        return [ana.Currents.get_current() for ana in self.get_analyses()]

    def get_values(self, string, f, runs=None, pbar=None, avrg=False, picklepath=None, flux_sort=False, plots=False, *args, **kwargs):
        runs = choose(runs, self.Runs)
        pbar = choose(pbar, 'redo' in kwargs and kwargs['redo'] or (True if picklepath is None else not all(file_exists(picklepath.format(run)) for run in runs)))
        self.info('Generating {} ...'.format(string), prnt=pbar)
        self.PBar.start(len(runs)) if pbar else do_nothing()
        values = []
        for ana in self.get_analyses(runs):
            values.append(f(ana, *args, **kwargs))
            self.PBar.update() if pbar else do_nothing()
        return values if plots else array(self.get_flux_average(array(values))) if avrg else array(values, dtype=object)[self.get_fluxes().argsort() if flux_sort else ...]

    def get_plots(self, string, f, runs=None, pbar=None, avrg=False, picklepath=None, *args, **kwargs):
        return self.get_values(string, f, runs, pbar, avrg, picklepath, False, True, *args, **kwargs)

    @staticmethod
    def get_mode(vs_time):
        return 'Time' if vs_time else 'Flux'

    def get_sm_std_devs(self, redo=False):
        pickle_path = self.make_pickle_path('Uniformity', 'STD', self.RunPlan, self.DUT.Number)

        def f():
            self.info('Getting STD of Signal Map ... ')
            return OrderedDict((key, ana.get_sm_std(redo=redo)) for key, ana in list(self.Analyses.items()))

        return do_pickle(pickle_path, f, redo=redo)

    def get_sm_std(self, redo=False, low=False, high=False):
        pickle_path = self.make_pickle_path('Uniformity', 'SMSTD', self.RunPlan, self.DUT.Number, suf='{}{}'.format(int(low), int(high)))
        runs = self.get_runs_below_flux(110) if low else self.get_runs_above_flux(2000) if high else self.Runs

        def f():
            return mean_sigma([v for run, v in list(self.get_sm_std_devs(redo).items()) if run in runs]) if runs else ufloat(0, 0)

        return do_pickle(pickle_path, f, redo=redo)

    def get_pulse_heights(self, *args, **kwargs):
        return []

    def get_pulse_height(self):
        return mean_sigma(self.get_pulse_heights())

    def get_efficiencies(self, suf='3', redo=False):
        return self.get_values('efficiencies', self.Analysis.get_efficiency, picklepath=self.make_simple_pickle_path(suf=suf, sub_dir='Efficiency', run='{}'), redo=redo)

    def get_rate_dependence(self, redo=False, values=None):
        values = choose(values, self.get_pulse_heights(redo=redo, pbar=False))
        return mean_sigma(values)[1] / mean(values), (max(values) - min(values)) / mean(values)

    def print_rate_dependence(self, values=None):
        s1, s2 = self.get_rate_dependence(values=values)
        print('Rel STD:    {:2.1f}'.format(s1.n * 100))
        print('Rel Spread: {:2.1f} \\pm {:0.1f}'.format(s2.n * 100, s2.s * 100))

    def get_runs_below_flux(self, flux):
        return [key for key, ana in list(self.Analyses.items()) if ana.Run.Flux <= flux]

    def get_runs_above_flux(self, flux):
        return [key for key, ana in list(self.Analyses.items()) if ana.Run.Flux >= flux]

    def get_signal_spread(self, peaks=False, redo=False, rel=False):
        groups = split(self.get_pulse_heights(redo=redo, err=False, peaks=peaks, flux_sort=True), self.get_flux_splits(show=False))
        values = array([(value - mean_sigma(grp)[0]) / (mean_sigma(grp)[0] if rel else 1) for grp in groups for value in grp if grp.size > 1])
        return values if values.size > 1 else None

    def get_repr_error(self, flux=None, peaks=False, redo=False):
        values = self.get_signal_spread(peaks, redo)
        return self.get_repr_error_old(flux, show=False) if values is None else mean_sigma(values, err=False)[1]

    def get_repr_error_old(self, flux, show=True, redo=False):

        pickle_path = self.make_pickle_path('Errors', 'Repr', self.RunPlan, self.DUT.Number, suf=flux)

        def f():
            runs = self.get_runs_below_flux(flux)
            if not runs:
                return 0
            values = self.get_pulse_heights(runs=runs, redo=redo, err=False)
            gr = Draw.make_tgrapherrors(self.get_fluxes(runs=runs), values, title='Pulse Heights Below {f} kHz/cm^{{2}}'.format(f=flux))
            gr.Fit('pol0', 'qs{s}'.format(s='' if show else '0'))
            format_histo(gr, x_tit='Flux [kHz/cm^{2}]', y_tit='Mean Pulse Height [au]', y_off=1.7)
            self.Draw(gr, 'ReprErrors', show, draw_opt='ap', lm=.14, prnt=show, stats=set_statbox(fit=True))
            if len(values) == 1:
                return .01  # take 1% if there is only one measurement below the given flux
            return mean_sigma(values, err=False)[1]

        return do_pickle(pickle_path, f, redo=redo)

    def get_uniformities(self, use_fwc=True, low_flux=False, high_flux=False, redo=False):
        runs = self.get_runs_below_flux(110) if low_flux else self.get_runs_above_flux(2000) if high_flux else self.Runs
        picklepath = self.make_simple_pickle_path('Uniformity', int(use_fwc), 'Signal', run='{}')
        return self.get_values('Getting uniformities', self.Analysis.get_uniformity, runs, picklepath=picklepath, redo=redo, use_fwc=use_fwc)

    def get_mean_uniformity(self, use_fwc=True, redo=False, low_flux=False, high_flux=False):
        values = self.get_uniformities(use_fwc, redo, low_flux, high_flux)
        return [mean_sigma(values[:, i][where(values[:, i] > 0)[0]]) for i in range(values[0].size)]

    @staticmethod
    def get_x_tit(vs_time):
        return 'Time [hh:mm]' if vs_time else 'Flux [kHz/cm^{2}]'

    @staticmethod
    def get_x_draw(vs_time=False):
        return {'logx': not vs_time, 'grid': vs_time}

    def get_tax_off(self, vs_time, rel_time=False):
        return None if not vs_time else self.StartTime if rel_time else 0

    @staticmethod
    def get_range(vs_time, x_range=None):
        return x_range if vs_time else Bins.FluxRange

    def get_x_args(self, vs_time=False, rel_time=False, x_range=None, draw_args=False):
        hist_kwargs = {'x_tit': self.get_x_tit(vs_time), 't_ax_off': self.get_tax_off(vs_time, rel_time), 'x_range': self.get_range(vs_time, x_range), 'x_off': None if vs_time else 1.2}
        return {**hist_kwargs, **self.get_x_draw(vs_time)} if draw_args else hist_kwargs

    def get_cmd_strings(self, cmd, kwargs):
        return '?'.join(['python analyse.py {} {} -tc {} -d -cmd {} -kw {}'.format(run, self.DUT.Number, self.TCString, cmd, kwargs) for run in self.Runs])
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    @staticmethod
    def draw_legend(graphs, x=.17):
        return Draw.legend(graphs, [g.GetTitle() for g in graphs], ['l' if i < 2 else 'p' for i in range(len(graphs))], x1=x, y1=.21, w=.2)

    def make_pulse_height_graph(self, bin_width=None, vs_time=False, first_last=True, redo=False, legend=True, corr=True, err=True, avrg=False, peaks=False):
        x, (ph0, ph) = self.get_x_var(vs_time, avrg=avrg), [self.get_pulse_heights(bin_width, redo, corr=corr, err=i, avrg=avrg, peaks=peaks) for i in [False, err]]
        g = Draw.make_tgrapherrors(x, ph0, title='stat. error', color=Draw.color(2, 1), markersize=1, lw=2)
        g_errors = Draw.make_tgrapherrors(x, ph, title='full error', marker=0, color=Draw.color(2, 0), markersize=0, lw=2)
        g1, g_last = [Draw.make_tgrapherrors([x[i].n], [ph[i].n], title='{} run'.format('last' if i else 'first'), marker=22 - i, color=2, markersize=1.5) for i in [0, -1]]
        graphs = [g_errors, g] + ([g1, g_last] if first_last and not avrg else [])
        mg = self.Draw.multigraph(graphs, 'Pulse Height', color=False, show=False)
        if legend:
            mg.GetListOfFunctions().Add(self.draw_legend(graphs))
        if vs_time:
            for ix, iph, flux in zip(x, ph0, self.get_fluxes()):
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(Draw.tlatex(ix.n, iph.n + iph.s * 1.2, '{:0.0f}'.format(flux.n), color=1, align=21, size=.02))
        return mg

    def draw_low_scale(self, avrg=True, yoff=.07):
        x, y = self.get_fluxes(avrg=avrg, corr=False), self.get_pulse_heights(err=False, avrg=avrg)
        y /= y[argmin(x)].n
        self.Draw.graph(x, y, y_tit='Scaled Pulse Height [mV]', y_range=[1 - yoff, 1 + yoff], **self.get_x_args(draw_args=True), lm=.12)
        self.print_rate_dependence(y)

    def draw_mean_scale(self, avrg=True, yoff=.07):
        x, y = self.get_fluxes(avrg=avrg, corr=False), self.get_pulse_heights(err=False, avrg=avrg)
        y /= mean(y).n
        self.Draw.graph(x, y, y_tit='Scaled Pulse Height [mV]', y_range=[1 - yoff, 1 + yoff], **self.get_x_args(draw_args=True), lm=.12)
        self.print_rate_dependence(y)

    def draw_mid_mean_scale(self, avrg=True, yoff=.07):
        x, y = self.get_fluxes(avrg=avrg, corr=False), self.get_pulse_heights(err=False, avrg=avrg)
        y /= mean(y[1:-1]).n
        self.Draw.graph(x, y, y_tit='Scaled Pulse Height [mV]', y_range=[1 - yoff, 1 + yoff], **self.get_x_args(draw_args=True), lm=.12)
        self.print_rate_dependence(y)

    def draw_scaled_pulse_heights(self, scale=1, binning=None, vs_time=False, show=True, y_range=None, redo=False, scale_to_low=False, avrg=False, peaks=False):
        """ Shows the scaled pulse heights of the single runs. """
        mg = self.make_pulse_height_graph(binning, vs_time, first_last=not vs_time, redo=redo, avrg=avrg, peaks=peaks)
        scale_multigraph(mg, scale, scale_to_low)
        self.Draw(mg, show=show, lm=.14, draw_opt='ap', **self.get_x_draw(vs_time), gridy=True, bm=.18)
        Draw.irradiation(make_irr_string(self.RunSelection.get_irradiation()))
        format_histo(mg, y_tit='Scaled Pulse Height', y_off=1.75, y_range=choose(y_range, [.95, 1.05]), ndivx=503, center_y=True, **self.get_x_args(vs_time))
        self.Draw.save_plots('ScaledPulseHeights{}'.format('Time' if vs_time else 'Flux'))
        return mg.GetListOfGraphs()[0]

    def draw_pulse_heights(self, bin_width=None, vs_time=False, show=True, show_first_last=True, y_range=None, corr=True, redo=False, err=True, avrg=False, fit=False, peaks=False):
        """ Shows the pulse heights of the runs. """
        mg = self.make_pulse_height_graph(bin_width, vs_time, show_first_last, redo, corr=corr, err=err, avrg=avrg, peaks=peaks)
        y_range = choose(y_range, ax_range(get_graph_y(mg.GetListOfGraphs()[0], err=False), fl=.5, fh=.2, rnd=True))
        format_histo(mg, color=None, y_tit='Signal Pulse Height [mV]', y_off=1.75, draw_first=True, y_range=y_range, **self.get_x_args(vs_time))
        if fit:
            mg.GetListOfGraphs()[0].Fit('pol0', 'qs')
        self.Draw(mg, 'PulseHeight{mod}'.format(mod=self.get_mode(vs_time)), show=show, lm=.14, draw_opt='ap', **self.get_x_draw(vs_time), stats=set_statbox(fit=fit, form='2.1f', stats=fit))
        return mg

    def draw_full_pulse_height(self, bin_width=10000, show=True, rel_t=True, redo=False, with_flux=True):
        """ Shows the pulse heights bins of all runs vs time. """
        def f():
            self.info('Getting pulse heights ...')
            # add a zero at the end for the time in between runs
            ph = concatenate([append(get_hist_vec(ana.draw_pulse_height(bin_width, corr=True, redo=redo, show=False, save=False)[0]), ufloat(0, 0)) for ana in self.get_analyses()])
            h0 = Draw.make_histo('Full Pulse Height', self.get_time_bins(bin_width))
            for i, v in enumerate(ph, 1):
                h0.SetBinContent(i, v.n)
                h0.SetBinError(i, v.s)
            return format_histo(h0, title='Full Pulse Height', x_tit='Time [hh:mm]', y_tit='Mean Pulse Height [mV]', fill_color=Draw.FillColor, y_off=.8)

        h1 = do_pickle(self.make_simple_pickle_path('Full', bin_width, 'PulseHeight'), f, redo=redo)
        format_histo(h1, y_range=[0, h1.GetMaximum() * 1.05], t_ax_off=self.get_tax_off(True, rel_t), stats=0)
        c = self.Draw(h1, show=show, draw_opt='hist', w=1.5, h=.75, lm=.065, gridy=True, rm=.1 if with_flux else None)
        if with_flux:
            h2 = self.draw_flux(rel_time=rel_t, show=False)
            c.cd()
            Draw.tpad('pr', margins=[.065, .1, c.GetBottomMargin(), .1], transparent=True, logy=True)
            x_range = [h2.GetXaxis().GetXmin(), h2.GetXaxis().GetXmax()]
            format_histo(h2, title=' ', fill_color=2, fill_style=3002, lw=1, y_range=[1, h2.GetMaximum() * 1.2], stats=0, y_off=1.05, x_range=x_range)
            self.Draw(h2, canvas=c, draw_opt='histy+same', rm=.1)
        self.Draw.save_plots('FullPulseHeight')
        return h1

    def draw_splits(self, m=2, show=True, normalise=False):
        x, y = self.get_x_var(), self.get_values('split pulse heights', DUTAnalysis.get_split_ph, m=m).T
        y /= mean(y, axis=1).reshape(m ** 2, 1) if normalise else 1
        graphs = [self.Draw.graph(x, iy, y_tit='Pulse Height [mV]', **self.get_x_args(), show=False) for iy in y]
        mg = self.Draw.multigraph(graphs, 'Region Pulse Heights', ['{}'.format(i) for i in range(y.shape[0])], show=show, draw_opt='alp', logx=True)
        format_histo(mg, **self.get_x_args())

    def draw_signal_legend(self):
        pass

    def draw_signal_distributions(self, bin_width=None, redo=False, logy=False, show=True):
        """Shows a stack of the signal distributions."""
        stack = THStack('hsd', 'Pulse Height Distributions')
        histos = self.get_values('signal distributions', self.Analysis.draw_signal_distribution, show=False, prnt=False, redo=redo, bin_width=bin_width)
        for i, h in enumerate(histos):
            format_histo(h, lw=2, color=self.Draw.get_color(len(histos)), fill_color=0, fill_style=4000, stats=0)
            h.Scale(1 / h.GetMaximum())
            stack.Add(h, 'hist')
        format_histo(stack, y_off=1.55, draw_first=True, x_tit='Pulse Height [au]', y_tit='Number of Entries')
        self.Draw(stack, 'SignalDistributions{}'.format('Log' if logy else ''), show=show, lm=.13, draw_opt='nostack', logy=logy, leg=self.make_flux_legend(histos))

    def draw_pulls_below_flux(self, bin_width=.5, flux=150, show=True):
        suffix = '{}_{}'.format(bin_width, flux)
        return do_pickle(self.make_simple_pickle_path('PhPullBel', suffix, 'Ph_fit'), self.draw_pulls, redo=show, runs=self.get_runs_below_flux(flux), bin_width=bin_width, show=show)

    def draw_pulls(self, runs=None, bin_width=.5, show=True):
        s = THStack('s_phd', 'Pulse Height Distributions')
        runs = choose(runs, self.Runs)
        histos = self.get_values('Generate ph pulls ...', self.Analysis.draw_ph_pull, runs, show=False, bin_width=bin_width, fit=False, save=False)
        for i, h in enumerate(histos):
            format_histo(h, fill_color=4000, stats=0, line_color=self.Draw.get_color(len(runs)))
            s.Add(h)
        self.Draw(s, show=show, draw_opt='', leg=self.make_flux_legend(histos, runs))
        ls = s.GetStack().Last()
        format_histo(s, 'Fit Results', stats=True, x_range=ax_range(0, 0, .5, .5, ls))
        self.Draw.save_plots('PulseHeightDistributions')
        return ls.GetMean(), ls.GetStdDev()

    def draw_signal_spread(self, rel=False, peaks=False, redo=False, show=True, save=True):
        values = self.get_signal_spread(peaks, redo, rel)
        if values is None:
            warning('Not enough data for signal spread ...')
            return
        t = 'Relative ' if rel else ''
        values *= 100 if rel else 1
        bins = make_bins(*ax_range(values, fl=.5, fh=.5), n=3 * sqrt(values.size)) if rel else [40, -2, 2]
        h = self.Draw.distribution(values, bins, '{}Signal Spread'.format(t), x_tit='{}Diff [{}]'.format(t, '%' if rel else 'mV'), y_tit='Number of Entries', show=False)
        self.Draw(h, 'SignalSpread', lm=.11, show=show, save=save)
        return values

    def draw_fwhm(self, arg=1, vs_time=False, use_fcw=True, show=True, redo=False):
        x, y = self.get_x_var(vs_time), self.get_uniformities(use_fcw, redo)[:, arg]
        ph_var = 'FWC' if use_fcw else 'Mean'
        var, unit, tit = [ph_var, 'FWHM', 'FWHM/{}'.format(ph_var)][arg], ' [mV]' if arg < 2 else '', [ph_var, 'Full Width Half Maximum', 'Uniformity'][arg]
        g = Draw.make_tgrapherrors(x, y, title=tit, y_tit='{}{}'.format(var, unit), y_off=1.3, **self.get_x_args(vs_time))
        self.Draw(g, 'FWHM', show, lm=.12, logx=not vs_time)
        return g

    def draw_means(self, vs_time=False, show=True, redo=False):
        return self.draw_fwhm(0, vs_time, False, show, redo)

    def draw_fcw(self, vs_time=False, show=True, redo=False):
        return self.draw_fwhm(0, vs_time, True, show, redo)

    def draw_uniformity(self, vs_time=False, show=True, redo=False):
        return self.draw_fwhm(2, vs_time, show=show, redo=redo)

    def draw_ph_slope(self, vs_time=False, show=True):
        y = [fit2u(ana.draw_pulse_height(show=False)[0].Fit('pol1', 'qs0'), par=1) * 60 for ana in self.get_analyses()]
        self.Draw.graph(self.get_x_var(vs_time), y, 'Pulse Height Slope', **self.get_x_args(vs_time, draw_args=True), y_tit='Slope [mV/min]', show=show)
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region MODEL
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
        self.Draw(f1, logx=True)
        format_histo(f1, y_range=[0, 2 * b], lw=2, x_range=[.1, 1e5])
        f2.Draw('same')
        f3.SetLineStyle(7)
        f3.Draw('same')
        Draw.add(f2, f3)

    def fit_pulse_height(self, logx=True):
        values = self.get_pulse_heights(pbar=False)
        g = Draw.make_tgrapherrors('gfph', 'Pulse Height Fit', x=self.get_fluxes(), y=values)

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
        self.Draw(g, logx=logx, lm=.12, stats=set_statbox(fit=True))
        l1 = Draw.horizontal_line(f.GetParameter(0), .1, 1e6, color=2, style=7, w=2)
        l2 = Draw.horizontal_line(f.GetParameter(0) - f.GetParameter(1) + f.GetParameter(3), .1, 1e6, color=4, style=7, w=2)
        leg = Draw.make_legend(w=.2)
        leg.AddEntry(l1, 'c', 'l')
        leg.AddEntry(l2, 'c - c_{1} + c_{2}', 'l')
        leg.Draw()
        f.SetLineStyle(7)
        f.Draw('same')
        g.GetListOfFunctions()[-1].Draw()  # redraw stats
        return f
    # endregion MODEL
    # ----------------------------------------

    # ----------------------------------------
    # region CURRENT
    def draw_currents(self, v_range=None, rel_time=False, averaging=True, with_flux=False, c_range=None, f_range=None, draw_opt='al', show=True):
        self.Currents.draw(rel_time=rel_time, v_range=v_range, averaging=averaging, with_flux=with_flux, c_range=c_range, f_range=f_range, show=show, draw_opt=draw_opt)

    def draw_iv(self, show=True):
        g = self.Currents.draw_iv(show=False)
        self.Draw(g, 'IV', draw_opt='ap', logy=True, lm=.12, show=show)

    def draw_current_flux(self, c_range=None, fit=True, show=True):
        currents = [ufloat(0, 0) if c is None else c for c in self.get_values('', self.Analysis.get_current, pbar=False)]
        g = self.Draw.graph(self.get_fluxes(rel_error=.01), currents, title='Leakage Current vs. Flux', show=show, lm=.13, draw_opt='ap', logx=True, logy=True)
        format_histo(g, x_tit='Flux [kHz/cm^{2}]', y_tit='Current [nA]', y_off=1.3, y_range=choose(c_range, [.1, max(currents).n * 2]), x_range=Bins.FluxRange)
        if fit:
            f = TF1('fcf', 'pol1', .1, 1e5)
            f.SetParLimits(0, .1, 5)
            f.SetParLimits(1, 1e-5, 5e-3)
            g.Fit('fcf', 'q')
        Draw.irradiation(make_irr_string(self.get_irradiation()))
        format_statbox(g, fit=fit, all_stat=True, w=.22)
        self.Draw.save_plots('FluxCurrent', show=show)
        return g
    # endregion CURRENT
    # ----------------------------------------

    # ----------------------------------------
    # region TELESCOPE
    def draw_beam_current(self, bin_width=60, rel_t=True, show=True):
        h = TH1F('hr1', 'Beam Current of Run Plan {r}'.format(r=self.RunPlan), *self.get_raw_time_bins(bin_width))
        values = [get_hist_vec(ana.Tel.draw_beam_current(bin_width, show=False, save=False)) for ana in list(self.Analyses.values())]
        values = concatenate([append(run_values, run_values[-1]) for run_values in values])  # add last flux value for time bin between runs
        for i, value in enumerate(values, 1):
            h.SetBinContent(i, value.n)
            h.SetBinError(i, value.s)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', y_off=.85, fill_color=Draw.FillColor, stats=0, markersize=.3, t_ax_off=self.StartTime if rel_t else 0)
        self.Draw(h, 'AllBeamRate', show=show, draw_opt='hist', x=1.5, y=.75, lm=.065)

    def draw_flux(self, bin_width=5, rel_time=True, show=True, redo=False):
        def f():
            h1 = TH1F('hff', 'Flux Profile', *self.get_raw_time_bins(bin_width))
            values = [get_hist_vec(ana.Tel.draw_flux(bin_width, show=False, prnt=False)) for ana in list(self.Analyses.values())]
            values = concatenate([append(run_values, ufloat(0, 0)) for run_values in values])  # add extra zero for time bin between runs
            for i, value in enumerate(values, 1):
                h1.SetBinContent(i, value.n)
                h1.SetBinError(i, value.s)
            return h1
        h = do_pickle(self.make_simple_pickle_path('FullFlux', bin_width, 'Flux'), f, redo=redo)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Flux [kHz/cm^{2}]', t_ax_off=self.get_tax_off(True, rel_time), fill_color=Draw.FillColor, y_range=Bins.FluxRange, stats=0)
        self.Draw(h, 'FluxEvo', w=1.5, h=.75, show=show, logy=True, draw_opt='hist')
        return h
    # endregion TELESCOPE
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL MAP
    def draw_signal_map(self, fid=False, res=.7, hitmap=False, redo=False, show=True):
        histos = self.get_values('signal maps', self.Analysis.draw_signal_map, show=False, prnt=False, hitmap=hitmap, res=res, redo=redo, fid=fid)
        for h in histos[1:]:
            histos[0].Add(h)
        format_histo(histos[0], title='Cumulative {} Map'.format('Hit' if hitmap else 'Signal'))
        self.Draw(histos[0], 'Cumulative{}Map'.format('Hit' if hitmap else 'Signal'), show=show, lm=.12, rm=.16, draw_opt='colz', leg=self.FirstAnalysis.Cut.get_fid())

    def draw_hitmap(self, fid=False, res=.7, redo=False, show=True):
        self.draw_signal_map(fid, res, True, redo, show)

    def draw_signal_maps(self, hitmap=False, redo=False):
        histos = self.get_plots('{} maps'.format('hit' if hitmap else 'signal'), self.Analysis.draw_signal_map, show=False, prnt=False, redo=redo, hitmap=hitmap)
        glob_max = round_up_to(max([h.GetMaximum() for h in histos]), 5) + 5
        glob_min = round_down_to(min([h.GetMinimum() for h in histos]), 5) - 5
        for i, h in enumerate(histos):
            format_histo(h, z_range=[glob_min, glob_max]) if not hitmap else do_nothing()
            self.Draw(h, '{n}Map{nr}'.format(nr=str(i).zfill(2), n='Hit' if hitmap else 'Signal'), show=False, draw_opt='colz', rm=.16, lm=.12, prnt=False)

    def draw_hitmaps(self, redo=False):
        self.draw_signal_maps(hitmap=True, redo=redo)

    def draw_signal_map_ratio(self, run1, run2, m=10, n=10, grid=True, show=True):
        h1, h2 = [self.Analyses[run].split_signal_map(m, n, show=False)[0] for run in [run1, run2]]
        xbins, ybins = self.Analyses[run1].split_signal_map(m, n, show=False)[1:]
        h = h1.Clone('hnew')
        for i in range((h.GetNbinsX() + 2) * (h.GetNbinsY() + 2)):
            v1, v2 = h1.GetBinContent(i), h2.GetBinContent(i)
            h.SetBinEntries(i, 1)
            h.SetBinContent(i, v1 / v2 if v2 else -1)
        format_histo(h, z_range=[0, 3], stats=0, z_tit='Pulse Height Ratio', title='Signal Map Ratio of Run {} & {}'.format(run1, run2))
        self.Draw(h, lm=.12, rm=.16, draw_opt='colzsame', show=show)
        Draw.grid(xbins, ybins, width=2) if grid else do_nothing()
        self.Draw.save_plots('SigMapRatio{}{}'.format(run1, run2))

    def draw_signal_spreads(self, vs_time=True, rel_time=True, show=True):
        spreads = [ana.get_signal_spread(prnt=False) for ana in self.get_analyses()]
        g = Draw.make_tgrapherrors(self.get_x_var(vs_time), spreads, title='Relative Spread')
        format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='Relative Spread [%]', y_off=1.2, t_ax_off=self.get_tax_off(vs_time, rel_time))
        self.Draw(g, 'RelativeSpread', lm=.12, logx=not vs_time, show=show)

    def draw_sm_std(self, vs_time=False, redo=False, show=True):
        y_values = list(self.get_sm_std_devs(redo).values())
        x_values = self.get_times() if vs_time else self.get_fluxes()
        g = Draw.make_tgrapherrors(x_values, y_values, title='STD of the Signal Map')
        format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='rel. STD', y_off=1.3, t_ax_off=0 if vs_time else None)
        self.Draw(g, 'STDSigMap', show, lm=.12, logx=not vs_time)
    # endregion SIGNAL MAP
    # ----------------------------------------

    # ----------------------------------------
    # region BEAM PROFILE
    def draw_beam_info(self, mode='x', fit_margin=.5, vs_time=True, rel_time=True, show=True):
        tits = ['Mean', 'Sigma']
        values = self.get_values('beam profile {}'.format(mode), self.Analysis.draw_beam_profile, show=False, fit=True, fit_range=fit_margin, mode=mode, prnt=False)
        values = [[fit2u(value, par=par) for value in values] for par in [1, 2]]
        graphs = [Draw.make_tgrapherrors(self.get_x_var(vs_time, rel_time), vals, title='{} of the Beam Profile in {}'.format(tit, mode.title())) for tit, vals in zip(tits, values)]
        c = Draw.canvas('Beam Infos {}'.format(mode.title()), x=1.5, y=.75, divide=2, show=show)
        for i, g in enumerate(graphs, 1):
            format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='{tit} [cm]'.format(tit=tits[i - 1]), y_off=1.8, t_ax_off=self.get_tax_off(vs_time, rel_time))
            self.Draw(g, 'BeamProfile{}{}{:1.0f}'.format(tits[i - 1], mode.title(), fit_margin * 100), lm=.125, show=False, logx=not vs_time)
            self.Draw(g, canvas=c.cd(i), lm=.125, show=False, logx=not vs_time)
        self.Draw.save_plots('BeamProfile{}{:1.0f}'.format(mode.title(), fit_margin * 100), show=show, canvas=c)
        return graphs

    def draw_xy_beam_info(self, vs_time=True, show=True, fitx=.5, fity=.5):
        graphs = concatenate([self.draw_beam_info(vs_time=vs_time, show=False, mode=m, fit_margin=f) for m, f in zip(['x', 'y'], [fitx, fity])])
        c = Draw.canvas('Beam Profiles', x=1.5, y=1.5, divide=(2, 2), show=show)
        for i, g in enumerate(graphs, 1):
            format_histo(g, y_off=1.3)
            self.Draw(g, logx=not vs_time, canvas=c.cd(i))
        self.Draw.save_plots('BeamProfileOverview', show=show, canvas=c)

    def draw_beam_profiles(self, mode='x', show=True):
        histos = self.get_values('beam profiles in {}'.format(mode), self.Analysis.draw_beam_profile, show=False, prnt=False, fit=False, mode=mode)
        leg = Draw.make_legend(nentries=self.NRuns)
        stack = THStack('sbp', 'AllBeamProfiles{mod}'.format(mod=mode.title()))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
            format_histo(h, lw=2, stats=0, normalise=True, line_color=self.Draw.get_color(self.NRuns), sumw2=False, fill_style=4000)
            stack.Add(h)
            leg.AddEntry(h, '{0:6.2f} kHz/cm^{{2}}'.format(flux.n), 'l')
        self.Draw(stack, 'AllBeamProfiles{mod}'.format(mod=mode.title()), draw_opt='nostack', show=show, leg=leg)
    # endregion BEAM PROFILE
    # ----------------------------------------

    # ----------------------------------------
    # region TRACKS
    def draw_chi2(self, mode=None, show=True):
        mod_str = '' if mode is None else mode
        histos = self.get_values('chi squares {}'.format(mod_str), self.Analysis.draw_chi2, show=False, prnt=False, mode=mode)
        cut_value = self.FirstAnalysis.Cut.get_chi2(choose(mod_str, 'x'))
        cuts = array([get_quantile(h, [.99, cut_value / 100]) for h in histos])  # .99 for drawing range
        stack = THStack('hx2', '#chi^{{2}}{}'.format(' in {}'.format(mod_str) if mod_str else ''))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
            format_histo(h, stats=0, color=self.Draw.get_color(self.NRuns), lw=2, normalise=True, sumw2=False, fill_style=4000, fill_color=4000)
            stack.Add(h)
        format_histo(stack, x_tit='#chi^{2}', y_tit='Fraction of Events [%]', y_off=1.5, draw_first=True, x_range=ax_range(0, cuts.max(), fh=.4))
        legend = self.make_flux_legend(histos)
        self.Draw(stack, '', show, lm=.15, draw_opt='nostack', leg=legend)
        if mod_str:
            line = Draw.vertical_line(cuts.min(), -1e9, 1e9, color=2, style=2)
            legend.AddEntry(line, 'cut: {}% q'.format(cut_value), 'l')
            histos[0].GetListOfFunctions().Add(line)
            histos[0].GetListOfFunctions().Add(legend)
        self.Draw.save_plots('AllChi2{mod}'.format(mod=mod_str.title()))
        return stack

    def draw_all_chi2s(self, show=True):
        stacks = [self.draw_chi2(mode, show=False) for mode in [None, 'x', 'y']]
        c = Draw.canvas('AllChiSquares', x=2, y=.66, divide=3, show=show)
        for i, s in enumerate(stacks, 1):
            self.Draw(s, draw_opt='nostack', canvas=c.cd(i), show=show)
        self.Draw.save_plots('AllChiSquares', show=show, canvas=c)

    def draw_angle(self, mode='x', show=True):
        histos = self.get_values('angle distribution {}'.format(mode), self.Analysis.draw_angle, mode=mode, show=False, prnt=False)
        legend = self.make_flux_legend(histos)
        stack = THStack('has', 'Track Angles in {mode}'.format(mode=mode.title()))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
            format_histo(h, stats=0, color=self.Draw.get_color(self.NRuns), normalise=True, sumw2=False, fill_color=4000, fill_style=4000)
            stack.Add(h)
        histos[0].GetListOfFunctions().Add(legend)
        format_histo(stack, x_tit='Angle [deg]', y_tit='Fraction of Events [%]', y_off=1.9, draw_first=True, x_range=[-4, 4])
        self.Draw(stack, 'AllTrackAngles{mod}'.format(mod=mode.title()), show, lm=.15, draw_opt='nostack', leg=legend)
        return stack

    def draw_both_angles(self, show=True):
        stacks = [self.draw_angle(mode, show=False) for mode in ['x', 'y']]
        c = Draw.canvas('AllAngles', x=1.5, y=.75, divide=2, show=show)
        for i, s in enumerate(stacks, 1):
            self.Draw(s, draw_opt='nostack', canvas=c.cd(i), show=show, lm=.15)
        self.Draw.save_plots('AllChiAngles', show=show, canvas=c)

    def draw_mean_angles(self, mode='x', vs_time=True, rel_time=True, show=True, redo=False):
        values = [value[0] for value in self.get_values('mean angles {}'.format(mode), self.Analysis.get_mean_angle, mode=mode, redo=redo)]
        g = Draw.make_tgrapherrors(self.get_x_var(vs_time, rel_time), values, title='Mean Track Angles in {}'.format(mode.title()))
        format_histo(g, y_tit='Mean Angle [deg]', y_off=1.3, **self.get_x_args(vs_time, rel_time))
        self.Draw(g, 'MeanTrackAngle'.format(mode.title()), show, logx=not vs_time)
    # endregion TRACKS
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE PLOTS
    def draw_run_currents(self):
        self.get_values('currents', self.Analysis.get_current)
        self.get_values('currents', self.Analysis.draw_current, show=False)

    def draw_chi2s(self):
        self.get_values('chi2s', self.Analysis.draw_chi2s, show=False, prnt=False)

    def draw_angles(self):
        self.get_values('angles', self.Analysis.draw_angles, show=False, prnt=False)

    def draw_occupancies(self):
        self.get_values('occupancies', self.Analysis.draw_occupancies, show=False, prnt=False, cluster=True)
    # end region GENERATE PLOTS
    # ----------------------------------------

    def draw_efficiencies(self, vs_time=False, show=True):
        x, y = self.get_x_var(vs_time), self.get_efficiencies()
        self.Draw.graph(x, y, 'Efficiencies', **self.get_x_args(vs_time, draw_args=True), y_tit='Effciency [%]', show=show, lm=.12, y_off=1.8)


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
        z.Currents.draw()
        input('Press any button to exit')
    if pargs.draw:
        z.draw_all(pargs.redo)
    if pargs.prnt:
        print(z.get_cmd_strings(pargs.command, pargs.kwargs))
