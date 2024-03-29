#! /usr/bin/env python
from ROOT import THStack, TF1, TMath
from numpy import sort, log, argmin, argmax
from uncertainties.umath import sqrt as usqrt  # noqa

import plotting.latex as tex
from helpers.utils import *
from plotting.draw import *
from src.analysis import Analysis
from src.currents import Currents
from src.dut_analysis import DUTAnalysis, Bins, reload_tree
from src.run_selection import RunPlan, RunSelection


class AnalysisCollection(Analysis):
    """ Analysis of the various runs of a single runplan. """

    StartTime = None
    PhUnit = '[mV]'
    PhTit = f'Pulse Height {PhUnit}'

    def __init__(self, name, dut_nr, testcampaign=None, load_tree=True, verbose=False):

        # RUN SELECTION
        self.Ensemble = RunPlan(name, testcampaign, dut_nr, verbose) if isfloat(name) else RunSelection(name, verbose)
        self.Runs = array(self.Ensemble.get_runs())
        self.NRuns = self.Runs.size
        self.RunPlan = str(self.Ensemble)
        self.DUT = self.Ensemble.DUT
        self.Type = self.Ensemble.Type
        self.Fluxes = self.Ensemble.get_fluxes()
        self.MinFluxRun, self.MaxFluxRun = self.get_high_low_rate_runs()

        super(AnalysisCollection, self).__init__(self.Ensemble.TCString, self.Ensemble.res_dir, sub_dir=self.Ensemble.save_dir, verbose=verbose)
        self.print_start(name, prnt=load_tree, dut=self.DUT.Name)
        self.print_start_info()

        # Loading the Trees and Time Vectors in Parallel
        self.LoadTree = load_tree

        # Loading the Single Analyses
        self.Analysis = self.load_dummy()  # dummy to get access to the methods
        self.Analyses = self.load_analyses()
        self.FirstAnalysis = self.Analyses[0]
        self.LastAnalysis = self.Analyses[-1]
        self.Bins = self.FirstAnalysis.Bins if load_tree else None
        AnalysisCollection.StartTime = self.FirstAnalysis.Run.StartTime if load_tree else time_stamp(self.FirstAnalysis.Run.LogStart)

        # Sub Classes
        self.Currents = Currents(self)

        self.print_finished() if verbose else print()

    def __str__(self):
        return f'{self.__class__.__name__} of {self.Ensemble} and {self.DUT}'

    def __repr__(self):
        return self.__str__()

    def __del__(self):
        if self.Verbose:
            print('\n good bye... ')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        from plotting.save import SaveDraw
        SaveDraw.File = None

    def i(self):
        print(f'{self.DUT} @ {tex.si(self.DUT.Bias, fmt="+.0f", unt="V")[0]}, Runplan {self.Ensemble.Name}, {tc2str(self.TCString, short=False)}.')

    @property
    def server_save_dir(self):
        return '' if 'RS' in self.Ensemble else Path('diamonds', str(self.DUT), self.TCString, str(self.Ensemble))

    # ----------------------------------------
    # region RATE DEPENDENCE PARAMETERS
    @save_pickle('p0Chi2')
    def pol0_chi2(self, _redo=False):
        return FitRes(self.draw_avrg_ph(redo=_redo, show=False, prnt=False).Fit('pol0', 'qs0')).get_chi2()

    @save_pickle('RelSpread')
    def rel_spread(self, _redo=False):
        x = self.get_pulse_heights(avrg=True)
        s = (max(x) - min(x)) / mean(uarr2n(x))
        return AsymVar(s.n, s.s1, min(s.n, s.s0) - s.s1) if type(s) is AsymVar else s

    def print_rate_pars(self, prec=1):
        print(f'Chi²:       {self.pol0_chi2():3.{prec}f}')
        print(f'Rel Spread: {self.rel_spread() * 100:3.{prec}f}%')
    # endregion RATE DEPENDENCE PARAMETERS
    # ----------------------------------------

    def draw_all(self, redo=False):
        pass

    def show_information(self):
        print_table(rows=self.get_values('', self.Analysis.show_information, pbar=False, prnt=False, ret_row=True), header=self.FirstAnalysis.info_header)

    def print_loaded(self):
        print('\033[1A\rRuns {0}-{1} were successfully loaded!{2}\n'.format(self.Runs[0], self.Runs[-1], 20 * ' '))

    def make_flux_legend(self, histos, runs=None, x1=.76, y2=.88, draw_style='l'):
        runs = choose(runs, self.Runs)
        legend = self.Draw.make_legend(x1, y2, nentries=len(runs))
        for h, run in zip(histos, runs):
            legend.AddEntry(h, flux2str(self.get_fluxes()[list(self.Runs).index(run)].n, prec=0), draw_style)
        return legend

    def set_verbose(self, status):
        super().set_verbose(status)
        for ana in self.get_analyses():
            ana.set_verbose(status)

    def save_all(self):
        t = info('creating all data ...')
        self.save_coll_plots()
        self.Draw.close_file()
        self.save_data()
        self.save_plots()
        print_elapsed_time(t, color='green')

    @quiet
    def save_plots(self):
        self.parallel(self.Analysis.save_plots, print_link=False)

    @quiet
    def save_data(self):
        for i, data in enumerate(self.parallel(self.Analysis.get_data)):
            self.Analyses[i].save_data(data)

    @quiet
    def save_coll_plots(self, prnt=True):
        self.draw_flux(show=False)
        self.draw_currents(show=False, fname='Currents')
        self.draw_pulse_heights(show=False)
        self.Draw.print_http('plots.html', force_print=prnt)

    @staticmethod
    @reload_tree
    def prep_f(ana, i, *args, f, **kwargs):
        PBAR.update(i) if PBAR.PBar and not PBAR.is_finished() else do_nothing()
        return f(ana, *args, **kwargs)

    @quiet
    def parallel(self, f, *args, runs=None, pbar=True, **kwargs):
        PBAR.start(self.NRuns if runs is None else len(runs)) if pbar else do_nothing()
        with Pool() as pool:
            res = pool.starmap(partial(self.prep_f, f=f, **kwargs), [(ana, i, *args) for i, ana in enumerate(self.get_analyses(runs))])
            PBAR.set_last()
            return res

    def test_anas(self, f=None, *args, **kwargs):
        for ana in self.Analyses:
            try:
                getattr(ana, choose(f, 'save_data'))(*args, **kwargs)
                print(ana)
            except Exception as err:
                warning(f'{ana}, {err}')

    # ----------------------------------------
    # region INIT
    @property
    def info_header(self):
        return ['Run', 'Flux [kHz/cm2]', 'Bias [V]', 'Start', 'Duration [hh:mm]']

    @property
    def info_rows(self):
        bias = [f'{bias:+8.0f}' for bias in self.Ensemble.get_biases()]
        times = [datetime.fromtimestamp(duration - 3600).strftime('%H:%M').rjust(15) for duration in self.Ensemble.get_durations()]
        return array([self.Runs, [f'{flux.n:14.1f}' for flux in self.Fluxes], bias, self.Ensemble.get_start_times(), times]).T

    def print_start_info(self):
        print_table(header=self.info_header, rows=self.info_rows, prnt=self.Verbose)

    def get_high_low_rate_runs(self):
        return self.Runs[argmin(self.Fluxes)], self.Runs[argmax(self.Fluxes)]

    @staticmethod
    def load_dummy():
        return DUTAnalysis

    def copy_raw_files(self, force=False):
        if self.LoadTree and not self.Ensemble.final_files_exist and not self.Ensemble.raw_files_exist or force:
            self.Ensemble.copy_raw_files()

    def load_analyses(self, r0=0):
        self.copy_raw_files()
        with Pool() as pool:
            res = pool.starmap(self.Analysis, [(run.Number, dut, run.TCString, self.LoadTree, self.Verbose, False) for run, dut in zip(self.Ensemble.Runs, self.Ensemble.get_dut_nrs()) if run > r0])
        if self.LoadTree:
            for r in res:
                r.reload_tree_()
                r.Cut.generate_fiducial()
        return res

    def reload_anas(self, r0, load_tree=None):
        self.LoadTree = choose(load_tree, self.LoadTree)
        self.Analyses = self.load_analyses(r0)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region BINNING
    def get_binning(self, bw=None, rel_time=False):
        b = concatenate([ana.Bins.get_raw_time(bw, t_from_event=True)[1] for ana in self.Analyses])
        return [b.size - 1, b - (self.FirstAnalysis.Run.StartTime if rel_time else 0)]

    def get_time_bins(self, bw=None, only_edges=False):
        b = self.fix_t_arrays([ana.Bins.get_time(bw)[1] for ana in self.Analyses])
        b = array([[b[0], b[-1]] for b in b]).flatten() if only_edges else concatenate(b)
        return [b.size - 1, b]

    def get_raw_time_bins(self, bw=None, only_edges=False, t_from_event=False):
        b = self.fix_t_arrays([ana.Bins.get_raw_time(bw, t_from_event=t_from_event)[1] for ana in self.Analyses])
        b = array([[b[0], b[-1]] for b in b]).flatten() if only_edges else concatenate(b)
        return [b.size - 1, b]

    def fix_t_arrays(self, t_arrays):
        """ Add the logged time between two runs if the start time of a run is less than the stop time of the previous. """
        for i in range(len(t_arrays) - 1):
            delta = t_arrays[i + 1][0] - t_arrays[i][-1]  # check the time difference from the first bin of the next array and the last of the current
            if delta < 0:
                t_arrays[i + 1] += self.get_break_time(i) - delta
        return t_arrays

    def get_break_time(self, ind):
        # because the log stop is not so reliable...
        return (self.Analyses[ind + 1].Run.LogStart - self.Analyses[ind].Run.LogStart - self.Analyses[ind].Run.Duration).total_seconds()
    # endregion BINNING
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_runs(self):
        return self.Runs

    def get_analyses(self, runs=None):
        return self.Analyses if runs is None else [ana for ana in self.Analyses if ana.Run.Number in runs]

    def get_hv_name(self):
        return self.Currents.Name

    @quiet
    def get_values(self, what, f, runs=None, pbar=True, avrg=False, picklepath=None, flux_sort=False, plots=False, **kwargs):
        runs = choose(runs, self.Runs)
        redo = 'redo' in kwargs and kwargs['redo'] or '_redo' in kwargs and kwargs['_redo']
        if picklepath is not None and all([file_exists(picklepath.format(run)) for run in runs]) and not redo:
            values = [load_pickle(picklepath.format(run)) for run in runs]
        else:
            self.info(f'Generating {what} ...', prnt=pbar)
            values = self.parallel(f, runs=runs, pbar=pbar, **kwargs)
        return values if plots else array(self.get_flux_average(array(values))) if avrg else array(values, dtype=object)[self.get_fluxes().argsort() if flux_sort else ...]

    def get_plots(self, string, f, runs=None, pbar=True, avrg=False, picklepath=None, **kwargs):
        return self.get_values(string, f, runs, pbar, avrg, picklepath, False, True, **kwargs)

    def get_fluxes(self, plane=None, corr=True, full_size=False, runs=None, avrg=False, pbar=None, rel=False, redo=False):
        picklepath = self.get_pickle_path(f'Flux', make_suffix(0, plane, corr, 1, full_size, self.DUT.Number), 'Telescope', dut='')
        pbar = False if not self.FirstAnalysis.has_branch('rate') else pbar
        values = self.get_values('fluxes', DUTAnalysis.get_flux, runs, pbar, avrg=avrg, picklepath=picklepath, plane=plane, corr=corr, full_size=full_size, redo=redo)
        return array([ufloat(v.n, v.n * .01) for v in values]) if rel else values

    def flux_strings(self, prec=0, runs=None):
        return flux2str(uarr2n(self.get_fluxes(runs)), prec)

    def bias_strings(self, root=True):
        return bias2str([ana.DUT.Bias for ana in self.Analyses], root)

    @save_pickle('Splits', sub_dir='Flux')
    def get_flux_splits(self, corr=True, w=.1, n=50, _redo=False):
        return self.draw_flux_splits(corr, w, n, show=False)

    def draw_flux_splits(self, corr=True, w=.1, n=50, **kwargs):
        x = sort([flux.n for flux in self.get_fluxes(corr=corr)])
        stats = set_statbox(entries=True, w=.2)
        h = self.Draw.distribution(x, log_bins(n, 1, 1e5), 'Flux Splits', **prep_kw(kwargs, **Draw.mode(2), **self.get_x_args(draw=True), file_name='FluxSplits', stats=stats))
        split_bins = histogram(x, self.find_flux_binning(h, w))[0]
        return cumsum(split_bins[where(split_bins > 0)])[:-1]

    def find_flux_binning(self, h, width=.1):
        maxima = find_maxima(h, 20, sigma=1, sort_x=True)[:, 0]
        b = [0.] + [v for i in maxima for v in [i / 10 ** width, i * 10 ** width]] + [1e5]
        return b if all(diff(b) > 0) else self.find_flux_binning(h, width - .01)

    def get_flux_average(self, values):
        values = values[self.get_fluxes().argsort()]  # sort by ascending fluxes
        return array([mean_sigma(lst)[0] for lst in split(values, self.get_flux_splits())])  # split into sub-lists of similar flux and take average

    def get_times(self, runs=None):
        return self.get_values('times', DUTAnalysis.get_time, runs, pbar=False)

    def get_events(self, redo=False):
        return self.get_values('events', DUTAnalysis.get_events, picklepath=self.make_simple_hdf5_path('', 'AllCuts', 'Events', '{}'), redo=redo)

    def get_n_events(self, flux_sort=False, redo=False):
        return self.get_values('n_events', DUTAnalysis.get_n_entries, picklepath=self.get_pickle_path(sub_dir='Entries'), _redo=redo, flux_sort=flux_sort)

    def get_x_var(self, vs_time=False, vs_irrad=False, avrg=False):
        return self.get_times() if vs_time else self.get_irradiations() if vs_irrad else self.get_fluxes(avrg=avrg)

    def get_irradiation(self):
        return self.FirstAnalysis.get_irradiation()

    def get_irradiations(self):
        return array([ana.get_irradiation() for ana in self.Analyses]) / 1e15

    def get_attenuator(self):
        return self.FirstAnalysis.get_attenuator()

    def get_hv_device(self):
        return self.FirstAnalysis.Currents.Name

    def get_currents(self):
        return [ana.Currents.get() for ana in self.Analyses]

    @staticmethod
    def get_mode(vs_time):
        return 'Time' if vs_time else 'Flux'

    def get_sm_stds(self, runs=None, redo=False):
        pickle_path = self.make_simple_pickle_path('Signal', f'AllCuts_{self.FirstAnalysis.Cut.get_chi2()}_None', sub_dir='SignalMaps', run='{}')
        return self.get_values('SM STD', self.Analysis.get_sm_std, runs, picklepath=pickle_path, redo=redo)

    def get_sm_std(self, redo=False, low=False, high=False):
        def f():
            runs = self.get_runs_below_flux(110) if low else self.get_runs_above_flux(2000) if high else self.Runs
            return mean_sigma(self.get_sm_stds(runs, redo))[0]
        return do_pickle(self.make_simple_pickle_path('SMSTD', f'{low:d}{high:d}', 'Uniformity'), f, redo=redo)

    def get_pulse_heights(self, *args, **kwargs):
        return array([])

    def get_pulse_height(self, *args, **kwargs):
        return mean_sigma(self.get_pulse_heights(*args, **kwargs))

    def get_efficiencies(self, suf='3', avrg=False, redo=False):
        return self.get_values('efficiencies', self.Analysis.get_efficiency, picklepath=self.get_pickle_path(suf=suf, sub_dir='Efficiency'), avrg=avrg, redo=redo)

    def get_runs_below_flux(self, flux):
        return self.Runs[self.get_fluxes() <= flux]

    def get_runs_above_flux(self, flux):
        return self.Runs[self.get_fluxes() >= flux]

    def get_signal_spread(self, peaks=False, redo=False, rel=False):
        groups = split(self.get_pulse_heights(redo=redo, err=False, peaks=peaks, flux_sort=True), self.get_flux_splits())
        values = array([(value - mean_sigma(grp)[0]) / (mean_sigma(grp)[0] if rel else 1) for grp in groups for value in grp if grp.size > 1])
        return values if values.size > 1 else None

    def get_repr_error(self, flux=None, peaks=False, redo=False):
        values = self.get_signal_spread(peaks, redo)
        return self.get_low_flux_std(flux) if values is None else mean_sigma(values, err=False)[1]

    def sys_error(self):
        return self.MainConfig.get_value('Parameters', 'systematic error', dtype=float)

    def calc_max_spread(self, redo=False):
        s = array([(max(g) - min(g)) / mean_sigma(g)[0].n for g in split(self.get_pulse_heights(redo=redo, err=False, flux_sort=True), self.get_flux_splits()) if g.size > 1])
        return max(s)

    def calc_rel_sys_error_old(self):
        x = self.get_pulse_heights(err=False)
        e_stat, (m, e_full) = mean_sigma([v.s for v in x])[0], mean_sigma(x)
        return usqrt(e_full ** 2 - e_stat ** 2) / m

    def subtract_grp_mean(self, err=True):
        x, y = sorted(self.get_fluxes()), self.get_pulse_heights(err=err, flux_sort=True)
        return array([[iy, ix - mean_sigma(g[:, 1])[0].n] for g in split(array([x, y]).T, self.get_flux_splits()) for iy, ix in g if g.shape[0] > 1]).T  # take only groups with more than 1 entry

    @save_pickle('RelSErr')
    def calc_rel_sys_error(self, _redo=False):
        x, y = self.subtract_grp_mean(err=False)
        f = Draw.make_tf1('chi', lambda i: FitRes(self.Draw.graph(x, add_err(y, i), show=False, save=False).Fit('pol0', 'qs')).get_chi2(), 0, 1)
        v = f.GetX(1)  # get the error size such that the chi2 is 1
        return ufloat(v, abs(v - f.GetX(1.1))) / self.get_pulse_height(err=False)[0]

    @save_pickle('LowFlux', sub_dir='Errors', suf_args=0)
    def get_low_flux_std(self, flux, _redo=False):
        x = self.get_pulse_heights(runs=self.get_runs_below_flux(flux), redo=_redo, err=False)
        return 0 if not x.size else .01 if x.size == 1 else mean_sigma(x, err=False)[1]

    def get_uniformities(self, use_fwc=True, low_flux=False, high_flux=False, avrg=False, redo=False):
        runs = self.get_runs_below_flux(110) if low_flux else self.get_runs_above_flux(2000) if high_flux else self.Runs
        picklepath = self.get_pickle_path('Uniformity', int(use_fwc), 'Signal')
        return self.get_values('uniformities', self.Analysis.get_uniformity, runs, picklepath=picklepath, _redo=redo, use_fwc=use_fwc, avrg=avrg)

    def get_mean_uniformity(self, use_fwc=True, redo=False, low_flux=False, high_flux=False):
        values = self.get_uniformities(use_fwc, redo, low_flux, high_flux)
        return [mean_sigma(values[:, i][where(values[:, i] > 0)[0]]) for i in range(values[0].size)]

    @staticmethod
    def get_x_tit(vs_time=False, vs_irrad=False):
        return 'Time [hh:mm]' if vs_time else 'Fluence [10^{15}n/cm^{2}]' if vs_irrad else 'Flux [kHz/cm^{2}]'

    @staticmethod
    def get_x_draw(vs_time=False):
        return {'logx': not vs_time, 'grid': vs_time}

    @staticmethod
    def get_tax_off(vs_time, rel_time=False, off=0):
        return None if not vs_time else AnalysisCollection.StartTime if rel_time else off

    @staticmethod
    def get_range(t, i, x_range=None):
        return x_range if t else Bins.FluenceRange if i else Bins.FluxRange

    @staticmethod
    def get_x_args(vs_time=False, rel_time=False, vs_irrad=False, draw=True, off=0, **kwargs):
        kwargs = prep_kw(kwargs, x_tit=AnalysisCollection.get_x_tit(vs_time, vs_irrad), t_ax_off=AnalysisCollection.get_tax_off(vs_time, rel_time, off),
                         x_range=AnalysisCollection.get_range(vs_time, vs_irrad), x_off=None if vs_time or vs_irrad else 1.1)
        return {**kwargs, **AnalysisCollection.get_x_draw(vs_time or vs_irrad)} if draw else kwargs

    def get_cmd_strings(self, cmd, kwargs):
        return '?'.join(['python analyse.py {} {} -tc {} -d -cmd {} -kw {}'.format(run, self.DUT.Number, self.TCString, cmd, kwargs) for run in self.Runs])

    def get_pickle_path(self, name='', suf='', sub_dir=None, dut=None, camp=None):
        return self.FirstAnalysis.make_simple_pickle_path(name, suf, sub_dir, run='{}', dut=dut, camp=camp)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def ph_spread(self, avrg=True, rel=True):
        x = self.get_pulse_heights(avrg=avrg)
        x = self.scale_ph(x, avrg) if rel else x
        return max(x) - min(x)

    @staticmethod
    def draw_legend(graphs, x=.17):
        return Draw.legend(graphs, [g.GetTitle() for g in graphs], ['l' if i < 2 else 'p' for i in range(len(graphs))], x1=x, y1=.21, w=.2)

    def scale_ph(self, x, avrg=False):
        """ scale the ph to the mean of the pulse heights in the 'min flux range' """
        f, (fmin, fmax) = self.get_fluxes(avrg=avrg), self.MainConfig.get_list('MAIN', 'min flux range')
        x0 = array(x)[where((f >= fmin) & (f <= fmax))]
        scale = mean(x0).n if is_ufloat(x0[0]) else mean(x0[:, 0])  # no error to keep the correct errors of single measurements
        return x / scale

    def make_pulse_height_graph(self, bw=None, t=False, irr=False, first_last=True, redo=False, legend=True, corr=True, err=True, avrg=False, peaks=False, scale=False, ls=.04, e_sys=None):
        x, (ph0, ph) = self.get_x_var(t, irr, avrg), [self.get_pulse_heights(bw, redo and not i, corr=corr, err=e, avrg=avrg, peaks=peaks, e_sys=e_sys) for i, e in enumerate([False, err])]
        ph0, ph = [self.scale_ph(iph, avrg) for iph in [ph0, ph]] if scale else [ph0, ph]
        g = Draw.make_tgraph(x, ph0, title='stat. error', color=Draw.color(2, 1), markersize=.7, y_tit=self.PhTit)
        g_errors = Draw.make_tgraph(x, ph, title='full error', marker=0, color=Draw.color(2, 0), markersize=0, y_tit=self.PhTit)
        g1, g_last = [Draw.make_tgraph([x[i].n], [ph[i].n], title='{} run'.format('last' if i else 'first'), marker=22 - i, color=2, markersize=1.5) for i in [0, -1]]
        graphs = ([g_errors] if err else []) + ([] if avrg else [g]) + ([g1, g_last] if first_last and not avrg and not t else [])
        mg = self.Draw.multigraph(graphs, 'Pulse Height', color=None, show=False)
        if legend and not avrg:
            mg.GetListOfFunctions().Add(self.Draw.legend(graphs, [g.GetTitle() for g in graphs], ['p' if 'run' in g.GetTitle() else 'l' for g in graphs], bottom=True, left=True))
        if t:
            for ix, iph, flux in zip(x, ph0, self.get_fluxes()):
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(Draw.tlatex(ix.n, iph.n + iph.s * 1.2, '{:0.0f}'.format(flux.n), color=1, align=21, size=ls))
        return mg

    def draw_low_scale(self, avrg=True, yoff=.07):
        x, y = self.get_fluxes(avrg=avrg, corr=False), self.get_pulse_heights(err=False, avrg=avrg)
        y /= y[argmin(x)].n
        self.Draw.graph(x, y, y_tit=f'Scaled {self.PhTit}', y_range=[1 - yoff, 1 + yoff], **self.get_x_args(draw=True), lm=.12)

    def draw_mean_scale(self, avrg=True, yoff=.07):
        x, y = self.get_fluxes(avrg=avrg, corr=False), self.get_pulse_heights(err=False, avrg=avrg)
        y /= mean(y).n
        self.Draw.graph(x, y, y_tit=f'Scaled {self.PhTit}', y_range=[1 - yoff, 1 + yoff], **self.get_x_args(draw=True), lm=.12)

    def draw_mid_mean_scale(self, avrg=True, yoff=.07):
        x, y = self.get_fluxes(avrg=avrg, corr=False), self.get_pulse_heights(err=False, avrg=avrg)
        y /= mean(y[1:-1]).n
        self.Draw.graph(x, y, y_tit=f'Scaled {self.PhTit}', y_range=[1 - yoff, 1 + yoff], **self.get_x_args(draw=True), lm=.12)

    def draw_scaled_pulse_heights(self, bw=None, t=False, redo=False, avrg=False, yr=.05, irr=True, **dkw):
        """ Shows the scaled pulse heights of the single runs. """
        x, y = self.get_x_var(t, avrg=avrg), self.scale_ph(self.get_pulse_heights(avrg=avrg, bw=bw, redo=redo), avrg)
        g = self.Draw.graph(x, y, 'RelPH', **prep_kw(dkw, y_tit=f'Relative Pulse Height', y_range=[1 - yr, 1 + yr], ndivx=503, show=False))
        leg = Draw.irradiation(irr2str(self.Ensemble.get_irradiation())) if irr else None
        return self.Draw.graph(g, **prep_kw(dkw, **Draw.mode(2 if t else 1), **self.get_x_args(t), gridy=True, leg=leg, file_name=fname('RelPH', avrg, t)))

    def draw_avrg_ph(self, bw=None, t=False, i=False, corr=True, err=True, fit=False, redo=False, **dkw):
        x, y = self.get_x_var(t, i, avrg=True), self.get_pulse_heights(bw, corr=corr, err=err, redo=redo, avrg=True, **dkw)
        g = self.Draw.graph(x, y, 'AvrPH', y_tit=self.PhTit)
        g.Fit('pol0', 'qs') if fit else do_nothing()
        return self.Draw.graph(g, **prep_kw(dkw, **self.get_x_args(t, vs_irrad=i, draw=True), stats=set_statbox(fit=fit, form='2.1f', stats=fit), file_name=fname('PH', True, t, i)))

    def draw_pulse_heights(self, bw=None, t=False, irr=False, show_first_last=True, legend=True, corr=True, redo=False, err=True, avrg=False, fit=False, peaks=False, ls=.04, e_sys=None, **dkw):
        """ Shows the pulse heights of the runs. """
        mg = self.make_pulse_height_graph(bw, t, irr, show_first_last, redo, legend, corr, err, avrg, peaks, ls=ls, e_sys=e_sys)
        mg.GetListOfGraphs()[0].Fit('pol0', f'qs') if fit else do_nothing()
        stats = set_statbox(fit=fit, form='2.1f', stats=fit)
        m = Draw.mode(2, y_off=.85, lm=.1) if t else Draw.mode(1, lm=.14, y_off=1.45)
        return self.Draw(mg, **prep_kw(dkw, **self.get_x_args(t, vs_irrad=irr, draw=True), file_name=fname('PH', avrg, t, irr), stats=stats, color=None, **m))

    @save_pickle('Full', sub_dir='PH', suf_args='all')
    def get_full_ph(self, bin_size=None, _redo=False):
        g = self.get_plots('ph trends', self.Analysis.get_pulse_height_trend, bin_size=bin_size, _redo=_redo, picklepath=self.get_pickle_path('Trend', make_suffix(self, bin_size, 1), 'PH'))
        ph = concatenate([append(graph_y(i), ufloat(0, 0)) for i in g])  # add a zero after each run for the bin in between
        return self.Draw.distribution(ph, self.get_time_bins(bin_size), 'Full Pulse Height', **self.get_x_args(True), y_tit=f'Mean {self.PhTit}', show=False)

    def draw_full_pulse_height(self, bin_size=None, rel_t=True, with_flux=True, redo=False, show=True, **dkw):
        """ Shows the pulse heights bins of all runs vs time. """
        if with_flux:
            self.draw_flux(rel_time=rel_t, fill_color=2, fill_style=3002, draw_opt='histy+', rm=.1, x_off=10, l_off_x=10, y_range=Bins.FluxRange, show=show)
        c = self.Draw.tpad(transparent=True, c=get_last_canvas()) if with_flux else None
        h = self.get_full_ph(bin_size, _redo=redo)
        return self.Draw.distribution(h, canvas=c, **prep_kw(dkw, y_range=[0, h.GetMaximum() * 1.05], **self.get_x_args(True, rel_t, draw=True), stats=0, **Draw.mode(2),
                                      rm=.1 if with_flux else None, fill_style=3002, file_name=f'FullPulseHeight{"Flux" if with_flux else""}', y_tit=self.PhTit))

    def draw_splits(self, m=2, show=True, normalise=False):
        x, y = self.get_x_var(), self.get_values('split pulse heights', DUTAnalysis.get_split_ph, m=m).T
        y /= mean(y, axis=1).reshape(m ** 2, 1) if normalise else 1
        graphs = [self.Draw.graph(x, iy, y_tit=self.PhTit, **self.get_x_args(), show=False) for iy in y]
        mg = self.Draw.multigraph(graphs, 'Region Pulse Heights', ['{}'.format(i) for i in range(y.shape[0])], show=show, draw_opt='alp', logx=True)
        format_histo(mg, **self.get_x_args())

    def draw_signal_legend(self):
        pass

    def draw_signal_distributions(self, bw=None, redo=False, **dkw):
        """Shows a stack of the signal distributions."""
        h = self.get_plots('signal distributions', self.Analysis.get_signal_distribution, picklepath=self.get_pickle_path('Disto', '1_0', 'PH'),  _redo=redo, bw=bw)
        return self.Draw.stack(h, 'PHDists', **prep_kw(dkw, leg_titles=self.flux_strings(), y_off=1, scale=True, file_name='PHDists'))

    def draw_ph_trends(self, bw=None, redo=False, **dkw):
        bw = choose(bw, Bins.Size)
        g = self.get_plots('ph trends', self.Analysis.get_pulse_height_trend, picklepath=self.get_pickle_path('Trend', f'1_20_{bw}', 'PH'), bw=bw, _redo=redo)
        g = [shift_graph(ig, ox=-graph_x(ig)[0].n) for ig in g]
        return self.Draw.multigraph(g, 'PHTrends', self.flux_strings(), **prep_kw(dkw, gridy=True, **self.get_x_args(vs_time=True, off=-1, draw=True), file_name='PHTrends'))

    def draw_pulls_below_flux(self, bw=.5, flux=150, show=True):
        suffix = '{}_{}'.format(bw, flux)
        return do_pickle(self.make_simple_pickle_path('PhPullBel', suffix, 'Ph_fit'), self.draw_pulls, redo=show, runs=self.get_runs_below_flux(flux), bin_width=bw, show=show)

    def draw_pulls(self, runs=None, bw=.5, show=True):
        s = THStack('s_phd', 'Pulse Height Distributions')
        runs = choose(runs, self.Runs)
        histos = self.get_values('Generate ph pulls ...', self.Analysis.draw_ph_pull, runs, show=False, bin_width=bw, fit=False, save=False)
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
        b = bins.make(*ax_range(values, fl=.5, fh=.5), nb=3 * sqrt(values.size)) if rel else [40, -2, 2]
        h = self.Draw.distribution(values, b, '{}Signal Spread'.format(t), x_tit='{}Diff [{}]'.format(t, '%' if rel else 'mV'), y_tit='Number of Entries', show=False)
        self.Draw(h, 'SignalSpread', lm=.11, show=show, save=save)
        return values

    def draw_fwhm(self, arg=1, t=False, vs_irrad=False, use_fwc=True, avrg=False, redo=False, **dkw):
        x, y = self.get_x_var(t, vs_irrad, avrg), self.get_uniformities(use_fwc, redo, avrg=avrg)[:, arg]
        ph_var = 'FWC' if use_fwc else 'Mean'
        var, unit, tit = [ph_var, 'FWHM', f'FWHM/{ph_var}'][arg], f' {self.PhUnit}' if arg < 2 else '', [ph_var, 'Full Width Half Maximum', 'Uniformity'][arg]
        g = self.Draw.graph(x[y != 0], y[y != 0], **prep_kw(dkw, title=tit, y_tit=f'{var}{unit}', **self.get_x_args(t, vs_irrad=vs_irrad, draw=True)))
        if vs_irrad and arg == 2:
            u_scvd = self.MainConfig.get_ufloat('Parameters', 'uniformity scCVD')
            Draw.make_tgraph([0], [u_scvd], color=634, marker=21).Draw('p')
            g_scvd = self.Draw.make_tgraph([0, 10], [u_scvd] * 2, color=634, fill_color=634, opacity=.2, marker=21)
            g_scvd.Draw('le3')
            self.Draw.legend([g, g_scvd], ['pCVD', 'scCVD'], ['p', 'plf'])
        self.Draw.save_plots(fname('Uni' if arg == 2 else var, avrg, t))
        return g

    def draw_means(self, vs_time=False, vs_irrad=False, avrg=False, redo=False, **kwargs):
        return self.draw_fwhm(0, vs_time, vs_irrad, False, avrg, redo, **kwargs)

    def draw_fwc(self, vs_time=False, vs_irrad=False, avrg=False, redo=False, **kwargs):
        return self.draw_fwhm(0, vs_time, vs_irrad, True, avrg, redo, **kwargs)

    def draw_uniformity(self, vs_time=False, vs_irrad=False, use_fwc=True, avrg=False, redo=False, **kwargs):
        return self.draw_fwhm(2, vs_time, vs_irrad, use_fwc, avrg, redo, **kwargs)

    def draw_ph_slope(self, vs_time=False, show=True):
        y = [fit2u(ana.draw_pulse_height(show=False)[0].Fit('pol1', 'qs0'), par=1) * 60 for ana in self.get_analyses()]
        self.Draw.graph(self.get_x_var(vs_time), y, 'Pulse Height Slope', **self.get_x_args(vs_time, draw=True), y_tit='Slope [mV/min]', show=show)

    def draw_expo(self, bw=None, fr=None, y2=None, **dkw):
        g = [ana.fit_expo(bw, show=False, save=False, shift=True, fr=fr, y2=y) for ana, y in zip(self.Analyses, [None] * self.NRuns if y2 is None else y2)]
        l_tits = [f'{t}, {flux2str(f.n)}' for t, f in zip(['down scan', 'up scan'], self.get_fluxes())]
        return self.Draw.multigraph(g, 'PHExpo', l_tits, **prep_kw(dkw, gridy=True, **self.get_x_args(vs_time=True, off=-1, draw=True), file_name='PHExpo', wleg=.3))

    def draw_tconst(self, bw=None, avrg=False, **dkw):
        x, y = self.get_x_var(avrg=avrg), self.get_values('expo fits', self.Analysis.expo_fit, bw=bw, picklepath=self.get_pickle_path('ExpoFit', 3, 'PH'))
        y = self.get_flux_average(y[:, -1]) if avrg else y[:, -1]
        return self.Draw.graph(x, y / 60, **prep_kw(dkw, **self.get_x_args(), y_tit='Time Constant [min]', file_name='TConst'))

    def draw_0fit(self, err=True, ey=0, **dkw):
        x, y = self.subtract_grp_mean(err)
        y = add_err(y, ey)
        g = self.Draw.graph(x, y, y_tit='Shifted Pulse Heigt [mV]')
        g.Fit('pol0', 'qs')
        return self.Draw(g, **prep_kw(dkw, **self.get_x_args(), stats=set_statbox(fit=True), file_name='ShiftedPHFit'))
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
        g = Draw.make_tgraph('gfph', 'Pulse Height Fit', x=self.get_fluxes(), y=values)

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
        format_histo(g, x_tit='Flux [kHz/cm^{2}]', y_tit=self.PhTit, y_off=1.4, x_range=self.Bins.FluxRange if logx else [0, 10100])
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
    def draw_currents(self, v_range=None, rel_time=False, averaging=True, with_flux=False, c_range=None, f_range=None, cunit='nA', **dkw):
        self.Currents.draw(rel_time=rel_time, v_range=v_range, averaging=averaging, with_flux=with_flux, c_range=c_range, f_range=f_range, cunit=cunit, **prep_kw(dkw, savename='Currents'))

    def draw_iv(self, show=True):
        g = self.Currents.draw_iv(show=False)
        self.Draw(g, 'IV', draw_opt='ap', logy=True, lm=.12, show=show)

    def draw_current_flux(self, c_range=None, fit=True, avrg=False, show=True):
        x, y = self.get_fluxes(avrg=avrg), [ufloat(0, 0) if c is None else c for c in self.get_values('', self.Analysis.get_current, pbar=False, avrg=avrg)]
        g = self.Draw.graph(x, y, title='Current vs. Flux', show=show, **self.get_x_args(draw=True))
        format_histo(g, y_tit='Current [nA]', y_range=choose(c_range, [0, max(y).n * 1.5]))
        g.Fit(Draw.make_f('fcf', 'pol1', .1, 5e5, pars=[.2, 1e-3], limits=[[.1, 5], [1e-6, 5e-3]]), f'q{"" if fit else "0"}')
        Draw.irradiation(irr2str(self.get_irradiation()))
        format_statbox(g, fit=fit, w=.3, form='.1e')
        self.Draw.save_plots('FluxCurrent', show=show)
        return g
    # endregion CURRENT
    # ----------------------------------------

    # ----------------------------------------
    # region TELESCOPE
    def draw_beam_current(self, bw=60, rel_t=True, show=True):
        h = TH1F('hr1', 'Beam Current of Run Plan {r}'.format(r=self.RunPlan), *self.get_raw_time_bins(bw))
        values = [hist_values(ana.Tel.draw_beam_current(bw, show=False, save=False)) for ana in self.Analyses]
        values = concatenate([append(run_values, run_values[-1]) for run_values in values])  # add last flux value for time bin between runs
        for i, value in enumerate(values, 1):
            h.SetBinContent(i, value.n)
            h.SetBinError(i, value.s)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', y_off=.85, fill_color=Draw.FillColor, stats=0, markersize=.3, t_ax_off=self.StartTime if rel_t else 0)
        self.Draw(h, 'AllBeamRate', show=show, draw_opt='hist', x=1.5, y=.75, lm=.065)

    @save_pickle('Prof', sub_dir='Flux', suf_args='[0]')
    def get_flux_prof(self, bw=5, _redo=False):
        x = [hist_values(ana.Tel.draw_flux(bw, show=False, prnt=False)) for ana in self.Analyses]
        x = concatenate([append(run_values, ufloat(0, 0)) for run_values in x])  # add extra zero for time bin between runs
        return self.Draw.distribution(x, self.get_raw_time_bins(bw), 'Flux Profile', **self.get_x_args(True), y_tit='Flux [kHz/cm^{2}]', show=False)

    def draw_flux(self, bw=5, rel_time=True, redo=False, **kwargs):
        h = self.get_flux_prof(bw, _redo=redo)
        return self.Draw.distribution(h, **prep_kw(kwargs, file_name='FluxProfile', **self.get_x_args(True, rel_time, draw=True), **Draw.mode(2), logy=True,
                                                   y_tit='Flux [kHz/cm^{2}]', stats=False))

    def draw_flux_ratio(self, avrg=False, **dkw):
        x, (r1, r2) = self.get_fluxes(avrg=avrg), [self.get_fluxes(pl, rel=True, avrg=avrg) for pl in [1, 2]]
        self.Draw.graph(x, r1 / r2, 'FluxRatio', **prep_kw(dkw, y_tit='Flux Plane1/Plane2', **self.get_x_args(draw=True), **Draw.mode(2, bm=.23), **dfname('FluxRatio', avrg)))
    # endregion TELESCOPE
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL MAP
    def draw_signal_map(self, fid=False, res=.7, redo=False, square=False, scale=False, **kwargs):
        pickle_path = self.get_pickle_path('SM', make_suffix(self.FirstAnalysis, res, None, fid, 0, square), 'Maps')
        histos = self.get_plots('signal maps', self.Analysis.get_signal_map, res=res, square=square, _redo=redo, fid=fid, picklepath=pickle_path)
        for h in histos[1:]:
            histos[0].Add(h)
        histos[0].Scale(1 / self.get_pulse_height()[0].n) if scale else do_nothing()
        rz = array([histos[0].GetMinimum(), histos[0].GetMaximum()]) * 1 / self.get_pulse_height()[0].n if scale else None
        h = self.Draw.prof2d(histos[0], title='Cumulative Pulse Height Map', **prep_kw(kwargs, centre=4, pal=53, z_range=rz, z_tit='Relative Pulse Height' if scale else None))
        self.Draw.save_plots('CumSignalMap2D', **kwargs)
        return h

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
        g = Draw.make_tgraph(self.get_x_var(vs_time), spreads, title='Relative Spread')
        format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='Relative Spread [%]', y_off=1.2, t_ax_off=self.get_tax_off(vs_time, rel_time))
        self.Draw(g, 'RelativeSpread', lm=.12, logx=not vs_time, show=show)

    def draw_sm_std(self, vs_time=False, redo=False, show=True):
        y_values = list(self.get_sm_stds(redo).values())
        x_values = self.get_times() if vs_time else self.get_fluxes()
        g = Draw.make_tgraph(x_values, y_values, title='STD of the Signal Map')
        format_histo(g, x_tit=self.get_x_tit(vs_time), y_tit='rel. STD', y_off=1.3, t_ax_off=0 if vs_time else None)
        self.Draw(g, 'STDSigMap', show, lm=.12, logx=not vs_time)
    # endregion SIGNAL MAP
    # ----------------------------------------

    # ----------------------------------------
    # region BEAM PROFILE
    def beam_pars(self, avrg=False, redo=False):
        x = self.get_values('beam pars', self.Analysis.beam_pars, picklepath=self.get_pickle_path('BeamPars', sub_dir='Tracks', dut=''), _redo=redo)
        d = concatenate(x.T)  # = [mx, my, sx, sy]
        return array([self.get_flux_average(i) for i in d] if avrg else d)

    def draw_beam_centres(self, avrg=False, redo=False, **dkw):
        x, y = self.get_x_var(avrg=avrg), self.beam_pars(avrg, redo)[:2]
        g = [self.Draw.graph(x, iy, show=False, y_tit='Beam Centre [mm]') for iy in y]
        return self.Draw.multigraph(g, 'BeamCentre', ['y direction', 'x direction'], **prep_kw(dkw, **self.get_x_args(), gridy=True, file_name='BeamCentre'))

    def draw_beam_width(self, avrg=False, redo=False, **dkw):
        x, y = self.get_x_var(avrg=avrg), self.beam_pars(avrg, redo)[2:]
        g = [self.Draw.graph(x, iy, show=False, y_tit='Beam Width [mm]') for iy in y]
        return self.Draw.multigraph(g, 'BeamWidth', ['x direction', 'y direction'], **prep_kw(dkw, **self.get_x_args(), gridy=True, file_name='BeamWidth'))

    def draw_beam_profiles(self, mode='x', show=True):
        histos = self.get_values('beam profiles in {}'.format(mode), self.Analysis.draw_beam_profile, show=False, prnt=False, fit=False, mode=mode)
        leg = Draw.make_legend(nentries=self.NRuns)
        s = THStack('sbp', 'AllBeamProfiles{mod}'.format(mod=mode.title()))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
            format_histo(h, lw=2, stats=0, normalise=True, line_color=self.Draw.get_color(self.NRuns), sumw2=False, fill_style=4000)
            s.Add(h)
            leg.AddEntry(h, '{0:6.2f} kHz/cm^{{2}}'.format(flux.n), 'l')
        self.Draw(s, 'AllBeamProfiles{mod}'.format(mod=mode.title()), draw_opt='nostack', show=show, leg=leg)
    # endregion BEAM PROFILE
    # ----------------------------------------

    # ----------------------------------------
    # region TRACKS
    def draw_chi2(self, mode=None, show=True):
        mod_str = '' if mode is None else mode
        histos = self.get_values('chi squares {}'.format(mod_str), self.Analysis.draw_chi2, show=False, prnt=False, mode=mode)
        cut_value = self.FirstAnalysis.Cut.get_chi2(choose(mod_str, 'x'))
        cuts = array([get_quantile(h, [.99, cut_value / 100]) for h in histos])  # .99 for drawing range
        s = THStack('hx2', '#chi^{{2}}{}'.format(' in {}'.format(mod_str) if mod_str else ''))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
            format_histo(h, stats=0, color=self.Draw.get_color(self.NRuns), lw=2, normalise=True, sumw2=False, fill_style=4000, fill_color=4000)
            s.Add(h)
        format_histo(s, x_tit='#chi^{2}', y_tit='Fraction of Events [%]', y_off=1.5, draw_first=True, x_range=ax_range(0, max(cuts), fh=.4))
        legend = self.make_flux_legend(histos)
        self.Draw(s, '', show, lm=.15, draw_opt='nostack', leg=legend)
        if mod_str:
            line = Draw.vertical_line(min(cuts), -1e9, 1e9, color=2, style=2)
            legend.AddEntry(line, 'cut: {}% q'.format(cut_value), 'l')
            histos[0].GetListOfFunctions().Add(line)
            histos[0].GetListOfFunctions().Add(legend)
        self.Draw.save_plots('AllChi2{mod}'.format(mod=mod_str.title()))
        return s

    def draw_all_chi2s(self, show=True):
        stacks = [self.draw_chi2(mode, show=False) for mode in [None, 'x', 'y']]
        c = Draw.canvas('AllChiSquares', x=2, y=.66, divide=3, show=show)
        for i, s in enumerate(stacks, 1):
            self.Draw(s, draw_opt='nostack', canvas=c.cd(i), show=show)
        self.Draw.save_plots('AllChiSquares', show=show, canvas=c)

    def draw_angle(self, mode='x', show=True):
        histos = self.get_values('angle distribution {}'.format(mode), self.Analysis.draw_angle, mode=mode, show=False, prnt=False)
        legend = self.make_flux_legend(histos)
        s = THStack('has', 'Track Angles in {mode}'.format(mode=mode.title()))
        for i, (h, flux) in enumerate(zip(histos, self.get_fluxes())):
            format_histo(h, stats=0, color=self.Draw.get_color(self.NRuns), normalise=True, sumw2=False, fill_color=4000, fill_style=4000)
            s.Add(h)
        histos[0].GetListOfFunctions().Add(legend)
        format_histo(s, x_tit='Angle [deg]', y_tit='Fraction of Events [%]', y_off=1.9, draw_first=True, x_range=[-4, 4])
        self.Draw(s, 'AllTrackAngles{mod}'.format(mod=mode.title()), show, lm=.15, draw_opt='nostack', leg=legend)
        return s

    def draw_both_angles(self, show=True):
        stacks = [self.draw_angle(mode, show=False) for mode in ['x', 'y']]
        c = Draw.canvas('AllAngles', x=1.5, y=.75, divide=2, show=show)
        for i, s in enumerate(stacks, 1):
            self.Draw(s, draw_opt='nostack', canvas=c.cd(i), show=show, lm=.15)
        self.Draw.save_plots('AllChiAngles', show=show, canvas=c)

    def draw_mean_angles(self, mode='x', vs_time=True, rel_time=True, show=True, redo=False):
        values = [value[0] for value in self.get_values('mean angles {}'.format(mode), self.Analysis.get_mean_angle, mode=mode, redo=redo)]
        g = Draw.make_tgraph(self.get_x_var(vs_time, rel_time), values, title='Mean Track Angles in {}'.format(mode.title()))
        format_histo(g, y_tit='Mean Angle [deg]', y_off=1.3, **self.get_x_args(vs_time, rel_time))
        self.Draw(g, 'MeanTrackAngle'.format(mode.title()), show, logx=not vs_time)
    # endregion TRACKS
    # ----------------------------------------

    # ----------------------------------------
    # region CUTS
    def draw_purity(self, avrg=False, redo=False, **dkw):
        x, y = self.get_fluxes(avrg=avrg), self.get_values('purity', self.Analysis.purity, picklepath=self.get_pickle_path('Purity', sub_dir='Cuts'), avrg=avrg, redo=redo)
        return self.Draw.graph(x, y * 100, **prep_kw(dkw, title='Purity', **self.get_x_args(draw=True), y_tit='Purity [%]', y_range=[0, 105], file_name='CutPurity', gridy=True))
    # endregion CUTS
    # ----------------------------------------

    # ----------------------------------------
    # region GENERATE PLOTS
    def save_full_plots(self, name, f, ftype='pdf', **kwargs):
        for i, plot in enumerate(self.get_plots(f.__name__, f, show=False, **{key: value for key, value in kwargs.items() if 'redo' in key})):
            self.Draw.save_full(plot, f'{name}-{i}', ftype=ftype, **kwargs)

    def draw_run_currents(self):
        self.get_values('currents', self.Analysis.get_current)
        self.get_values('currents', self.Analysis.draw_current, show=False)

    def draw_chi2s(self):
        self.get_values('chi2s', self.Analysis.draw_chi2s, show=False, prnt=False)

    def draw_angles(self):
        self.get_values('angles', self.Analysis.draw_angles, show=False, prnt=False)

    def draw_occupancies(self):
        self.get_values('occupancies', self.Analysis.draw_occupancies, show=False, prnt=False, cluster=True)
    # endregion GENERATE PLOTS
    # ----------------------------------------

    def draw_efficiencies(self, avrg=False, t=False, **dkw):
        x, y = self.get_x_var(t, avrg=avrg), self.get_efficiencies()
        self.Draw.graph(x, y, 'Efficiencies', **prep_kw(dkw, **self.get_x_args(t, draw=True), y_tit='Effciency [%]'))

    def draw_eff_vs_current(self, **dkw):
        x, y = self.get_currents(), self.get_efficiencies()
        self.Draw.graph(x, y, 'Eff vs Current', **prep_kw(dkw, x_tit='Current [nA]', y_tit='Effciency [%]', file_name='EffCurrent'))

    def draw_gain(self, fit=True, irr=False, **dkw):
        v = self.get_pulse_heights()
        x, y = self.get_x_var(vs_time=not irr, vs_irrad=irr), add_perr(v[1::2] / v[::2], .005)
        g = self.Draw.graph(x[::2], y, y_tit='Gain Ratio')
        g.Fit('pol0', 'qs') if fit else do_nothing()
        x_args = self.get_x_args(not irr, draw=True, vs_irrad=irr, **{} if irr else {'x_tit': 'Time [YY/MM]', 'tform': '%y/%m'})
        return self.Draw(g, **prep_kw(dkw, file_name=fname('GainRatio', i=irr), **x_args, stats=set_statbox(fit=fit, form='.3f')))


def fname(n, avrg=False, t=False, i=False):
    return f'{n}{"Time" if t else ""}{"Avr" if avrg else ""}{"Irr" if i else ""}'


def dfname(n, avrg=False, t=False):
    return {'file_name': f'{n}{"Time" if t else ""}{"Avr" if avrg else ""}'}


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
