# --------------------------------------------------------
#       Analysis collection child for pad detectors
# revised on Nov 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from pad.pulser_collection import PulserCollection
from src.analysis_collection import *
from pad.analysis import PadAnalysis, in1d


class PadCollection(AnalysisCollection):
    """ Analysis of the various runs of a single runplan. """

    def __init__(self, run_plan, dut_nr, test_campaign=None, load_tree=True, verbose=False):
        AnalysisCollection.__init__(self, run_plan, dut_nr, test_campaign, load_tree, verbose)
        # Sub Classes
        self.Pulser = PulserCollection(self)

    # ----------------------------------------
    # region INIT
    @staticmethod
    def load_dummy():
        return PadAnalysis
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region RESULTS
    def draw_all(self, redo=False):
        t0 = self.info('Generate all plots ... ')
        self.draw_pulse_heights(show=False, redo=redo)
        self.Pulser.draw_pulse_heights(show=False, redo=redo)
        self.draw_pedestals(show=False, redo=redo)
        self.draw_pedestals(show=False, sigma=True, redo=redo)
        self.Pulser.draw_pedestals(show=False, redo=redo)
        self.Pulser.draw_pedestals(show=False, redo=redo, sigma=True)
        if 'voltage' not in self.Type.lower():
            self.draw_scaled_pulse_heights(show=False)
            self.draw_scaled_pulse_heights(show=False, vs_time=True)
        self.draw_currents(show=False, draw_opt='al')
        self.draw_flux(show=False)
        self.draw_ph_currents(show=False)
        self.draw_signal_distributions(show=False, redo=redo)
        self.draw_signal_maps(redo=redo)
        self.draw_hitmaps(redo=redo)
        self.draw_run_currents()
        self.draw_chi2s()
        self.draw_angles()
        self.draw_occupancies()
        self.get_values('timing', PadAnalysis.draw_timing)
        print_elapsed_time(t0)

    @quiet
    def save_coll_plots(self):
        super(PadCollection, self).save_coll_plots()
        self.draw_pedestals(show=False)
        self.draw_pedestals(show=False, sigma=True)
        self.Pulser.draw_pulse_heights(show=False)
        self.Pulser.draw_sigmas(show=False)
    # endregion RESULTS
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_pulse_heights(self, bin_width=None, redo=False, runs=None, corr=True, err=True, pbar=True, avrg=False, peaks=False, flux_sort=False):
        picklepath = None if peaks else self.get_pickle_path('Fit', int(corr), 'PH')
        pbar = False if peaks else pbar
        x = self.get_values('pulse heights', self.Analysis.get_pulse_height, runs, pbar, avrg, picklepath, bin_size=bin_width, redo=redo, corr=corr, peaks=peaks, flux_sort=flux_sort)
        return array([ufloat(ph.n, sqrt(ph.s ** 2 + ((self.get_sys_error() if err else 0) * ph.n) ** 2)) for ph in x])

    def get_pedestals(self, runs=None, sigma=False, flux_sort=False, avrg=False, redo=False):
        picklepath = self.FirstAnalysis.Pedestal.make_simple_pickle_path(suf='AllCuts_ab2', run='{}')
        return self.get_values('pedestals', self.Analysis.get_pedestal, runs, par=[1, 2][sigma], redo=redo, flux_sort=flux_sort, picklepath=picklepath, avrg=avrg)

    def get_peak_fluxes(self, corr=True, avrg=False, rel=False):
        values = self.get_values('peak Flux', f=self.Analysis.get_peak_flux, pbar=False, prnt=False, avrg=avrg, corr=corr)
        return array([ufloat(v.n, v.n * .01) for v in values]) if rel else values

    def get_additional_peaks(self, start=None, end=None):
        picklepath = self.make_simple_pickle_path('NAdd', '' if start is None else '{}_{}'.format(start, end), 'Peaks', '{}')
        return self.get_values('number of additional peaks', f=self.Analysis.get_n_peaks, picklepath=picklepath, start_bunch=start, end_bunch=end)

    def get_pulser_pulse_heights(self, avrg=False, redo=False):
        return self.Pulser.get_pulse_heights(avrg=avrg, redo=redo)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL/PEDESTAL
    def draw_ph_currents(self, yr=.06, show=True):
        ph, pul, t = self.get_pulse_heights(err=False), self.Pulser.get_pulse_heights(), self.get_times()
        graphs, cur = [self.Draw.make_tgrapherrors(t, y / mean(y), color=col) for y, col in [(ph, 2), (pul, 4)]], self.Currents.draw(show=False)
        mg = self.Draw.multigraph(graphs, 'signal', bm=0, color=False, show=False, y_tit='Pulse Height [mV]')
        c = Draw.canvas('Pulse Height Current', w=2, h=1.2, transp=True, show=show)
        self.Draw.Legends.draw(c.cd(), all_pads=False)
        self.Draw(mg, info_leg=False, canvas=Draw.tpad(pos=[0, .6, 1, 1], gridx=True, gridy=True, transparent=True, bm=0, tm=0, lm=.08))
        Draw.legend(graphs, ['signal', 'pulser'], 'p', y2=.97, scale=2, w=.12, fix=True)
        self.Draw(cur, gridy=True, info_leg=False, canvas=Draw.tpad(pos=[0, .1, 1, .6], gridx=True, gridy=True, transparent=True, c=c, bm=.2, tm=0, lm=.08))
        x_range = ax_range(cur.GetX()[0], cur.GetX()[cur.GetN() - 1], .1, .1)
        format_histo(mg, y_range=[1 - yr, 1 + yr], **self.get_x_args(True, x_range=x_range), lab_size=.08, tit_size=.085, center_y=True, y_off=.5, ndivy=504)
        format_histo(cur, **self.get_x_args(True, x_range=x_range), lab_size=.08, tit_size=.085, center_y=True, y_off=.5)

    def draw_pedestals(self, vs_time=False, sigma=False, redo=False, avrg=False, **kwargs):
        x, y = self.get_x_var(vs_time, avrg=avrg), self.get_pedestals(redo=redo, sigma=sigma, avrg=avrg)
        y_tit, tit = f'{["Pedestal", "Noise"][sigma]}', f'{["Pedestal", "Noise"][sigma]} {["Flux", "Time"][vs_time]} {"Avrg" if avrg else ""}'
        kwargs = prep_kw(kwargs, title=tit, y_tit=f'{y_tit} [mV]', **self.get_x_args(vs_time, draw=True), color=810, file_name=tit.replace(' ', ''), markersize=.8)
        return self.Draw.graph(x, y, **kwargs)

    def draw_noise(self, vs_time=False, redo=False, show=True):
        return self.draw_pedestals(vs_time, True, redo, show)

    def draw_ped_spread(self, redo=False, show=True):
        values = self.get_pedestals(redo=redo, flux_sort=True)
        values = [v - mean(lst) for lst in split(values, self.get_flux_splits(show=False)) for v in lst if len(lst) > 1]  # subtract flux group mean
        self.Draw.distribution(values, Bins.make(-.5, .5, n=20), 'Relative Pedestal Spread', x_tit='Relative Pedestal', file_name='PedestalSpread', show=show)

    def draw_noise_spread(self, redo=False, show=True):
        x = self.get_pedestals(sigma=True, redo=redo)
        return self.Draw.distribution([(v - mean(x)).n for v in x], Bins.make(-.3, .3, n=20), 'Relative Noise Spread', file_name='NoiseSpread', x_tit='Relative Noise', show=show)

    def draw_snrs(self, vs_time=False, show=True):
        x, y = self.get_x_var(vs_time), self.get_values('', self.Analysis.calc_snr, pbar=False)
        self.Draw.graph(x, y, title='Signal to Noise Rations', y_tit='SNR', **self.get_x_args(vs_time), show=show, y_range=ax_range(y, rnd=True), logx=not vs_time)

    def compare_signal_vs_peak_height(self, i0=0, i1=-1, ym=.055, cft=False, show=True, redo=False):
        """draws the pulse height ratio of two runs vs the peak time. """
        def f():
            func = self.Analysis.draw_signal_vs_cft if cft else self.Analysis.draw_ph_peaktime
            (x0, y0), y1 = get_hist_vecs(func(self.Analyses[i0], show=False)), get_hist_vec(func(self.Analyses[i1], show=False))
            return x0, y1 / where(y0 == 0, 1e10, y0)  # make ratio of zero entries 0
        x, y = do_pickle(self.make_simple_pickle_path('SigPeakRatio', int(cft), sub_dir='Peaks', dut='{}{}'.format(i0, i1)), f, redo=redo)
        flux0, flux1 = [make_flux_string(self.Analyses[i].get_flux().n) for i in [i0, i1]]
        self.Draw.graph(x, y, title='Signal Ratio of {} and {}'.format(flux0, flux1), x_tit='Signal Peak Time [ns]', y_tit='Signal Ratio', y_range=array([-ym, ym]) + 1, gridy=True, show=show)
        return mean_sigma(y[y > .1])[0]

    def compare_signal_vs_cft(self, i0=0, i1=-1, ym=.055, show=True, redo=False):
        return self.compare_signal_vs_peak_height(i0, i1, ym, True, show, redo)

    def compare_all_sig_vs_peakheight(self, ym=.05, show=True):
        y = [self.compare_signal_vs_peak_height(0, i, show=False) for i in range(1, self.NRuns)]
        self.Draw.graph(self.get_fluxes()[1:], y, title='Signal Ratio Vs Flux', y_tit='Signal Ratio', y_range=array([-ym, ym]) + 1, **self.get_x_args(False), gridy=True, show=show)

    def find_sm_correlation(self, sm1, sm2):
        return self.FirstAnalysis.find_best_sm_correlation([sm1, sm2])
    # endregion SIGNAL/PEDESTAL
    # ----------------------------------------

    # ----------------------------------------
    # region CUTS
    def draw_bucket_ratio(self, fit=True, avrg=False, redo=False, **kwargs):
        x, y = self.get_fluxes(), 100 * self.get_values('bucket ratios', PadAnalysis.get_bucket_ratio, picklepath=self.get_pickle_path('Ratio', '1', 'Bucket'), all_cuts=True, _redo=redo, avrg=avrg)
        g = self.Draw.graph(x, y, 'Bucket Ratio', y_tit='Fraction of Bucket Events [%]', **prep_kw(kwargs, y_range=[0, max(y).n * 1.5], markersize=.7, **self.get_x_args(draw=True)))
        if fit:
            g.Fit(self.Draw.make_f(None, 'pol1', 0, 4e7, pars=[0, 1e-3], fix=0), 'qs')
            format_statbox(g, fit=True and 'stats' not in kwargs, form='.2e')
        return g

    def draw_bucket_ph(self, show=True, redo=False):
        pickle_path = self.make_simple_pickle_path('Fit', '{}_1_b2_nobucket'.format(Bins.Size), 'Ph_fit', '{}')
        x, y = self.get_fluxes(), self.get_values('bucket ph', self.Analysis.get_bucket_pulse_heights, picklepath=pickle_path, redo=redo).T
        graphs = [Draw.make_tgrapherrors(x, iy) for iy in y]
        mg = self.Draw.multigraph(graphs, 'Bucket Pulse Heights', ['no bucket', 'w/o thresh', 'with thresh'], y_tit='Pulse Height [mV]', show=show, logx=True)
        format_histo(mg, **self.get_x_args(), y_range=ax_range(concatenate(y), 0, .3, .3))

    def get_bucket_tp_ratios(self, avrg=False, redo=False):
        return self.get_values('bucket tp ratios', PadAnalysis.get_bucket_tp_ratio, picklepath=self.get_pickle_path('TPRatio', '1', 'Bucket'), _redo=redo, all_cuts=True, avrg=avrg)

    def draw_bucket_tp_ratio(self, avrg=False, redo=False, **kwargs):
        x, y = self.get_fluxes(), self.get_bucket_tp_ratios(avrg, redo)
        return self.Draw.graph(x, y, 'Bucket Trigger Phase Ratio', y_tit='Percentage of Trigger Phase', **self.get_x_args(draw=True), **kwargs)

    @save_pickle('TPRatio')
    def get_bucket_tpr(self, _redo=False):
        x = concatenate([ana.get_tree_vec('trigger_phase', ana.Cut.get_bucket(all_cuts=True)) for ana in self.Analyses])
        return calc_eff(values=in1d(x, self.FirstAnalysis.Cut.get_trigger_phases()))

    def create_bucket_estimate(self):
        print(f'bucket scale = {FitRes(self.draw_bucket_ratio(show=False).GetListOfFunctions()[0])[1]:.2e}')
        print(f'bucket tpr = {eff2u(self.get_bucket_tpr() / 100):.4f}')
    # endregion CUTS
    # ----------------------------------------

    # ----------------------------------------
    # region TIMING
    def draw_peak_timings(self, vs_time=False, redo=False, show=True):
        """ Shows the means of the signal peak distributions. """
        x, y = self.get_x_var(vs_time), self.get_values('peak timings', self.Analysis.get_peak_timing, redo=redo, picklepath=self.make_simple_pickle_path('PeakVals', sub_dir='Timing', run='{}'))
        self.Draw.graph(x, y, 'Peak Times', **self.get_x_args(vs_time), y_tit='Peak Time [ns]', logx=not vs_time, show=show)
    # endregion TIMING
    # ----------------------------------------

    # ----------------------------------------
    # region PEAKS
    def draw_n_peaks(self, logy=True, show=True):
        fluxes, n_events = self.get_fluxes(), self.get_n_events()
        n_peaks = array([mean(ana.Peaks.draw_additional(show=False)) for ana in self.get_analyses()]) / n_events
        g = self.Draw.graph(fluxes, n_peaks / min(n_peaks), 'Number of Additional Peaks', **self.get_x_args(False), y_tit='Number of Additional Peaks', logy=logy, logx=True, show=show)
        g.Fit('pol1', 'qs')
        format_statbox(g, fit=True, center_x=True)

    def draw_ph_vs_peaks(self, show=True):
        peak_heights = array([ana.Peaks.find_additional(scale=True, show=False) for ana in self.get_analyses()]) / self.get_fluxes()
        peak_heights /= mean(peak_heights)
        pulse_heights = self.get_pulse_heights()
        g = Draw.make_tgrapherrors('gpr', 'Pulse Height vs Normalised Peak Height', y=pulse_heights, x=peak_heights)
        format_histo(g, y_tit='Pulse Height', x_tit='Normalised Peak Height')
        self.Draw(g, show=show, lm=.12)

    def draw_fluxes(self, show=True):
        g = [self.Draw.make_tgrapherrors(arange(1, self.NRuns + 1), v) for v in [self.get_fluxes(pl, 1, 1) for pl in [1, 2, None]] + [self.get_peak_fluxes()]]
        self.Draw.multigraph(g, 'Flux Comparison', ['Plane 1', 'Plane 2', 'Mean 1 & 2', 'Peaks'], x_tit='Run', y_tit=self.get_x_tit(), show=show)

    def draw_flux_ratios(self, show=True):
        peak_flux = self.get_peak_fluxes(rel=True)
        g = [self.Draw.make_tgrapherrors(self.get_fluxes(), v / peak_flux) for v in [self.get_fluxes(pl, 1, 1, rel=True) for pl in [1, 2, None]]]
        mg = self.Draw.multigraph(g, 'Plane/Peak Flux Ratios', ['Plane 1', 'Plane 2', 'Mean 1 & 2'], y_tit='Plane/Peak Flux Ratio', show=show, **self.get_x_draw())
        format_histo(mg, **self.get_x_args(x_tit='Mean Plane Flux [kHz/cm^{2}]'))

    def compare_fluxes(self, plane=None, logy=True, corr=True, avrg=False, y_range=None, show=True):
        x, y = self.get_peak_fluxes(corr, avrg, rel=True), self.get_fluxes(plane, corr, full_size=True, avrg=avrg, rel=True)
        tit, xtit, ytit = 'FAST-OR Flux vs Peak Flux', 'Peak Flux [kHz/cm^{2}]', f'{"Mean Plane" if plane is None else f"Plane {plane}"} Flux [kHz/cm^{{2}}]'
        g = self.Draw.graph(x, y, tit, **self.get_x_args(x_tit=xtit, draw=True), y_range=choose(y_range, Bins.FluxRange), y_tit=ytit, logy=logy, show=show)
        g.Fit('pol1', 'qs')
        format_statbox(g, fit=True, center_x=True, form='.2f')

    def compare_tu_fluxes(self, logy=True, show=True):
        x, y = self.get_fluxes(), [ana.Run.Flux for ana in self.get_analyses()]
        g = self.Draw.graph(x, y, 'TU Flux vs Logged Flux', **self.get_x_args(), y_range=Bins.FluxRange, y_tit='Logged Flux [kHz/cm^{2}]', logx=True, logy=logy, show=show)
        g.Fit('pol1', 'qs')
        format_statbox(g, fit=True, center_x=True)

    def draw_bunch_systematics(self, bunch=0, show=True):
        """draws the number of additional peaks in bunch <n> relative to the average number of additional peaks per bunch vs. flux"""
        x, y, y0 = self.get_fluxes(), self.get_additional_peaks(bunch, bunch + 1), self.get_additional_peaks() / self.FirstAnalysis.Peaks.NBunches
        self.Draw.graph(x, y / y0, 'Bunch Systematics',  **self.get_x_args(), show=show, logx=True, y_tit='N Peaks in Bunch {} / Average Peaks per Bunch'.format(bunch))
    # endregion PEAKS
    # ----------------------------------------
