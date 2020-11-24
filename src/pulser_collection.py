# --------------------------------------------------------
#       Pulser sub class for AnalysisCollection
# created by Michael Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.sub_ana_collection import SubCollection
from helpers.draw import *
from src.pulser import PulserAnalysis


class PulserCollection(SubCollection):

    def __init__(self, ana_collection):
        super().__init__(ana_collection)
        if self.Ana.FirstAnalysis.Tree.Hash():
            self.Analysis = self.Ana.FirstAnalysis.Pulser
        self.Type = '?' if 'pulser' not in self.Ana.FirstAnalysis.Run.Info else self.Ana.FirstAnalysis.Run.Info['pulser']

    # ----------------------------------------
    # region GET
    def get_analyses(self, runs=None):
        return [ana.Pulser for run, ana in self.Analyses.items() if run in choose(runs, self.Ana.Runs)]

    def get_pulse_heights(self, corr=True, beam_on=True, sigma=False, redo=False):
        pickle_path = self.Analysis.make_simple_pickle_path('HistoFit', '{}_{}'.format(int(corr), 'Beam{}'.format(int(beam_on))), run='{}')
        return self.get_values('pulser pulse heights', PulserAnalysis.get_pulse_height, corr=corr, beam_on=beam_on, redo=redo, par=[1, 2][sigma], picklepath=pickle_path)

    def get_sigmas(self, corr=True, beam_on=True, redo=False):
        return self.get_pulse_heights(corr, beam_on, sigma=True, redo=redo)

    def get_pedestals(self, sigma=False, beam_on=True, redo=False):
        pickle_path = self.Analysis.make_simple_pickle_path('Pedestal', str(int(beam_on)), run='{}')
        return self.get_values('pulser pedestals', PulserAnalysis.get_pedestal, sigma=sigma, beam_on=beam_on, picklepath=pickle_path, redo=redo)

    def get_noise(self):
        return self.get_pedestals(sigma=True)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def draw_pulse_height(self, bin_size=30000, redo=False, show=True, rel_t=True):
        """Draw time evolution of the pulse height"""
        f, picklepath = PulserAnalysis.draw_pulse_height, self.Analysis.make_simple_pickle_path('PHT', 'PulserBeamOn{}'.format(self.Analysis.Bins(bin_size)), run='{}')
        x, y = get_hist_vecs(array(self.get_plots('pulser pulse heights', f, show=False, prnt=False, picklepath=picklepath, redo=redo, bin_size=bin_size), object)[:, 0], err=False)
        self.Draw.profile(x, y, self.get_time_bins(bin_size), 'Pulser Pulse Height Trend', **self.get_x_args(True, rel_t), y_tit='Pulser Pulse Height [mV]', w=2, show=show, stats=0,
                          y_range=ax_range(y, 0, .3, .3))

    def draw_pulse_heights(self, sigma=False, vs_time=False, scaled=False, corr=True, beam_on=True, show_flux=True, legend=True, fit=False, redo=False, show=True):
        x, y = self.get_x(vs_time), self.get_pulse_heights(corr, beam_on, sigma, redo)
        y /= mean(y) if scaled else 1
        marker, colors, ms = [20, 22, 23], [self.Draw.get_color(10, 6)] + [int(Draw.Colors[0])] * 2, [1, 2, 2]
        graphs = [self.Draw.make_tgrapherrors(ix, iy, markersize=ms[i], color=colors[i], marker=marker[i]) for i, (ix, iy) in enumerate([(x, y), ([x[0].n], [y[0].n]), ([x[-1].n], [y[-1].n])])]
        mg = self.Draw.multigraph(graphs[:(1 if vs_time else 3)], 'Pulser Pulse Height', ['data', 'first', 'last'] if legend else None, **self.get_x_draw(vs_time), color=False, show=show, lm=.12)
        format_histo(mg, **self.get_x_args(vs_time), y_range=ax_range(y, 0, .5, 1), y_tit='Pulser Pulse Height [mV]')
        Draw.info('Pulser Type: {}'.format(self.Type))
        if vs_time and show_flux:
            for ix, iy, flux in zip(x, y, self.get_fluxes()):
                mg.GetListOfGraphs()[0].GetListOfFunctions().Add(Draw.tlatex(ix.n, iy.n + iy.s * 1.2, '{:1.0f}'.format(flux.n), color=1, align=21, size=.02))
        if fit:
            graphs[0].Fit('pol0', 'q')
            format_statbox(graphs[0], fit=True, x2=.41)
        return mg

    def draw_scaled_pulse_heights(self, sigma=False, vs_time=False, corr=True, beam_on=True, show_flux=True, legend=True, fit=False, redo=False, show=True):
        return self.draw_pulse_heights(sigma, vs_time, True, corr, beam_on, show_flux, legend, fit, redo, show)

    def draw_sigmas(self, vs_time=False, corr=True, beam_on=True, show_flux=True, legend=True, fit=False, redo=False, show=True):
        return self.draw_pulse_heights(True, vs_time, False, corr, beam_on, show_flux, legend, fit, redo, show)

    def draw_distributions(self, show=True, corr=True):
        histos = self.get_plots('pulser distributions', PulserAnalysis.draw_distribution, show=False, corr=corr, picklepath=self.Analysis.make_simple_pickle_path('Disto', '{}12'.format(int(corr))))
        self.Draw.stack(histos, 'Pulser Distributions', self.Ana.get_flux_strings(), scale=True, show=show)

    def compare_pulse_heights(self, sigma=False, scaled=False, show=True):
        graphs = [self.draw_pulse_heights(sigma, scaled=scaled, show=False, beam_on=ib, corr=ic).GetListOfGraphs()[0] for ib, ic in [(1, 1), (1, 0), (0, 1), (0, 0)]]
        self.Draw.multigraph(graphs, 'Pulse Pulse Heights', ['BeamOn w/ corr.', 'BeamOn w/0 corr.', 'BeamOff w/ corr.', 'BeamOff w/0 corr.'], show=show, logx=True)
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region PEDESTAL
    def draw_pedestals(self, vs_time=False, rel_time=False, sigma=False, beam_on=True, redo=False, show=True):
        x, y = self.get_x(vs_time, rel_time), self.get_pedestals(sigma, beam_on, redo)
        tit = 'Pulser {}'.format('Noise' if sigma else 'Pedestal')
        return self.Draw.graph(x, y, tit, **self.get_x_args(vs_time, rel_time, draw_args=True), y_tit='{} [mV]'.format(tit), show=show)

    def compare_pedestals(self, sigma=False, show=True):
        graphs = [self.Ana.draw_pedestals(sigma=sigma, show=False)] + [self.draw_pedestals(sigma=sigma, beam_on=i, show=False) for i in [1, 0]]
        self.Draw.multigraph(graphs, 'Pedestal Comparison', ['Signal', 'Pulser', 'Beam Off'], logx=True, show=show)
    # endregion PEDESTAL
    # ----------------------------------------

    # ----------------------------------------
    # region RATE
    def draw_rate(self, bin_size=50, cut=None, rel_t=True, show=True):
        x, y = get_hist_vecs(self.get_plots('pulser rate', PulserAnalysis.draw_rate, show=False, prnt=False, bin_size=bin_size, cut=cut), err=False)
        self.Draw.profile(x, y, self.Ana.get_raw_time_bins(bin_size), 'Pulser Rate', **self.get_x_args(True, rel_t), y_tit='Pulser Rate [%]', w=2, show=show, stats=0,  y_range=ax_range(y, 0, .3, .3))

    def draw_rates(self, vs_time=False, redo=False, show=True):
        x, y = self.get_x(vs_time), self.get_values('pulser rates', PulserAnalysis.get_rate, prnt=False, redo=redo, picklepath=self.Analysis.make_simple_pickle_path('Rate', run='{}'))
        self.Draw.graph(x, y, 'Pulser Rates', y_tit='Pulser Rate [%]', **self.get_x_args(vs_time, draw_args=True), show=show)

    def draw_fractions(self, vs_time=False, show=True):
        x, y = self.get_x(vs_time), self.get_values('', PulserAnalysis.get_real_fraction, pbar=False)
        self.Draw.graph(x, y, 'Real Pulser Fraction', y_tit='Pulser Fraction [%]', **self.get_x_args(vs_time, draw_args=True), show=show)

    def make_rate_plots(self):
        self.get_values('pulser rates', PulserAnalysis.draw_rate, show=False, prnt=False)
    # endregion RATE
    # ----------------------------------------

    def draw_peak_times(self, sigma=False, vs_time=False, show=True, redo=False):
        x, y = self.get_x(vs_time), self.get_values('pulser peak times', PulserAnalysis.get_peak_time, redo=redo, sigma=sigma, picklepath=self.Analysis.make_simple_pickle_path('PeakTime', run='{}'))
        self.Draw.graph(x, y, 'Pulser Peak Times', y_tit='Pulser Peak Time [ns]', **self.get_x_args(vs_time, draw_args=True), show=show, lm=.13, y_off=1.6)
