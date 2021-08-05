#!/usr/bin/env python
# --------------------------------------------------------
#       Main class for Rate Pixel Analysis
# created some time in 2016 by D. Sanz (sandiego@phys.ethz.ch), maintained by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------


from numpy import insert, sum

from pixel.calibration import Calibration
from pixel.cut import PixCut
from pixel.run import PixelRun
from src.dut import Plane
from src.dut_analysis import *


class PixAnalysis(DUTAnalysis):
    def __init__(self, run_number, dut, test_campaign=None, load_tree=True, verbose=False, prnt=True):

        DUTAnalysis.__init__(self, run_number, dut, test_campaign, load_tree, verbose, prnt)

        # Main
        self.N = self.DUT.Number + self.Run.NTelPlanes - 1

        if self.Tree.Hash():
            self.Calibration = Calibration(self)

        self.print_finished(prnt=prnt)

    # ----------------------------------------
    # region INIT
    @staticmethod
    def init_run(run_number, testcampaign, load_tree, verbose):
        return PixelRun(run_number, testcampaign, load_tree, verbose)

    def init_cut(self):
        return PixCut(self)

    def update_config(self):
        self.Config.read(join(self.Dir, 'config', self.TCString, 'PixelConfig.ini'))

    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_ph_var(self, plane=None):
        return f'cluster_charge[{choose(plane, self.N)}]'

    @save_pickle('Fit', sub_dir='PH', suf_args='all')
    def get_pulse_height(self, bin_size=None, cut=None, _redo=False):
        return self.draw_pulse_height(bin_size, cut, show=False)[1][0]

    def get_vcal(self, redo=False):
        h = self.draw_vcal_distribution(show=False, redo=redo)
        return ufloat(h.GetMean(), h.GetMeanError())

    def get_nhits(self, cut=None):
        return self.get_tree_vec(f'n_hits[{self.N}]', self.Cut(cut), dtype='u2')

    def get_track_vars(self, mm=True, local=False):
        return self.Cut.get_track_vars(self.DUT.Number - 1, mm, local)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region OCCUPANCY
    def draw_occupancy(self, roc=None, name=None, cluster=True, tel_coods=False, cut='', show=True):
        """ draw hitmap or cluster map """
        return self.Tel.draw_occupancy(choose(roc, self.N), choose(name, self.DUT.Name, roc), cluster, tel_coods, cut, show)

    def draw_occupancy_trend(self, cut=None, fid=False, bin_size=None, **kwargs):
        cut = self.Cut.generate_custom(exclude='fiducial' if not fid else []) if cut is None else self.Cut(cut)
        x, y, zz = self.get_tree_vec([self.get_t_var()] + self.Tel.get_hit_vars(self.N), cut)
        h = self.Draw.histo_3d(x, y, zz, self.Bins.get_time(bin_size) + self.Bins.get_pixel(), y_tit='col', z_tit='row')
        x, px, py = get_3d_profiles(h, 'zy')
        y0, y1 = [[ufloat(ip.GetMean(), ip.GetMeanError()) for ip in p] for p in [px, py]]
        g = [self.Draw.graph(x, y, x_tit='Time [hh:mm]', y_tit='Mean Pixel Pos', show=False) for y in [y0, y1]]
        return self.Draw.multigraph(g, 'Hit Position Trend', ['Column', 'Row'], t_ax_off=0, **kwargs)
    # endregion OCCUPANCY
    # ----------------------------------------

    # ----------------------------------------
    # region DISTRIBUTIONS
    def draw_adc_distribution(self, cut=None, col=None, row=None, pix=None, **kwargs):
        x = self.get_tree_vec('adc', self.Cut(cut) + self.Cut.generate_masks(col, row, pix, exclude=False).Value + self.Cut.get_plane(), 'i2')
        return self.Draw.distribution(x, Bins.get_adc(), 'ADC Distribution', **prep_kw(kwargs, logy=True))

    @save_pickle('VcalDisto', suf_args='all')
    def get_vcal_disto(self, cut=None, col=None, row=None, pix=None, vcal=True, _redo=False):
        cut = self.Cut(cut) + self.Cut.generate_masks(col, row, pix, exclude=False).Value + self.Cut.get_ncluster()
        n, v = self.get_nhits(cut), self.Calibration.get_vcals(*self.get_tree_vec(['col', 'row', 'adc'], cut + self.Cut.get_plane(), dtype='i2'))
        v = sum(insert(v, cumsum(n).astype('i').repeat(max(n) - n), 0).reshape(n.size, max(n)), axis=1)  # fill arrays with zeros where there are less than max hits
        v *= (1 if vcal else Bins.Vcal2El)
        return self.Draw.distribution(v, title='Pulse Height Distribution', x_tit=f'Pulse Height [{"VCAL" if vcal else "e"}]', show=False)

    def draw_vcal_distribution(self, cut=None, col=None, row=None, pix=None, vcal=True, redo=False, **kwargs):
        h = self.get_vcal_disto(cut, col, row, pix, vcal, _redo=redo)
        return self.Draw.distribution(h, **kwargs, filename=f'PHDisto{"V" if vcal else "E"}')

    @save_pickle('PH', suf_args='all')
    def get_signal_disto(self, roc=None, cut=None, vcal=False, _redo=False):
        x = self.get_tree_vec(self.get_ph_var(roc), self.Cut(cut)) / (Bins.Vcal2El if vcal else 1)
        return self.Draw.distribution(x, title='Pulse Height Distribution', x_tit=f'Pulse Height [{"VCAL" if vcal else "e"}]', show=False)

    def draw_threshold(self, x=1500, y0=0, y1=1, show=True):
        if show:
            return self.Draw.y_axis(x, y0, y1, f'threshold #approx {x}e', off=.3, line=True, opt='-L')

    def draw_signal_distribution(self, roc=None, cut=None, vcal=False, redo=False, draw_thresh=False, **kwargs):
        h = self.get_signal_disto(roc, cut, vcal, _redo=redo)
        t = self.draw_threshold(1500, 0, h.GetMaximum(), draw_thresh)
        return self.Draw.distribution(h, **prep_kw(kwargs, x_range=ax_range(10, 10, fl=.2, fh=.5, h=h), leg=t))

    def draw_ncluster_disto(self, n=1, cut=None, redo=False, **kwargs):
        return self.draw_signal_distribution(cut=self.Cut.make(f'{n}cl', self.Cut(cut) + self.Cut.get_ncluster(n)), redo=redo, **kwargs)

    def draw_nhit_disto(self, n=1, cut=None, redo=False, **kwargs):
        return self.draw_signal_distribution(cut=self.Cut.make(f'{n}hit', self.Cut(cut) + self.Cut.get_ncluster(1) + self.Cut.get_nhit(n)), redo=redo, **kwargs)

    def draw_nhit_distos(self, nmax=4, cut=None, redo=False, **kwargs):
        h = [self.draw_nhit_disto(n, cut, redo, show=False) for n in range(1, nmax + 1)]
        return self.Draw.stack(h, 'NHit Distributions', ['1 hit'] + [f'{n} hits' for n in range(2, nmax + 1)], **kwargs)
    # endregion DISTRIBUTIONS
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def draw_pulse_height(self, bin_size=30000, cut=None, **kwargs):
        """ Pulse height analysis vs event for a given cut. If no cut is provided it will take all. """
        (x, y), bins = self.get_tree_vec([self.get_t_var(), self.get_ph_var()], self.Cut(cut)), self.Bins.get_time(bin_size)
        h = self.Draw.profile(x, y, bins, **prep_kw(kwargs, x_tit='Time [hh:mm]', y_tit='Pulse Height [e]', y_off=1.8, lm=.17, graph=True, stats=set_statbox(fit=True), t_ax_off=0))
        fit = FitRes(h.Fit('pol0', 'qs'))
        self.Draw.save_plots(f'PulseHeight{bin_size}')
        return h, fit
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region 2D DISTRIBUTIONS
    def draw_adc_map(self, cut=None, **kwargs):
        x, y, zz = self.get_tree_vec(['col', 'row', 'adc'], self.Cut(cut) + self.Cut.get_plane())
        return self.Draw.prof2d(x, y, zz, Bins.get_pixel(), **prep_kw(kwargs, x_tit='col', y_tit='row', z_tit='ADC'))

    def draw_threshold_map(self, vcal=True, cols=None, rows=None, pix=None, **kwargs):
        x, y, zz = self.Calibration.get_thresholds(cols, rows, pix, vcal).T
        return self.Draw.prof2d(x, y, zz, Bins.get_pixel(), 'Artificial Threshold Map', **prep_kw(kwargs, x_tit='col', y_tit='row', z_tit=f'0Treshold [{"VCAL" if vcal else "ke"}]'))

    def draw_adc_fixed_vcal_map(self, vcal=200, **kwargs):
        cols, rows = self.Cut.get_fid_lines()
        x, y, zz = array([[col, row, self.Calibration.get_adc(col, row, vcal)] for col in cols for row in rows]).T
        return self.Draw.prof2d(x, y, zz, Bins.get_pixel(), f'ADC Map (VCAL={vcal}', **prep_kw(kwargs, x_tit='col', y_tit='row', z_tit='ADC'))

    def draw_sig_map_disto(self, res=None, cut=None, fid=True, x_range=None, redo=False, normalise=False, ret_value=False, ph_bins=None, show=True, save=True):
        super(PixAnalysis, self).draw_sig_map_disto(res, cut, fid, x_range, redo, normalise, ret_value, ph_bins=self.Bins.get_ph(), show=show, save=save)
    # endregion 2D DISTRIBUTIONS
    # ----------------------------------------

    # ----------------------------------------
    # region EFFICIENCY
    def get_efficiency_cut(self, trig_phase=True):
        return self.Cut.generate_custom(exclude=None if trig_phase else 'trigger_phase', prnt=False)

    def get_eff_var(self, plane=None):
        return f'n_hits[{choose(plane, self.N)}] > 0'

    def get_hit_efficiency(self, plane=None, cut=None):
        return calc_eff(values=self.get_tree_vec(self.get_eff_var(plane), choose(cut, self.get_efficiency_cut()), dtype=bool))

    def draw_eff_vs_chi2(self, **kwargs):
        x, e = self.get_tree_vec(['chi2_tracks', self.get_eff_var()], self.get_efficiency_cut())
        self.Draw.efficiency(x, e, find_bins(x, lfac=0, lq=0), title='Efficiency vs Chi2', **prep_kw(kwargs, x_tit='#chi^{2}'))

    def draw_hit_efficiency(self, cut=None, bin_size=None, rel_time=False, **kwargs):
        (x, e), bins = self.get_tree_vec([self.get_t_var(), self.get_eff_var()], choose(cut, self.get_efficiency_cut())), self.Bins.get_time(bin_size, cut)
        g = self.Draw.efficiency(x, e, bins, 'Hit Efficiency', **prep_kw(kwargs, **self.get_t_args(rel_time), y_range=[-5, 115], gridy=True, draw_opt='apz'))
        fit = FitRes(g.Fit('pol0', 'qs'))
        self.Draw.stats(fit, width=.35, y2=.35, names=['Efficiency'])
        self.Draw.preliminary()
        self.Draw.save_plots('HitEfficiency', **kwargs)
        return fit if fit.Parameter(0) is not None else 0

    def draw_efficiency_map(self, res=None, fid=False, **kwargs):
        x, y, zz = self.get_tree_vec(self.get_track_vars() + [self.get_eff_var()], self.get_efficiency_cut() if fid else self.Cut.generate_custom('fiducial'))
        tit, (xtit, ytit), ztit = 'Efficiency Map', [f'Track Position {i} [mm]' for i in ['X', 'Y']], 'Efficiency [%]'
        self.Draw.prof2d(x, y, zz * 100, Bins.get_global(res), tit, **prep_kw(kwargs, x_tit=xtit, y_tit=ytit, z_tit=ztit))
        self.Draw.preliminary()
        self.draw_fid_cut()
        self.Draw.save_plots('Efficiency Map')

    def get_fiducial_cell(self, n):
        x1, x2, y1, y2 = self.Cut.CutConfig['fiducial']
        nx = int(round((x2 - x1) / Plane.PX))
        return round(x1 + Plane.PX * (n % nx), 4), round(y1 + Plane.PY * (n / nx), 4)

    def draw_cell_efficiency(self, nbins=None, **dkw):
        x, y, e = self.get_tree_vec(self.get_track_vars() + [self.get_eff_var()], self.get_efficiency_cut())
        bins = None if nbins is None else [nbins, 0, Plane.PX, nbins, 0, Plane.PY]
        self.Draw.prof2d(x % Plane.PX, y % Plane.PY, e * 100, bins, 'Cell Efficiency', **prep_kw(dkw, x_tit='Track X [mm]', y_tit='Track Y [mm]', z_tit='Efficiency [%]'))

    def draw_efficiency_vs_trigphase(self, **kwargs):
        x, e = self.get_tree_vec(['trigger_phase[1]', self.get_eff_var()], self.get_efficiency_cut(trig_phase=False))
        return self.Draw.efficiency(x, e, make_bins(-.5, 10), 'Trigger Phase Efficiency', **prep_kw(kwargs, x_tit='Trigger Phase', x_range=[-1, 10], draw_opt='bap'))
    # endregion EFFICIENCY
    # ----------------------------------------

    # ----------------------------------------
    # region TRIGGER PHASE
    def draw_trigger_phase(self, cut=None, **kwargs):
        self.Tel.draw_trigger_phase(dut=True, cut=cut, **kwargs)

    def draw_trigger_phase_offset(self, cut=None, **dkw):
        x, y = self.get_tree_vec([f'trigger_phase[{i}]' for i in [1, 0]], choose(cut, self.Cut.generate_custom('trigger_phase', prnt=False)))
        return self.Draw.distribution(x - y, make_bins(-9.5, 10), **prep_kw(dkw, ndivx=20, x_tit='#Delta Trigger Phase', stats=set_entries()))

    def draw_tphase_offset_trend(self, bin_width=None, cut=None, **dkw):
        t, y0, y1 = self.get_tree_vec([self.get_t_var()] + [f'trigger_phase[{i}]' for i in [1, 0]], choose(cut, self.Cut.generate_custom('trigger_phase', prnt=False)))
        return self.Draw.profile(t, y0 - y1, self.Bins.get_time(bin_width, cut), 'Trigger Phase vs Time', **prep_kw(dkw, graph=True, y_tit='Trigger Phase', y_range=[-9, 9], **self.get_t_args()))
    # endregion TRIGGER PHASE
    # ----------------------------------------

    # ----------------------------------------
    # region ALIGNMENT
    def get_alignment(self):
        from pixel.alignment import PixAlignment
        return PixAlignment

    def draw_correlation(self, tel_plane=None, offset=0, bin_size=200):
        self.get_alignment()(self.Run.Converter, tel_plane).draw_correlation(offset, bin_size)
    # endregion ALIGNMENT
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_dut_hits(self, dut2=None, cut=None, **dkw):
        duts = [self.get_next_dut(dut2), self.DUT]
        x, y = self.get_tree_vec([f'n_hits[{dut.Number + self.Run.NTelPlanes - 1}]' for dut in duts], self.Cut(cut))
        x_tit, y_tit = [f'Number of Hits in {dut.Name}' for dut in duts]
        return self.Draw.histo_2d(x, y, [w for i in [x, y] for w in make_bins(*find_range(i, 0, 1.5))], **prep_kw(dkw, x_tit=x_tit, y_tit=y_tit, logz=True))

    def draw_hit_pie(self, dut2=None):
        duts = [self.get_next_dut(dut2), self.DUT]
        x, y = self.get_tree_vec([f'n_hits[{dut.Number + self.Run.NTelPlanes - 1}]' for dut in duts])
        labels = ['No Hits'] + [f'{dut.Name} Hit' for dut in duts] + ['Both Hits']
        e = [count_nonzero(i) for i in [(x == 0) & (y == 0), (x > 0) & (y == 0), (x == 0) & (y > 0), (x > 0) & (y > 0)]]
        self.Draw.pie(labels, e, offset=.05, h=.04, r=.2, text_size=.025, angle3d=70, angle_off=250, label_format='%txt (%perc)')

    def draw_cluster_size(self, cut='', show=True):
        return self.Tel.draw_cluster_size(self.N, self.DUT.Name, cut, show)

    def draw_residual(self, mode=None, cut=None, **dkw):
        return self.Tracks.draw_residual(self.N, mode=mode, cut=cut, **dkw)
    # endregion DRAW
    # ----------------------------------------


if __name__ == '__main__':
    pargs = init_argparser(run=139, tc='201810', dut=1, has_verbose=True, tree=True)
    z = PixAnalysis(pargs.run, pargs.dut, pargs.testcampaign, pargs.tree)
