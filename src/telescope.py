#!/usr/bin/env python
# --------------------------------------------------------
#       analysis class for the telescope
# revised on Oct 4th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TH1F, TCut
from numpy import arange, concatenate, mean, sqrt, full, ceil
from src.binning import Bins
from helpers.draw import Draw, format_histo, set_statbox, set_entries, format_statbox, get_hist_vec, fill_hist, warning, time_stamp, do_nothing, update_canvas, mean_sigma, ufloat, do_pickle, calc_eff
from src.sub_analysis import SubAnalysis, choose


class Telescope(SubAnalysis):
    """ Class for the analysis of the telescope specific stuff of a single run. """

    def __init__(self, analysis):
        super().__init__(analysis, sub_dir=analysis.Run.Number, pickle_dir='Telescope', dut=False)

        self.StartTime = self.Run.StartTime if self.Tree else time_stamp(self.Run.LogStart)

    # ----------------------------------------
    # region GET
    def get_flux(self, rel_error=0, show=False):
        return self.calculate_flux(prnt=False, show=show) if self.Tree.Hash and self.has_branch('rate') else self.Run.get_flux(rel_error)

    def get_rate_var(self, plane, flux=False):
        if not self.has_branch('rate'):
            return warning('The "rate" branch does not exist in this tree')
        area = self.Run.get_unmasked_area()[plane] if plane in self.Run.get_unmasked_area() else self.Run.Plane.Area / 100
        # rate[0] is scintillator --> add + 1 to trigger planes
        return 'rate[{p}]{a}'.format(p=plane + 1, a='/{}'.format(area) if flux else '')

    def get_flux_var(self, plane):
        return self.get_rate_var(plane, flux=True)

    @staticmethod
    def get_col_var(plane, cluster=True, tel_coods=False):
        return ('cluster_xpos_local[{}] * 10' if tel_coods else 'cluster_col[{}]' if cluster else 'col').format(plane)

    @staticmethod
    def get_row_var(plane, cluster=True, tel_coods=False):
        return ('cluster_ypos_local[{}] * 10' if tel_coods else 'cluster_row[{}]' if cluster else 'row').format(plane)

    def get_hit_vars(self, plane, cluster=True, tel_coods=False):
        return [self.get_col_var(plane, cluster, tel_coods), self.get_row_var(plane, cluster, tel_coods)]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region HITS
    def draw_cluster_size(self, roc, name=None, cut='', show=True):
        values = self.get_tree_vec(var='cluster_size[{}]'.format(roc), cut=self.Cut(cut))
        h = self.Draw.distribution(values, Bins.make(0, 50), 'Cluster Size {d}'.format(d='ROC {n}'.format(n=roc) if name is None else name), logy=True, stats=set_entries())
        format_histo(h, x_tit='Cluster Size', y_off=1.3, fill_color=Draw.FillColor, x_range=[0, h.FindLastBinAbove(5)])
        self.Draw.save_plots('ClusterSize{}'.format(roc), show=show)
        return h

    def draw_n_clusters(self, roc=0, name=None, cut='', y_range=None, x_range=None, show=True):
        values = self.get_tree_vec(var='n_clusters[{}]'.format(roc), cut=self.Cut(cut))
        h = self.Draw.distribution(values, Bins.make(0, 50), 'Number of Clusters {d}'.format(d='ROC {n}'.format(n=roc) if name is None else name), logy=True, stats=set_entries())
        format_histo(h, x_tit='Number of Clusters', y_off=1.3, fill_color=Draw.FillColor, x_range=choose(x_range, [0, h.FindLastBinAbove(2) + 1]), y_range=y_range)
        self.Draw.save_plots('NClusters{}'.format(roc), show=show)
        return h

    def draw_event(self, event, show=True, grid=True):
        x, y, p = self.get_tree_vec(self.get_hit_vars(0, cluster=False) + ['plane'], nentries=1, firstentry=event)
        c = self.Draw.canvas(w=1.5, h=1.5, divide=(int(ceil(sqrt(self.NRocs))), int(ceil(sqrt(self.NRocs)))), show=show)
        for i in range(self.NRocs):
            self.Draw.histo_2d(x[p == i], y[p == i], Bins.get_pixel(), 'Hits', draw_opt='col', x_tit='col', y_tit='row', rm=.03, stats=0, show=show, canvas=c.cd(i + 1))
            self.draw_pixel_grid() if grid else do_nothing()

    @staticmethod
    def draw_pixel_grid():
        n, x, n, y = Bins.get_pixel()
        Draw.grid(x, y, color=921)

    def draw_occupancy(self, plane, name=None, cluster=True, tel_coods=False, cut='', show=True, prnt=True):
        name = 'ROC {i}'.format(i=plane) if name is None else name
        cut_string = self.Cut(cut) + TCut('' if cluster else 'plane == {}'.format(plane))
        x, y = self.get_tree_vec(var=self.get_hit_vars(plane, cluster, tel_coods), cut=cut_string)
        bins = Bins.get_native_global() if tel_coods else Bins.get_pixel()
        h = self.Draw.histo_2d(x, y, bins, '{h} Occupancy {n}'.format(n=name, h='Cluster' if cluster else 'Hit'), show=show, draw_opt='colz', z_off=1.4)
        format_histo(h, x_tit='x [mm]' if tel_coods else 'col', y_tit='y [mm]' if tel_coods else 'row', y_off=1.2)
        self.Draw.save_plots('{}Map{}'.format('Cluster' if cluster else 'Hit', plane), prnt=prnt)
        return h

    def draw_occupancies(self, planes=None, cut='', cluster=True, show=True, prnt=True):
        histos = [self.draw_occupancy(plane, cluster=cluster, cut=cut, show=False, prnt=False) for plane in (range(self.NRocs) if planes is None else planes)]
        c = self.Draw.canvas('Hitmaps', w=1.5, h=1.5, divide=(2, 2), show=show)
        for i, h in enumerate(histos, 1):
            h.SetStats(0)
            pad = c.cd(i)
            pad.SetBottomMargin(.15)
            h.Draw('colz')
        self.Draw.save_plots('HitMaps', show=show, prnt=prnt)

    def draw_hit_efficiency(self, plane=0, cut=None, y_range=None, show=True):
        e, x = self.get_tree_vec(var=['n_clusters[{}]'.format(plane)] + [self.get_t_var()], cut=choose(cut, self.Cut.get('pulser')))
        self.Draw.efficiency(x, e > 0, self.Bins.get_time(), x_tit='Time[hh:mm]', stats=set_entries(), y_range=y_range, show=show, t_ax_off=0)
        self.Draw.info('Efficiency: {:1.2f}%'.format(calc_eff(values=e > 0)[0]))
    # endregion HITS
    # ----------------------------------------

    # ----------------------------------------
    # region TIME
    def draw_trigger_phase(self, dut=False, cut=None, show=True):
        values = self.get_tree_vec(var='trigger_phase[{}]'.format(1 if dut else 0), cut=self.Cut.generate_custom(exclude=['trigger_phase']) if cut is None else TCut(cut))
        self.Draw.distribution(values, Bins.make(0, 10), 'Trigger Phase', x_tit='Trigger Phase', y_off=1.95, fill_color=Draw.FillColor, show=show, lm=.145, stats=set_entries())

    def draw_trigger_phase_trend(self, dut=False, bin_width=None, cut=None, show=True):
        values, t = self.get_tree_vec(var=['trigger_phase[{}]'.format(1 if dut else 0), self.get_t_var()], cut=self.Cut.generate_custom(exclude=['trigger_phase']) if cut is None else TCut(cut))
        p = self.Draw.profile(t, values, self.Bins.get_time(bin_width, cut), '{} Trigger Phase vs Time'.format('DUT' if dut else 'TEL'), show=show, lm=.16, stats=set_entries())
        format_histo(p, x_tit='Time [hh:mm]', y_tit='Trigger Phase', y_off=1.8, fill_color=Draw.FillColor, t_ax_off=self.StartTime)

    def draw_time(self, show=True, corr=False):
        t = self.Run.Time / 1000 if corr else self.get_tree_vec(var=self.get_t_var())
        t = t[t >= 0]
        g = self.Draw.graph(t - t[0], arange(t.size), title='Events vs Time', draw_opt='al', lm=.13, tm=.11, show=show)
        fit = g.Fit('pol1', 'qs')
        self.info('Average data taking rate: {r:5.1f} Hz'.format(r=fit.Parameter(1)))
        format_histo(g, title='Time vs Events', x_tit='Event Number', y_tit='Time [s]', y_off=1.5)
    # endregion TIME
    # ----------------------------------------

    # ----------------------------------------
    # region RATE
    def draw_beam_current(self, bin_width=30, cut='', rel_t=True, prof=True, show=True, save=True):
        if not self.has_branch('beam_current'):
            return warning('Branch "beam_current" does not exist!')
        values, t = self.get_tree_vec(var=['beam_current', self.get_t_var()], cut=TCut('beam_current < 2500') + TCut(cut))
        if prof:
            h = self.Draw.profile(t, values, self.Bins.get_raw_time(bin_width), 'Beam Current [mA]', w=1.5, h=.75, lm=.08, draw_opt='hist', fill_color=Draw.FillColor)
        else:
            h = self.Draw.graph(concatenate([t, [t[-1]]]), concatenate([values, [0]]), w=1.5, h=.75, title='Beam Current [mA]', lm=.08, draw_opt='afp', fill_color=Draw.FillColor)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', markersize=.4, t_ax_off=self.StartTime if rel_t else 0, x_range=None if prof else [h.GetX()[0], h.GetX()[t.size]])
        format_statbox(h, all_stat=True)
        self.Draw.save_plots('BeamCurrent{}'.format(h.ClassName()[1]), show=show, save=save)
        return h

    def draw_plane_rate(self, plane=1, flux=False, rel_t=True, show=True):
        """ Draws the single plane rates versus time. The first entry of the vector corresponds to the scintillator rate """
        rate, t = self.get_tree_vec(var=[self.get_rate_var(plane, flux), self.get_t_var()], cut='beam_current < 10000 && rate[{}]<1e9'.format(plane + 1))
        g = self.Draw.graph(concatenate([t, [t[-1]]]), concatenate([rate, [0]]), title='Rate of Plane {n}'.format(n=plane), draw_opt='afp', lm=.08, w=1.5, h=.75, show=show)
        format_histo(g, x_tit='Time [hh:mm]', y_tit='{} [Hz]'.format('Flux' if flux else 'Rate'), fill_color=Draw.FillColor, markersize=.4, t_ax_off=self.StartTime if rel_t else 0)
        update_canvas()

    def draw_flux(self, bin_width=5, cut='', rel_time=True, show=True, prnt=True):
        cut = TCut('beam_current < 10000 && rate[{0}] < 1e9 && rate[{1}] < 1e9 && rate[{0}] && rate[{1}]'.format(*self.Run.TriggerPlanes + 1)) + TCut(cut)
        if self.has_branch('rate'):
            flux1, flux2, t = self.get_tree_vec(var=[self.get_flux_var(p) for p in self.Run.TriggerPlanes] + [self.get_t_var()], cut=cut)
            flux = mean([flux1, flux2], axis=0)[1:] / 1000
        else:
            t, flux = self.Run.Time / 1000, full(self.Run.NEvents - 1, self.get_flux().n)
        p = self.Draw.profile(t[1:], flux, self.Bins.get_raw_time(bin_width=bin_width), 'Flux Profile', fill_color=Draw.FillColor, draw_opt='hist', lm=.08, w=1.5, h=.75, show=show)
        format_histo(p, x_tit='Time [hh:mm]', y_tit='Flux [kHz/cm^{2}]', markersize=1, t_ax_off=self.StartTime if rel_time else 0, stats=0, y_range=[0, p.GetMaximum() * 1.2])
        self.Draw.save_plots('FluxProfile', prnt=prnt, show=show)
        return p

    def draw_bc_vs_rate(self, cut='', show=True):
        cut = TCut('beam_current < 10000 && rate[{0}] < 1e9 && rate[{1}] < 1e9 && rate[{0}] && rate[{1}]'.format(*self.Run.TriggerPlanes + 1)) + TCut(cut)
        flux1, flux2, bc = self.get_tree_vec(var=[self.get_flux_var(p) for p in self.Run.TriggerPlanes] + ['beam_current'], cut=cut)
        h = self.Draw.histo_2d(bc, mean([flux1, flux2], axis=0)[1:] / 1000, title='Correlation between Beam Current and Flux', show=show)
        format_histo(h, x_tit='Beam Current [mA]', y_tit='Flux [kHz/cm^{2}]', y_off=1.3, stats=0)

    def calculate_flux(self, show=False, prnt=True):
        def f():
            h = self.draw_flux(cut=self.Cut.generate_custom(include=['beam stops', 'event range'], prnt=prnt), show=False, prnt=prnt)
            values = get_hist_vec(h, err=False)
            m, s = mean_sigma(values[values > 0], err=False)
            h = TH1F('hfl', 'Flux Distribution', int(sqrt(h.GetNbinsX()) * 2), m - 3 * s, m + 4 * s)
            fill_hist(h, values)
            max_val = h.GetBinCenter(h.GetMaximumBin())
            fit = h.Fit('gaus', 'qs{}'.format('' if show else 0), '', max_val * .9, max_val * 1.1)
            format_histo(h, 'Fit Result', y_tit='Number of Entries', x_tit='Flux [kHz/cm^{2}]', fill_color=Draw.FillColor, y_off=1.3)
            self.Draw(h, lm=.13, show=show, prnt=prnt, stats=set_statbox(fit=True, all_stat=True))
            fm, fs = fit.Parameter(1), fit.Parameter(2)
            m, s = (fm, fs) if fs < fm / 2. and fit.Ndf() and fit.Chi2() / fit.Ndf() < 10 else (m, s)
            return ufloat(m, s + .05 * m)

        return do_pickle(self.make_simple_pickle_path('Flux'), f, redo=show)
    # endregion RATE
    # ----------------------------------------
