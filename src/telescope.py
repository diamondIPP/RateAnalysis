from __future__ import print_function

from analysis import *
from ROOT import TH2F, TH1F, TCut, TF1, TGraph, TH1I, TProfile, TMultiGraph, TH2I, THStack
from numpy import log, split, polyfit, polyval, genfromtxt, arctan, sin, cos, tan, rad2deg

from binning import *
from selector import run_selector
from helpers.fit import Langau


class TelecopeAnalysis(Analysis):
    """ Class for the analysis of the telescope specific stuff of a single run. """

    # TODO make subclass and not inherited class...
    def __init__(self, run=None, test_campaign=None, tree=True, t_vec=None, verbose=False):
        """
        :param run: inits a new Run instance if number is provided or takes the provided Run instance
        :param test_campaign:   if None is provided: uses default from main config
        :param verbose:
        in order to save tel plots use argument both_dias=True in save_plots or save_tel_histo
        """
        Analysis.__init__(self, test_campaign, verbose)

        # Run
        self.Run = run_selector(run, self.TCString, tree, t_vec, verbose)

        # basics
        self.TelSaveDir = str(self.Run.Number).zfill(3)
        self.Tree = self.Run.Tree
        self.StartTime = self.Run.StartTime if self.Tree else time_stamp(self.Run.LogStart)

        if self.Tree:
            self.Cut = Cut(self)
            self.NRocs = self.Run.NPlanes
            self.StartEvent = self.Cut.get_min_event()
            self.EndEvent = self.Cut.get_max_event()
            self.Bins = Bins(self.Run, cut=self.Cut)

    def get_t_var(self):
        return 'time / 1000.' if self.Run.TimeOffset is None else '(time - {}) / 1000.'.format(self.Run.TimeOffset)

    def get_flux(self, corr=False, rel_error=0, show=False):
        return self._get_flux(prnt=False, show=show) if self.Tree and self.has_branch('rate') else self.Run.get_flux(rel_error)

    def get_time(self):
        return self.Run.get_time()

    def get_tree_vecs(self, strings, cut=None, dtypes=None):
        n = self.Tree.Draw(':'.join(strings), self.Cut(cut), 'goff')
        dtypes = [None] * len(strings) if dtypes is None else dtypes
        return [self.Run.get_root_vec(n, i, dtype) for i, dtype in enumerate(dtypes)]

    # ----------------------------------------
    # region TRACKS
    def _draw_cluster_size(self, roc, name=None, cut='', show=True):
        h = TH1I('h_cs', 'Cluster Size {d}'.format(d='ROC {n}'.format(n=roc) if name is None else name), 10, 0, 10)
        self.Tree.Draw('cluster_size[{d}]>>h_cs'.format(d=roc), TCut(cut), 'goff')
        self.format_statbox(entries=True)
        format_histo(h, x_tit='Cluster Size', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.save_tel_histo(h, 'ClusterSize', show, logy=True)
        return h

    def draw_n_clusters(self, roc=0, name=None, cut='', show=True):
        h = TH1I('h_cs', 'Number of Clusters {d}'.format(d='ROC {n}'.format(n=roc) if name is None else name), 10, 0, 10)
        self.Tree.Draw('n_clusters[{d}]>>h_cs'.format(d=roc), TCut(cut), 'goff')
        self.format_statbox(entries=True)
        format_histo(h, x_tit='Number of Clusters', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.save_tel_histo(h, 'NCluster', show, logy=True)
        return h

    def draw_event(self, event, plane, show=True, grid=True):
        cut = 'plane == {r}'.format(r=plane)
        h = TH2F('h_ed{i}'.format(i=plane), 'Event Hits for Plane {r}'.format(r=plane), *self.Bins.get_pixel())
        self.Tree.Draw('row:col>>h_ed{i}'.format(i=plane), cut, 'goff', 1, event)
        format_histo(h, x_tit='col', y_tit='row', y_off=1.3, stats=0)
        self.object(h, draw_opt='col', show=show)
        self.draw_pixel_grid() if grid else do_nothing()
        self.save_tel_plots('EventDisplay{e}_{p}'.format(e=event, p=plane))

    def draw_pixel_grid(self):
        n, x, n, y = self.Bins.get_pixel()
        self.grid(x, y, color=921)

    def get_events(self, cut=None, redo=False):
        return self.Run.get_root_vec(dtype='i4', var='Entry$', cut=self.Cut(cut))

    def get_plane_hits(self):
        t = self.info('getting plane hits...', next_line=False)
        self.Tree.SetEstimate(self.Tree.GetEntries(self.Cut.CutStrings['tracks'].GetTitle()) * self.NRocs)
        n = self.Tree.Draw('cluster_xpos_local:cluster_ypos_local', self.Cut.CutStrings['tracks'], 'goff')
        x, y = [array(split(vec, arange(self.NRocs, n, self.NRocs))).T * 10 for vec in self.Run.get_root_vecs(n, 2)]
        z_positions = genfromtxt(join(self.Run.Converter.TrackingDir, 'ALIGNMENT', 'telescope{}.dat'.format(self.Run.Converter.TelescopeID)), skip_header=1, usecols=[1, 5])
        z_positions = z_positions[where(z_positions[:, 0] > -1)][:, 1] * 10
        self.add_to_info(t)
        return x, y, z_positions

    def calculate_residuals(self, x, y, z_positions, steps=0):
        dx, dy, da = None, None, None
        for _ in range(steps + 1):
            dx, dy = self.calc_residual_step(x, y, z_positions)
            da = self.align(x, y, dx, dy)
        return dx, dy, da

    @staticmethod
    def align(x, y, dx, dy):
        bins = [arange(-4, 4.01, .2), arange(-1, 1.02, .03), arange(-4, 4.15, .15), arange(-1, 1.02, .02)]
        histos = [[TH2F('h{}{}'.format(i, j), 'rotres', bins[i].size - 1, bins[i], bins[i + 1].size - 1, bins[i + 1]) for i in arange(0, 4, 2)] for j in range(x.shape[0])]
        angles = []
        for h, ix, iy, idx, idy in zip(histos, x, y, dx, dy):
            h[0].FillN(ix.size, ix.astype('d'), idy.astype('d'), full(ix.size, 1, dtype='d'))
            h[1].FillN(ix.size, iy.astype('d'), idx.astype('d'), full(ix.size, 1, dtype='d'))
            [format_histo(h[i], x_tit='X [mm]', y_tit='dY [mm]', y_off=1.3, pal=53) for i in range(2)]
            fits = [h[i].Fit('pol1', 'qs0') for i in range(2)]
            angles.append([fits[i].Parameter(1) for i in range(2)])
        angles = arctan(((array(angles)[:, 0] - array(angles)[:, 1]) / 2))  # calculate the mean of the two angles (one has the opposite sign of the other)
        rot = array([[[cos(a), -sin(a)], [sin(a), cos(a)]] for a in angles])
        for i, irot in enumerate(rot):
            x[i], y[i] = irot.dot(array([x[i], y[i]]))
        # self.format_statbox(entries=True, x=.82)
        # self.draw_histo(histos[0][0], draw_opt='colz', rm=.16)
        x += mean(dx, axis=1).reshape((x.shape[0], 1))
        y += mean(dy, axis=1).reshape((x.shape[0], 1))
        return angles

    @staticmethod
    def calc_residual_step(x, y, z_positions):
        x_fits, y_fits = [polyfit(z_positions, vec, 1) for vec in [x, y]]
        return [polyval(fits, z_positions.reshape((x.shape[0], 1))) - vec for fits, vec in [(x_fits, x), (y_fits, y)]]

    def draw_raw_residuals(self, roc=None, steps=0, show=True):
        x, y, z_positions = self.get_plane_hits()
        dx, dy, da = self.calculate_residuals(x, y, z_positions, steps)
        x_bins, y_bins = arange(-1, 1.03, .03),  arange(-1, 1.02, .02)
        histos = [TH2F('h{}'.format(i), 'Raw Residuals ROC {}'.format(i), x_bins.size - 1, x_bins, y_bins.size - 1, y_bins) for i in range(dx.shape[0])]
        self.format_statbox(entries=True, x=.82)
        for h, ix, iy in zip(histos, dx, dy):
            h.FillN(ix.size, ix.astype('d'), iy.astype('d'), full(ix.size, 1, dtype='d'))
            format_histo(h, x_tit='dX [mm]', y_tit='dY [mm]', y_off=1.3, pal=53)
        if roc is not None:
            self.object(histos[roc], draw_opt='colz', lm=.12, rm=.16, show=show, logz=True)
        else:
            c = self.canvas('crr', 'Raw Residuals', 1.5, 1.5, show=show, divide=(2, 2))
            for i, h in enumerate(histos, 1):
                self.object(h, draw_opt='colz', lm=.12, rm=.16, show=show, canvas=c.cd(i), logz=True)
            self.save_plots('RawResisudal', canvas=c, show=show)
        return dx, dy, da

    def draw_mean_residuals(self, roc=0, mode='x', steps=5):
        x, y, z_positions = self.get_plane_hits()
        values = array([[abs(mean(vec[roc])) for vec in self.calculate_residuals(x, y, z_positions)][:2] for _ in range(steps)])
        g = self.make_tgrapherrors('gmr', 'Mean Residual in {} for Each Step'.format(mode.title()), x=arange(steps), y=values[:, 0 if mode == 'x' else 1])
        format_histo(g, x_tit='Step #', y_tit='Mean Residual [mm]', y_off=1.5, x_range=[-.5, steps - .5], ndivx=505)
        self.object(g, draw_opt='ap', lm=.13)
        return array([[sum(values[:, 0])], [sum(values[:, 1])]])

    def draw_residual_angles(self, roc=0, steps=10):
        x, y, z_positions = self.get_plane_hits()
        angles = array([abs(rad2deg(tan(self.calculate_residuals(x, y, z_positions)[-1][roc]))) for _ in xrange(steps)])
        g = self.make_tgrapherrors('gmra', 'Residual Angles for Each Step', x=arange(steps), y=angles)
        format_histo(g, x_tit='Step #', y_tit='Residual Angle [#circ]', y_off=1.5, x_range=[-.5, steps - .5], ndivx=505)
        self.object(g, draw_opt='ap', lm=.13)

    # endregion TRACKS
    # ----------------------------------------

    # ----------------------------------------
    # region PIXEL HITS
    def _draw_occupancy(self, plane, name=None, cluster=True, tel_coods=False, cut='', show=True, prnt=True):
        name = 'ROC {i}'.format(i=plane) if name is None else name
        bins = self.Bins.get_native_global(mm=True) if tel_coods else self.Bins.get_pixel()
        set_root_warnings(False)
        h = TH2F('h_hm{i}'.format(i=plane), '{h} Occupancy {n}'.format(n=name, h='Hit' if not cluster else 'Cluster'), *bins)
        cut_string = self.Cut() if cut is None else TCut(cut)
        cut_string += 'plane == {0}'.format(plane) if not cluster else ''
        draw_string = 'cluster_row[{i}]:cluster_col[{i}]' if cluster else 'row:col'
        draw_string = 'cluster_ypos_local[{i}] * 10:cluster_xpos_local[{i}] * 10' if tel_coods else draw_string
        self.format_statbox(entries=True, x=.83)
        self.Tree.Draw('{ds}>>h_hm{i}'.format(ds=draw_string.format(i=plane), i=plane), cut_string, 'goff')
        format_histo(h, x_tit='x [mm]' if tel_coods else 'col', y_tit='y [mm]' if tel_coods else 'row', y_off=1.2)
        self.save_tel_histo(h, 'HitMap{0}'.format(plane), show, draw_opt='colz', rm=.15, prnt=prnt)
        return h

    def draw_occupancies(self, planes=None, cut='', cluster=True, show=True, prnt=True):
        planes = range(4) if planes is None else list(planes)
        histos = [self._draw_occupancy(plane, cluster=cluster, cut=cut, show=False, prnt=False) for plane in planes]
        set_root_output(show)
        c = self.canvas('c_hm', 'Hitmaps', x=1.5, y=1.5, divide=(2, 2), show=show)
        for i, h in enumerate(histos, 1):
            h.SetStats(0)
            pad = c.cd(i)
            pad.SetBottomMargin(.15)
            h.Draw('colz')
        self.save_tel_plots('HitMaps', sub_dir=self.TelSaveDir, show=show, prnt=prnt)

    def draw_tracking_map(self, at_dut=1, res=.7, cut='', show=True, prnt=True):
        set_root_output(False)
        h = TH2I('htm{}'.format(at_dut), 'Tracking Map', *self.Bins.get_global(res, mm=True))
        x_var, y_var = [self.Cut.get_track_var(at_dut - 1, v) for v in ['x', 'y']]
        self.Tree.Draw('{y} * 10:{x} * 10 >> htm{}'.format(at_dut, x=x_var, y=y_var), TCut(cut), 'goff')
        format_histo(h, x_tit='Track Position X [mm]', y_tit='Track Position Y [mm]', y_off=1.4, z_off=1.5, z_tit='Number of Hits', ncont=50, ndivy=510, ndivx=510, stats=0)
        self.save_tel_histo(h, 'TrackMap{}'.format(at_dut), lm=.12, rm=.16, draw_opt='colz', show=show, x_fac=1.15, prnt=prnt)
        return h

    def draw_beam_profile(self, at_dut=1, mode='x', fit=True, fit_range=.8, res=.7, show=True, prnt=True):
        h = self.draw_tracking_map(at_dut, res, show=False, prnt=prnt)
        p = h.ProjectionX() if mode.lower() == 'x' else h.ProjectionY()
        format_histo(p, title='Profile {}'.format(mode.title()), name='pbp{}'.format(self.Run.Number), y_off=1.3, y_tit='Number of Hits', fill_color=self.FillColor)
        self.format_statbox(all_stat=True)
        self.object(p, lm=.13, show=show)
        if fit:
            fit = self.fit_beam_profile(p, fit_range)
        self.save_plots('BeamProfile{}{}'.format(at_dut, mode.title()), both_dias=True, prnt=prnt)
        return FitRes(fit) if fit else p

    @staticmethod
    def fit_beam_profile(p, fit_range):
        x_peak = p.GetBinCenter(p.GetMaximumBin())
        x_min, x_max = [p.GetBinCenter(i) for i in [p.FindFirstBinAbove(0), p.FindLastBinAbove(0)]]
        return p.Fit('gaus', 'qs', '', x_peak - (x_peak - x_min) * fit_range, x_peak + (x_max - x_peak) * fit_range)
    # endregion PIXEL HITS
    # ----------------------------------------

    def _draw_trigger_phase(self, dut=False, cut=None, show=True):
        cut_string = self.Cut.generate_custom(exclude=['trigger_phase']) if cut is None else TCut(cut)
        h = TH1I('h_tp', 'Trigger Phase', 10, 0, 10)
        self.Tree.Draw('trigger_phase[{r}]>>h_tp'.format(r=1 if dut else 0), cut_string, 'goff')
        self.format_statbox(entries=True)
        format_histo(h, x_tit='Trigger Phase', y_tit='Number of Entries', y_off=1.95, fill_color=self.FillColor)
        self.save_tel_histo(h, '{m}TriggerPhase'.format(m='DUT' if dut else 'Tel'), show, lm=.145)

    def _draw_trigger_phase_time(self, dut=False, bin_width=None, cut=None, show=True):
        h = TProfile('htpt', 'Trigger Phase Offset vs Time - {}'.format('DUT' if dut else 'TEL'), *self.Bins.get(bin_width, vs_time=True))
        self.Tree.Draw('trigger_phase[{}]:{}>>htpt'.format(1 if dut else 0, self.get_t_var()), self.Cut.generate_custom(exclude='trigger_phase') if cut is None else cut, 'goff')
        self.format_statbox(entries=True, y=0.88)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Trigger Phase', y_off=1.8, fill_color=self.FillColor, t_ax_off=self.Run.StartTime)
        self.save_histo(h, 'TPTime', show, lm=.16)

    def draw_pix_map(self, n=1, start=None, plane=1):
        start_event = self.StartEvent if start is None else start
        h = TH2F('h', 'Pixel Map', 52, 0, 51, 80, 0, 79)
        self.Tree.GetEntry(start_event)
        for pln, col, row, adc in zip(self.Tree.plane, self.Tree.col, self.Tree.row, self.Tree.adc):
            if pln == plane:
                h.SetBinContent(col + 1, row + 1, -adc)
        c = TCanvas('c', 'Pixel Map', 1000, 1000)
        c.SetBottomMargin(.15)
        c.SetRightMargin(.14)
        h.SetStats(0)
        h.GetZaxis().SetTitle('adc [au]')
        h.GetZaxis().SetTitleOffset(1.3)
        format_histo(h, x_tit='col', y_tit='row')
        h.Draw('colz')
        self.Objects.append([c, h])
        self.save_plots('PixMapPlane{pln}{evts}'.format(pln=plane, evts=n), sub_dir=self.TelSaveDir)

    # ----------------------------------------
    # region TIME
    def draw_time(self, show=True, corr=False):
        n = self.Tree.Draw(self.get_t_var(), '', 'goff')
        t = self.Run.get_root_vec(n) if not corr else self.Run.Time / 1000.
        t -= t[0]
        # gr = self.make_tgrapherrors('g_t', 'Time vs Events', x=arange(t.size, dtype='d'), y=t)
        gr = TGraph(t.size, arange(t.size, dtype='d'), t)
        fit = gr.Fit('pol1', 'qs')
        self.info('Average data taking rate: {r:5.1f} Hz'.format(r=1 / fit.Parameter(1)))
        format_histo(gr, 'g_t', 'Time vs Events', x_tit='Event Number', y_tit='Time [s]', y_off=1.5)
        self.object(gr, show=show, draw_opt='al', lm=.13, rm=.08)

    def get_event_at_time(self, seconds, rel=False):
        return self.Run.get_event_at_time(seconds, rel)
    # endregion TIME
    # ----------------------------------------

    # ----------------------------------------
    # region RATE
    def draw_beam_current_prof(self, bin_width=30, cut='', rel_t=True, show=True, save=True):
        if not self.has_branch('beam_current'):
            return warning('Branch "beam_current" does not exist!')
        set_root_output(False)
        p = TProfile('pbc', 'Beam Current Profile', *self.Bins.get_raw_time(bin_width))
        self.Tree.Draw('beam_current / 1000.:{}>>pbc'.format(self.get_t_var()), TCut('beam_current < 2500') + TCut(cut), 'goff')
        self.format_statbox(all_stat=True, y=.4)
        format_histo(p, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', fill_color=self.FillColor, markersize=.4, t_ax_off=self.Run.StartTime if rel_t else 0, y_range=[0, p.GetMaximum() * 1.1])
        c = self.object(p, draw_opt='hist', lm=.08, w=1.5, h=.75, ind=None, show=show)
        self.save_histo(p, 'BeamCurrentProf', canvas=c, draw_opt='esame', save=save)
        return p

    def draw_beam_current(self, cut='', rel_t=True, show=True, save=True):
        if not self.has_branch('beam_current'):
            return warning('Branch "beam_current" does not exist!')
        n = self.Tree.Draw('beam_current/1000.:{}'.format(self.get_t_var()), TCut('beam_current < 2500') + TCut(cut), 'goff')
        current, t = self.Run.get_root_vecs(n, 2)
        g = self.make_tgrapherrors('gbc', 'Beam Current', x=concatenate([t, [t[-1]]]), y=concatenate([current, [0]]))
        format_histo(g, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', fill_color=self.FillColor, markersize=.4, t_ax_off=self.Run.StartTime if rel_t else 0,
                     x_range=[g.GetX()[0], g.GetX()[n]])
        self.save_tel_histo(g, 'BeamCurrent', draw_opt='afp', lm=.08, x_fac=1.5, y_fac=.75, ind=None, show=show, save=save)
        return g

    def draw_rate(self, plane=1, flux=False, rel_t=True, show=True):
        """ Draws the single plane rates versus time. The first entry of the vector corresponds to the scintillator rate """
        if not self.has_branch('rate'):
            warning('The "rate" branch does not exist in this tree')
            return
        area = self.Run.get_unmasked_area()[plane] if plane in self.Run.get_unmasked_area() else .01 * .015 * 4160
        n = self.Tree.Draw('rate[{p}] {a}:{t}'.format(p=plane, a='/{}'.format(area) if flux else '', t=self.get_t_var()), 'beam_current < 10000 && rate[{}]<1e9'.format(plane), 'goff')
        rate = [self.Tree.GetV1()[i] for i in range(n)]
        t = [self.Tree.GetV2()[i] for i in range(n)]
        g = self.make_tgrapherrors('gpr', 'Rate of Plane {n}'.format(n=plane), x=t + [t[-1]], y=rate + [0])
        format_histo(g, x_tit='Time [hh:mm]', y_tit='Rate [Hz]', fill_color=self.FillColor, markersize=.4, t_ax_off=self.Run.StartTime if rel_t else 0)
        self.save_tel_histo(g, 'Plane{n}Rate'.format(n=plane), draw_opt='afp', lm=.08, x_fac=1.5, y_fac=.75, ind=None, show=show)

    def draw_flux(self, bin_width=5, cut='', rel_t=True, show=True, prnt=True):
        set_root_warnings(OFF)
        p = TProfile('pf', 'Flux Profile', *self.Bins.get_raw_time(bin_width=bin_width))
        p1, p2 = self.Run.TriggerPlanes
        a1, a2 = self.Run.get_unmasked_area().values()
        cut = TCut('beam_current < 10000 && rate[{0}] < 1e9 && rate[{1}] < 1e9 && rate[{0}] && rate[{1}]'.format(p1 + 1, p2 + 1)) + TCut(cut)
        # rate[0] is scintillator
        self.Tree.Draw('(rate[{p1}] / {a1} + rate[{p2}] / {a2}) / 2000 : {t}>>pf'.format(p1=p1 + 1, p2=p2 + 1, a1=a1, a2=a2, t=self.get_t_var()), cut, 'goff', self.Run.NEvents, 1)
        y_range = [0, p.GetMaximum() * 1.2]
        format_histo(p, x_tit='Time [hh:mm]', y_tit='Flux [kHz/cm^{2}]', fill_color=Draw.FillColor, markersize=1, t_ax_off=self.Run.StartTime if rel_t else 0, stats=0, y_range=y_range)
        self.Draw(p, 'FluxProfile', draw_opt='hist', lm=.08, w=1.5, h=.75, show=show, prnt=prnt)
        return p

    def draw_bc_vs_rate(self, cut='', show=True):
        g1 = self.draw_flux(cut=cut, show=False)
        g2 = self.draw_beam_current(cut=cut, show=False)
        fluxes = [g1.GetY()[i] for i in range(g1.GetN())]
        beam_currents = [g2.GetY()[i] for i in range(g2.GetN())]
        xbins = [int(max(beam_currents) + 10 - sorted(beam_currents)[3]), sorted(beam_currents)[3], max(beam_currents) + 10]
        ybins = [int(sqrt(g1.GetN()) * 4), sorted(fluxes)[3], max(fluxes)]
        print(xbins + ybins)
        h = TH2F('hbcr', 'Correlation between Beam Current and Flux', *(xbins + ybins))
        for flux, beam_cur in zip(fluxes, beam_currents):
            h.Fill(beam_cur, flux)
        format_histo(h, x_tit='Beam Current [mA]', y_tit='Flux [kHz/cm^{2}]', y_off=1.3, stats=0)
        self.save_tel_histo(h, 'BeamCurrentFlux', lm=.13, rm=.18, ind=None, show=show, draw_opt='colz')

    def calculate_flux(self, show=False, prnt=True):

        pickle_path = self.make_simple_pickle_path(sub_dir='Flux', dut='')

        def f():
            format_statbox(fit=True, entries=6)
            h = self.draw_flux(cut=self.Cut.generate_custom(include=['beam_interruptions', 'event_range'], prnt=prnt), show=False, prnt=prnt)
            values = [h.GetBinContent(i) for i in range(h.GetNbinsX()) if h.GetBinContent(i) and h.GetBinContent(i) < 1e6]
            m, s = mean_sigma(values)
            h = TH1F('hfl', 'Flux Distribution', int(sqrt(h.GetNbinsX()) * 2), m - 3 * s, m + 4 * s)
            for val in values:
                h.Fill(val)
            max_val = h.GetBinCenter(h.GetMaximumBin())
            fit = h.Fit('gaus', 'qs{}'.format('' if show else 0), '', max_val * .9, max_val * 1.1)
            format_histo(h, 'Fit Result', y_tit='Number of Entries', x_tit='Flux [kHz/cm^{2}]', fill_color=Draw.FillColor, y_off=1.3)
            self.Draw(h, lm=.13, show=show, prnt=prnt)
            m, s = fit.Parameter(1), fit.Parameter(2)
            m, s = (m, s) if s < m / 2. and fit.Ndf() and fit.Chi2() / fit.Ndf() < 10 else mean_sigma(values)
            return make_ufloat((m, s + .05 * m))

        return do_pickle(pickle_path, f, redo=show)
    # endregion RATE
    # ----------------------------------------

    def save_tree(self, cut=None):
        f = TFile('test.root', 'RECREATE')
        t = self.Tree.CloneTree(0)
        n = self.Tree.Draw('Entry$', self.Cut(cut), 'goff')
        good_events = self.Run.get_root_vec(n, dtype='i4')
        self.PBar.start(n)
        for i, ev in enumerate(good_events):
            self.Tree.GetEntry(ev)
            t.Fill()
            self.PBar.update(i)
        f.cd()
        t.Write()
        macro = self.Run.RootFile.Get('region_information')
        if macro:
            macro.Write()
        f.Write()
        f.Close()
        self.info('successfully saved tree with only cut events.')

    def fit_langau(self, h=None, nconv=30, show=True, chi_thresh=8, fit_range=None):
        h = self.draw_signal_distribution(show=show) if h is None and hasattr(self, 'draw_signal_distribution') else h
        fit = Langau(h, nconv, fit_range)
        fit.get_parameters()
        fit(show=show)
        get_last_canvas().Modified()
        get_last_canvas().Update()
        if fit.get_chi2() > chi_thresh and nconv < 80:
            self.Count += 5
            self.info('Chi2 too large ({c:2.2f}) -> increasing number of convolutions by 5'.format(c=fit.get_chi2()))
            fit = self.fit_langau(h, nconv + self.Count, chi_thresh=chi_thresh, show=show)
        print('MPV: {:1.1f}'.format(fit.get_mpv()))
        self.Count = 0
        self.add(fit)
        return fit


if __name__ == '__main__':

    pargs = init_argparser(run=23, tc='201908', tree=True, has_verbose=True)
    z = TelecopeAnalysis(pargs.run, pargs.testcampaign, pargs.tree, pargs.verbose)
