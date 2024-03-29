#!/usr/bin/env python
# --------------------------------------------------------
#       analysis class for the telescope tracks
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TCut, TF1, TMultiGraph, TMath
from numpy import log, genfromtxt, rad2deg, polyfit, polyval, tan, delete, deg2rad, cumsum
from scipy.stats import norm

from helpers.utils import arctan, save_pickle, update_pbar, quiet, PBAR
from plotting.draw import *
from plotting.fit import Gauss
from src.dut import Plane
from src.sub_analysis import SubAnalysis
from functools import partial


class Tracks(SubAnalysis):

    M = ['x', 'y']

    def __init__(self, analysis):
        super().__init__(analysis, sub_dir=analysis.Run.Number, pickle_dir='Tracks', dut=False)

    # ----------------------------------------
    # region GET
    def get_efficiency(self):
        return calc_eff(values=self.get_tree_vec(var='n_tracks'))

    def get_angles(self, mode='x'):
        a = self.get_tree_vec(var='angle_{}'.format(mode))
        return a[a != -999]

    @save_pickle('TrackAngle', suf_args='all')
    def get_mean_angle(self, mode='x', _redo=False):
        return mean_sigma(self.get_angles(mode))

    @save_pickle('Res', suf_args='all')
    def get_residual(self, roc=0, mode='x', cut='', _redo=False):
        return self.draw_residual(roc, mode, cut, ret_res=True, show=False)

    @save_pickle('URes', suf_args='all')
    def get_unbiased_residual(self, roc=0, mode='x', cut='', _redo=False):
        return self.draw_unbiased_residual(roc, mode, cut, fit=True, show=False)

    def get_residuals(self, mode='x', cut='', unbias=False, redo=False):
        return array([(self.get_unbiased_residual(roc, mode, cut, _redo=redo) if unbias else self.get_residual(roc, mode, cut, _redo=redo)) for roc in range(self.Run.NTelPlanes)])

    @update_pbar
    def get_resolution(self, mode='x', cut='', unbias=False, redo=False):
        return self.draw_resolution(mode, cut=cut, unbias=unbias, show=False, redo=redo)

    @update_pbar
    def get_mean_resolution(self, cut='', unbias=True, redo=False):
        return mean([self.get_resolution(m, cut, unbias, redo) for m in ['x', 'y']])

    @save_pickle('Res', suf_args='all')
    def get_chi2_residual(self, roc, chi2, mode='x', _redo=False):
        self.Cut.set_chi2(chi2)
        values = self.get_tree_vec(var='residuals_{m}[{r}]*1e4'.format(m=mode, r=roc), cut=self.Cut())
        return mean_sigma(values)[1]

    def get_chi2_residuals(self, roc, chi2s=None, mode='x'):
        chi2s = choose(chi2s, arange(10, 101, 10))
        return [self.get_chi2_residual(roc, chi2, mode) for chi2 in chi2s]

    def get_z_positions(self, e=0):
        x = genfromtxt(join(self.Run.Converter.TrackingDir, 'data', 'alignments.txt'), usecols=[0, 2, 6])
        return array([ufloat(ix, e) if e else ix for ix in x[(x[:, 0] == self.Run.Converter.TelescopeID) & (x[:, 1] > -1)][:, 2] * 10])  # [mm]

    @staticmethod
    def get_vars(local=False):
        return [f'cluster_{p}pos_{"local" if local else "tel"}' for p in ['x', 'y']]

    @staticmethod
    def get_res_var(mode=None):
        return f'residuals{"_{}".format(mode.lower()) if mode else ""}'

    @staticmethod
    def ax_tits(pixel=False):
        return Plane.AxTits if pixel else {f'{i.lower()}_tit': f'Track Position {i} [mm]' for i in ['X', 'Y']}

    def get_plane_hits(self, local=True, add_cut=''):
        t = self.info('getting plane hits...', endl=False)
        self.Tree.SetEstimate(self.Cut.get_size('tracks', excluded=False) * self.NRocs)
        ntel, loc = self.Run.NTelPlanes, array(self.get_tree_vec(self.get_vars(local), self.Cut['tracks'] + add_cut)).T * 10
        n = self.get_tree_vec('total_clusters', self.Cut['tracks'] + add_cut, dtype='i')
        self.add_to_info(t)
        return loc.T.reshape((2, -1, ntel)) if ntel == self.NRocs else loc[append(0, cumsum(n)[:-1]).repeat(ntel) + tile(arange(ntel, dtype='i'), n.size)].T.reshape((2, n.size, ntel))
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw(self, n=0, angle=None):
        add = TCut('n_hits_total == 4') + ('' if angle is None else 'angle_x < {0} + .05 & angle_x > {0} - .05'.format(angle))
        cut = self.Cut.generate_custom(include='tracks', add=add, name='t1{}'.format(choose(angle, '')))
        (x, y), z_ = 10 * array(self.get_tree_vec(self.get_vars(), cut, nentries=1, firstentry=self.Ana.get_events(cut)[n])), self.get_z_positions()
        self.Draw.graph(z_, x, y_range=ax_range(min(x), min(x) + .3, .3, .3), y_tit='x').Fit('pol1', 'q')
        self.Draw.graph(z_, y, y_range=ax_range(min(y), min(y) + .3, .3, .3), y_tit='y').Fit('pol1', 'q')

    def draw_n(self, show=True):
        p = self.Draw.profile(*self.get_tree_vec(var=[self.get_t_var(), 'n_tracks']), self.Bins.get_raw_time(10), 'Number of Tracks vs. Time', x_tit='Time [hh:mm}', y_tit='Number of Tracks',
                              t_ax_off=0, fill_color=Draw.FillColor, show=False)
        self.Draw(p, 'NTracksTime', draw_opt='hist', show=show, stats=set_statbox(all_stat=True, center_y=True))

    def draw_residual_vs_chi2(self, roc, mode='x', step_size=10, y_range=None, show=True):
        chi2s = arange(10, 101, step_size)
        residuals = self.get_chi2_residuals(roc, chi2s, mode)
        g = Draw.make_tgraph(x=chi2s, y=residuals, title='Tracking Resolution in {} for Plane {}'.format(mode.title(), roc))
        format_histo(g, y_tit='Residual Standard Deviation [#mum]', y_off=1.4, y_range=y_range)
        self.Draw(g, 'TrackRes', draw_opt='alp', show=show)
        return g

    def draw_residuals_vs_chi2(self, show=True):
        mg = TMultiGraph('mgtr', 'Tracking Resolution')
        leg = Draw.make_legend(y2=.41, nentries=4)
        for roc, mode in zip([1, 1, 2, 2], ['x', 'y', 'x', 'y']):
            g = self.draw_residual_vs_chi2(roc, mode, show=False)
            format_histo(g, color=self.Draw.get_color(4))
            mg.Add(g, 'pl')
            leg.AddEntry(g, 'ROC {} in {}'.format(roc, mode.title()), 'pl')
        y_range = [0, max(g.GetY()[i] for g in mg.GetListOfGraphs() for i in range(g.GetN())) * 1.1]
        format_histo(mg, x_tit='#chi^{2} [quantile]', y_tit='Residual Standard Deviation [#mum]', y_off=1.5, y_range=y_range, draw_first=True)
        self.Draw(mg, 'EventOffsets', show, draw_opt='ap', leg=leg, lm=.13)

    def draw_rel_signal_vs_angle(self, ph0=1, amax=4, amin=None, **dkw):
        a0, a1 = sorted([choose(amin, -amax), amax])
        f = Draw.make_tf1('a', lambda x: ph0 / cos(deg2rad(x)), a0, a1)
        self.Draw(f, **prep_kw(dkw, file_name='SigAngle', x_tit='Track Angle [deg]', y_tit='Relative Signal', lm=.144, y_off=1.5, y_range=[1, 1.01]))
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region BEAM
    def draw_beam_profile(self, mode='x', fit=True, res=.7, **dkw):
        p = getattr(self.Ana.draw_hitmap(res, show=False, prnt=False), f'Projection{mode.title()}')()
        self.fit_beam_profile(p) if fit else do_nothing()
        format_histo(p, title=f'Profile {mode.title()}', name=self.Draw.get_name('p'), y_tit='Number of Tracks')
        return self.Draw(p, **prep_kw(dkw, file_name=f'BeamProfile{mode.title()}', stats=set_statbox(all_stat=True, fit=fit), draw_opt=''))

    @staticmethod
    def fit_beam_profile(p):
        """fit the beam profile only at the center, because the edges are cut of by the pixel mask"""
        return Gauss(p, thresh=.05, fl=0, fh=0).fit(fl=-.2, fh=-.2)[1:]

    @save_pickle('BeamPars', suf_args='all')
    def beam_pars(self, res=None, _redo=False):
        h = self.Ana.draw_hitmap(res, show=False, prnt=False)
        return array([self.fit_beam_profile(p) for p in [h.ProjectionX(), h.ProjectionY()]])

    def draw_pixel_pdf(self, m=0, **dkw):
        s, p = self.beam_pars()[m][1].n * 1000, self.Run.Plane.P[m] * 1000  # mm -> um
        f = self.Draw.make_f('ppdf', 'gaus', 0, p, pars=[1 / (2 * s ** 2), s - p / 2, s])
        self.Draw(f, **prep_kw(dkw, x_tit='Pixel Size [#mum]', y_tit='Probability', tm=.05, y_range=[0, f(s - p)]))
    # endregion BEAM
    # ----------------------------------------

    # ----------------------------------------
    # region ANGLE
    def draw_angle(self, i=0, cut='', show_cut=False, **dkw):
        """ Shows the angle distribution of the tracks. """
        x, m = self.get_tree_vec(f'angle_{self.M[i]}', self.Cut(cut)), self.M[i].title()
        h = self.Draw.distribution(x[x > -900], **prep_kw(dkw, rf=1, lf=1, title=f'Track Angle in {m}', x_tit=f'Track Angle {m} [deg]', stats=set_statbox(all_stat=True, form='.2f')))
        self.draw_angle_cut(i) if show_cut else do_nothing()
        self.Draw.save_plots(f'TrackAngle{m}', **dkw)
        return h

    def draw_angle_cut(self, i=0):
        (m, s), n = self.Cut.calc_angle(i), self.Cut.get_config('track angle sigma', dtype=float)
        b = [Draw.box(m + f * s * n, -10, f * 100, 1e7, line_color=2, width=2, fillcolor=2, style=7, opacity=.2) for f in [-1, 1]][1:]
        self.Draw.legend(b, [f'cut ({n} #sigma)'], 'lf', y2=.773, margin=.45, w=.3)

    def draw_angles(self, **dkw):
        h = [self.draw_angle(mode, show=False, prnt=False, normalise=True) for mode in range(2)]
        self.Draw.stack(h, 'TrackAngles', [f'Angle in {m.title()}' for m in self.M], **prep_kw(dkw))
    # endregion ANGLE
    # ----------------------------------------

    # ----------------------------------------
    # region CHI2
    def chi2_var(self, m=None):
        return f'chi2_{"tracks" if m is None else self.M[m]}'

    def chi2_vars(self):
        return [self.chi2_var(i) for i in range(2)]

    def chi2_tit(self, m=None):
        return 'Tracks' if m is None else self.M[m].title()
    
    def draw_chi2(self, m=None, bin_size=None, fit=False, cut='', x_range=None, show_cut=False, **dkw):
        x, t = self.get_tree_vec(self.chi2_var(m), cut=self.Cut(cut)), self.chi2_tit(m)
        x_range = choose(x_range, [0, quantile(x[(x > -900) & (x < 100)], .99)])
        h = self.Draw.distribution(x[x > -900], self.Bins.get_chi2(bin_size), f'Chisquare {t}', x_tit='#chi^{2}', x_range=x_range, **dkw)
        fit_chi2(h, m, show=fit)
        self.draw_chi2_cut(m) if show_cut else do_nothing()
        format_statbox(h, entries=True, fit=fit)
        self.Draw.save_plots(f'Chi2{t}', **dkw)
        return h

    def draw_chi2_cut(self, m):
        b = Draw.box(0, -10, self.Cut.calc_chi2(m), 1e7, line_color=2, width=2, fillcolor=2, style=7, opacity=.2)
        self.Draw.legend([b], [f'cut ({self.Cut.get_chi2(m):d}#kern[.3]{{%}} qu.)'], 'lf', y2=.817, margin=.45, w=.3)

    def draw_chi2s(self, show=True, prnt=True):
        self.draw_chi2(fit=True, show=show, prnt=prnt)
        x_range = [0, get_last_canvas().GetUxmax()]
        self.draw_chi2(0, show_cut=True, show=show, x_range=x_range, prnt=prnt)
        self.draw_chi2(1, show_cut=True, show=show, x_range=x_range, prnt=prnt)

    def draw_xy_chi2(self, **dkw):
        h = [self.draw_chi2(m, show=False, normalise=True, **rm_key(dkw, 'show')) for m in Tracks.M]
        self.Draw.stack(h, 'Chi-squares', Tracks.M, **prep_kw(dkw))
        
    def get_chi2s(self, e=18, g=None):
        x, y = self.get_z_positions(), self.get_plane_hits(local=False, add_cut='cluster_size == 1')
        e, g = e * Plane.P / sqrt(12), choose(g, sqrt(5))
        return [polyfit(x, y[i].T, deg=1, full=True, w=1 / array([1e-5, e[i] / g, e[i] / g * 1.2, e[i]]))[1] / 2 for i in range(2)]

    @save_pickle('ChiQ', suf_args='all')
    def get_chi2_cut(self, q=None, e: Any = 15, g=None, _redo=False):
        return quantile(self.get_chi2s(e[0] if is_iter(e) else e, g)[0], q=choose(q, self.Cut.get_chi2() / 100))

    def get_chi2_x(self, e=18, g=None):
        return self.get_chi2s(e, g)[0]

    def get_chi2_y(self, e=18, g=None):
        return self.get_chi2s(e, g)[1]

    def draw_fit_chi2s(self, e=18, g=None, fit=False, show_cut=False, **dkw):
        h = [self.Draw.distribution(x, normalise=True, x0=0, x_tit='#chi^{2} / DOF', show=False, file_name=f'chi_{["x", "y"][i]}') for i, x in enumerate(self.get_chi2s(e, g))]
        s = self.Draw.stack(h, 'FitChi2', None, **dkw)
        self.Draw.make_tf1('fd', lambda x: TMath.GammaDist(x, 1, 0, 2) * h[1].GetBinWidth(1), 0, 20, npx=500).Draw('same') if fit else do_nothing()
        b = self.Draw.box(self.get_chi2_cut(e=e, g=g), -1, 100, 1, line_color=2, width=2, fillcolor=2, style=7, opacity=.2) if show_cut else do_nothing()
        self.Draw.legend(h + [b], ['x', 'y', 'cut'], styles=['l', 'l', 'lf'], margin=.45, w=.2)
        format_histo(s, **prep_kw(dkw, lw=2))
        self.Draw.save_plots('FitChi2')
    # endregion CHI2
    # ----------------------------------------

    # ----------------------------------------
    # region RESIDUALS
    def draw_residual(self, roc, mode='x', cut='', fit=False, ret_res=False, **dkw):
        x = self.get_tree_vec(f'{self.get_res_var(mode)}[{roc}]', self.Cut(cut)) * 1e4  # convert to [um]
        tit = f'{mode.title() if mode else ""} Residuals for Plane {roc}'
        h = self.Draw.distribution(x, **prep_kw(dkw, show=False, title=tit, x_tit='Distance [#mum]', normalise=True))
        f = self.fit_residual(h)
        self.Draw(h, **prep_kw(dkw, file_name=f'{mode.title() if mode else ""}ResidualRoc{roc}', leg=f.Fit, y_off=1.5, lm=.14, stats=set_statbox(fit=fit, all_stat=True, fit_opt=10)))
        return f[2] if ret_res else h

    def get_residual_fit(self, plane=0, m='x', cut=None):
        """ :return FWHM fit of the residual with a Gaussian."""
        return fit_fwhm(self.draw_residual(plane, m, cut, show=False))[1:]

    def draw_xy_residual(self, roc, cut='', f=.5, show_cut=False, **dkw):
        x, y = array(self.get_tree_vec([f'{self.get_res_var(m)}[{roc}]' for m in ['x', 'y']], self.Cut(cut))) * 1e4  # convert to [um]
        ((mx, sx), (my, sy)), n = self.Cut.calc_res_sigmas() * 1e4, self.Cut.get_config('rhit sigma', dtype=float)
        rcut = Draw.ellipse(sx * n, sy * n, mx, my, show=False) if show_cut else None
        return self.Draw.histo_2d(x, y, bins.find(x, f, f, nbins=2) + bins.find(y, f, f, nbins=2), **prep_kw(dkw, x_tit='Residual in X [#mum]', y_tit='Residual in Y [#mum]', leg=rcut,
                                                                                                             file_name=f'XYResidual{roc}'))

    def draw_unbiased_residual(self, roc=0, mode='x', cut='', fit=False, **dkw):
        """ fit the track without the plane under test and calculate residuals. """
        x, y = self.get_z_positions()[:self.Run.NTelPlanes], self.get_plane_hits(local=False, add_cut=cut)[0 if mode == 'x' else 1].T
        fits = polyfit(delete(x, roc), delete(y, roc, axis=0), deg=1)
        v = (polyval(fits, x[roc]) - y[roc]) * 1e3  # to mm -> um
        tit = 'Unbiased Residuals in {} for Plane {}'.format(mode.title(), roc)
        h = self.Draw.distribution(v, bins.make(-1000, 1000, 2), tit, **prep_kw(dkw, y_off=1.5, x_tit='Distance [#mum]', normalise=True, lm=.14, stats=set_statbox(fit=fit, all_stat=True, fit_opt=10)))
        res = mean_sigma(v, err=0)[1] if 'chi2' in self.Cut.get_name(cut) else self.fit_residual(h, show=fit)[2]
        return res if fit else h

    def draw_raw_residuals(self, roc=None, steps=0, show=True):
        x, y, z_positions = self.get_plane_hits()
        dx, dy, da = self.calculate_residuals(x, y, z_positions, steps)
        x_bins, y_bins = arange(-1, 1.03, .03),  arange(-1, 1.02, .02)
        histos = [TH2F('h{}'.format(i), 'Raw Residuals ROC {}'.format(i), x_bins.size - 1, x_bins, y_bins.size - 1, y_bins) for i in range(dx.shape[0])]
        for h, ix, iy in zip(histos, dx, dy):
            fill_hist(h, ix, iy)
            format_histo(h, x_tit='dX [mm]', y_tit='dY [mm]', y_off=1.3, pal=53)
        if roc is not None:
            self.Draw(histos[roc], draw_opt='colz', lm=.12, rm=.16, show=show, logz=True, stats=set_entries())
        else:
            c = self.Draw.canvas('Raw Residuals', w=1.5, h=1.5, show=show, divide=(2, 2))
            for i, h in enumerate(histos, 1):
                self.Draw(h, draw_opt='colz', lm=.12, rm=.16, show=show, canvas=c.cd(i), logz=True, stats=set_entries())
        return dx, dy, da

    def draw_mean_residuals(self, roc=0, mode='x', steps=10):
        x, y, z_positions = self.get_plane_hits()
        values = array([[abs(mean(vec[roc])) for vec in self.calculate_residuals(x, y, z_positions)][:2] for _ in range(steps)])
        g = self.Draw.graph(arange(steps), values[:, {'x': 0, 'y': 1}[mode]], title='Mean Residual in {} for Each Step'.format(mode.title()), draw_opt='ap', lm=.13, logy=True)
        format_histo(g, x_tit='Step #', y_tit='Mean Residual [mm]', y_off=1.5, x_range=[-.5, steps - .5], ndivx=505)
        return array([[sum(values[:, 0])], [sum(values[:, 1])]])

    def draw_residual_angles(self, roc=0, steps=10):
        x, y, z_positions = self.get_plane_hits()
        angles = array([abs(rad2deg(tan(self.calculate_residuals(x, y, z_positions)[-1][roc]))) for _ in range(steps)])
        g = self.Draw.graph(arange(steps), angles, title='Residual Angles', draw_opt='ap', lm=.13)
        format_histo(g, x_tit='Step #', y_tit='Residual Angle [#circ]', y_off=1.5, x_range=[-.5, steps - .5], ndivx=505)

    def calculate_residuals(self, x, y, z_positions, steps=0):
        dx, dy, da = None, None, None
        for _ in range(steps + 1):
            dx, dy = self.calc_residual_step(x, y, z_positions)
            da = self.align(x, y, dx, dy)
        return dx, dy, da

    @staticmethod
    def calc_residual_step(x, y, z_positions):
        x_fits, y_fits = [polyfit(z_positions, vec, 1) for vec in [x, y]]
        return [polyval(fits, z_positions.reshape((x.shape[0], 1))) - vec for fits, vec in [(x_fits, x), (y_fits, y)]]

    @staticmethod
    def align(x, y, dx, dy):
        b = [arange(-4, 4.01, .2), arange(-1, 1.02, .03), arange(-4, 4.15, .15), arange(-1, 1.02, .02)]
        histos = [[TH2F('h{}{}'.format(i, j), 'rotres', b[i].size - 1, b[i], b[i + 1].size - 1, b[i + 1]) for i in arange(0, 4, 2)] for j in range(x.shape[0])]
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
        x += mean(dx, axis=1).reshape((x.shape[0], 1))
        y += mean(dy, axis=1).reshape((x.shape[0], 1))
        return angles

    @staticmethod
    def fit_residual(h, r=.2, show=True):
        ymax = FitRes(h.Fit('gaus', 'qs0', '', -1000, 1000))[0].n
        fit_range = [h.GetBinCenter(f(ymax * r)) for f in [h.FindFirstBinAbove, h.FindLastBinAbove]]
        return Gauss(h, fit_range, npx=500, fl=0, fh=0).fit(draw=show)

    @staticmethod
    def fit_residual_old(h, show=True):
        fit = Draw.make_f('f', 'gaus(0) + gaus(3)', npx=500)
        sigma = (get_fwhm(h) / (2 * sqrt(2 * log(2)))).n
        fit.SetParLimits(2, sigma / 2, 5 * sigma)
        fit.SetParLimits(5, sigma / 2, 5 * sigma)
        fit.SetParameters(h.GetMaximum() / 10, 0, sigma * 5, h.GetMaximum(), 0, sigma)
        fit.SetParNames('c1', 'mean1', '#sigma1', 'c2', 'mean2', '#sigma2')
        fit.SetNpx(500)
        h.Fit(fit, 'q{}'.format('' if show else 0))
        f2 = TF1('f2', 'gaus')
        f2.SetParameters(fit.GetParameters())
        f2.SetLineStyle(2)
        if show:
            f2.Draw('same')
        h.GetListOfFunctions().Add(f2)
        return min(abs(fit.GetParameter(i)) for i in [2, 5])

    def calc_res(self, x1=0, x2=0, x3=1, mode='x'):
        z_ = sorted(self.get_z_positions())
        x = array([0, x1, x2, x3]) * (self.Run.Plane.PX if mode == 'x' else self.Run.Plane.PY)
        g = self.Draw.graph(z_, [ufloat(i, .15/sqrt(12)) for i in x])
        g.Fit('pol1')
        f = polyfit(z_, x, 1)
        print(polyval(f, 0) * 1000)
        return f
    # endregion RESIDUALS
    # ----------------------------------------

    # ----------------------------------------
    # region RESOLUTION
    def draw_resolution(self, mode='x', cut='', n=1e5, unbias=False, redo=False, **dkw):
        z_ = self.get_z_positions()[:self.Run.NTelPlanes]
        r = uarr2n(self.get_residuals(mode=mode, cut=cut, unbias=unbias, redo=redo))
        x_range = -20, max(z_) + 20
        self.Draw.graph(z_, [ufloat(0, ex) for ex in r], **prep_kw(dkw, x_tit='z [mm]', y_tit='{} [#mum]'.format(mode.lower()), y_range=[-195, 195], x_range=x_range))
        x = array([norm.rvs(0, ir, size=int(n)) for ir in r])
        fits = array(polyfit(z_, x, deg=1))
        p = linspace(*x_range, 100)
        z0 = polyval(fits, p.reshape(p.size, 1))
        ex = array([ufloat(m, s) for m, s in [mean_sigma(iz, err=False) for iz in z0]])
        g = self.Draw.make_tgraph(p, ex, fill_color=634, opacity=.5)
        g.Draw('e3')
        xm, ym = mean(z_), mean_sigma(polyval(fits, mean(z_)))[1].n
        Draw.arrow(xm, xm, -ym, ym, width=2, opt='<|>', size=.02)
        Draw.tlatex(xm, ym + 10, f'{ym:2.0f}#kern[.1]{{#mum}}', align=21)
        return min(e.s for e in ex)

    @quiet
    def draw_resolution_vs_chi2(self, m=0, s=10, **dkw):
        x = arange(max(s, 6), 100.1, s, dtype='i')
        PBAR.start(x.size, counter=True)
        r = self.get_mean_resolution if m is None else partial(self.get_resolution, mode=self.M[m])
        y = [r(cut=self.Cut.make(f'chi{q}', self.Cut.generate_chi2s(q)), unbias=True) for q in x]
        return self.Draw.graph(x, y, **prep_kw(dkw, x_tit='Quantile [%]', y_tit='Resolution [#mum]', y_range=[0, max(y) * 1.2]))

    def draw_xy_res_vs_chi2(self, s=10, **dkw):
        g = [self.draw_resolution_vs_chi2(m, s=s, show=False) for m in range(2)]
        self.Draw.multigraph(g, 'TelRes', self.M, **prep_kw(dkw))
    # endregion RESOLUTION
    # ----------------------------------------


def get_chi2(k=2, s=1):
    return Draw.make_tf1('fd', lambda x: TMath.GammaDist(x, k / 2, 0, 2) * s, 0, 20, npx=500)


def fit_chi2(h, mode, show=True):
    fit = TF1('f', '[0]*TMath::GammaDist(x, {ndf}/2, 0, 2)'.format(ndf=4 if mode == 'tracks' else 2))
    h.Fit(fit, 'qs{}'.format('' if show else 0))
    return fit
