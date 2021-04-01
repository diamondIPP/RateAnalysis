#!/usr/bin/env python
# --------------------------------------------------------
#       analysis class for the telescope tracks
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TCut, TF1, TMultiGraph, THStack
from numpy import log, genfromtxt, rad2deg, polyfit, polyval, tan, delete, deg2rad
from src.sub_analysis import SubAnalysis
from helpers.draw import *
from scipy.stats import norm


class Tracks(SubAnalysis):

    def __init__(self, analysis):
        super().__init__(analysis, sub_dir=analysis.Run.Number, pickle_dir='Tracks', dut=False)

    # ----------------------------------------
    # region GET
    def get_efficiency(self):
        return calc_eff(values=self.get_tree_vec(var='n_tracks'))

    def get_angles(self, mode='x'):
        a = self.get_tree_vec(var='angle_{}'.format(mode))
        return a[a != -999]

    def get_mean_angle(self, mode='x', redo=False):
        def f():
            return mean_sigma(self.get_angles(mode))
        return do_pickle(self.make_simple_pickle_path('TrackAngle{}'.format(mode)), f, redo=redo)

    def get_residual(self, roc=0, mode='x', cut='', redo=False):
        pickle_path = self.make_simple_pickle_path('Res{}{}'.format(mode.title(), roc), self.Cut.get_name(cut))
        return do_pickle(pickle_path, self.draw_residual, redo=redo, roc=roc, cut=cut, ret_res=True, show=False, mode=mode)

    def get_unbiased_residual(self, roc=0, mode='x', cut='', redo=False):
        pickle_path = self.make_simple_pickle_path('URes{}{}'.format(mode.title(), roc), self.Cut.get_name(cut))
        return do_pickle(pickle_path, self.draw_unbiased_residual, redo=redo, roc=roc, cut=cut, fit=True, show=False, mode=mode)

    def get_residuals(self, mode='x', cut='', unbias=False, redo=False):
        return array([(self.get_unbiased_residual(roc, mode, cut, redo) if unbias else self.get_residual(roc, mode, cut, redo)) for roc in range(self.NRocs)])

    def get_resolution(self, mode='x', cut=''):
        return self.draw_resolution(mode, cut=cut, show=False)

    def get_chi2_residual(self, roc, chi2, mode='x', redo=False):
        def f():
            self.Cut.set_chi2(chi2)
            values = self.get_tree_vec(var='residuals_{m}[{r}]*1e4'.format(m=mode, r=roc), cut=self.Cut())
            return mean_sigma(values)[1]
        return do_pickle(self.make_simple_pickle_path('Res{}'.format(mode.title()), chi2, dut=roc), f, redo=redo)

    def get_chi2_residuals(self, roc, chi2s=None, mode='x'):
        chi2s = choose(chi2s, arange(10, 101, 10))
        return [self.get_chi2_residual(roc, chi2, mode) for chi2 in chi2s]

    def get_z_positions(self, e=0):
        x = genfromtxt(join(self.Run.Converter.TrackingDir, 'data', 'alignments.txt'), usecols=[0, 2, 6])
        return array([ufloat(ix, e) if e else ix for ix in x[(x[:, 0] == self.Run.Converter.TelescopeID) & (x[:, 1] > -1)][:, 2] * 10])  # [mm]

    @staticmethod
    def get_vars(local=False):
        return ['cluster_{}pos_{}'.format(n, 'local' if local else 'tel') for n in ['x', 'y']]

    def get_plane_hits(self, local=True, add_cut=''):
        t = self.info('getting plane hits...', endl=False)
        self.Tree.SetEstimate(self.Cut.get_size('tracks', excluded=False) * self.NRocs)
        x, y = self.get_tree_vec(self.get_vars(local), self.Cut['tracks'] + add_cut)
        x, y = [array(split(vec, arange(self.NRocs, x.size, self.NRocs))).T * 10 for vec in [x, y]]
        self.add_to_info(t)
        return x, y, self.get_z_positions()
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw(self, n=0, angle=None):
        add = TCut('total_hits == 4') + ('' if angle is None else 'angle_x < {0} + .05 & angle_x > {0} - .05'.format(angle))
        cut = self.Cut.generate_custom(include='tracks', add=add, name='t1{}'.format(choose(angle, '')))
        (x, y), z_ = 10 * array(self.get_tree_vec(self.get_vars(), cut, nentries=1, firstentry=self.Ana.get_events(cut)[n])), self.get_z_positions()
        self.Draw.graph(z_, x, y_range=ax_range(min(x), min(x) + .3, .3, .3), y_tit='x').Fit('pol1', 'q')
        self.Draw.graph(z_, y, y_range=ax_range(min(y), min(y) + .3, .3, .3), y_tit='y').Fit('pol1', 'q')

    def draw_n(self, show=True):
        p = self.Draw.profile(*self.get_tree_vec(var=[self.get_t_var(), 'n_tracks']), self.Bins.get_raw_time(10), 'Number of Tracks vs. Time', x_tit='Time [hh:mm}', y_tit='Number of Tracks',
                              t_ax_off=0, fill_color=Draw.FillColor, show=False)
        self.Draw(p, 'NTracksTime', draw_opt='hist', show=show, stats=set_statbox(all_stat=True, center_y=True))

    def draw_chi2(self, mode=None, fit=False, x_range=None, cut='', normalise=None, show=True, save=True, prnt=True, show_cut=False):
        mode = 'tracks' if mode is None else mode
        h = Draw.make_histo('Chisquare {}'.format(mode.title()), [500, 0, 100])
        self.Tree.Draw('chi2_{}>>{}'.format(mode, h.GetName()), TCut('n_tracks > 0') + TCut(cut), 'goff')
        y_tit = '{} of Entries'.format('Number' if normalise is None else 'Percentage')
        format_histo(h, x_tit='#chi^{2}', y_tit=y_tit, y_off=2, x_range=choose(x_range, [0, get_quantile(h, .99)]), normalise=normalise, fill_color=Draw.FillColor)
        fit_chi2(h, mode, show=fit)
        self.Draw(h, show=show, prnt=prnt, lm=.13, stats=set_statbox(entries=True, fit=fit))
        self.draw_chi2_cut(mode, show=show_cut)
        self.Draw.save_plots('Chi2{0}'.format(mode.title()), show=show, save=save, prnt=prnt)
        return h

    def draw_chi2_cut(self, mode, show=True):
        if show:
            chi2 = self.Cut.calc_chi2(mode)
            line = Draw.vertical_line(chi2, -100, 1e6, style=7, w=2, color=2)
            legend = Draw.make_legend(.75, y2=.83, nentries=1, margin=.35)
            legend.AddEntry(line, 'cut ({}%)'.format(self.Cut.get_chi2(mode)), 'l')
            legend.Draw()

    def draw_chi2s(self, show=True, prnt=True):
        self.draw_chi2(fit=True, show=show, prnt=prnt)
        x_range = [0, get_last_canvas().GetUxmax()]
        self.draw_chi2('x', show_cut=True, show=show, x_range=x_range, prnt=prnt)
        self.draw_chi2('y', show_cut=True, show=show, x_range=x_range, prnt=prnt)

    def draw_angle(self, mode='x', cut=None, show_cut=False, normalise=None, show=True, prnt=True):
        """ Shows the angle distribution of the tracks. """
        h = Draw.make_histo('Track Angle Distribution in {}'.format(mode.title()), self.Bins.get_angle())
        self.Tree.Draw('angle_{}>>{}'.format(mode, h.GetName()), self.Cut('angle_{}>-900'.format(mode) if cut is None else cut), 'goff')
        y_tit = '{} of Entries'.format('Number' if normalise is None else 'Percentage')
        format_histo(h, x_tit='Track Angle {} [deg]'.format(mode.title()), y_tit=y_tit, y_off=2, normalise=normalise, fill_color=Draw.FillColor)
        self.Draw(h, show=show, lm=.14, prnt=prnt, stats=None)
        self.draw_angle_cut(show=show_cut)
        self.Draw.save_plots('TrackAngle{mod}'.format(mod=mode.upper()), prnt=prnt)
        return h

    def draw_angle_cut(self, show=True):
        if show:
            xmax = -self.Cut.get_track_angle()
            line = Draw.vertical_line(-xmax, -100, 1e6, style=7, w=2, color=2)
            Draw.vertical_line(xmax, -100, 1e6, style=7, w=2, color=2)
            legend = Draw.make_legend(.65, y2=.73, nentries=1, margin=.35, scale=1.3)
            legend.AddEntry(line, 'cut ({} deg)'.format(xmax), 'l')
            legend.Draw()

    def draw_angles(self, show=True, prnt=True):
        histos = [self.draw_angle(mode, show=False, prnt=False) for mode in ['x', 'y']]
        leg = Draw.make_legend(nentries=2, w=.25)
        stack = THStack('has', 'Track Angles')
        for h in histos:
            format_histo(h, stats=False, color=self.Draw.get_color(2))
            leg.AddEntry(h, 'Angle in {}'.format(h.GetTitle()[-1]), 'l')
            stack.Add(h)
        self.Draw(stack, 'TrackAngles', lm=.14, leg=leg, draw_opt='nostack', show=show, prnt=prnt)

    def draw_residual_vs_chi2(self, roc, mode='x', step_size=10, y_range=None, show=True):
        chi2s = arange(10, 101, step_size)
        residuals = self.get_chi2_residuals(roc, chi2s, mode)
        g = Draw.make_tgrapherrors(x=chi2s, y=residuals, title='Tracking Resolution in {} for Plane {}'.format(mode.title(), roc))
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

    def draw_rel_signal_vs_angle(self, a_max=4):
        f = Draw.make_tf1('a', lambda x: 1 / cos(deg2rad(x)), -a_max, a_max)
        format_histo(f, x_tit='Track Angle [deg]', y_tit='Relative Signal', y_off=1.6, y_range=[1, 1.01])
        self.Draw(f, lm=.121)

    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region RESIDUALS
    def draw_residual(self, roc, mode='x', cut='', x_range=None, fit=False, ret_res=False, show=True):
        x = self.get_tree_vec(var='residuals{}[{}]'.format('_{}'.format(mode.lower()) if mode else '', roc), cut=self.Cut(cut)) * 1e4  # convert to [um]
        tit = '{} Residuals for Plane {}'.format(mode.title(), roc)
        h = self.Draw.distribution(x, make_bins(-1000, 1000, 2), tit, y_off=2.0, x_tit='Distance [#mum]', x_range=x_range, show=False, normalise=True)
        res = self.fit_residual(h, show=fit)
        self.Draw(h, '{}ResidualsRoc{}'.format(mode.title(), roc), show, lm=.14, stats=set_statbox(fit=fit, all_stat=True))
        return res if ret_res else h

    def draw_unbiased_residual(self, roc=0, mode='x', cut='', x_range=None, fit=False, show=True):
        """ fit the track without the plane under test and calculate residuals. """
        x, y, z_ = self.get_plane_hits(local=False, add_cut=cut)
        var = x if mode == 'x' else y
        fits = polyfit(delete(z_, roc), delete(var, roc, axis=0), deg=1)
        v = (polyval(fits, z_[roc]) - var[roc]) * 1e3  # to mm -> um
        tit = 'Unbiased Residuals in {} for Plane {}'.format(mode.title(), roc)
        h = self.Draw.distribution(v, make_bins(-1000, 1000, 2), tit, y_off=2.0, x_tit='Distance [#mum]', x_range=x_range, show=show, normalise=True, lm=.14)
        res = mean_sigma(v, err=0)[1] if 'chi2' in self.Cut.get_name(cut) else self.fit_residual(h, show=fit)
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

    def draw_resolution(self, mode='x', cut='', n=1e5, y_range=None, unbias=False, show=True):
        z_ = self.get_z_positions()
        r = self.get_residuals(mode=mode, cut=cut, unbias=unbias)
        x_range = -20, max(z_) + 20
        self.Draw.graph(z_, [ufloat(0, ex) for ex in r], x_tit='z [mm]', y_tit='{} [#mum]'.format(mode.lower()), y_range=choose(y_range, [-200, 200]), x_range=x_range, show=show)
        x = array([norm.rvs(0, ir, size=int(n)) for ir in r])
        fits = array(polyfit(z_, x, deg=1))
        p = linspace(*x_range, 100)
        z0 = polyval(fits, p.reshape(p.size, 1))
        ex = array([ufloat(m, s) for m, s in [mean_sigma(iz, err=False) for iz in z0]])
        g = self.Draw.make_tgrapherrors(p, ex, fill_color=634, opacity=.5)
        g.Draw('e3')
        xm, ym = mean(z_), mean_sigma(polyval(fits, mean(z_)))[1].n
        Draw.arrow(xm, xm, -ym, ym, width=2, opt='<|>', size=.02)
        Draw.tlatex(xm, ym + 10, '{:2.0f}#mum'.format(ym), align=21)
        return min(e.s for e in ex)

    @staticmethod
    def calc_residual_step(x, y, z_positions):
        x_fits, y_fits = [polyfit(z_positions, vec, 1) for vec in [x, y]]
        return [polyval(fits, z_positions.reshape((x.shape[0], 1))) - vec for fits, vec in [(x_fits, x), (y_fits, y)]]

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
        # format_statbox(entries=True, x=.82)
        # self.draw_histo(histos[0][0], draw_opt='colz', rm=.16)
        x += mean(dx, axis=1).reshape((x.shape[0], 1))
        y += mean(dy, axis=1).reshape((x.shape[0], 1))
        return angles

    @staticmethod
    def fit_residual(h, show=True):
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


def fit_chi2(h, mode, show=True):
    fit = TF1('f', '[0]*TMath::GammaDist(x, {ndf}/2, 0, 2)'.format(ndf=4 if mode == 'tracks' else 2))
    h.Fit(fit, 'qs{}'.format('' if show else 0))
    return fit
