# --------------------------------------------------------
#       sub-analysis class for efficiency of the pixel detectors
# created on December 17th 2021 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from helpers.utils import save_pickle, deepcopy, file_exists, do_nothing, update_pbar, array
from plotting.draw import FitRes, find_bins, choose, prep_kw, calc_eff, quantile, arange, get_graph_x, set_statbox
from plotting.fit import Erf, Draw
from src.binning import Bins, make_bins, Plane
from src.cut import Cut
from pixel.analysis import PixAnalysis


class Efficiency(PixAnalysis):

    def __init__(self, parent):  # noqa
        self.__dict__.update(parent.__dict__)
        self.PickleSubDir = 'Efficiency'

        self.UseRhit = False
        self.Cut = deepcopy(self.Cut).remove('rhit' if self.UseRhit else None, ret=True)

    def get_values(self, cut=None):
        return self.get_tree_vec(self.get_var(), self.Cut(cut), dtype='?')

    def get_var(self, sigma=None, plane=None, percent=False):
        var = f'(n_clusters[{choose(plane, self.N)}] > 0){" * 100" if percent else ""}'
        return Cut.to_string(self.Cut.generate_rhit(sigma) + Cut.make('', var)) if self.UseRhit else var

    @update_pbar
    @save_pickle(suf_args='all', field='UseRhit')
    def get(self, cut=None, _redo=False):
        return calc_eff(values=self.get_values(cut))

    def draw(self, cut=None, bin_size=None, rel_time=False, **dkw):
        x, e = self.get_tree_vec([self.get_t_var(), self.get_var()], choose(cut, self.Cut()))
        g = self.Draw.efficiency(x, e, self.Bins.get_time(choose(bin_size, 1000 if x.size // 20 < 1000 or x.size / 1000 < 20 else x.size // 20)), show=False)  # min bin size of 1000 max 20 points
        fit, tit = FitRes(g.Fit('pol0', 'qs')), 'Hit Efficiency'
        self.Draw(g, **prep_kw(dkw, title=tit, **self.get_t_args(rel_time), y_range=[-5, 115], gridy=True, draw_opt='apz', stats=set_statbox(fit=True, form='.2f', center_y=True), file_name='HitEff'))
        return fit if fit.Parameter(0) is not None else 0

    @save_pickle('Map', suf_args='all')
    def get_map(self, res=None, fid=False, cut=None, local=False, _redo=False):
        x, y, zz = self.get_tree_vec(self.get_track_vars(local=local) + [self.get_var()], self.Cut(cut) if cut or fid else self.Cut.exclude('fiducial'))
        tit, (xtit, ytit), ztit = 'Efficiency Map', [f'Track Position {i} [mm]' for i in ['X', 'Y']], 'Efficiency [%]'
        return self.Draw.prof2d(x, y, zz * 100, Bins.get_global(res), tit, x_tit=xtit, y_tit=ytit, z_tit=ztit, leg=self.Cut.get_fid(), show=False)

    def draw_map(self, res=None, fid=False, cut=None, local=False, redo=False, **dkw):
        self.Draw(self.get_map(res, fid, cut, local, _redo=redo), **prep_kw(dkw, file_name='EfficiencyMap'))

    def get_mod_vars(self, mx=1, my=1, ox=0, oy=0, cut=None, max_angle=None, expand=True, **kwargs):
        return super(Efficiency, self).get_mod_vars(mx, my, ox, oy, zvar=self.get_var(percent=True), cut=cut, expand=expand)

    def draw_in_cell(self, nbins=None, ox=0, oy=0, cut=None, max_angle=None, **dkw):
        """ in 3D cell"""
        return self.draw_in(self.DUT.PX, self.DUT.PY, ox, oy, nbins, cut, max_angle, **prep_kw(dkw, title='In Cell Efficiency', z_tit='Efficiency [%]'))

    def draw_in_pixel(self, nbins=None, ox=0, oy=0, cut=None, max_angle=None, **dkw):
        """ in pixel of ROC"""
        return self.draw_in(Plane.PX, Plane.PY, ox, oy, nbins, cut, max_angle, **prep_kw(dkw, title='In Pixel Efficiency', z_tit='Efficiency [%]'))

    def draw_vs_chi2(self, **dkw):
        x, e = self.get_tree_vec(['chi2_tracks', self.get_var()], self.Cut.exclude('chi2_x', 'chi2_y'))
        self.Draw.efficiency(x, e, find_bins(x, lfac=0, lq=0), title='Efficiency vs Chi2', **prep_kw(dkw, x_tit='Track #chi^{2}'))

    def draw_vs_chi2_cut(self, step=5, **dkw):
        cx, cy, e = self.get_tree_vec(['chi2_x', 'chi2_y', self.get_var()], self.Cut.exclude('chi2_x', 'chi2_y'))
        x = arange(100, step=step) + step
        qx, qy = [quantile(c, x / 100) for c in [cx, cy]]
        y = [calc_eff(values=e[(cx < qx[i]) & (cy < qy[i])]) for i in range(qx.size)]
        return self.Draw.graph(x, y, 'Efficiency vs Chi2 Cut', **prep_kw(dkw, x_tit='Cut Quantile [%]', y_tit='Efficiency [%]'))

    def draw_vs_trigphase(self, **kwargs):
        x, e = self.get_tree_vec(['trigger_phase[1]', self.get_var()], self.Cut.exclude('trigger_phase'))
        return self.Draw.efficiency(x, e, make_bins(-.5, 10), **prep_kw(kwargs, title='Trigger Phase Efficiency', x_tit='Trigger Phase', x_range=[-1, 10], draw_opt='bap'))

    def draw_vs_cuts(self, cuts=None, short=False, redo=False, **dkw):
        cuts = choose(cuts, self.Cut.get_consecutive(short))
        self.PBar.start(len(cuts), counter=True) if redo or not file_exists(self.make_simple_pickle_path(suf=f'{[*cuts.values()][-1].GetName()}_{self.UseRhit}')) else do_nothing()
        x, y = arange(len(cuts)), array([self.get(cut, _redo=redo) for cut in cuts.values()])
        return self.Draw.graph(x, y, title='Efficiency for Consecutive Cuts', y_tit='Efficiency [%]', **prep_kw(dkw, draw_opt='ap', gridy=True, x_range=[-1, len(y)], bin_labels=cuts.keys()))

    def draw_vs_angle(self, **dkw):
        x, y, e = self.get_tree_vec(['angle_x', 'angle_y', self.get_var()], self.Cut.exclude('track angle x', 'track angle y'))
        self.Draw.prof2d(x, y, e * 100, **prep_kw(dkw, x_tit='Angle X', y_tit='Angle Y', z_tit='Efficiency [%]'))

    def _find_alignment(self, p, w, show=False):
        g = self.Draw.make_graph_from_profile(p)
        (x0, x1), m = get_graph_x(g, err=False)[[0, -1]], p.GetMean()
        self.Draw(g, x_tit='Coordinate', y_tit='Efficiency [%]') if show else do_nothing()
        f0, f1 = Erf(g, [x0, m]).fit(show=show).Fit, Erf(g, [m, x1]).fit(show=show).Fit
        return self.Draw.make_tf1(None, lambda x: f0(x) - f1(x + w), x0, m).GetX(0)

    @save_pickle('Align', suf_args='[0]')
    def find_alignment(self, res=.5, show=False, _redo=False):
        e = self.get_map(res)
        return self._find_alignment(e.ProfileX(), self.DUT.W, show), self._find_alignment(e.ProfileY(), self.DUT.H, show)

    def draw_pdf(self, k=9, n=10, **dkw):
        (med, e0, e1), m = calc_eff(k, n) / 100, (k + 1) / (n + 2)
        f_str = 'TMath::Factorial([1] + 1) / (TMath::Factorial([0]) * TMath::Factorial([1] - [0])) * x ** [0] * (1 - x) ** ([1] - [0])'
        f = self.Draw(self.Draw.make_f('epdf', f_str, pars=[k, n], npx=1000), 'Efficiency PDF', **prep_kw(dkw, x_tit='Efficiency', y_tit='Probability'))
        ff = self.Draw(self.Draw.make_f(None, f_str, med - e0, med + e1, pars=[k, n], lw=0, npx=1000, fill_style=1001, fill_color=1, opacity=.3), draw_opt='same')
        ls = [Draw.vertical_line(m, 0, f(m), style=7), Draw.vertical_line(med, 0, f(med))]
        _ = [Draw.arrow(med - e0, med, *[f(med) / 4] * 2, size=.01), Draw.arrow(med + e1, med, *[f(med) / 4] * 2, size=.01), f.Draw('same')]
        Draw.legend(ls + [ff], ['mean', 'mode', 'CL = 1#sigma'], ['l', 'l', 'f'], left=True)
