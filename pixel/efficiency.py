# --------------------------------------------------------
#       sub-analysis class for efficiency of the pixel detectors
# created on December 17th 2021 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from helpers.utils import save_pickle, deepcopy, file_exists, do_nothing, update_pbar, array
from plotting.draw import FitRes, find_bins, choose, prep_kw, calc_eff, quantile, arange, get_graph_x
from plotting.fit import Erf
from src.binning import Bins, Plane, make_bins
from src.cut import Cut
from src.sub_analysis import SubAnalysis


class Efficiency(SubAnalysis):

    def __init__(self, pix_analysis):
        super().__init__(pix_analysis, pickle_dir='Efficiency')
        self.N = self.Ana.N
        self.UseRhit = False
        self.Cut = deepcopy(self.Ana.Cut).remove('rhit' if self.UseRhit else None, ret=True)

    def get_values(self, cut=None):
        return self.get_tree_vec(self.get_var(), self.Cut(cut), dtype='?')

    def get_var(self, sigma=None, plane=None):
        var = f'n_clusters[{choose(plane, self.N)}] > 0'
        return Cut.to_string(self.Cut.generate_rhit(sigma) + Cut.make('', var)) if self.UseRhit else var

    @update_pbar
    @save_pickle(suf_args='all', field='UseRhit')
    def get(self, cut=None, _redo=False):
        return calc_eff(values=self.get_values(cut))

    def draw(self, cut=None, bin_size=None, rel_time=False, **kwargs):
        (x, e), bins = self.get_tree_vec([self.get_t_var(), self.get_var()], choose(cut, self.Cut())), self.Bins.get_time(bin_size, cut)
        g = self.Draw.efficiency(x, e, bins, **prep_kw(kwargs, title='Hit Efficiency', **self.get_t_args(rel_time), y_range=[-5, 115], gridy=True, draw_opt='apz'))
        fit = FitRes(g.Fit('pol0', 'qs'))
        self.Draw.stats(fit, width=.35, y2=.35, names=['Efficiency'])
        self.Draw.preliminary()
        self.Draw.save_plots('HitEfficiency', **kwargs)
        return fit if fit.Parameter(0) is not None else 0

    @save_pickle('Map', suf_args='all')
    def get_map(self, res=None, fid=False, cut=None, local=False, _redo=False):
        x, y, zz = self.get_tree_vec(self.Ana.get_track_vars(local=local) + [self.get_var()], self.Cut(cut) if cut or fid else self.Cut.exclude('fiducial'))
        tit, (xtit, ytit), ztit = 'Efficiency Map', [f'Track Position {i} [mm]' for i in ['X', 'Y']], 'Efficiency [%]'
        return self.Draw.prof2d(x, y, zz * 100, Bins.get_global(res), tit, x_tit=xtit, y_tit=ytit, z_tit=ztit, leg=self.Cut.get_fid(), show=False)

    def draw_map(self, res=None, fid=False, cut=None, local=False, redo=False, **dkw):
        self.Draw(self.get_map(res, fid, cut, local, _redo=redo), **prep_kw(dkw, file_name='Efficiency Map'))

    def get_mod_vars(self, mx, my, ox=0, oy=0, cut=None, expand=True):
        return self.Ana.get_mod_vars(mx, my, ox, oy, self.get_var(), self.Cut(cut), expand)

    def draw_in(self, mx=1, my=1, ox=0, oy=0, nbins=None, cut=None, **dkw):
        x, y, e = self.get_mod_vars(mx, my, ox, oy, cut, expand=True)
        return self.Ana.draw_in(x, y, e * 100, mx * Plane.PX * 1e3, my * Plane.PY * 1e3, nbins, **prep_kw(dkw, title='In Cell Effciency', z_tit='Efficiency [%]'))

    def draw_in_cell(self, nbins=None, ox=0, oy=0, cut=None, **dkw):
        """ in 3D cell"""
        return self.draw_in(self.DUT.GX, self.DUT.GY, ox, oy, nbins, cut, **dkw)

    def draw_in_pixel(self, nbins=None, ox=0, oy=0, cut=None, **dkw):
        """ in pixel of ROC"""
        return self.draw_in(1, 1, ox, oy, nbins, cut, **prep_kw(dkw, title='In Pixel Efficiency'))

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
