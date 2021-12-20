# --------------------------------------------------------
#       sub-analysis class for efficiency of the pixel detectors
# created on December 17th 2021 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.sub_analysis import SubAnalysis
from helpers.utils import save_pickle, deepcopy
from src.binning import Bins, Plane, make_bins
from plotting.draw import FitRes, find_bins, choose, prep_kw, calc_eff
from src.cut import Cut


class Efficiency(SubAnalysis):

    def __init__(self, pix_analysis):
        super().__init__(pix_analysis, pickle_dir='Efficiency')
        self.N = self.Ana.N
        self.Cut = deepcopy(self.Ana.Cut).remove('rhit', ret=True)

    def get_values(self, cut=None, rhit=True):
        return self.get_tree_vec(self.get_var(rhit=rhit), self.Cut(cut), dtype='?')

    def get_var(self, plane=None, rhit=True, sigma=None):
        var = f'n_clusters[{choose(plane, self.N)}] > 0'
        return Cut.to_string(self.Cut.generate_rhit(sigma) + Cut.make('', var)) if rhit else var

    @save_pickle(sub_dir='Efficiency', suf_args='all')
    def get(self, cut=None, rhit=True, _redo=False):
        return calc_eff(values=self.get_values(cut, rhit))

    def draw_eff_vs_chi2(self, **kwargs):
        x, e = self.get_tree_vec(['chi2_tracks', self.get_var()], self.Cut())
        self.Draw.efficiency(x, e, find_bins(x, lfac=0, lq=0), title='Efficiency vs Chi2', **prep_kw(kwargs, x_tit='#chi^{2}'))

    def draw_hit_efficiency(self, cut=None, bin_size=None, rel_time=False, **kwargs):
        (x, e), bins = self.get_tree_vec([self.get_t_var(), self.get_var()], choose(cut, self.Cut())), self.Bins.get_time(bin_size, cut)
        g = self.Draw.efficiency(x, e, bins, **prep_kw(kwargs, title='Hit Efficiency', **self.get_t_args(rel_time), y_range=[-5, 115], gridy=True, draw_opt='apz'))
        fit = FitRes(g.Fit('pol0', 'qs'))
        self.Draw.stats(fit, width=.35, y2=.35, names=['Efficiency'])
        self.Draw.preliminary()
        self.Draw.save_plots('HitEfficiency', **kwargs)
        return fit if fit.Parameter(0) is not None else 0

    def draw_map(self, res=None, fid=False, cut=None, **kwargs):
        x, y, zz = self.get_tree_vec(self.Ana.get_track_vars() + [self.get_var()], self.Cut(cut) if cut or fid else self.Cut.exclude('fiducial'))
        tit, (xtit, ytit), ztit = 'Efficiency Map', [f'Track Position {i} [mm]' for i in ['X', 'Y']], 'Efficiency [%]'
        self.Draw.prof2d(x, y, zz * 100, Bins.get_global(res), tit, **prep_kw(kwargs, x_tit=xtit, y_tit=ytit, z_tit=ztit))
        self.Draw.preliminary()
        self.Ana.draw_fid_cut()
        self.Draw.save_plots('Efficiency Map')

    def get_fiducial_cell(self, n):
        x1, x2, y1, y2 = self.Cut.CutConfig['fiducial']
        nx = int(round((x2 - x1) / Plane.PX))
        return round(x1 + Plane.PX * (n % nx), 4), round(y1 + Plane.PY * (n / nx), 4)

    def draw_cell_efficiency(self, nbins=None, **dkw):
        x, y, e = self.get_tree_vec(self.Ana.get_track_vars() + [self.get_var()], self.Cut())
        bins = None if nbins is None else [nbins, 0, Plane.PX, nbins, 0, Plane.PY]
        self.Draw.prof2d(x % Plane.PX, y % Plane.PY, e * 100, bins, 'Cell Efficiency', **prep_kw(dkw, x_tit='Track X [mm]', y_tit='Track Y [mm]', z_tit='Efficiency [%]'))

    def draw_efficiency_vs_trigphase(self, **kwargs):
        x, e = self.get_tree_vec(['trigger_phase[1]', self.get_var()], self.Cut.exclude('trigger_phase'))
        return self.Draw.efficiency(x, e, make_bins(-.5, 10), **prep_kw(kwargs, title='Trigger Phase Efficiency', x_tit='Trigger Phase', x_range=[-1, 10], draw_opt='bap'))
