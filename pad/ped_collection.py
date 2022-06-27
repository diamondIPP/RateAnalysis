# --------------------------------------------------------
#       Pulser sub class for AnalysisCollection
# created by Michael Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.sub_ana_collection import SubCollection, fname
from plotting.draw import *
import plotting.latex as tex
from pad.pedestal import PedestalAnalysis


class PedCollection(SubCollection):

    PhTit = 'Pedestal [mV]'
    NoiseTit = 'Noise [mV]'
    ErrSys = SubCollection.MainConfig.get_value('MAIN', 'systematic ped error', default=0.)

    def __init__(self, ana_collection):
        super().__init__(ana_collection)
        if self.Ana.FirstAnalysis.Tree.Hash():
            self.Analysis = self.Ana.FirstAnalysis.Pulser

    # ----------------------------------------
    # region GET
    def get_analyses(self, runs=None):
        return [ana.Pedestal for ana in self.Analyses if ana.Run.Number in choose(runs, self.Ana.Runs)]

    def get_pickle_path(self, name='', suf='', sub_dir=None, dut=None, camp=None):
        return super().get_pickle_path(name, suf, 'Pedestal', dut, camp)

    def noise(self, avrg=False, flux_sort=False, runs=None, redo=False):
        return self.get_values('noise', PedestalAnalysis.get_noise, runs, None, avrg, self.get_pickle_path(), flux_sort, redo=redo)

    def mean(self, avrg=False, flux_sort=False, runs=None, err=True, redo=False):
        x = self.get_values('mean', PedestalAnalysis.get_mean, runs, None, False, self.get_pickle_path(), flux_sort, redo=redo)
        x = add_err(x, PedCollection.ErrSys if err else 0)  # add sys error to all runs
        return self.Ana.get_flux_average(x) if avrg else x
    # endregion GET
    # ----------------------------------------

    def draw(self, t=False, avrg=False, redo=False, **dkw):
        x, y = self.get_x(t, avrg=avrg), self.mean(avrg, redo=redo)
        self.Draw.graph(x, y, **prep_kw(dkw, title='Ped', y_tit=self.PhTit, **self.get_x_args(t), color=810, file_name=fname('Ped', avrg, t)))

    def draw_noise(self, t=False, avrg=False, redo=False, **dkw):
        x, y = self.get_x(t, avrg=avrg), self.noise(avrg, redo=redo)
        self.Draw.graph(x, y, **prep_kw(dkw, title='Noise', y_tit=self.NoiseTit, **self.get_x_args(t), color=810, file_name=fname('Noise', avrg, t)))

    def draw_dists(self, **dkw):
        h = self.get_plots('pedestal distributions', PedestalAnalysis.get_dist, picklepath=self.get_pickle_path('Dist'), prnt=False)
        return self.Draw.stack(h, 'PedDists', self.flux_strings(), **prep_kw(dkw, scale=True))

    def draw_trends(self, **dkw):
        g = self.get_plots('pedestal trends', PedestalAnalysis.get_trend, picklepath=self.get_pickle_path('Trend'))
        g = [shift_graph(ig, ox=-get_graph_x(ig)[0].n) for ig in g]
        return self.Draw.multigraph(g, 'PedTrends', self.flux_strings(), **prep_kw(dkw, gridy=True, **self.get_x_args(vs_time=True, off=-1, draw=True), file_name='PedTrends'))

    def print(self, prec=2):
        header = tex.bold(f'Flux {tex.unit("khzcm")}', f'Pedestal {tex.unit("mV")}', f'Noise {tex.unit("mV")}')
        m, s, f = self.mean(), self.noise(), uarr2n(self.get_fluxes())
        print(tex.table(header, rows=[tex.num(f[i], fmt='.0f') + tex.num(m[i], s[i], fmt=f'.{prec}f') for i in range(m.size)]))
