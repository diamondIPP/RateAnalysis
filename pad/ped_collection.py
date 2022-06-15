# --------------------------------------------------------
#       Pulser sub class for AnalysisCollection
# created by Michael Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.sub_ana_collection import SubCollection
from plotting.draw import *
import plotting.latex as tex
from pad.pedestal import PedestalAnalysis


class PedCollection(SubCollection):

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

    def noise(self, redo=False):
        return self.get_values('noise', PedestalAnalysis.get_noise, picklepath=self.get_pickle_path(), redo=redo)

    def mean(self, redo=False):
        return self.get_values('mean', PedestalAnalysis.get_mean, picklepath=self.get_pickle_path(), redo=redo)
    # endregion GET
    # ----------------------------------------

    def draw_dists(self, **dkw):
        h = self.get_plots('pedestal distributions', PedestalAnalysis.get_dist, picklepath=self.get_pickle_path('Dist'), prnt=False)
        return self.Draw.stack(h, 'PedDists', self.flux_strings(), **prep_kw(dkw, scale=True))

    def print(self, prec=2):
        header = tex.bold(f'Flux {tex.unit("khzcm")}', f'Pedestal {tex.unit("mV")}', f'Noise {tex.unit("mV")}')
        m, s, f = self.mean(), self.noise(), uarr2n(self.get_fluxes())
        print(tex.table(header, rows=[tex.num(f[i], fmt='.0f') + tex.num(m[i], s[i], fmt=f'.{prec}f') for i in range(m.size)]))
