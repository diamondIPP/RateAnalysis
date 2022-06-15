# --------------------------------------------------------
#       Pulser sub class for AnalysisCollection
# created by Michael Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.sub_ana_collection import SubCollection
from plotting.draw import *
from helpers.utils import flux2str
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
    # endregion GET
    # ----------------------------------------

    def draw_dists(self, **dkw):
        h = self.get_plots('pedestal distributions', PedestalAnalysis.get_dist, picklepath=self.get_pickle_path('Dist'), prnt=False)
        return self.Draw.stack(h, 'PedDists', flux2str(uarr2n(self.get_fluxes())), **prep_kw(dkw, scale=True))
