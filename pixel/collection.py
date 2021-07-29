#!/usr/bin/env python
# --------------------------------------------------------
#       Analyses of high rate pixel analyses
# created on April 5th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.analysis_collection import AnalysisCollection
from pixel.analysis import PixAnalysis, init_argparser, array, prep_kw


class PixCollection(AnalysisCollection):

    def __init__(self, run_plan, dut_nr, test_campaign=None, load_tree=True, verbose=False):
        AnalysisCollection.__init__(self, run_plan, dut_nr, test_campaign, load_tree, verbose)

    @staticmethod
    def load_dummy():
        return PixAnalysis

    def get_hit_efficiecies(self):
        return array(f[0] for f in self.get_values('hit eff', self.Analysis.draw_hit_efficiency, show=False))

    def draw_hit_efficiencies(self, **kwargs):
        x, y = self.get_fluxes(), self.get_hit_efficiecies()
        self.Draw.graph(x, y, 'Hit Efficiencies', **prep_kw(kwargs, y_tit='Hit Efficiency [%]', **self.get_x_args(draw=True), draw_opt='alp'))

    def generate_threshold_pickle(self):
        pass


if __name__ == '__main__':

    p = init_argparser(run=5, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)
    p.add_argument('-r', '--runs', action='store_true')
    p.add_argument('-d', '--draw', action='store_true')
    p.add_argument('-rd', '--redo', action='store_true')
    pargs = p.parse_args()

    z = PixCollection(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
    z.print_loaded()
    if pargs.runs:
        z.Currents.draw_indep_graphs()
        input('Press any button to exit')
    if pargs.draw:
        z.draw_all(pargs.redo)
