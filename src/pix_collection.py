#!/usr/bin/env python
# --------------------------------------------------------
#       Analyses of high rate pixel analyses
# created on April 5th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis_collection import AnalysisCollection
from pix_analysis import PixAnalysis, format_histo, init_argparser


class PixCollection(AnalysisCollection):

    def __init__(self, run_plan, dut_nr, test_campaign=None, load_tree=True, verbose=False):
        AnalysisCollection.__init__(self, run_plan, dut_nr, test_campaign, load_tree, verbose)

    def load_analysis(self, run_number):
        return PixAnalysis(run_number, self.DUTNumber, self.TCString, self.Threads[run_number].Tuple, self.Threads[run_number].Time, self.Verbose, prnt=False)

    @staticmethod
    def load_dummy():
        return PixAnalysis

    def draw_hit_efficiencies(self, show=True):
        fits = [ana.draw_hit_efficiency(show=False) for ana in self.Analyses.itervalues()]
        y, ey = [fit.Parameter(0) for fit in fits], [fit.ParError(0) for fit in fits]
        x, ex = self.get_fluxes().values(), [flux * .1 for flux in self.get_fluxes().itervalues()]
        g = self.make_tgrapherrors('ghes', 'Hit Efficiencies', y=y, x=x, ex=ex, ey=ey)
        format_histo(g, x_tit='Flux [kHz/cm^{2}]', y_tit='Hit Efficiency [%]', y_off=1.5)
        self.save_histo(g, 'HitEfficiencies', lm=.12, draw_opt='alp', show=show, logx=True, bm=.17)

    def generate_threshold_pickle(self):
        pass


if __name__ == '__main__':

    p = init_argparser(run=5, dia=1, tree=True, verbose=True, collection=True, return_parser=True)
    p.add_argument('-r', '--runs', action='store_true')
    p.add_argument('-d', '--draw', action='store_true')
    p.add_argument('-rd', '--redo', action='store_true')
    pargs = p.parse_args()

    z = PixCollection(pargs.runplan, pargs.dia, pargs.testcampaign, pargs.tree, pargs.verbose)
    z.print_loaded()
    if pargs.runs:
        z.Currents.draw_indep_graphs()
        raw_input('Press any button to exit')
    if pargs.draw:
        z.draw_all(pargs.redo)
