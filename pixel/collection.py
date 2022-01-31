#!/usr/bin/env python
# --------------------------------------------------------
#       Analyses of high rate pixel analyses
# created on April 5th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.analysis_collection import AnalysisCollection, fname
from pixel.analysis import PixAnalysis, init_argparser, array, prep_kw, ufloat, sqrt, quiet


class PixCollection(AnalysisCollection):

    PhTit = 'Pulse Height [vcal]'

    def __init__(self, run_plan, dut_nr, test_campaign=None, load_tree=True, verbose=False):
        AnalysisCollection.__init__(self, run_plan, dut_nr, test_campaign, load_tree, verbose)

        self.N = self.FirstAnalysis.N

    @staticmethod
    def load_dummy():
        return PixAnalysis

    @quiet
    def save_coll_plots(self):
        self.draw_efficiencies(show=False)
        self.draw_cluster_sizes(show=False)
        super(PixCollection, self).save_coll_plots()

    # ----------------------------------------
    # region GET
    def get_pulse_heights(self, bin_width=None, redo=False, runs=None, avrg=False, flux_sort=False, err=False, **kw):
        picklepath = self.get_pickle_path('Fit', sub_dir='PH')
        x = self.get_values('pulse heights', self.Analysis.get_pulse_height, runs, avrg=avrg, picklepath=picklepath, bin_size=bin_width, redo=redo, flux_sort=flux_sort)
        return array([ufloat(ph.n, sqrt(ph.s ** 2 + ((self.get_sys_error() if err else 0) * ph.n) ** 2)) for ph in x])

    def get_efficiencies(self, suf=None, avrg=False, redo=False):
        return super(PixCollection, self).get_efficiencies('0', avrg, redo)

    def get_cluster_sizes(self, avrg=False, redo=False):
        return self.get_values('cluster sizes', self.Analysis.get_cluster_size, None, 1, avrg, self.get_pickle_path('CS', '', 'Telescope', self.N), redo=redo)
    # endregion GET
    # ----------------------------------------

    def draw_efficiencies(self, avrg=False, t=False, redo=False, **dkw):
        x, y = self.get_x_var(t, avrg=avrg), self.get_efficiencies(avrg=avrg, redo=redo)
        return self.Draw.graph(x, y, 'Hit Efficiencies', **prep_kw(dkw, y_tit='Hit Efficiency [%]', **self.get_x_args(t, draw=True), draw_opt='alp', file_name=fname('Efficiencies', avrg, t)))

    def draw_cluster_sizes(self, avrg=False, t=False, redo=False, **dkw):
        x, y = self.get_x_var(t, avrg=avrg), self.get_cluster_sizes(avrg, redo)
        return self.Draw.graph(x, y[:, 0], 'Cluster Size', **prep_kw(dkw, y_tit='Cluster Size', **self.get_x_args(t, draw=True), draw_opt='alp', file_name=f'ClusterSize{"Avr" if avrg else ""}'))


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
