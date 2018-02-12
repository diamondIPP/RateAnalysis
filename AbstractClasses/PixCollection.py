#!/usr/bin/env python
# --------------------------------------------------------
#       Collection of high rate pixel analyses
# created on April 5th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from AnalysisCollection import AnalysisCollection
from Elementary import Elementary
from RunSelection import RunSelection
from time import time
from argparse import ArgumentParser
from Utils import get_resolution, load_root_files
from PixAnalysis import PixAnalysis
from Run import Run


class PixCollection(AnalysisCollection):

    def __init__(self, run_selection, threads=None, verbose=False):
        AnalysisCollection.__init__(self, run_selection, threads=threads, verbose=verbose)

    def add_analyses(self):
        """ Creates and adds Analysis objects with run numbers in runs. """
        for run in self.Runs:
            run_class = Run(run, tree=self.Threads[run].Tuple, verbose=self.verbose)
            analysis = PixAnalysis(run_class, self.selection.SelectedDiamondNr)
            self.collection[analysis.run.RunNumber] = analysis
            self.current_run_number = analysis.run.RunNumber

    def generate_slope_pickle(self):
        pass

    def generate_threshold_pickle(self):
        pass


if __name__ == "__main__":
    st = time()
    p = ArgumentParser()
    p.add_argument('runplan', nargs='?', default=16.1)
    p.add_argument('dia', nargs='?', default=1, type=int)
    p.add_argument('-tc', '--testcampaign', nargs='?', default='')
    p.add_argument('-t', '--tree', action='store_false')
    p.add_argument('-v', '--verbose', action='store_false')
    args = p.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    a = Elementary(tc, True, get_resolution())
    sel = RunSelection(testcampaign=tc)
    sel.select_runs_from_runplan(args.runplan, ch=args.dia)
    a.print_banner('STARTING PIX-ANALYSIS COLLECTION FOR RUNPLAN {0}'.format(args.runplan))
    a.print_testcampaign()

    t = load_root_files(sel, args.tree)
    z = PixCollection(sel, threads=t, verbose=args.verbose)
    z.print_loaded()
    z.print_elapsed_time(st, 'Instantiation')
