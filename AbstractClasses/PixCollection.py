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
from Utils import get_resolution
from PixAnalysis import PixAnalysis


class PixCollection(AnalysisCollection):

    def __init__(self, run_selection, load_tree=True, verbose=False):
        AnalysisCollection.__init__(self, run_selection, load_tree=load_tree, verbose=verbose)

    def add_analyses(self, load_tree):
        """ Creates and adds Analysis objects with run numbers in runs. """
        for run in self.runs:
            analysis = PixAnalysis(run, self.selection.SelectedDiamondNr, load_tree=load_tree, verbose=self.verbose)
            self.collection[analysis.run.RunNumber] = analysis
            self.current_run_number = analysis.run.RunNumber

    def generate_slope_pickle(self):
        pass

    def generate_threshold_pickle(self):
        pass


if __name__ == "__main__":
    st = time()
    main_parser = ArgumentParser()
    main_parser.add_argument('runplan', nargs='?', default=16.1)
    main_parser.add_argument('dia', nargs='?', default=1, type=int)
    main_parser.add_argument('-tc', '--testcampaign', nargs='?', default='')
    main_parser.add_argument('-t', '--tree', default=True, action='store_false')
    main_parser.add_argument('-v', '--verbose', default=True, action='store_false')
    args = main_parser.parse_args()
    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    a = Elementary(tc, True, get_resolution())
    sel = RunSelection(testcampaign=tc)
    sel.select_runs_from_runplan(args.runplan, ch=args.dia)
    a.print_banner('STARTING PIX-ANALYSIS COLLECTION FOR RUNPLAN {0}'.format(args.runplan))
    a.print_testcampaign()

    z = PixCollection(sel, load_tree=args.tree, verbose=args.verbose)
    z.print_loaded()
    z.print_elapsed_time(st, 'Instantiation')
