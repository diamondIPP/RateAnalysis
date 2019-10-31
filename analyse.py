#!/usr/bin/env python
# --------------------------------------------------------
#       general script to choose the correct analysis for a given run
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from sys import path
from os.path import realpath, dirname, join
path.append(join(dirname(realpath(__file__)), 'src'))
from utils import *
from run import Run

args = init_argparser(run=139, tc='201810', dia=1, tree=True, verbose=True, collection=True)

if not args.collection:
    from pad_analysis import PadAnalysis
    from pix_analysis import PixAnalysis
    dummy = Run(args.runplan, args.testcampaign, tree=False)
    if dummy.get_type() == 'pad':
        z = PadAnalysis(args.runplan, args.dia, args.testcampaign, args.tree, None, args.verbose)
    elif dummy.get_type() == 'pixel':
        z = PixAnalysis(args.runplan, args.dia, args.testcampaign, args.tree, None, args.verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')
else:
    from analysis_collection import AnalysisCollection
    z = AnalysisCollection(args.runplan, args.dia, args.testcampaign, args.tree, args.verbose)
