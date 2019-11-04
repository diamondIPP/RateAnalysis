#!/usr/bin/env python
# --------------------------------------------------------
#       general script to choose the correct analysis for a given run
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from sys import path
from os.path import realpath, dirname, join
path.append(join(dirname(realpath(__file__)), 'src'))
from draw import *

args = init_argparser(run=139, tc='201810', dut=1, tree=True, has_verbose=True, has_collection=True)
# pargs = init_argparser(run_number=23, tc='201908', dut=1, tree=True, has_verbose=True, has_collection=True)

if not args.collection:
    from pad_analysis import PadAnalysis
    from pix_analysis import PixAnalysis
    from run import Run
    dummy = Run(args.runplan, args.testcampaign, tree=False)
    if dummy.get_type() == 'pad':
        z = PadAnalysis(args.runplan, args.dut, args.testcampaign, args.tree, None, args.verbose)
    elif dummy.get_type() == 'pixel':
        z = PixAnalysis(args.runplan, args.dut, args.testcampaign, args.tree, None, args.verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')
else:
    from pad_collection import PadCollection
    from pix_collection import PixCollection
    from run_selection import RunSelection
    dummy = RunSelection(args.testcampaign, args.runplan, args.dut, verbose=False)
    if dummy.get_selected_type() == 'pad':
        z = PadCollection(args.runplan, args.dut, args.testcampaign, args.tree, args.verbose)
    elif dummy.get_selected_type() == 'pixel':
        z = PixCollection(args.runplan, args.dut, args.testcampaign, args.tree, args.verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')

