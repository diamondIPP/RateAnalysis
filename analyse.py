#!/usr/bin/env python
# --------------------------------------------------------
#       general script to choose the correct analysis for a given run
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from sys import path
from os.path import realpath, dirname, join
path.append(join(dirname(realpath(__file__)), 'src'))
from draw import *
from selector import analysis_selector, collection_selector

args = init_argparser(run=139, tc='201810', dut=1, tree=True, has_verbose=True, has_collection=True)
# pargs = init_argparser(run_number=23, tc='201908', dut=1, tree=True, has_verbose=True, has_collection=True)

if not args.collection:
    z = analysis_selector(args.runplan, args.dut, args.testcampaign, args.tree, args.verbose)
else:
    z = collection_selector(args.runplan, args.dut, args.testcampaign, args.tree, args.verbose)
