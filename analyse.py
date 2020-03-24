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

# pargs = init_argparser(run_number=23, tc='201908', dut=1, tree=True, has_verbose=True, has_collection=True)
aparser = init_argparser(run=88, tc=None, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)

aparser.add_argument('-d', '--draw', action='store_true', help='make all plots')
aparser.add_argument('-rd', '--redo', action='store_true', help='redo all plots')
pargs = aparser.parse_args()

if not pargs.collection:
    z = analysis_selector(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
else:
    z = collection_selector(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
    if pargs.draw:
        z.draw_all(redo=pargs.redo)
