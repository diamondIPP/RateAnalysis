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
from json import loads


aparser = init_argparser(run=178, tc=None, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)

aparser.add_argument('-d', '--draw', action='store_true', help='make all plots')
aparser.add_argument('-rd', '--redo', action='store_true', help='redo all plots')
aparser.add_argument('-cmd', '--command', nargs='?', help='method to be executed')
aparser.add_argument('-kw', '--kwargs', nargs='?', help='key word arguments as dict {"show": 1}', default='{}')
pargs = aparser.parse_args()

if not pargs.collection:
    z = analysis_selector(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
    try:
        p = z.Peaks if pargs.tree else None
        w = z.Waveform if pargs.tree else None
        t = z.Timing if pargs.tree else None
    except AttributeError:
        pass
    if pargs.draw:
        get_attribute(z, pargs.command)(**loads(pargs.kwargs))
else:
    z = collection_selector(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
    if pargs.draw:
        z.draw_all(redo=pargs.redo)
