#!/usr/bin/env python
# --------------------------------------------------------
#       general script to choose the correct analysis for a given run
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from helpers.draw import *  # import everything so that methods are available in ipython


def analysis_selector(run, dut, tc, tree, verbose=False):
    from src.run import Run
    dummy = Run(int(run), tc, tree=False)
    if dummy.get_type() == 'pad':
        from src.pad_analysis import PadAnalysis
        return PadAnalysis(int(run), dut, tc, tree, None, verbose)
    elif dummy.get_type() == 'pixel':
        from src.pix_analysis import PixAnalysis
        return PixAnalysis(int(run), dut, tc, tree, None, verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')


def collection_selector(rp, dut, tc, tree, verbose=False):
    from src.run_selection import RunSelection
    dummy = RunSelection(tc, rp, dut, verbose=False)
    if dummy.get_selected_type() == 'pad':
        from src.pad_collection import PadCollection
        return PadCollection(rp, dut, tc, tree, verbose)
    elif dummy.get_selected_type() == 'pixel':
        from src.pix_collection import PixCollection
        return PixCollection(rp, dut, tc, tree, verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')


if __name__ == "__main__":

    aparser = init_argparser(run=171, tc=None, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)

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
            c = z.Run.Converter
        except AttributeError:
            pass
        if pargs.draw:
            get_attribute(z, pargs.command)(**loads(pargs.kwargs))
    else:
        z = collection_selector(pargs.runplan, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
        if pargs.draw:
            z.draw_all(redo=pargs.redo)
