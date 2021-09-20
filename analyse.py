#!/usr/bin/env python
# --------------------------------------------------------
#       general script to choose the correct analysis for a given run
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from helpers.draw import *  # import everything so that methods are available in ipython
from src.binning import Bins
_ = Bins.Verbose


def analysis_selector(run, dut, tc, tree, verbose=False, prnt=True):
    from src.run import Run
    dummy = Run(int(run), tc, load_tree=False)
    if dummy.get_type() == 'pad':
        from pad.analysis import PadAnalysis
        return PadAnalysis(int(run), dut, tc, tree, verbose, prnt)
    elif dummy.get_type() == 'pixel':
        from pixel.analysis import PixAnalysis
        return PixAnalysis(int(run), dut, tc, tree, verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')


def collection_selector(name, dut, tc, tree, verbose=False):
    from src.run_selection import RunPlan, RunSelection
    dummy = RunPlan(name, tc, dut, verbose) if isfloat(name) else RunSelection(name, verbose)
    if dummy.DUTType == 'pad':
        from pad.collection import PadCollection
        if 'voltage' in dummy.Type:
            from src.voltage_scan import make_volage_scan
            return make_volage_scan(PadCollection)(name, dut, tc, tree, verbose)
        return PadCollection(name, dut, tc, tree, verbose)
    elif dummy.DUTType == 'pixel':
        from pixel.collection import PixCollection
        if 'voltage' in dummy.Type:
            from src.voltage_scan import make_volage_scan
            return make_volage_scan(PixCollection)(name, dut, tc, tree, verbose)
        return PixCollection(name, dut, tc, tree, verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')


if __name__ == '__main__':

    aparser = init_argparser(run=392, tc=None, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)

    aparser.add_argument('-d', '--draw', action='store_true', help='make all plots')
    aparser.add_argument('-rd', '--redo', action='store_true', help='redo all plots')
    aparser.add_argument('-rm', '--remove', action='store_true', help='remove metadata')
    aparser.add_argument('-rc', '--reconvert', action='store_true', help='remove current file and reconvert the run')
    aparser.add_argument('-cmd', '--command', nargs='?', help='method to be executed')
    aparser.add_argument('-kw', '--kwargs', nargs='?', help='key word arguments as dict {"show": 1}', default='{}')
    pargs = aparser.parse_args()
    pargs.tree = False if pargs.remove or pargs.reconvert else pargs.tree  # don't load tree if only metadata should be removed

    from src.analysis import Analysis
    this_tc = Analysis.find_testcampaign(pargs.testcampaign)

    if not pargs.collection and isint(pargs.runplan):
        z = analysis_selector(pargs.runplan, pargs.dut, this_tc, pargs.tree, pargs.verbose)
        try:
            if pargs.tree:
                if z.Run.Type == 'pixel':
                    cal = z.Calibration
                else:
                    p = z.Peaks
                    pul = z.Pulser
                    ped = z.Pedestal
                    w = z.Waveform
                tr = z.Tracks
                t = z.Tel
            c = z.Run.Converter
            cut = z.Cut
        except AttributeError:
            pass
        if pargs.draw:
            get_attribute(z, pargs.command)(**loads(pargs.kwargs))
    elif isfloat(pargs.runplan) or pargs.collection:
        z = collection_selector(pargs.runplan, pargs.dut, this_tc, pargs.tree, pargs.verbose)
        if pargs.draw:
            z.draw_all(redo=pargs.redo)
    else:
        from src.runplan_selection import DiaScans
        z = DiaScans(pargs.runplan, pargs.verbose)

    if pargs.remove:
        z.remove_metadata(all_subdirs=True)
        critical(f'Removed all pickles for {z}')

    if pargs.reconvert:
        z.Run.Converter.reconvert()
        z = analysis_selector(pargs.runplan, pargs.dut, this_tc, True, pargs.verbose)
