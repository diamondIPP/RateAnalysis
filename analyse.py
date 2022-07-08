#!/usr/bin/env python
# --------------------------------------------------------
#       general script to choose the correct analysis for a given run
# created on Oct 15th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from plotting.draw import *  # import everything so that methods are available in ipython
import plotting.latex as tex  # noqa
from src.binning import Bins  # noqa
from numpy import *
from helpers.utils import *


def analysis_selector(run_, dut_, tc, tree, verbose=False, prnt=True):
    from src.run import Run
    dummy = Run(int(run_), tc, load_tree=False)
    if dummy.get_type() == 'pad':
        from pad.analysis import PadAnalysis
        return PadAnalysis(int(run_), dut_, tc, tree, verbose, prnt)
    elif dummy.get_type() == 'pixel':
        from pixel.analysis import PixAnalysis
        return PixAnalysis(int(run_), dut_, tc, tree, verbose, prnt)
    else:
        critical('wrong run type: has to be in [pad, pixel]')


def collection_selector(tag, dut_, tc, tree, verbose=False):
    from src.run_selection import RunPlan, RunSelection
    dummy = RunPlan(tag, tc, dut_, verbose) if isfloat(tag) else RunSelection(tag, verbose)
    if dummy.DUTType == 'pad':
        from pad.collection import PadCollection
        if 'voltage' in dummy.Type:
            from src.voltage_scan import PadVScan
            return PadVScan(tag, dut_, tc, tree, verbose)
        if 'angle' in dummy.Type:
            from src.angle_scan import PadAScan
            return PadAScan(tag, dut_, tc, tree, verbose)
        return PadCollection(tag, dut_, tc, tree, verbose)
    elif dummy.DUTType == 'pixel':
        from pixel.collection import PixCollection
        if 'voltage' in dummy.Type:
            from src.voltage_scan import PixVScan
            return PixVScan(tag, dut_, tc, tree, verbose)
        if 'angle' in dummy.Type:
            from src.angle_scan import PixAScan
            return PixAScan(tag, dut_, tc, tree, verbose)
        return PixCollection(tag, dut_, tc, tree, verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')


if __name__ == '__main__':

    aparser = init_argparser(run='392', tc=None, dut=1, tree=True, has_verbose=True, has_collection=True, return_parser=True)

    aparser.add_argument('-d', '--draw', action='store_true', help='make all plots')
    aparser.add_argument('-rd', '--redo', action='store_true', help='redo all plots')
    aparser.add_argument('-rm', '--remove', action='store_true', help='remove metadata')
    aparser.add_argument('-rc', '--reconvert', action='store_true', help='remove current file and reconvert the run')
    aparser.add_argument('-cmd', '--command', nargs='?', help='method to be executed')
    aparser.add_argument('-kw', '--kwargs', nargs='?', help='key word arguments as dict {"show": 1}', default='{}')
    pargs = aparser.parse_args()

    from src.analysis import Analysis
    from src.run import Run
    from src.run_selection import Ensemble

    this_tc = Analysis.find_testcampaign(pargs.testcampaign)
    run = Run(None, this_tc, load_tree=False, verbose=False)
    runs = list(run.load_run_info_file().keys())
    e = Ensemble()

    run_plans = list(load_json(e.Dir.joinpath(run.MainConfig.get('SELECTION', f'run plan file')))[run.TCString].keys())
    run_plans += list(load_json(e.Dir.joinpath(run.MainConfig.get('SELECTION', f'run selection file'))))

    if pargs.runplan in runs and not pargs.collection:
        if pargs.reconvert:
            analysis_selector(pargs.runplan, pargs.dut, this_tc, False, False, False).Run.Converter.reconvert()
        if pargs.remove:
            z = analysis_selector(pargs.runplan, pargs.dut, this_tc, False, False, False)
            z.remove_metadata(all_subdirs=True)
            info(f'removed all pickles for {z}', blank_lines=1)
        z = analysis_selector(pargs.runplan, pargs.dut, this_tc, pargs.tree, pargs.verbose)
        try:
            if pargs.tree:
                if z.Run.Type == 'pixel':
                    cal = z.Calibration
                    e = z.Efficiency
                else:
                    p = z.Peaks
                    pul = z.Pulser
                    ped = z.Pedestal
                    w = z.Waveform
                    tm = z.Timing
                tr = z.Tracks
                t = z.Tel
            run = z.Run
            dut = z.DUT
            c = z.Run.Converter
            cut = z.Cut
        except AttributeError:
            pass
        if pargs.draw:
            from json import loads
            get_attribute(z, pargs.command)(**loads(pargs.kwargs))
    elif rp2str(pargs.runplan) in run_plans or pargs.collection:
        if pargs.remove:
            z = collection_selector(pargs.runplan, pargs.dut, this_tc, False, False)
            z.remove_metadata(all_subdirs=True)
            z.remove_tc_metadata()
            info(f'removed all pickles for {z}', blank_lines=1)
        z = collection_selector(pargs.runplan, pargs.dut, this_tc, pargs.tree, pargs.verbose)
        try:
            e = z.Ensemble
            p = z.Pedestal
            pul = z.Pulser
            dut = z.DUT
        except AttributeError:
            pass
        if pargs.draw:
            z.draw_all(redo=pargs.redo)
    else:
        from src.runplan_selection import DiaScans
        z = DiaScans(pargs.runplan, pargs.verbose)
    d = z.Draw
