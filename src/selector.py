# --------------------------------------------------------
#       methods to select which run or collection to take
# created on Nov 7th 2019 by M. Reichmann
# --------------------------------------------------------
from utils import critical


def run_selector(run, tc, tree, t_vec, verbose):
    from pixel_run import PixelRun
    from pad_run import PadRun
    from run import Run
    dummy = Run(run, tc, tree=False)
    if dummy.get_type() == 'pad':
        return PadRun(run, tc, tree, t_vec, verbose)
    elif dummy.get_type() == 'pixel':
        return PixelRun(run, tc, tree, t_vec, verbose)
    critical('wrong run type: has to be in [pad, pixel]')


def analysis_selector(run, dut, tc, tree, verbose=False):
    from run import Run
    from pad_analysis import PadAnalysis
    from pix_analysis import PixAnalysis
    dummy = Run(int(run), tc, tree=False)
    if dummy.get_type() == 'pad':
        return PadAnalysis(int(run), dut, tc, tree, None, verbose)
    elif dummy.get_type() == 'pixel':
        return PixAnalysis(int(run), dut, tc, tree, None, verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')


def collection_selector(rp, dut, tc, tree, verbose=False):
    from pad_collection import PadCollection
    from pix_collection import PixCollection
    from run_selection import RunSelection
    dummy = RunSelection(tc, rp, dut, verbose=False)
    if dummy.get_selected_type() == 'pad':
        return PadCollection(rp, dut, tc, tree, verbose)
    elif dummy.get_selected_type() == 'pixel':
        return PixCollection(rp, dut, tc, tree, verbose)
    else:
        critical('wrong run type: has to be in [pad, pixel]')
