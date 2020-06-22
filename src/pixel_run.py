#!/usr/bin/env python
# --------------------------------------------------------
#       pad run class
# created on Oct 21st 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from run import Run, join
from utils import ensure_dir, init_argparser


class PixelRun(Run):
    """ Run class containing all the information for a single run. """

    def __init__(self, number=None, test_campaign=None, tree=True, t_vec=None, verbose=False):
        """
        :param number: if None is provided it creates a dummy run
        :param test_campaign: if None is provided ...
        :param tree: root_tree object, if None is given it will start the converter
        :param t_vec: time sequence of the run, if None is provide it will generate a corrected one
        :param verbose: turn on more output
        """

        # Settings
        self.Type = 'pixel'

        Run.__init__(self, number, test_campaign, tree, t_vec, verbose)

    def load_rootfile_dirname(self):
        return ensure_dir(join(self.TCDir, 'root', self.Type))


if __name__ == '__main__':

    args = init_argparser(23, tc='201809', tree=True)
    z = PixelRun(args.run, tree=args.tree, test_campaign=args.testcampaign)
