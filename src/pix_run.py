#!/usr/bin/env python
# --------------------------------------------------------
#       pixel run class
# created on Oct 21st 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from run import Run
from helpers.utils import init_argparser


class PixelRun(Run):
    """ Run class containing all the information for a single run. """

    def __init__(self, number=None, testcampaign=None, load_tree=True, verbose=False):
        """
        :param number: if None is provided it creates a dummy run
        :param testcampaign: if None is provided ...
        :param load_tree: root_tree object, if None is given it will start the converter
        :param verbose: turn on more output
        """

        # Settings
        self.Type = 'pixel'

        Run.__init__(self, number, testcampaign, load_tree, verbose)


if __name__ == '__main__':

    args = init_argparser(489, tc='201610', tree=True)
    z = PixelRun(args.run, testcampaign=args.testcampaign, load_tree=args.tree)
