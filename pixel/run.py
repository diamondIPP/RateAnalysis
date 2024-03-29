#!/usr/bin/env python
# --------------------------------------------------------
#       pixel run class
# created on Oct 21st 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.run import Run
from src.dut import PixelDUT
from helpers.utils import init_argparser


class PixelRun(Run):
    """ Run class child with pixel extension. """

    def __init__(self, number=None, testcampaign=None, load_tree=True, verbose=False):
        """
        :param number: if None is provided it creates a dummy run
        :param testcampaign: if None is provided ...
        :param load_tree: root_tree object, if None is given it will start the converter
        :param verbose: turn on more output
        """

        Run.__init__(self, number, testcampaign, load_tree, verbose)
        # Settings
        self.Type = 'pixel'

    @property
    def dut(self):
        return PixelDUT


if __name__ == '__main__':

    # e.g.: (138, 201510), (489, 201610), (293, 201508)
    args = init_argparser(489, tc='201610', tree=True)
    z = PixelRun(args.run, testcampaign=args.testcampaign, load_tree=args.tree)
