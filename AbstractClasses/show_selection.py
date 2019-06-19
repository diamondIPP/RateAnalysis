#!/usr/bin/env python
# --------------------------------------------------------
#       small script to show the information about the run plans of a selection
# created on June 19th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from DiamondRateScans import DiaScans
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('selection', nargs='?', help='Name of the selection of runplans', type=str, default='test')
args = parser.parse_args()

z = DiaScans(args.selection)
z.show_selection()