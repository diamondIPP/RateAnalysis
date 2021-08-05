#!/usr/bin/env python
# --------------------------------------------------------
#       small script to merge up to three ROOT trees
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from argparse import ArgumentParser
from src.run import Run, loads, join, info, PBar, TFile, zeros
from src.analysis import Analysis

parser = ArgumentParser()
parser.add_argument('runs', help='list of runs: [112,113]')
parser.add_argument('-tc', nargs='?', default=None, help='Beam Test')
args = parser.parse_args()

tc = Analysis.find_testcampaign(args.tc)
runs = [Run(run_nr, tc, verbose=False) for run_nr in loads(args.runs)]
trees = [run.Tree for run in runs]
run = runs[0]


# NewFile = TFile(root_files[0].GetName(), 'RECREATE')
new_file = TFile(join(run.RootFileDir, f'merged-{"-".join(str(r.Number) for r in runs)}.root'), 'RECREATE')

info('cloning tree...')
new_tree = trees[0].CloneTree()
new_tree.ResetBranchAddresses()
event_number = zeros(1, dtype='i')

at_entry = run.NEvents

pbar = PBar(sum(r.NEvents for r in runs[1:]))
for tree in trees[1:]:
    new_tree.SetBranchAddress('event_number', event_number)
    new_tree.CopyAddresses(tree)
    for ev in tree:
        event_number[0] = at_entry
        at_entry += 1
        new_tree.Fill()
        pbar.update()
    tree.ResetBranchAddresses()

new_file.cd()
new_tree.Write()
macro = run.RootFile.Get('region_information')
if macro:
    macro.Write()
new_file.Write()
