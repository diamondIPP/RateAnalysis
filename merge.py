#!/usr/bin/env python
# --------------------------------------------------------
#       small script to merge up to three ROOT trees
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TFile
from argparse import ArgumentParser
from glob import glob
from array import array
from sys import stdout

def get_root_files(arg):
    if arg is None:
        return None
    try:
        return TFile(glob('*{r}*'.format(r=str(arg).zfill(3)))[0])
    except IndexError:
        raise ValueError('No root file found for run {r}'.format(r=arg))

parser = ArgumentParser()
parser.add_argument('r1', type=int)
parser.add_argument('r2', type=int)
parser.add_argument('r3', nargs='?', default=None, type=int)

args = parser.parse_args()
root_files = [get_root_files(i) for i in [args.r1, args.r2, args.r3] if get_root_files(i) is not None]

# NewFile = TFile(root_files[0].GetName(), 'RECREATE')
new_file = TFile('test.root', 'RECREATE')

trees = [f.Get('tree') for f in root_files]
print 'cloning tree...'
new_tree = trees[0].CloneTree()
new_tree.ResetBranchAddresses()

event_number = array('i', [0])

at_entry = trees[0].GetEntries()

for tree in trees[1:]:
    new_tree.SetBranchAddress('event_number', event_number)
    new_tree.CopyAddresses(tree)
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        event_number[0] = at_entry
        at_entry += 1
        print '\r', event_number[0],
        stdout.flush()
        new_tree.Fill()
    tree.ResetBranchAddresses()
print

new_file.cd()
new_tree.Write()
macro = root_files[0].Get('region_information')
if macro:
    macro.Write()
new_file.Write()
