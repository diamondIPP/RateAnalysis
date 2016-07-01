#!/usr/bin/python

from os import system
runs = {'201508': {106: [3], 145: [0], 462: [3]},
        '201510': {171: [0], 216: [0], 392: [0, 3], 434: [0, 3]}}

for tc, vals in sorted(runs.iteritems()):
    for run, chs in sorted(vals.iteritems()):
        for ch in chs:
            print tc, run, ch
            cmd = 'python readtree.py ~/sdvlp/eudaq-drs4/conf/test{tc}00{run}.root {ch}'.format(tc=tc[2:], run=run, ch=ch)
            system(cmd)
