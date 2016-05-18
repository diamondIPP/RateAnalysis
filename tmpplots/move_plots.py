#!/usr/bin/python

from glob import glob
from shutil import copy
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', default='S129')
args = parser.parse_args()

names = glob('../Results2015*/{dia}/runpl*/*/PulseHeight*Fl*'.format(dia=args.d))
for name in names:
    sub_names = name.split('/')
    rp = '.rp{0}'.format(sub_names[-3][-2:])
    save_name = '{0}/'.format(args.d)
    if sub_names[-1].endswith('.eps'):
        save_name += sub_names[-1]
    elif sub_names[-1].endswith('.png'):
        save_name += sub_names[-1][:-4] + rp + sub_names[-1][-4:]
    else:
        save_name += sub_names[-1][:-5] + rp + sub_names[-1][-5:]
    copy(name, save_name)
