#!/usr/bin/python

import os
from glob import glob
from shutil import copy
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', default='S129')
args = parser.parse_args()
dias = ['S129', 'poly-B2', 'poly-D']

for dia in dias:
    names = glob('../Results2015*/{dia}/runpl*/*/SignalDist*'.format(dia=dia))
    for name in names:
        sub_names = name.split('/')
        file_name = sub_names[-1]
        if '.' in file_name:
            extension = file_name.split('.')[-1]
            file_name = file_name.split('.')[0]
        elif file_name.endswith('root'):
            extension = 'root'
            file_name = file_name[:-4]
        else:
            extension = sub_names[-1][-3:]
            file_name = file_name[:-3]
        print name,file_name,extension
        rp = 'rp{0}'.format(sub_names[-3][-2:])
        save_name = '{dia}/{ext}/'.format(dia=dia,ext=extension)
        folder = save_name
        save_name+='{dia}_'.format(dia=dia.replace('-',''))
        try:
            os.makedirs(folder)
        except:
            pass
        if sub_names[-3][-2:] not in ['23','15','51','11','05','03','91']:
            continue
        tc = '_2015' + sub_names[-5][-2:] + '_'
        if sub_names[-1][-3:] in ['pdf', 'eps']:
            #print  sub_names[-1][-8:]
            save_name += sub_names[-1][:-8] + tc + sub_names[-1][-8:]
        elif sub_names[-1][-3:] in ['png']:
            #print  sub_names[-1][-4:]
            save_name += sub_names[-1][:-4] + tc + rp + sub_names[-1][-4:]
        else:
            #print sub_names[-1][-5:]
            save_name += sub_names[-1][:-4] + tc + rp + '.'+sub_names[-1][-4:]
        print save_name
        copy(name, save_name)
