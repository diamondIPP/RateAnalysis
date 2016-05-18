#!/usr/bin/python
from ROOT import TFile
from sys import argv

print argv
f = TFile(argv[1])
t = f.Get('tree')
