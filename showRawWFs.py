#!/usr/bin/env python

from ROOT import TFile, TCanvas, TH2F, gStyle, TH1F
from argparse import ArgumentParser
from numpy import mean
from array import array

parser = ArgumentParser()
parser.add_argument('run', nargs=1)
args = parser.parse_args()
filename = args.run[0]
print 'analysing file:', filename

f = TFile(filename)
for key in f.GetListOfKeys():
    print key
tree = f.Get('tree')
n = tree.GetEntries()
print 'the tree has {0} entries'.format(n)
h1 = TH2F('wf', 'Waveform', 1024, 0, 1023, 1000, -500, 500)
h3 = TH2F('wfp', 'Waveform', 1024, 0, 1023, 1000, -500, 500)
h1.SetStats(0)
h3.SetStats(0)
gStyle.SetPalette(55)
tree.Draw('wf0:Iteration$>>wf', '!(pulser)', 'goff')
tree.Draw('wf0:Iteration$>>wfp', 'pulser', 'goff')
c = TCanvas('c', 'WaveForm', 900, 300)
c.SetRightMargin(.045)
h1.GetXaxis().SetTitle('DRS4 Units')
h3.GetXaxis().SetTitle('DRS4 Units')
h1.Draw('col')
c2 = TCanvas('c2', 'WaveForm', 900, 300)
c2.SetRightMargin(.045)
h3.Draw('col')
signal_range = [int(raw_input('enter min signal range: ')), int(raw_input('enter max signal range: '))]
mean_sig = 0
for i in xrange(10):
    tree.GetEntry(n/2 + i)
    mean_sig += mean(tree.wf0)
print 'mean signal is:', mean_sig / 10
pol = 1 if mean_sig > 0 else -1
h2 = TH1F('pv', 'Peak Values', 40, signal_range[0], signal_range[1])
for i in xrange(n/2, n):
    tree.GetEntry(i)
    lst = list(tree.wf0)
    vec = lst[signal_range[0]:signal_range[1]]
    extrema = lst.index(min(vec)) if mean_sig < 0 else lst.index(max(vec))
    h2.Fill(extrema)
c1 = TCanvas('c1', 'WaveForm', 900, 300)
c1.SetRightMargin(.045)
h2.Draw()
c.SaveAs('RawWaveforms.png')
c1.SaveAs('PeakValues.png')
raw_input('Press any key to end: ')