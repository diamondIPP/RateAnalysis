from ROOT import TCanvas, TF2, Double, TMath, gRandom, TFile, TTree
import numpy as np
from array import array
from datetime import datetime

def LateralShape(x,par):
    result = 0.1
    norm = 1.
    sigma = 0.04
    for i in xrange(len(par)/2):
        result += norm*TMath.Gaus(x[0], par[2*i], sigma)*TMath.Gaus(x[1], par[2*i+1], sigma)
    return result

file = TFile('track_info.root','RECREATE')
track_info_tree = TTree('track_info', 'MC track_info')
track_x = array('f',[0])
track_y = array('f',[0])
integral50 = array('f',[0])
calibflag = array('i',[0])
track_info_tree.Branch('track_x', track_x, 'track_x/F')
track_info_tree.Branch('track_y', track_y, 'track_y/F')
track_info_tree.Branch('integral50', integral50, 'integral50/F')
track_info_tree.Branch('calibflag', calibflag, 'calibflag/I')

hits = 330000

MPV = 104
sigma = 10.5
center_x = -0.025
center_y = 0.2
sigma_x = 0.1
sigma_y = 0.12
xmin = -0.18
xmax = 0.13
ymin = 0.03
ymax = 0.35
canvas = TCanvas('canvas', 'canvas')
canvas.cd()

f_lateral = TF2('f_lateral', LateralShape, xmin, xmax, ymin, ymax, 12)
f_lateral.SetNpx(80)
f_lateral.SetNpy(80)
# 6 gaus centers: x1    y1     x2   y2     x3   y3    x4    y4     x5    y5   x6   y6
par = np.array([-0.06, 0.27, 0.02, 0.27, 0.02, 0.2, -0.06, 0.2, -0.06, 0.13, 0.02, 0.13])
f_lateral.SetParameters(par)
f_lateral.Draw('surf1')
raw_input('blah')
a = Double()
b = Double()

integral50_max = 5000
i = 0
j = 0
today = datetime.today()
seed = int((today-datetime(today.year, today.month, today.day , 0, 0, 0, 0)).total_seconds() % 1800 *1e6)
gRandom.SetSeed(seed)
while i < hits and j < 2*hits:
    f_lateral.GetRandom2(a,b)
    track_x[0] = a
    track_y[0] = b
    integral50[0] = gRandom.Landau(MPV, sigma)
    if xmin < track_x[0] < xmax and ymin < track_y[0] < ymax and integral50[0] < integral50_max:
        track_info_tree.Fill()
        i += 1
    j += 1
if not j<2*hits:
    print "abbruch.."

file.Write()