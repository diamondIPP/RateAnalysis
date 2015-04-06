from ROOT import TFile, TTree, gRandom
from array import array
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt


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

hits = 350000

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
integral50_max = 5000

calibflag[0] = 0

i = 0
j = 0
today = datetime.today()
seed = int((today-datetime(today.year, today.month, today.day , 0, 0, 0, 0)).total_seconds() % 1800 *1e6)
gRandom.SetSeed(seed)
while i < hits and j < 2*hits:
    track_x[0] = gRandom.Gaus(center_x, sigma_x)
    track_y[0] = gRandom.Gaus(center_y, sigma_y)
    integral50[0] = gRandom.Landau(MPV, sigma)
    if xmin < track_x[0] < xmax and ymin < track_y[0] < ymax and integral50[0] < integral50_max:
        track_info_tree.Fill()
        i += 1
    j += 1
if not j<2*hits:
    print "abbruch.."

file.Write()