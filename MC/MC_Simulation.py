from ROOT import TCanvas, TF2, Double, TMath, gRandom, TFile, TTree, TH2D
import numpy as np
import ROOT
from array import array
from datetime import datetime

SignalMode = 'Landau' # 'Landau' or 'Gaus'
HitDistributionMode = 'Import' # 'Manual' or 'Import' or 'Uniform'

today = datetime.today()
seed = int((today-datetime(today.year, today.month, today.day , 0, 0, 0, 0)).total_seconds() % 1800 *1e6)
gRandom.SetSeed(seed)

if HitDistributionMode is 'Import':
    CountHistoFile = ROOT.TFile('MCInputs/364counthisto.root')
    counthisto = CountHistoFile.Get('counthisto')


def LateralShape(x,par):
    result = 0.1
    norm = 1.
    sigma = 0.04
    for i in xrange(len(par)/2):
        result += norm*TMath.Gaus(x[0], par[2*i], sigma)*TMath.Gaus(x[1], par[2*i+1], sigma)
    return result

def SignalShape(x,par):
    '''

    :param x:   x[0]: x position
                x[1]: y position
    :param par: par[0]: number of peaks
                par[1]: mean bkg signal
                par[2]: peak height in percent of bkg
                par[3]: peak1 x position
                par[4]: peak1 y position
                par[5]: peak1 sigma in x
                par[6]: peak1 sigma in y
                ...
                par[3+4*n]: peakn x position
                par[4+4*n]: peakn y position
                par[5+4*n]: peakn sigma in x
                par[6+4*n]: peakn sigma in y
    :return:
    '''
    norm = par[1]*par[2]
    result = par[1]
    if par[1] == 0:
        norm = 1
    for i in xrange(int(par[0])):
        result += norm*TMath.Gaus(x[0], par[3+4*i], par[5+4*i])*TMath.Gaus(x[1], par[4+4*i], par[6+4*i])
    return result

def CreateRandomPeaks(xmin, xmax, ymin, ymax, bkg = 120, peak_height = 0.1, npeaks = None):
    if npeaks is None:
        npeaks = int(round(gRandom.Uniform(0,15)))
    parameters = np.zeros(3+4*npeaks)
    parameters[0] = npeaks
    parameters[1] = bkg
    parameters[2] = peak_height
    for i in xrange(npeaks):
        parameters[3+4*i] = gRandom.Uniform(xmin, xmax)
        parameters[4+4*i] = gRandom.Uniform(ymin, ymax)
        parameters[5+4*i] = gRandom.Uniform(0.02, 0.07)
        parameters[6+4*i] = gRandom.Uniform(0.02, 0.07)
    return parameters

# create track_info ROOT file
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

hits = 300000

MPV = 80
sigma = 11
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

if HitDistributionMode is 'Manual':
    f_lateral = TF2('f_lateral', LateralShape, xmin, xmax, ymin, ymax, 12)
    f_lateral.SetNpx(80)
    f_lateral.SetNpy(80)
    # 6 gaus centers: x1    y1     x2   y2     x3   y3    x4    y4     x5    y5   x6   y6
    par = np.array([-0.06, 0.27, 0.02, 0.27, 0.02, 0.2, -0.06, 0.2, -0.06, 0.13, 0.02, 0.13])
    f_lateral.SetParameters(par)
    # f_lateral.Draw('surf1')
    # raw_input('blah')
a = Double()
b = Double()

if SignalMode == 'Landau':
    SignalParameters = CreateRandomPeaks(xmin, xmax, ymin, ymax, peak_height=0.5, bkg=MPV)
elif SignalMode == 'Gaus':
    SignalParameters = CreateRandomPeaks(xmin, xmax, ymin, ymax, peak_height=0.5, bkg=100)
else:
    assert(False), "wrong SignalMode, choose `Gaus` or `Landau`"
f_signal = TF2('f_signal', SignalShape, xmin, xmax, ymin, ymax, len(SignalParameters))
f_signal.SetNpx(40)
f_signal.SetNpy(40)
f_signal.SetParameters(SignalParameters)
f_signal.Draw('surf1')
ROOT.gPad.Print('RealSignalDistribution.png')
f_signal.SaveAs('RealSignalDistribution.root')
answer = raw_input('for data creation, type `yes`: ')

if HitDistributionMode is 'Manual':
    HitsTemplate = f_lateral
elif HitDistributionMode is 'Import':
    HitsTemplate = counthisto
elif HitDistributionMode is not 'Uniform':
    assert(False), 'Wrong HitDistributionMode; HitDistributionMode has to be either `Manual`, `Import` or `Uniform`.'

if answer == 'yes':
    integral50_max = 5000
    i = 0
    j = 0
    while i < hits and j < 2*hits:
        if HitDistributionMode is 'Uniform':
            track_x[0] = gRandom.Uniform(xmin,xmax)
            track_y[0] = gRandom.Uniform(ymin,ymax)
        else:
            HitsTemplate.GetRandom2(a,b)
            track_x[0] = gRandom.Gaus(a, 0.004) # 20mu track resolution
            track_y[0] = gRandom.Gaus(b, 0.004)
        if SignalMode is 'Landau':
            integral50[0] = gRandom.Landau(f_signal(track_x[0], track_y[0]), sigma)
        else:
            integral50[0] = gRandom.Gaus(f_signal(track_x[0], track_y[0]), 0.6*f_signal(track_x[0], track_y[0])-33)

        if xmin < track_x[0] < xmax and ymin < track_y[0] < ymax and integral50[0] < integral50_max:
            track_info_tree.Fill()
            i += 1
        j += 1
    if not j<2*hits:
        print "abbruch.."

    file.Write()

    print "Toydata containing {:.0f} peaks generated.".format(SignalParameters[0])
    if SignalMode == 'Landau':
        for i in xrange(int(SignalParameters[0])):
            x = SignalParameters[3+4*i]
            y = SignalParameters[4+4*i]
            print "Peak {0:.0f} at position: ({1:.3f}/{2:.3f}) with Laundau Response MPV: {3:.2f} Sigma: {4:.1f}".format(i+1, x, y, f_signal(x, y),sigma)
    else:
        integral50[0] = gRandom.Gaus(f_signal(track_x[0], track_y[0]), 0.6*f_signal(track_x[0], track_y[0])-33)
        for i in xrange(int(SignalParameters[0])):
            x = SignalParameters[3+4*i]
            y = SignalParameters[4+4*i]
            print "Peak {0:.0f} at position: ({1:.3f}/{2:.3f}) with Gaussian Response Mean: {3:.2f} Sigma: {4:.1f}".format(i+1, x, y, f_signal(x, y), 0.6*f_signal(x, y)-33)
