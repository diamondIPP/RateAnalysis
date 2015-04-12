from RunClass import Run
from ROOT import TFile, TMath, gRandom, TCanvas, TTree, TF2, Double, gPad
from array import array
from datetime import datetime
import os
import numpy as np

class MCRun(Run):

    def __init__(self, validate = True, run_number = None):
        Run.__init__(self, validate=validate, run_number=run_number)
        self.HitDistributionMode = ''
        self.SignalMode = ''
        self.NumberOfHits = 300000
        self.SetHitDistributionMode('Manual')
        self.SetSignalMode('Landau')
        self.IsMonteCarlo = True
        self.Data = {
            'track_x': [],
            'track_y': [],
            'integral50': []
        }
        self.DataIsMade = False

    def SetHitDistributionMode(self, mode):
        if mode is 'Import' or mode is 'Manual' or mode is 'Uniform':
            self.HitDistributionMode = mode
        else:
            assert(False), 'Wrong Hit Distribution Mode - mode has to be `Import` or `Manual` or `Uniform`'
        if mode is 'Import':
            self.CountHistoFile = TFile('MCInputs/364counthisto.root')
            self.counthisto = self.CountHistoFile.Get('counthisto')

    def SetSignalMode(self, mode):
        if mode is 'Landau' or mode is 'Gaus':
            self.SignalMode = mode
        else:
            assert(False), 'Wrong Signal Mode - mode has to be `Landau` or `Gaus`'

    def SetNumberOfHits(self, hits):
        self.NumberOfHits = hits

    def ResetMC(self):
        self.DataIsMade = False
        self.Data = {
            'track_x': [],
            'track_y': [],
            'integral50': []
        }

    def ShowMCConfig(self):
        print '\nHit Distribution Mode is Set to: '+self.HitDistributionMode
        print 'Signal Mode is Set to: '+self.SignalMode
        print 'Number of Signals to produce: ', self.NumberOfHits, '\n'

    def Simulate(self, save=True, draw=True):
        '''

        :param save: if True: saves the root file as well as the true signal distribution
        :param draw: if True draws the signal distribution
        :param ask:
        :return:
        '''

        # Settings:
        MPV = 80 # background for Landau MPV distribution
        sigma = 11 # Landau scaling sigma
        MCRunPath = 'runs/MCrun_{0}/'.format(self.run_number)
        xmin = self.diamond.Position['xmin'] + 0.01
        xmax = self.diamond.Position['xmax'] - 0.01
        ymin = self.diamond.Position['ymin'] + 0.01
        ymax = self.diamond.Position['ymax'] - 0.01
        # xmin = -0.18
        # xmax = 0.13
        # ymin = 0.03
        # ymax = 0.35
        print self.diamond.Specifications['Name']
        print xmin
        print xmax
        print ymin
        print ymax

        def ManualHitDistribution(x, par):
            '''
            Probability density function for hit distribution based on
            6 2d gaussian in a rectangular pattern. The PDF is NOT
            normalized to 1
            :param x: array, x[0]: x position, x[1]: y position
            :param par: parameter array;
            :return:
            '''
            result = 0.1 # bkg
            norm = 1.
            sigma = 0.04
            for i in xrange(len(par)/2):
                result += norm*TMath.Gaus(x[0], par[2*i], sigma)*TMath.Gaus(x[1], par[2*i+1], sigma)
            return result

        def CreateRandomPeaksConfig(xmin, xmax, ymin, ymax, bkg = 120, peak_height = 0.1, npeaks = None):
            '''
            Creates the input parameters for SignalShape, which describes the peak distribution in the
            Signal distribution
            :param xmin: window parameter
            :param xmax: window parameter
            :param ymin: window parameter
            :param ymax: window parameter
            :param bkg: background of signal response
            :param peak_height: height of one peak in % of bkg, if bkg==0: the peaks heights are set to 1
            :param npeaks: number of peaks in window - if none it picks random between 0 and 15
            :return: array containing the parameters
            '''
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

        def SignalShape(x, par):
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

        # Set seed for random number generator:
        today = datetime.today()
        seed = int((today-datetime(today.year, today.month, today.day , 0, 0, 0, 0)).total_seconds() % 1800 *1e6)
        gRandom.SetSeed(seed)

        if save:
            pass
        # create track_info ROOT file
        if not os.path.exists(MCRunPath):
            os.makedirs(MCRunPath)
        file = TFile(MCRunPath+'track_info.root','RECREATE')
        self.track_info = TTree('track_info', 'MC track_info')
        track_x = array('f',[0])
        track_y = array('f',[0])
        integral50 = array('f',[0])
        calibflag = array('i',[0])
        self.track_info.Branch('track_x', track_x, 'track_x/F')
        self.track_info.Branch('track_y', track_y, 'track_y/F')
        self.track_info.Branch('integral50', integral50, 'integral50/F')
        self.track_info.Branch('calibflag', calibflag, 'calibflag/I')


        if self.HitDistributionMode is 'Manual':
            f_lateral = TF2('f_lateral', ManualHitDistribution, xmin, xmax, ymin, ymax, 12)
            f_lateral.SetNpx(80)
            f_lateral.SetNpy(80)
            # 6 gaus centers: x1    y1     x2   y2     x3   y3    x4    y4     x5    y5   x6   y6
            par = np.array([-0.06, 0.27, 0.02, 0.27, 0.02, 0.2, -0.06, 0.2, -0.06, 0.13, 0.02, 0.13])
            f_lateral.SetParameters(par)
        a = Double()
        b = Double()

        if self.SignalMode is 'Landau':
            SignalParameters = CreateRandomPeaksConfig(xmin, xmax, ymin, ymax, peak_height=0.5, bkg=MPV)
        else:
            SignalParameters = CreateRandomPeaksConfig(xmin, xmax, ymin, ymax, peak_height=0.5, bkg=100)
        f_signal = TF2('f_signal', SignalShape, xmin, xmax, ymin, ymax, len(SignalParameters))
        f_signal.SetNpx(40)
        f_signal.SetNpy(40)
        f_signal.SetParameters(SignalParameters)
        if draw:
            canvas = TCanvas('canvas', 'canvas')
            canvas.cd()
            f_signal.Draw('surf1')
            gPad.Print(MCRunPath+'RealSignalDistribution.png')
            answer = raw_input('for data creation, type `yes`: ')
        else:
            answer = 'yes'

        # Set the Hit distribution for Manual or Import
        if self.HitDistributionMode is 'Manual':
            HitsTemplate = f_lateral
        elif self.HitDistributionMode is 'Import':
            HitsTemplate = self.counthisto

        # Get Toy Data
        if answer == 'yes':
            integral50_max = 5000 # Maximum of Signal response allowed (data: 500 ?)
            i = 0
            j = 0
            while i < self.NumberOfHits and j < 2*self.NumberOfHits:
                # Get x and y
                if self.HitDistributionMode is 'Uniform':
                    track_x[0] = gRandom.Uniform(xmin,xmax)
                    track_y[0] = gRandom.Uniform(ymin,ymax)
                else:
                    HitsTemplate.GetRandom2(a,b)
                    track_x[0] = gRandom.Gaus(a, 0.002) # 20mu track resolution (first guess)
                    track_y[0] = gRandom.Gaus(b, 0.002)

                # Get Signal at x and y
                if self.SignalMode is 'Landau':
                    integral50[0] = gRandom.Landau(f_signal(track_x[0], track_y[0]), sigma)
                else:
                    integral50[0] = gRandom.Gaus(f_signal(track_x[0], track_y[0]), 0.6*f_signal(track_x[0], track_y[0])-33)

                # check if found values fulfill requirements
                if xmin < track_x[0] < xmax and ymin < track_y[0] < ymax and integral50[0] < integral50_max:
                    self.track_info.Fill()
                    self.Data['track_x'].append(track_x[0])
                    self.Data['track_y'].append(track_y[0])
                    self.Data['integral50'].append(integral50[0])
                    i += 1
                j += 1
            if not j<2*self.NumberOfHits: # if too many times requirements were not fulfilled
                assert(False), "Bad MC Parameters"

            # Save root file and true Signal Shape:
            if save:
                file.Write()
                f_signal.SaveAs(MCRunPath+'RealSignalDistribution.root')

            # Print Settings of created data:
            print "Toydata containing {:.0f} peaks generated.".format(SignalParameters[0])
            if self.SignalMode is 'Landau':
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

            self.DataIsMade = True
            print 'Data madeeee'
