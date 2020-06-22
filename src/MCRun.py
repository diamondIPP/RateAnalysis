from RunClass import Run
from ROOT import TFile, TMath, gRandom, TCanvas, TTree, TF2, Double, gPad, gStyle
from array import array
from datetime import datetime
import os
import types as t
import numpy as np
import ConfigParser

class MCRun(Run):

    def __init__(self, validate = True, run_number = None, verbose=False):
        Run.__init__(self, validate=validate, number=run_number, verbose=verbose)
        self.MCAttributes = {
            'NPeaks': None,
            'PeakHeight': 0.5,
            'PeakSigmaX_min': 0.02,
            'PeakSigmaX_max': 0.07,
            'PeakSigmaY_min': 0.02,
            'PeakSigmaY_max': 0.07,
            'SpecialDistribution': "No",
            'NumberOfHits': 99,
            'HitDistributionMode': '',
            'TrackResolution': 0.,
            'SignalMode': '',
            'Landau_MPV_bkg': 80,
            'Landau-Sigma': 10,
            'integral50_max': 500,
            'MCRunPath': 'runs/MC/MC_{0}/', # {0}: placeholder for run number
            'DrawRealDistribution': False,
            'Save': True
        }
        self.NumberOfHits = 300000
        self.IsMonteCarlo = True
        self.Data = {
            'track_x': [],
            'track_y': [],
            'integral50': []
        }
        self.DataIsMade = False
        self.ShowAndWait = False
        self.LoadMCConfig()

    def LoadMCConfig(self):
        parser = ConfigParser.ConfigParser()
        parser.read('config/MonteCarloConfig.cfg')
        SignalMode = parser.get('SIGNAL-DISTRIBUTION','SignalMode')

        # LOAD SIGNAL DISTRIBUTION SETTINGS
        npeaks = parser.get('SIGNAL-DISTRIBUTION','NPeaks')
        if npeaks == 'None':
            npeaks = None
        else:
            npeaks = int(npeaks)
        SpecialDistribution     = parser.get('SIGNAL-DISTRIBUTION','SpecialDistribution')
        assert(SpecialDistribution in ["No", "False", "4Peaks", "4peaks", "Central", "central", "L", "3Peaks", "Testpattern", "testpattern", "Test", "test"]), "Bad MC Config file. SpecialDistribution: "+SpecialDistribution+" in SIGNAL-DSITRIBUTION unknown"
        PeakHeight          = parser.getfloat('SIGNAL-DISTRIBUTION','PeakHeight')
        Landau_MPV_bkg      = parser.getint('SIGNAL-DISTRIBUTION','Landau_MPV_bkg')
        Landau_Sigma        = parser.getint('SIGNAL-DISTRIBUTION','Landau_Sigma')
        integral50_max      = parser.getint('SIGNAL-DISTRIBUTION','integral50_max')
        PeakSigmaX_min      = parser.getfloat('SIGNAL-DISTRIBUTION','PeakSigmaX_min')
        PeakSigmaX_max      = parser.getfloat('SIGNAL-DISTRIBUTION','PeakSigmaX_max')
        PeakSigmaY_min      = parser.getfloat('SIGNAL-DISTRIBUTION','PeakSigmaY_min')
        PeakSigmaY_max      = parser.getfloat('SIGNAL-DISTRIBUTION','PeakSigmaY_max')

        # LOAD HIT DISTRIBUTION SETTINGS
        NumberOfHits = int(parser.get('HIT-DISTRIBUTION','NumberOfHits'))
        TrackResolution = float(parser.get('HIT-DISTRIBUTION','TrackResolution'))
        HitDistributionMode = parser.get('HIT-DISTRIBUTION','HitDistributionMode')

        # LOAD SAVE SETTINGS
        MCRunPath = parser.get('SAVE','MCRunPath')
        Save = parser.getboolean('SAVE', 'Save')

        # LOAD DISPLAY SETTINGS
        ShowAndWait = parser.getboolean('DISPLAY','ShowAndWait')
        DrawRealDistribution = parser.getboolean('DISPLAY','DrawRealDistribution')


        # APPLY SETTINGS:
        self.MCAttributes['NPeaks'] = npeaks
        self.MCAttributes['PeakHeight'] = PeakHeight
        self.MCAttributes['SpecialDistribution'] = SpecialDistribution
        self.MCAttributes['NumberOfHits'] = NumberOfHits
        self.NumberOfHits = NumberOfHits
        self.SetHitDistributionMode(HitDistributionMode)
        self.MCAttributes['SignalMode'] = SignalMode
        self.MCAttributes['Landau_MPV_bkg'] = Landau_MPV_bkg
        self.MCAttributes['Landau_Sigma'] = Landau_Sigma
        self.MCAttributes['integral50_max'] = integral50_max
        self.MCAttributes['MCRunPath'] = MCRunPath
        self.MCAttributes['TrackResolution'] = TrackResolution
        self.MCAttributes['PeakSigmaX_min'] = PeakSigmaX_min
        self.MCAttributes['PeakSigmaX_max'] = PeakSigmaX_max
        self.MCAttributes['PeakSigmaY_min'] = PeakSigmaY_min
        self.MCAttributes['PeakSigmaY_max'] = PeakSigmaY_max
        self.MCAttributes['DrawRealDistribution'] = DrawRealDistribution
        self.MCAttributes['Save'] = Save

        self.ShowAndWait = ShowAndWait

    def SetHitDistributionMode(self, mode):
        if mode == 'Import' or mode == 'Manual' or mode == 'Uniform':
            self.MCAttributes['HitDistributionMode'] = mode
        else:
            assert(False), 'Wrong Hit Distribution Mode - mode has to be `Import` or `Manual` or `Uniform`'
        if mode == 'Import':
            self.CountHistoFile = TFile('MCInputs/364counthisto.root')
            self.counthisto = self.CountHistoFile.Get('counthisto')

    def SetSignalMode(self, mode):
        if mode == 'Landau' or mode == 'Gaus':
            self.MCAttributes['SignalMode'] = mode
        else:
            assert(False), 'Wrong Signal Mode - mode has to be `Landau` or `Gaus`'

    def SetNumberOfHits(self, hits): # Not nice
        self.NumberOfHits = hits
        self.MCAttributes['NumberOfHits'] = hits


    def SetOnlyCentralPeak(self, bool):
        if bool:
            self.MCAttributes['SpecialDistribution'] = "Central"
        else:
            self.MCAttributes['SpecialDistribution'] = "No"

    def ResetMC(self):
        self.DataIsMade = False
        self.Data = {
            'track_x': [],
            'track_y': [],
            'integral50': []
        }

    def ShowMCConfig(self):
        print '\nHit Distribution Mode is Set to: '+self.MCAttributes['HitDistributionMode']
        print 'Signal Mode is Set to: '+self.MCAttributes['SignalMode']
        print 'Number of Signals to produce: ', self.NumberOfHits
        if self.MCAttributes['NPeaks'] != None:
            print self.MCAttributes['NPeaks']," Peaks generated. Each Peak with amplitude: ", self.MCAttributes['PeakHeight'], '\n'
        else:
            print '\n'

    def Simulate(self, save=None, draw=None):
        '''

        :param save: if True: saves the root file as well as the true signal distribution
        :param draw: if True draws the signal distribution
        :param ask:
        :return:
        '''

        # Settings:
        MPV = self.MCAttributes['Landau_MPV_bkg'] # background for Landau MPV distribution
        sigma = self.MCAttributes['Landau_Sigma'] # Landau scaling sigma
        NPeaks = self.MCAttributes['NPeaks']
        MCRunPath = self.MCAttributes['MCRunPath'].format(self.RunNumber)
        xmin = self.diamond.Position['xmin'] + 0.01
        xmax = self.diamond.Position['xmax'] - 0.01
        ymin = self.diamond.Position['ymin'] + 0.01
        ymax = self.diamond.Position['ymax'] - 0.01
        center_x = (xmax+xmin)/2.
        center_y = (ymax+ymin)/2.
        if draw != None:
            assert(type(draw) == t.BooleanType), "draw has to be boolean type"
            self.MCAttributes['DrawRealDistribution'] = draw
        if save != None:
            assert(type(save) == t.BooleanType), "save argument has to be of type boolean"
            self.MCAttributes['Save'] = save

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

        def CreateRandomPeaksConfig(xmin, xmax, ymin, ymax, bkg = 120, peak_height = 0.5, npeaks = NPeaks):
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
            if npeaks == None:
                npeaks = int(round(gRandom.Uniform(0,15)))

            if self.MCAttributes['SpecialDistribution'] in ["Central", "central"]:
                npeaks = 1
                dxy = 0.02
                parameters = np.zeros(7)
                parameters[0] = npeaks
                parameters[1] = bkg
                parameters[2] = peak_height
                parameters[3] = gRandom.Uniform(center_x-dxy, center_x+dxy)
                parameters[4] = gRandom.Uniform(center_y-dxy, center_y+dxy)
                parameters[5] = gRandom.Uniform(self.MCAttributes['PeakSigmaX_min'], self.MCAttributes['PeakSigmaX_max'])
                parameters[6] = gRandom.Uniform(self.MCAttributes['PeakSigmaY_min'], self.MCAttributes['PeakSigmaY_max'])
            elif self.MCAttributes['SpecialDistribution'] in ["4Peaks", "4peaks", "L", "3Peaks"]:
                npeaks = 4
                dxy = 0.02
                parameters = np.zeros(3+4*npeaks)
                parameters[0] = npeaks
                parameters[1] = bkg
                parameters[2] = peak_height
                # peaks:
                peaknr = 0
                for x in [center_x-0.07, center_x+0.07]:
                    for y in [center_y-0.07, center_y+0.07]:
                        parameters[3+4*peaknr] = gRandom.Uniform(x-dxy, x+dxy)
                        parameters[4+4*peaknr] = gRandom.Uniform(y-dxy, y+dxy)
                        parameters[5+4*peaknr] = gRandom.Uniform(self.MCAttributes['PeakSigmaX_min'], self.MCAttributes['PeakSigmaX_max'])
                        parameters[6+4*peaknr] = gRandom.Uniform(self.MCAttributes['PeakSigmaY_min'], self.MCAttributes['PeakSigmaY_max'])
                        peaknr += 1
                if self.MCAttributes['SpecialDistribution'] in ["L", "3Peaks"]:
                    npeaks = 3
                    parameters[0] = npeaks
                    parameters = parameters[:3+4*npeaks]
            elif self.MCAttributes['SpecialDistribution'] in ["Testpattern", "testpattern", "Test", "test"]:
                npeaks = 8
                Center_x = gRandom.Uniform(center_x-0.01, center_x+0.01)
                Center_y = gRandom.Uniform(center_y-0.01, center_y+0.01)
                parameters = np.zeros(3+4*npeaks)
                parameters[0]  = npeaks
                parameters[1]  = bkg
                parameters[2]  = peak_height
                parameters[3]  = Center_x - 0.1
                parameters[4]  = Center_y + 0.08
                parameters[5]  = 0.02
                parameters[6]  = 0.02
                parameters[7]  = Center_x - 0.04
                parameters[8]  = Center_y + 0.08
                parameters[9]  = 0.02
                parameters[10] = 0.02
                parameters[11]  = Center_x - 0.1
                parameters[12]  = Center_y
                parameters[13]  = 0.025
                parameters[14] = 0.025
                parameters[15]  = Center_x
                parameters[16]  = Center_y
                parameters[17]  = 0.02
                parameters[18] = 0.02
                parameters[19]  = Center_x + 0.1
                parameters[20]  = Center_y
                parameters[21]  = 0.03
                parameters[22] = 0.03
                parameters[23]  = Center_x - 0.04
                parameters[24]  = Center_y - 0.06
                parameters[25]  = 0.015
                parameters[26] = 0.015
                parameters[27]  = Center_x - 0.12
                parameters[28]  = Center_y - 0.12
                parameters[29]  = 0.03
                parameters[30] = 0.03
                parameters[31]  = Center_x + 0.08
                parameters[32]  = Center_y - 0.12
                parameters[33]  = 0.04
                parameters[34] = 0.04
            else:
                parameters = np.zeros(3+4*npeaks)
                parameters[0] = npeaks
                parameters[1] = bkg
                parameters[2] = peak_height
                for i in xrange(npeaks):
                    parameters[3+4*i] = gRandom.Uniform(xmin, xmax)
                    parameters[4+4*i] = gRandom.Uniform(ymin, ymax)
                    parameters[5+4*i] = gRandom.Uniform(0.02, 0.07)
                    parameters[6+4*i] = gRandom.Uniform(0.02, 0.07)
            self.MCAttributes['NPeaks'] = npeaks
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

        # create track_info ROOT file
        if not os.path.exists(MCRunPath):
            os.makedirs(MCRunPath)
        if self.MCAttributes['Save']:
            file = TFile(MCRunPath+'track_info.root','RECREATE')
        self.track_info = TTree('track_info', 'MC track_info')
        track_x = array('f',[0])
        track_y = array('f',[0])
        integral50 = array('f',[0])
        calibflag = array('i',[0])
        calib_offset = array('i',[0])
        time_stamp = array('f', [0])
        self.track_info.Branch('track_x', track_x, 'track_x/F')
        self.track_info.Branch('track_y', track_y, 'track_y/F')
        self.track_info.Branch('integral50', integral50, 'integral50/F')
        self.track_info.Branch('calibflag', calibflag, 'calibflag/I')
        self.track_info.Branch('calib_offset', calib_offset, 'calib_offset/I')
        self.track_info.Branch('time_stamp', time_stamp, 'time_stamp/F')


        # Create Manual Hit Distribution:
        if self.MCAttributes['HitDistributionMode'] == 'Manual':
            dx = 0.08
            dy = 0.07
            f_lateral = TF2('f_lateral', ManualHitDistribution, xmin, xmax, ymin, ymax, 12)
            f_lateral.SetNpx(80)
            f_lateral.SetNpy(80)
            # 6 gaus centers:
            par = np.array([center_x-dx/2.,     # x1
                            center_y+dy,        # y1
                            center_x+dx/2.,     # x2
                            center_y+dy,        # y2
                            center_x+dx/2.,     # x3
                            center_y,           # y3
                            center_x-dx/2.,     # x4
                            center_y,           # y4
                            center_x-dx/2.,     # x5
                            center_y-dy,        # y5
                            center_x+dx/2.,     # x6
                            center_y-dy         # y6
                            ])
            f_lateral.SetParameters(par)
        a = Double()
        b = Double()

        # Generate Signal Distribution:
        if self.MCAttributes['SignalMode'] == 'Landau':
            self.SignalParameters = CreateRandomPeaksConfig(xmin, xmax, ymin, ymax, peak_height=self.MCAttributes['PeakHeight'], bkg=MPV)
        else:
            self.SignalParameters = CreateRandomPeaksConfig(xmin, xmax, ymin, ymax, peak_height=self.MCAttributes['PeakHeight'], bkg=100)
        f_signal = TF2('f_signal', SignalShape, xmin, xmax, ymin, ymax, len(self.SignalParameters))
        f_signal.SetNpx(40) # Resolution
        f_signal.SetNpy(40)
        f_signal.SetParameters(self.SignalParameters)
        if self.MCAttributes['DrawRealDistribution']:
            canvas = TCanvas('canvas', 'canvas')
            canvas.cd()
            gStyle.SetPalette(55) # a Rain Bow palette == used.
            gStyle.SetNumberContours(8)
            f_signal.Draw('surf1')
            gPad.Print(MCRunPath+'RealSignalDistribution.png')
            if self.ShowAndWait:
                answer = raw_input('for data creation, type `yes`: ')
            else:
                answer = 'yes'
        else:
            answer = 'yes'

        # Set the Hit distribution for Manual or Import
        if self.MCAttributes['HitDistributionMode'] == 'Manual':
            HitsTemplate = f_lateral
        elif self.MCAttributes['HitDistributionMode'] == 'Import':
            HitsTemplate = self.counthisto

        mc_start_timestamp = 42. # arbitrary timestamp for the first MC event
        MCEventDeltaTime = 30.*60./300000. # 1800s/300000Hits = 0.006s/Hit
        tmp_time = mc_start_timestamp
        # Generate Toy Data:
        if answer == 'yes':
            if self.verbose:
                self.ShowMCConfig()
            self.verbose_print('Creating Toy Data with {0} Hits'.format(self.NumberOfHits))
            integral50_max = self.MCAttributes['integral50_max'] # Maximum of Signal response allowed (data: 500 ?)
            i = 0
            j = 0
            while i < self.NumberOfHits and j < 2*self.NumberOfHits:
                # Get x and y
                if self.MCAttributes['HitDistributionMode'] == 'Uniform':
                    track_x[0] = gRandom.Uniform(xmin,xmax)
                    track_y[0] = gRandom.Uniform(ymin,ymax)
                else:
                    HitsTemplate.GetRandom2(a,b)
                    track_x[0] = gRandom.Gaus(a, self.MCAttributes['TrackResolution']) # 0.002 = 20mu track resolution (first guess)
                    track_y[0] = gRandom.Gaus(b, self.MCAttributes['TrackResolution'])

                # Get Signal at x and y
                if self.MCAttributes['SignalMode'] == 'Landau':
                    integral50[0] = gRandom.Landau(f_signal(track_x[0], track_y[0]), sigma)
                else:
                    integral50[0] = gRandom.Gaus(f_signal(track_x[0], track_y[0]), 0.6*f_signal(track_x[0], track_y[0])-33)

                # check if found values fulfill requirements
                if xmin < track_x[0] < xmax and ymin < track_y[0] < ymax and integral50[0] < integral50_max:
                    tmp_time += MCEventDeltaTime
                    time_stamp[0] = tmp_time
                    self.track_info.Fill()
                    self.Data['track_x'].append(track_x[0])
                    self.Data['track_y'].append(track_y[0])
                    self.Data['integral50'].append(integral50[0])
                    i += 1
                j += 1
            if not j<2*self.NumberOfHits: # if too many times requirements were not fulfilled
                assert(False), "Bad MC Parameters"

            # Save root file and true Signal Shape:
            if self.MCAttributes['Save']:
                file.Write()
                f_signal.SaveAs(MCRunPath+'RealSignalDistribution.root')
                self.verbose_print(MCRunPath + 'track_info.root has been written')
                self.TrackingPadAnalysis['ROOTFile'] = MCRunPath+'track_info.root'

            # Print Settings of created data:
            if self.verbose:
                print "\nToydata containing {:.0f} peaks generated.".format(self.SignalParameters[0])
                if self.MCAttributes['SignalMode'] == 'Landau':
                    for i in xrange(int(self.SignalParameters[0])):
                        x = self.SignalParameters[3+4*i]
                        y = self.SignalParameters[4+4*i]
                        print "Peak {0:.0f} at position: ({1:.3f}/{2:.3f}) with Laundau Response MPV: {3:.2f} Sigma: {4:.1f}".format(i+1, x, y, f_signal(x, y),sigma)
                else:
                    integral50[0] = gRandom.Gaus(f_signal(track_x[0], track_y[0]), 0.6*f_signal(track_x[0], track_y[0])-33)
                    for i in xrange(int(self.SignalParameters[0])):
                        x = self.SignalParameters[3+4*i]
                        y = self.SignalParameters[4+4*i]
                        print "Peak {0:.0f} at position: ({1:.3f}/{2:.3f}) with Gaussian Response Mean: {3:.2f} Sigma: {4:.1f}".format(i+1, x, y, f_signal(x, y), 0.6*f_signal(x, y)-33)

            self.DataIsMade = True
            self.verbose_print('Monte Carlo Toy Data created.\n')
