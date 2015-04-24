from AbstractClasses.AnalysisCollection import AnalysisCollection
from AbstractClasses.AnalysisClass import Analysis
from AbstractClasses.MCRun import MCRun
from AbstractClasses.Elementary import Elementary
from AbstractClasses.ConfigClass import Config
import ROOT
from ROOT import TFile, TTree
import os
import shutil
import gc
from array import array
from datetime import datetime
'''
To improve:
    import config file
    write config file for later imports
    prohibit memory overfill
    include in framework ?
'''

class MCPerformance(Elementary):

    def __init__(self, verbose = False):
        Elementary.__init__(self, verbose=verbose)
        # Settings:
        self.tries = 20 # number of repetitions for a certain signal height
        self.min_percent = 5   # Quantiles
        self.max_percent = 99  # Quantiles

    def DoSignalHeightScan(self, heights=None, hits_per_height=300000):
        gc.disable()
        starttime = datetime.today()


        # ROOT Logfile:
        path = "MC/Performance_Results/"+str(starttime)
        os.makedirs(path)
        rootfile = TFile(path+'/MCPerformanceLog.root','RECREATE')
        LogTree = TTree('LogTree', 'MC Log Tree')
        RealSignalAmplitude = array('f',[0])
        Repetition = array('i', [0])
        TrueNPeaks = array('i', [0])
        Ninjas = array('i',[0])
        Ghosts = array('i',[0])
        RecSA_Quantiles = array('f',[0])
        RecSA_MinMax = array('f', [0])

        LogTree.Branch('RealSignalAmplitude', RealSignalAmplitude, 'RealSignalAmplitude/F')
        LogTree.Branch('Repetition', Repetition, 'Repetition/I')
        LogTree.Branch('TrueNPeaks', TrueNPeaks, 'TrueNPeaks/I')
        LogTree.Branch('Ninjas', Ninjas, 'Ninjas/I')
        LogTree.Branch('Ghosts', Ghosts, 'Ghosts/I')
        LogTree.Branch('RecSA_Quantiles', RecSA_Quantiles, 'RecSA_Quantiles/F')
        LogTree.Branch('RecSA_MinMax', RecSA_MinMax, 'RecSA_MinMax/F')


        # copy Config files:
        shutil.copy("Configuration/MonteCarloConfig.cfg", path+"/MonteCarloConfig.cfg")
        shutil.copy("Configuration/AnalysisConfig.cfg", path+"/AnalysisConfig.cfg")


        if heights is None:
            heights = [0.001, 0.05, 0.08, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.5, 0.8, 1.0]

        # infofile:
        infofile = open(path+"/info.txt", "w")
        infofile.write("DoSignalHeightScan\n\n")
        infofile.write("Number of Repetitions for each Amplitude: "+str(self.tries)+"\n")
        infofile.write("Number of different Amplitudes:           "+str(len(heights))+"\n")
        infofile.write("Hits per Amplitude:                       "+str(hits_per_height)+"\n")

        success_prob = []
        ghost_prob = []
        cycle_nr = 0
        cycles = self.tries*len(heights)
        for height in heights: # add more statistics for each height, not just one try..
            fails = 0
            tot_ghosts = 0
            peaks_generated = 0
            for repetition in xrange(self.tries):
                cycle_nr += 1
                print "\n{0}th repetition with Signal height set to: {1}\n".format(repetition, height)
                run_object = MCRun(validate=False,verbose=self.verbose,run_number=364)
                print "newAnalysis = Analysis(run_object)"
                newAnalysis = Analysis(run_object)
                run_object.MCAttributes['PeakHeight'] = height
                run_object.SetNumberOfHits(hits_per_height)
                print "newAnalysis.FindMaxima()"
                newAnalysis.FindMaxima()
                print "newAnalysis.FindMinima()"
                newAnalysis.FindMinima()
                npeaks = newAnalysis.ExtremeAnalysis.ExtremaResults['TrueNPeaks']
                ninjas = newAnalysis.ExtremeAnalysis.ExtremaResults['Ninjas']
                ghosts = newAnalysis.ExtremeAnalysis.ExtremaResults['Ghosts']
                maxima = newAnalysis.ExtremeAnalysis.ExtremaResults['FoundMaxima']
                minima = newAnalysis.ExtremeAnalysis.ExtremaResults['FoundMinima']
                # Reconstruct Signal Amplitude:
                if len(maxima)*len(minima)>0:
                    maxbin = newAnalysis.ExtremeAnalysis.Pad.GetBinByCoordinates(*(maxima[0]))
                    maxbin.FitLandau()
                    minbin = newAnalysis.ExtremeAnalysis.Pad.GetBinByCoordinates(*(minima[0]))
                    minbin.FitLandau()
                    rec_sa_minmax = maxbin.Fit['MPV']/minbin.Fit['MPV']-1.
                else:
                    rec_sa_minmax = -99
                q = array('d', [1.*self.min_percent/100., 1.*self.max_percent/100.])
                y = array('d', [0,0])
                newAnalysis.ExtremeAnalysis.CreateMeanSignalHistogram()
                newAnalysis.ExtremeAnalysis.MeanSignalHisto.GetQuantiles(2, y, q)
                rec_sa_quantiles = y[1]/y[0]-1.

                # Fill ROOT file:
                RealSignalAmplitude[0] = height
                Repetition[0]          = repetition
                TrueNPeaks[0]          = npeaks
                Ninjas[0]              = ninjas
                Ghosts[0]              = ghosts
                RecSA_Quantiles[0]     = rec_sa_quantiles
                RecSA_MinMax[0]        = rec_sa_minmax
                LogTree.Fill()

                assert(npeaks>0), 'no peak in MC created'
                peaks_generated += npeaks
                fails += ninjas
                tot_ghosts += ghosts
                # self.AddAnalysis(newAnalysis)
                del newAnalysis
                del run_object
                elapsed_time = datetime.today() - starttime
                estimated_time = elapsed_time/cycle_nr*cycles
                remaining_time = estimated_time-elapsed_time
                print "\n\nAPPROXIMATED TIME LEFT: "+str(remaining_time)+"\n"
            success = 1.*(peaks_generated-fails)/peaks_generated
            ghost = 4.*ghosts/self.tries
            success_prob.append(success)
            ghost_prob.append(ghost)

        print "Write ROOT-file"
        rootfile.Write()
        print "ROOT File written. Write infofile"
        infofile.write("\nTotal Time elapsed: "+str(datetime.today() - starttime))
        infofile.close()
        print "infofile written"

        # canvas = ROOT.TCanvas('canvas', 'canvas') # HERE IT CRASHES DUE TO MEMORY PROBLEMS
        # canvas.cd()
        # graph1 = ROOT.TGraph()
        # graph1.SetNameTitle('graph1', 'success')
        # graph1.SaveAs(path+"/SuccessGraph.root")
        # graph2 = ROOT.TGraph()
        # graph2.SetNameTitle('graph2', 'ghosts')
        # graph2.SaveAs(path+"/GhostsGraph.root")
        # for i in xrange(len(heights)):
        #     graph1.SetPoint(i, heights[i], success_prob[i])
        #     graph2.SetPoint(i, heights[i], ghost_prob[i])
        # graph1.Draw('ALP*')
        # graph2.Draw('SAME LP*')
        # self.SavePlots("PerformanceResult", "png", path+"/")
        answer = raw_input('Wanna crash?')
        if answer == 'yes':
            gc.collect()

    def DoGhostScan(self):
        pass
