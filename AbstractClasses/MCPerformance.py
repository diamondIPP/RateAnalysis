from AbstractClasses.AnalysisCollection import AnalysisCollection
from AbstractClasses.AnalysisClass import Analysis
from AbstractClasses.MCRun import MCRun
from AbstractClasses.ConfigClass import Config
import ROOT
import gc
from datetime import datetime


class MCPerformance(AnalysisCollection):

    def DoSignalHeightScan(self, heights=None, hits_per_height=30000):
        gc.disable()
        tries = 15
        if heights is None:
            heights = [0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.5, 0.8, 1.0]
        success_prob = []
        ghost_prob = []
        cycle_nr = 0
        cycles = tries*len(heights)
        starttime = datetime.today()
        for height in heights: # add more statistics for each height, not just one try..
            fails = 0
            tot_ghosts = 0
            peaks_generated = 0
            for repetition in xrange(tries):
                cycle_nr += 1
                print "\n{0}th repetition with Signal height set to: {1}\n".format(repetition, height)
                self.run_object = MCRun(validate=False,verbose=True,run_number=364)
                newAnalysis = Analysis(self.run_object)
                self.run_object.MCAttributes['PeakHeight'] = height
                self.run_object.MCAttributes['NumberOfHits'] = hits_per_height
                newAnalysis.FindMaxima()
                npeaks = newAnalysis.MaximaAnalysis.MCResults['TrueNPeaks']
                ninjas = newAnalysis.MaximaAnalysis.MCResults['Ninjas']
                ghosts = newAnalysis.MaximaAnalysis.MCResults['Ghosts']
                assert(npeaks>0), 'no peak in MC created'
                peaks_generated += npeaks
                fails += ninjas
                tot_ghosts += ghosts
                self.AddAnalysis(newAnalysis)
                elapsed_time = datetime.today() - starttime
                estimated_time = elapsed_time/cycle_nr*cycles
                remaining_time = estimated_time-elapsed_time
                print "\n\nAPPROXIMATED TIME LEFT: "+str(remaining_time)+"\n"
            success = 1.*(peaks_generated-fails)/peaks_generated
            ghost = 1.*ghosts/tries
            success_prob.append(success)
            ghost_prob.append(ghost)


        canvas = ROOT.TCanvas('canvas', 'canvas')
        canvas.cd()
        graph1 = ROOT.TGraph()
        graph1.SetNameTitle('graph1', 'success')
        graph2 = ROOT.TGraph()
        graph2.SetNameTitle('graph2', 'ghosts')
        for i in xrange(len(heights)):
            graph1.SetPoint(i, heights[i], success_prob[i])
            graph2.SetPoint(i, heights[i], ghost_prob[i])
        graph1.Draw('ALP*')
        graph2.Draw('SAME LP*')
        answer = raw_input('Wanna crash?')
        if answer == 'yes':
            gc.collect()
