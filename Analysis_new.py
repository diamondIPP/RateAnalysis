from AbstractClasses.newAnalysis import Analysis
from AbstractClasses.RunClass import Run
from AbstractClasses.AnalysisCollection import AnalysisCollection
import ROOT


runnumbers1 = [415, 419, 420, 422, 426, 427, 428, 430, 432, 434, 436]
runnumbers2 = [342, 344, 346, 348, 350, 352, 354, 356, 358, 359]

runs1 = {}
runs2 = {}

collection = AnalysisCollection()
collection2 = AnalysisCollection()

for runnr1 in runnumbers1:
    print "col1", runnr1
    runs1[runnr1] = Run(runnr1, 3)

    analysis = Analysis(runs1[runnr1])
    analysis.RemoveBeamInterruptions()


    #print "beam interruptions: ", analysis.GetBeamInterruptions()

    #analysis.MakePreAnalysis()
    #canvas = ROOT.TCanvas("c1", "c1")
    #analysis.GetLandau(channel=0, canvas=canvas, color=ROOT.kBlue)
    #analysis.GetLandau(channel=3)

    # analysis.CreateHitMap(0)
    # analysis.CreateSignalMaps()

    collection.AddAnalysis(analysis)

for runnr2 in runnumbers2:
    print "col2"
    runs2[runnr2] = Run(runnr2, 3)

    analysis = Analysis(runs2[runnr2])
    analysis.RemoveBeamInterruptions()

    collection2.AddAnalysis(analysis)


print "CP1"
collection2.ShowSignalVSRate()
print "CP2"
collection.ShowSignalVSRate(canvas=collection2.ratecanvas)



raw_input("end")