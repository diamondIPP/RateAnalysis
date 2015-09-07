from AllClasses import *
from ROOT import TCanvas

runplanType     = "rate_scan"
diamond         = "II6-B2"
rate_analyses = {}



rateScanCanvas = TCanvas("rateScanCanvas", "rateScanCanvas")


sel = RunSelection()
runplans =sel.runplan[runplanType].keys()

for runplanNr in runplans:
    sel.UnselectAll()
    sel.SelectRunsFromRunPlan(number=runplanNr, type_=runplanType)      # select all runs from a certain run plan
    sel.UnSelectUnlessDiamond(diamondname=diamond)                      # unselect runs, which dont have the diamond in it's channels
    sel.SelectDiamondRuns(diamondname=diamond, only_selected_runs=True) # activate only a certain diamond for analysis

    print "Selected runs to analyze:"
    sel.ShowSelectedRuns(commentlength=0)

    savePath = "Results/{runplan}/{diamond}/".format(runplan=runplanType+str(runplanNr), diamond=diamond)
    sel.SetSaveDirectory(savePath)

    rate_analyses[runplanNr] = AnalysisCollection(sel)

    rate_analyses[runplanNr].MakePreAnalysis(savePlot=True)
    rate_analyses[runplanNr].ShowPulserRates()



    rate_analyses[runplanNr].ShowSignalVSRate(canvas=rateScanCanvas)
    sel.SavePlots(savename="SignalVSRate"+str(runplanNr)+".root", canvas=rateScanCanvas)

sel.UnselectAll()
sel.SelectRunsInRange(439,459)      # select all runs from a certain run plan
sel.SetChannels()

print "Selected runs to analyze:"
sel.ShowSelectedRuns(commentlength=0)

savePath = "Results/new1/"
sel.SetSaveDirectory(savePath)

rate_analyses[10] = AnalysisCollection(sel)

rate_analyses[10].MakePreAnalysis(savePlot=True)
rate_analyses[10].ShowPulserRates()
rate_analyses[10].ShowSignalVSRate(canvas=rateScanCanvas)
sel.SavePlots(savename="SignalVSRate"+str(10)+".root", canvas=rateScanCanvas)

sel.UnselectAll()
sel.SelectRunsInRange(463,478)      # select all runs from a certain run plan
sel.SetChannels()

print "Selected runs to analyze:"
sel.ShowSelectedRuns(commentlength=0)

savePath = "Results/new2/"
sel.SetSaveDirectory(savePath)

rate_analyses[11] = AnalysisCollection(sel)

rate_analyses[11].MakePreAnalysis(savePlot=True)
rate_analyses[11].ShowPulserRates()
rate_analyses[11].ShowSignalVSRate(canvas=rateScanCanvas)
sel.SavePlots(savename="SignalVSRate"+str(11)+".root", canvas=rateScanCanvas)