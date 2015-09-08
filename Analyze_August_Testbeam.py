from AllClasses import *
from ROOT import TCanvas

runplanType     = "rate_scan"
diamond1         = "II6-B2"
diamond2         = "II6-97"
rate_analyses = {}



rateScanCanvas = TCanvas("rateScanCanvas", "rateScanCanvas")


sel = RunSelection()
runplans =sel.runplan[runplanType].keys()
runplans = ["1"]

for runplanNr in runplans:
    sel.UnselectAll()
    #sel.SelectRunsFromRunPlan(number=runplanNr, type_=runplanType)      # select all runs from a certain run plan
    #sel.UnSelectUnlessDiamond(diamondname=diamond)                      # unselect runs, which dont have the diamond in it's channels
    sel.SelectRunsInRange(75, 77)
    sel.SelectDiamondRuns(diamondname=diamond1, only_selected_runs=True) # activate only a certain diamond for analysis
    sel.SelectDiamondRuns(diamondname=diamond2, only_selected_runs=True)

    print "Selected runs to analyze:"
    sel.ShowSelectedRuns(commentlength=0)

    rate_analyses[runplanNr] = AnalysisCollection(sel)

    sel.SetSaveDirectory("Results/{runplan}/PreAnalysis/".format(runplan=runplanType+str(runplanNr)))
    rate_analyses[runplanNr].MakePreAnalysises(savePlot=True)
    sel.SetSaveDirectory("Results/{runplan}/PulserRates/".format(runplan=runplanType+str(runplanNr)))
    rate_analyses[runplanNr].ShowPulserRates()

    for run in rate_analyses[runplanNr].collection.keys():
        sel.SetSaveDirectory("Results/{runplan}/SignalMaps/".format(runplan=runplanType+str(runplanNr)))
        rate_analyses[runplanNr].collection[run].ShowSignalMaps(saveplots=True)

        sel.SetSaveDirectory("Results/{runplan}/FFT/".format(runplan=runplanType+str(runplanNr)))
        rate_analyses[runplanNr].collection[run].ShowFFT(drawoption="mc", cut=True, savePlots=True)

        sel.SetSaveDirectory("Results/{runplan}/Cuts/".format(runplan=runplanType+str(runplanNr)))
        rate_analyses[runplanNr].collection[run].ShowSignalHisto(savePlots=True)
        rate_analyses[runplanNr].collection[run]._ShowHisto(signaldef="abs(median[{channel}])", infoid="Median", savePlots=True, binning=200, xmin=-50, xmax=80)
        rate_analyses[runplanNr].collection[run]._ShowHisto(signaldef="sig_spread", infoid="Spread", savePlots=True, binning=500, xmin=0, xmax=300)
        rate_analyses[runplanNr].collection[run].ShowPeakPosition()

    sel.SetSaveDirectory("Results/{runplan}/".format(runplan=runplanType+str(runplanNr)))
    rate_analyses[runplanNr].ShowSignalVSRate(canvas=rateScanCanvas)
    sel.SavePlots(savename="SignalVSRate"+str(runplanNr)+".root", canvas=rateScanCanvas)
    sel.SavePlots(savename="SignalVSRate"+str(runplanNr)+".png", canvas=rateScanCanvas)


#
#
# sel.UnselectAll()
# sel.SelectRunsInRange(439,459)      # select all runs from a certain run plan
# sel.SetChannels()
#
# print "Selected runs to analyze:"
# sel.ShowSelectedRuns(commentlength=0)
#
# savePath = "Results/new1/"
# sel.SetSaveDirectory(savePath)
#
# rate_analyses[10] = AnalysisCollection(sel)
#
# rate_analyses[10].MakePreAnalysis(savePlot=True)
# rate_analyses[10].ShowPulserRates()
# rate_analyses[10].ShowSignalVSRate(canvas=rateScanCanvas)
# sel.SavePlots(savename="SignalVSRate"+str(10)+".root", canvas=rateScanCanvas)
#
# sel.UnselectAll()
# sel.SelectRunsInRange(463,478)      # select all runs from a certain run plan
# sel.SetChannels()
#
# print "Selected runs to analyze:"
# sel.ShowSelectedRuns(commentlength=0)
#
# savePath = "Results/new2/"
# sel.SetSaveDirectory(savePath)
#
# rate_analyses[11] = AnalysisCollection(sel)
#
# rate_analyses[11].MakePreAnalysis(savePlot=True)
# rate_analyses[11].ShowPulserRates()
# rate_analyses[11].ShowSignalVSRate(canvas=rateScanCanvas)
# sel.SavePlots(savename="SignalVSRate"+str(11)+".root", canvas=rateScanCanvas)