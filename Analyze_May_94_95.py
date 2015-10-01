from AllClasses import *
from ROOT import TCanvas

Analysis.SetTestCampaign(201505)

runplanType     = "rate_scan"
diamond1         = "II-6-94"
diamond2         = "II-6-95"
rate_analyses = {}
mask = {
    "5": "none", # polyA_II-6_94_signal_mask_2015_05_30_0219.dat
    "8": "II-IV-95_II-IV-96_signal_mask_2015_05_31_1512.dat"
}


rateScanCanvas = TCanvas("rateScanCanvas", "rateScanCanvas")


myselection = RunSelection()
runplans = map(str, [5,8])

for runplanNr in runplans:
    # --- change the selected run plan and diamond ---
    myselection.UnselectAll()
    myselection.SelectRunsFromRunPlan(number=runplanNr, type_=runplanType)      # select all runs from a certain run plan
    # activate only diamond with names diamond1 or diamond2:
    myselection.SelectDiamondRuns(diamondname=diamond1, only_selected_runs=True)
    myselection.SelectDiamondRuns(diamondname=diamond2, only_selected_runs=True)

    # show the runs in the terminal:
    print "Selected runs to analyze:"
    myselection.ShowSelectedRuns(commentlength=0)

    # make AnalysisCollection object of selected runs:
    rate_analyses[runplanNr] = AnalysisCollection(myselection, maskfilename=mask[runplanNr])

    # PREANALYSIS plots:
    # myselection.SetSaveDirectory("Results/{runplan}/PreAnalysis/".format(runplan=runplanType+str(runplanNr)))
    # rate_analyses[runplanNr].MakePreAnalysises(savePlot=True)
    # myselection.SetSaveDirectory("Results/{runplan}/PulserRates/".format(runplan=runplanType+str(runplanNr)))
    # rate_analyses[runplanNr].ShowPulserRates()

    # do some stuff for each run individually:
    for run in rate_analyses[runplanNr].collection.keys():
        # myselection.SetSaveDirectory("Results/{runplan}/SignalMaps/".format(runplan=runplanType+str(runplanNr)))
        # rate_analyses[runplanNr].collection[run].ShowSignalMaps(saveplots=True)

        # myselection.SetSaveDirectory("Results/{runplan}/FFT/".format(runplan=runplanType+str(runplanNr)))
        # rate_analyses[runplanNr].collection[run].ShowFFT(drawoption="mc", cut=True, savePlots=True)

        # myselection.SetSaveDirectory("Results/{runplan}/Cuts/".format(runplan=runplanType+str(runplanNr)))
        # rate_analyses[runplanNr].collection[run].ShowSignalHisto(savePlots=True, gridx=True, xmin=-50, xmax=200, binning=250)
        # rate_analyses[runplanNr].collection[run]._ShowHisto(signaldef="abs(median[{channel}])", infoid="Median", drawruninfo=True, logy=True, gridx=True, savePlots=True, binning=200, xmin=-50, xmax=80)
        # rate_analyses[runplanNr].collection[run]._ShowHisto(signaldef="sig_spread[{channel}]", infoid="Spread", drawruninfo=True, savePlots=True, logy=True, gridx=True, binning=150, xmin=0, xmax=150)
        # rate_analyses[runplanNr].collection[run].ShowPeakPosition()

        myselection.SetSaveDirectory("Results/{runplan}/Overview/".format(runplan=runplanType+str(runplanNr)))
        rate_analyses[runplanNr].collection[run]._ShowPreAnalysisOverview(savePlot=True)


    # create Signal_VS_Rate plot:
    myselection.SetSaveDirectory("Results/{runplan}/".format(runplan=runplanType+str(runplanNr)))
    rate_analyses[runplanNr].ShowSignalVSRate(canvas=rateScanCanvas)
    myselection.SavePlots(savename="SignalVSRate"+str(runplanNr)+".root", canvas=rateScanCanvas)
    myselection.SavePlots(savename="SignalVSRate"+str(runplanNr)+".png", canvas=rateScanCanvas)
