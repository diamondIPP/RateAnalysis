import numpy as np
import ROOT
import sys # to get arguments
import os # for folder and file handling
from AbstractClasses.AnalysisClass import *
from AbstractClasses.RunClass import Run
from Runinfos.RunInfo import RunInfo
from AbstractClasses.ConfigClass import Pad2DHistConfig
#from Configuration.initialize_ROOT import initialize_ROOT


if __name__ == "__main__":
    #initialize_ROOT()
    # get run numbers to analyze from arguments (type: list):
    run_numbers = map(int, [i for i in sys.argv if i.isdigit() == True])

    # config
    show_plots = True

    minimum_statistics = 5 # don't draw bins which contain less than minimum_statistics hits
    run = Run(False)

    MaxNrOfRuns = 9999
    if(len(run_numbers) == 0):
        MaxNrOfRuns = int(raw_input("No runs specified. Set Maximal number of runs to analyze (0-999) press Enter to select all Runs: "))
        if MaxNrOfRuns == '':
            MaxNrOfRuns = 9999
            print "All runs selected..."
        assert (0 <= MaxNrOfRuns < 10000), "Invalid maximal run quantity"
        all_run_numbers =  RunInfo.runs.keys()
        run_numbers = []
        for i in all_run_numbers:
            run.SetRun(i)
            if run.current_run['data_type'] == 0: # 0 = data, 1 = pedestrial, 2 = voltagescan, 3 =  long run, 4 = other
                run_numbers += [i]
        run_numbers = run_numbers[:MaxNrOfRuns]
    print "Starting analysis with ",len(run_numbers)," runs:"
    run.ValidateRuns(run_numbers)
    print run_numbers


    collection = AnalysisCollectionClass()
    for run_number in run_numbers:

        if run.SetRun(run_number): # run has to exist and is setted

            if True or run.diamond.Specifications['Irradiation'] == 'no':
                print run_number

                #
                newAnalysis = Analysis(run)
                newAnalysis.DoAnalysis(minimum_statistics)
                #newAnalysis.CreatePlots(True,'2D_Signal_dist')
                #newAnalysis.CreateMeanSignalHistogram(True)
                newAnalysis.CreateBoth(True)
                #print newAnalysis.Pad.ListOfBins[newAnalysis.Pad.GetBinNumber(0.05,0.05)].GetEntries()
                newAnalysis.Pad.ShowBinXYSignalHisto(0.08,0.3)
                newAnalysis.Pad.ShowBinXYSignalHisto(-0.04,0.3)
                newAnalysis.CreateHitsDistribution()

                collection.AddAnalysis(newAnalysis)
        else:
            print "Analysis of run ",run_number, " failed"


    collection.CreateFWHMPlot()



    if show_plots: raw_input("Press ENTER to quit:")





