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
    run = Run()

    MaxNrOfRuns = 9999
    if(len(run_numbers) == 0):
        MaxNrOfRuns = int(raw_input("No runs specified. Set Maximal number of runs to analyze (0-999) press Enter to select all Runs: "))
        assert (0 <= MaxNrOfRuns < 1000), "Invalid maximal run quantity"
        all_run_numbers =  RunInfo.runs.keys()
        run_numbers = []
        for i in all_run_numbers:
            run.SetRun(i)
            if run.current_run['data_type'] == 0: # 0 = data, 1 = pedestrial, 2 = voltagescan, 3 =  long run, 4 = other
                run_numbers += [i]
    print "Starting analysis with ",len(run_numbers[:MaxNrOfRuns])," runs:"
    print run_numbers[:MaxNrOfRuns]


    collection = AnalysisCollectionClass()
    for run_number in run_numbers[:MaxNrOfRuns]:

        run.SetRun(run_number)

        print run_number

        new2DHistConfig = Pad2DHistConfig(50,-0.2,0.15,50,0.03,0.4)
        #
        newAnalysis = Analysis(run,new2DHistConfig)
        newAnalysis.DoAnalysis(minimum_statistics)
        newAnalysis.CreatePlots(True,'2D_Signal_dist')
        newAnalysis.CreateMeanSignalHistogram(True)
        newAnalysis.CreateBoth(True)
        collection.AddAnalysis(newAnalysis)


    collection.CreateFWHMPlot()



    if show_plots: raw_input("Press ENTER to quit:")





