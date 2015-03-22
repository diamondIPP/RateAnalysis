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
    assert (len(run_numbers) > 0), "No run specified. Please pass at least one run number as argument."

    print run_numbers

    # config
    show_plots = True

    minimum_statistics = 5 # don't draw bins which contain less than minimum_statistics hits
    run = Run()

    run_numbers =  RunInfo.runs.keys()

    collection = AnalysisCollectionClass()
    for run_number in run_numbers[:20]:

        run.SetRun(run_number)

        if run.current_run['data_type'] == 0: # 0 = data, 1 = pedestrial, 2 = voltagescan, 3 =  long run, 4 = other

            print run_number


            new2DHistConfig = Pad2DHistConfig(50,-0.2,0.15,50,0.03,0.4)
            #
            newAnalysis = Analysis(run,new2DHistConfig)
            newAnalysis.DoAnalysis(minimum_statistics)
            #newAnalysis.CreatePlots(True,'2D_Signal_dist')
            # newAnalysis.CreateMeanSignalHistogram(True)
            newAnalysis.CreateBoth(True)
            collection.AddAnalysis(newAnalysis)


    collection.CreateFWHMPlot()



    if show_plots: raw_input("Press ENTER to quit:")





