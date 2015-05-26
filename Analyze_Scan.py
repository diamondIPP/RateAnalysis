import numpy as np
import ROOT
import sys # to get arguments
import os # for folder and file handling
from AbstractClasses.AnalysisClass import Analysis
from AbstractClasses.AnalysisCollection import AnalysisCollection
from AbstractClasses.RunClass import Run
from AbstractClasses.MCRun import MCRun
from AbstractClasses.RunSelection import RunSelection
from Runinfos.RunInfo import RunInfo
from array import array
from AbstractClasses.ConfigClass import *
from AbstractClasses.MCPerformance import MCPerformance
#from Configuration.initialize_ROOT import initialize_ROOT
from AbstractClasses.Elementary import Elementary


if __name__ == "__main__":

    # config
    Diamond = "IIa-2"
    Bias = -500
    show_plots = True
    minimum_statistics = 100 # don't draw bins which contain less than minimum_statistics hits
    binning = 100

    if 'v' in sys.argv or '-v' in sys.argv or 'V' in sys.argv:
        verbose = True
    else:
        verbose = False

    if '-mc' in sys.argv or 'mc' in sys.argv or 'MC' in sys.argv:
        run = MCRun(validate=False,verbose=verbose)
    else:
        run = Run(validate=False,verbose=verbose)

    # Select Run Numbers:
    selection = RunSelection(verbose=verbose)
    selection.SelectDataType(0)
    selection.UnSelectUnlessDiamond(Diamond)
    selection.UnSelectUnlessBias(Bias)
    selection.ExcludeRun(854)
    run_numbers = selection.GetSelectedRuns()

    print "Starting analysis with ",len(run_numbers)," runs:"
    run.ValidateRuns(run_numbers)

    # SAVE PATH
    savepath = "Results/Scan/"+Diamond+"_"+str(Bias)+"/"

    collection = AnalysisCollection()
    for run_number in run_numbers:

        if run.SetRun(run_number): # run has to exist and is setted

            if True or run.diamond.Specifications['Irradiation'] == 'no':
                print "Starting ",run_number

                #
                newAnalysis = Analysis(run,Config(binning),verbose=verbose)
                newAnalysis.SetSaveDirectory(savepath+str(run_number)+"/")
                newAnalysis.DoAnalysis(minimum_statistics)
                newAnalysis.CreateBoth(saveplots=False)

                newAnalysis.FindMaxima(show=True, binning=binning, minimum_bincontent=minimum_statistics)
                newAnalysis.FindMinima(show=True, binning=binning, minimum_bincontent=minimum_statistics)

                newAnalysis.combined_canvas.cd(1)
                if run.IsMonteCarlo:
                    print "Run is MONTE CARLO"
                    if run.SignalParameters[0] > 0:
                        height = run.SignalParameters[4]
                        newAnalysis.ExtremeAnalysis.Pad.GetSignalInRow(height, show=True)
                    newAnalysis.combined_canvas.cd(1)
                    newAnalysis.ExtremeAnalysis.Pad.MaximaSearch.real_peaks.SetMarkerColor(ROOT.kBlue)
                    newAnalysis.ExtremeAnalysis.Pad.MaximaSearch.real_peaks.SetLineColor(ROOT.kBlue)
                    newAnalysis.ExtremeAnalysis.Pad.MaximaSearch.real_peaks.Draw('SAME P0')
                # newAnalysis.ExtremeAnalysis.Pad.MaximaSearch.found_extrema.SetMarkerColor(ROOT.kGreen+2)
                # newAnalysis.ExtremeAnalysis.Pad.MaximaSearch.found_extrema.Draw('SAME P0')
                minima = newAnalysis.ExtremeAnalysis.ExtremaResults['FoundMinima']
                for i in xrange(len(minima)):
                    text = ROOT.TText()
                    text.SetTextColor(ROOT.kBlue-4)
                    text.DrawText(minima[i][0]-0.01, minima[i][1]-0.005, 'low')
                maxima = newAnalysis.ExtremeAnalysis.ExtremaResults['FoundMaxima']
                for i in xrange(len(maxima)):
                    text = ROOT.TText()
                    text.SetTextColor(ROOT.kRed)
                    text.DrawText(maxima[i][0]-0.02, maxima[i][1]-0.005, 'high')

                ROOT.gStyle.SetPalette(53)
                ROOT.gStyle.SetNumberContours(999)
                newAnalysis.combined_canvas.Update()
                newAnalysis.SavePlots("Combined_Plot.png")


                newAnalysis.GetSignalHeight()
                newAnalysis.ShowTotalSignalHistogram(save=True, showfit=False)
                newAnalysis.SignalTimeEvolution(Mode = "Mean", show = True, time_spacing=3)
                newAnalysis.RateTimeEvolution(show = True, time_spacing=3)


                collection.AddAnalysis(newAnalysis)



        else:
            print "Analysis of run ",run_number, " failed"

    Elementary.SetSaveDirectory(savepath)
    collection.SignalHeightScan()
    collection.PeakComparison()

    canvas = ROOT.gROOT.GetListOfCanvases().FindObject("combined_canvasmaxima")
    if canvas:
        print "combined_canvasmaxima found"
        collection.PeakSignalEvolution(OnThisCanvas=canvas, BinRateEvolution=True)
    else:
        collection.PeakSignalEvolution(BinRateEvolution=True)

    if show_plots: raw_input("Press ENTER to quit:")

    del collection

    # print "\ngDirectory List:\n"
    # ROOT.gDirectory.GetList().ls()
    # print "\nListOfFiles:\n"
    # ROOT.gROOT.GetListOfFiles().ls()
    # print "\nListOfCanvases:\n"
    # ROOT.gROOT.GetListOfCanvases().ls()



# if len(maxima)*len(minima)>0:
#     maxbin = newAnalysis.ExtremeAnalysis.Pad.GetBinByCoordinates(*(maxima[0]))
#     maxbin.FitLandau()
#     minbin = newAnalysis.ExtremeAnalysis.Pad.GetBinByCoordinates(*(minima[0]))
#     minbin.FitLandau()
#     print '\nApproximated Signal Amplitude: {0:0.0f}% - (high/low approximation)\n'.format(100.*(maxbin.Fit['MPV']/minbin.Fit['MPV']-1.))
#
# min_percent = 5
# max_percent = 99
# q = array('d', [1.*min_percent/100., 1.*max_percent/100.])
# y = array('d', [0,0])
# newAnalysis.ExtremeAnalysis.MeanSignalHisto.GetQuantiles(2, y, q)
# print '\nApproximated Signal Amplitude: {0:0.0f}% - ({1:0.0f}%/{2:0.0f}% Quantiles approximation)\n'.format(100.*(y[1]/y[0]-1.), max_percent, min_percent)
# # newAnalysis.ExtremeAnalysis.Pad.MinimaSearch.found_extrema.SetMarkerColor(ROOT.kBlue-4)
# # newAnalysis.ExtremeAnalysis.Pad.MinimaSearch.found_extrema.Draw('SAME P0')