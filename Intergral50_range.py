from ROOT import TFile, gROOT
from AbstractClasses.RunClass import Run
from AbstractClasses.RunSelection import RunSelection


selection = RunSelection()
selection.select_all_runs()
# selection.SelectDiamondRuns("IIa-2")
# selection.UnSelectUnlessDataType(0)
# selection.UnSelectUnlessBias(-500)

runs = selection.get_selected_runs()

run = Run(validate=True)
mins = []
maxs = []
i = 0
for runnummer in runs:

    if run.set_run(runnummer):

        file = TFile(run.TrackingPadAnalysis["ROOTFile"])
        track_info = file.Get("track_info")

        track_info.Draw("-integral50")
        histo = gROOT.FindObject("htemp")
        if not histo:
            print "histo not found"
        else:
            mins.append(histo.GetXaxis().GetXmin())
            maxs.append(histo.GetXaxis().GetXmax())
        i += 1
        print "{0:0.2f}% done.".format(1.*i/len(runs))
    else:
        i += 1
        print "{0:0.2f}% done.".format(1.*i/len(runs))

mins.sort()
maxs.sort()
print "mins:"
print mins
print "maxs:"
print maxs
