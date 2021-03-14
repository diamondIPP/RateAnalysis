import ROOT
import numpy as np
from ROOT import TFile, TGraph, TGraphErrors, TF1, TCanvas

'''
To improve:
    add import of MCPerformanceAnalysis configuration
    include in framework ?
'''

# Settings:
hits = 300000
tries = 20
peaks = "Test"
binning = 200
min_bincontent = 50
pathExtension = ""

PerfResults = "MC/Performance_Results/"
# foldername = "2015-05-04 16-24-09.965268/"
foldername = "_" + str(min_bincontent) + "_" + str(binning) + "_" + str(hits) + "_" + pathExtension + "/"
filename = "MCPerformanceLog.root"
filepath = PerfResults + foldername + filename

file = TFile(filepath)

LogTree = file.Get("LogTree")

success_graph = TGraphErrors()
ghost_graph = TGraph()
minimas_graph = TGraph()
RecSA_MinMax_graph = TGraphErrors()
RecSA_Quantiles_graph = TGraphErrors()

success = []
real_amplitude = []
All_RecSA_MinMax = []
All_RecSA_Quantiles = []
tmp_success = np.zeros(tries)
tmp_RecSA_Q = np.zeros(tries)
tmp_RecSA_M = np.zeros(tries)
tmp_Ghosts = np.zeros(tries)
tmp_Minimas = np.zeros(tries)

count = 0
for i in range(LogTree.GetEntries()):
    # read the ROOT TTree
    LogTree.GetEntry(i)
    TrueNPeaks = LogTree.TrueNPeaks
    Ninjas = LogTree.Ninjas
    RealSignalAmplitude = LogTree.RealSignalAmplitude
    RecSA_Quantiles = LogTree.RecSA_Quantiles
    RecSA_MinMax = LogTree.RecSA_MinMax
    Ghosts = LogTree.Ghosts
    Minimas = LogTree.Minimas
    Repetition = LogTree.Repetition

    tmp_success[Repetition] = 1. * (TrueNPeaks - Ninjas) / TrueNPeaks
    tmp_RecSA_M[Repetition] = RecSA_MinMax
    tmp_RecSA_Q[Repetition] = RecSA_Quantiles
    tmp_Ghosts[Repetition] = Ghosts
    tmp_Minimas[Repetition] = Minimas

    if Repetition == tries - 1:
        mean_success = tmp_success.mean()
        mean_RecSA_Q = tmp_RecSA_Q.mean()
        mean_RecSA_M = tmp_RecSA_M.mean()
        # success.append(mean_success)
        # real_amplitude.append(RealSignalAmplitude)
        success_graph.SetPoint(count, RealSignalAmplitude, mean_success)
        success_graph.SetPointError(count, 0, tmp_success.std() / np.sqrt(tries))

        RecSA_MinMax_graph.SetPoint(count, RealSignalAmplitude, mean_RecSA_M)
        RecSA_MinMax_graph.SetPointError(count, 0, tmp_RecSA_M.std())

        RecSA_Quantiles_graph.SetPoint(count, RealSignalAmplitude, mean_RecSA_Q)
        RecSA_Quantiles_graph.SetPointError(count, 0, tmp_RecSA_Q.std())

        ghost_graph.SetPoint(count, RealSignalAmplitude, tmp_Ghosts.mean())
        minimas_graph.SetPoint(count, RealSignalAmplitude, tmp_Minimas.mean())

        count += 1

canvas = TCanvas("canvas", "canvas")
canvas.SetGrid()

pad = canvas.GetPad(0)
success_graph.SetNameTitle("success", "MC Performance Analysis Result ({0} Hits)".format(hits))
success_graph.GetXaxis().SetTitle("Relative Real Signal Amplitude")
success_graph.GetYaxis().SetTitle("Peak Finding Efficiency | Reconstructed Signal Amplitude")
success_graph.GetYaxis().SetRangeUser(0, 1.1)
success_graph.Draw("ALP*")

ghost_graph.SetNameTitle("ghost_graph", "Ghost Peaks")
ghost_graph.Draw("SAME LP*")

func = TF1("func", "x", 0, 1)
func.Draw("SAME")
RecSA_Quantiles_graph.SetMarkerColor(ROOT.kRed)
RecSA_Quantiles_graph.Draw("SAME P*")
PngSaveName = PerfResults + foldername + ("MC_PerformanceAnalysis_{0}_" + peaks + "_{1}rep_" + str(binning) + "_" + str(min_bincontent) + ".png").format(hits, tries)
RootSaveName = PerfResults + foldername + ("MC_PerformanceAnalysis_{0}_" + peaks + "_{1}rep_" + str(binning) + "_" + str(min_bincontent) + ".root").format(hits, tries)
ROOT.gPad.Print(PngSaveName)
ROOT.gPad.Print(RootSaveName)
# RecSA_MinMax_graph.SetMarkerColor(ROOT.kBlue)
# RecSA_MinMax_graph.Draw("SAME P*")

ghostcanvas = TCanvas("ghostcanvas", "ghostcanvas")
ghostcanvas.cd()
ghost_graph.GetXaxis().SetTitle("Relative Real Signal Amplitude")
ghost_graph.GetYaxis().SetTitle("NGhostPeaks/NRepetitions")
ghost_graph.Draw("ALP*")
ghostcanvas.Update()
ROOT.gPad.Print(PerfResults + foldername + "ghosts.png")

minimascanvas = TCanvas("minimascanvas", "minimascanvas")
minimascanvas.cd()
minimas_graph.SetNameTitle("minimas_graph", "Minimas found")
minimas_graph.GetXaxis().SetTitle("Relative Real Signal Amplitude")
minimas_graph.GetYaxis().SetTitle("NMinimas/NRepetitions")
minimas_graph.Draw("ALP*")
minimascanvas.Update()
ROOT.gPad.Print(PerfResults + foldername + "minimas.png")

input("finish")
