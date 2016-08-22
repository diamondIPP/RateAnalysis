from ROOT import TFile, TCanvas, TProfile, TH2D, kBlue, kRed, TH1D, gStyle, TF1, kWhite, gROOT, kBlack, kOrange, TGraphErrors, TGraph, kGreen, gPad, kError
from optparse import OptionParser
from ConfigParser import ConfigParser
from array import array
from numpy import mean, std
from math import sqrt
import os
import json
from glob import glob

__author__ = "Diego Alejandro"

def plotAll(s1, s2, s3, s4, cluster=False):
    if not cluster:
        s1.plot_rate_results('-600V_irradiated',kRed+1, False,3)
        s2.plot_rate_results('-1000V_irradiated',kOrange-1, False,3)
        s3.plot_rate_results('-600V_unirradiated',kBlue+1, False,3)
        s4.plot_rate_results('-1000V_unirradiated',kGreen+3, False,3)
        maximum = max(s1.graph_cuts.GetMaximum(),s2.graph_cuts.GetMaximum(),s3.graph_cuts.GetMaximum(),s4.graph_cuts.GetMaximum())
        s1.graph_cuts.SetMaximum(maximum)
        s2.graph_cuts.SetMaximum(maximum)
        s3.graph_cuts.SetMaximum(maximum)
        s4.graph_cuts.SetMaximum(maximum)
        canvas = TCanvas('canvas', 'canvas', 2100, 1500)
        canvas.SetLeftMargin(0.1)
        canvas.SetRightMargin(0.1)
        canvas.SetLogx(1)
        canvas.cd()
        s4.graph_cuts.Draw('apl')
        s3.graph_cuts.Draw('pl')
        s2.graph_cuts.Draw('pl')
        s1.graph_cuts.Draw('pl')
        canvas.BuildLegend(0.65, 0.7, 0.9, 0.9)
        s1.startPoint_cuts.Draw('p')
        s1.finishPoint_cuts.Draw('p')
        s2.startPoint_cuts.Draw('p')
        s2.finishPoint_cuts.Draw('p')
        s3.startPoint_cuts.Draw('p')
        s3.finishPoint_cuts.Draw('p')
        s4.startPoint_cuts.Draw('p')
        s4.finishPoint_cuts.Draw('p')
    
        canvas.SetTitle('canvas')
        canvas.SetName('canvas')
        canvas.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_Cl_sizes.root'.format(dir='.'))
        canvas.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_Cl_sizes.png'.format(dir='.'))
    else:
        s1.SetName('600V_irradiated')
        s2.SetName('1000V_irradiated')
        s3.SetName('600V_unirradiated')
        s4.SetName('1000V_unirradiated')
        s1.SetTitle('600V_irradiated')
        s2.SetTitle('1000V_irradiated')
        s3.SetTitle('600V_unirradiated')
        s4.SetTitle('1000V_unirradiated')
        maximum = max(s1.GetMaximum(),s2.GetMaximum(),s3.GetMaximum(),s4.GetMaximum())
        s1.SetMaximum(maximum)
        s2.SetMaximum(maximum)
        s3.SetMaximum(maximum)
        s4.SetMaximum(maximum)
        canvas = TCanvas('canvas', 'canvas', 2100, 1500)
        canvas.SetLeftMargin(0.1)
        canvas.SetRightMargin(0.1)
        canvas.SetLogx(1)
        canvas.cd()
        s4.Draw('apl')
        s3.Draw('pl')
        s2.Draw('pl')
        s1.Draw('pl')
        canvas.BuildLegend(0.65, 0.7, 0.9, 0.9)
        canvas.SaveAs('{dir}/Cluster_Vs_Flux_with_cuts.root'.format(dir='.'))
        canvas.SaveAs('{dir}/Cluster_Vs_Flux_with_cuts.png'.format(dir='.'))

def doClusterAll(d4c1_1, d4c2_1, d4c3_1, d4c4_1, d5c1_1, d5c2_1, d5c3_1, d5c4_1, d4c1_2, d4c2_2, d4c3_2, d4c4_2, d5c1_2, d5c2_2, d5c3_2, d5c4_2):
    graphClus4_1 = TGraphErrors(7)
    graphClus4_1.SetNameTitle('clusters_4_1', 'clusters_4_1')
    RateAnalysisExtraction.set_1D_options(d4c1_1, graphClus4_1, 'clusters_4_1', 'Flux (kHz/cm^{2})', 'mean Pix per cluster',kBlue+1, 20, 3, 0, 2.5)

    graphClus5_1 = TGraphErrors(7)
    graphClus5_1.SetNameTitle('clusters_5_1', 'clusters_5_1')
    RateAnalysisExtraction.set_1D_options(d5c1_1, graphClus5_1, 'clusters_5_1', 'Flux (kHz/cm^{2})', 'mean Pix per cluster',kOrange-1, 20, 3, 0, 2.5)

    graphClus4_2 = TGraphErrors(7)
    graphClus4_2.SetNameTitle('clusters_4_2', 'clusters_4_2')
    RateAnalysisExtraction.set_1D_options(d4c1_2, graphClus4_2, 'clusters_4_2', 'Flux (kHz/cm^{2})', 'mean Pix per cluster',kRed+1, 20, 3, 0, 2.5)

    graphClus5_2 = TGraphErrors(7)
    graphClus5_2.SetNameTitle('clusters_5_2', 'clusters_5_2')
    RateAnalysisExtraction.set_1D_options(d5c1_2, graphClus5_2, 'clusters_5_2', 'Flux (kHz/cm^{2})', 'mean Pix per cluster',kGreen+3, 20, 3, 0, 2.5)

    doClusterOne(d4c1_1, d4c2_1, d4c3_1, d4c4_1, graphClus4_1, 1)
    doClusterOne(d5c1_1, d5c2_1, d5c3_1, d5c4_1, graphClus5_1, 1)
    doClusterOne(d4c1_2, d4c2_2, d4c3_2, d4c4_2, graphClus4_2, 2)
    doClusterOne(d5c1_2, d5c2_2, d5c3_2, d5c4_2, graphClus5_2, 2)
    plotAll(graphClus4_1, graphClus4_2, graphClus5_1, graphClus5_2, True)

def doClusterOne(d4c1_1, d4c2_1, d4c3_1, d4c4_1, graphClus4_1, type=1):
    histos = {'2k': TH1D('2k', '2k', 4, 0.5, 4.5), '20k': TH1D('20k', '20k', 4, 0.5, 4.5), '80k': TH1D('80k', '80k', 4, 0.5, 4.5), '200k': TH1D('200k', '200k', 4, 0.5, 4.5), '600k': TH1D('600k', '600k', 4, 0.5, 4.5), '2M': TH1D('2M', '2M', 4, 0.5, 4.5), '5M': TH1D('5M', '5M', 4, 0.5, 4.5)}
    for run in d4c1_1.runs:
        if run in d4c1_1.run2k:
            temp = histos['2k'].GetBinContent(1)
            histos['2k'].SetBinContent(1, temp + d4c1_1.projs_cuts[run].GetEntries())
        elif run in d4c1_1.run20k:
            temp = histos['20k'].GetBinContent(1)
            histos['20k'].SetBinContent(1, temp + d4c1_1.projs_cuts[run].GetEntries())
        elif run in d4c1_1.run80k:
            temp = histos['80k'].GetBinContent(1)
            histos['80k'].SetBinContent(1, temp + d4c1_1.projs_cuts[run].GetEntries())
        elif run in d4c1_1.run200k:
            temp = histos['200k'].GetBinContent(1)
            histos['200k'].SetBinContent(1, temp + d4c1_1.projs_cuts[run].GetEntries())
        elif run in d4c1_1.run600k:
            temp = histos['600k'].GetBinContent(1)
            histos['600k'].SetBinContent(1, temp + d4c1_1.projs_cuts[run].GetEntries())
        elif run in d4c1_1.run20k:
            temp = histos['2M'].GetBinContent(1)
            histos['2M'].SetBinContent(1, temp + d4c1_1.projs_cuts[run].GetEntries())
        else:
            temp = histos['5M'].GetBinContent(1)
            histos['5M'].SetBinContent(1, temp + d4c1_1.projs_cuts[run].GetEntries())
    for run in d4c2_1.runs:
        if run in d4c2_1.run2k:
            temp = histos['2k'].GetBinContent(2)
            histos['2k'].SetBinContent(2, temp + d4c2_1.projs_cuts[run].GetEntries())
        elif run in d4c2_1.run20k:
            temp = histos['20k'].GetBinContent(2)
            histos['20k'].SetBinContent(2, temp + d4c2_1.projs_cuts[run].GetEntries())
        elif run in d4c2_1.run80k:
            temp = histos['80k'].GetBinContent(2)
            histos['80k'].SetBinContent(2, temp + d4c2_1.projs_cuts[run].GetEntries())
        elif run in d4c2_1.run200k:
            temp = histos['200k'].GetBinContent(2)
            histos['200k'].SetBinContent(2, temp + d4c2_1.projs_cuts[run].GetEntries())
        elif run in d4c2_1.run600k:
            temp = histos['600k'].GetBinContent(2)
            histos['600k'].SetBinContent(2, temp + d4c2_1.projs_cuts[run].GetEntries())
        elif run in d4c2_1.run20k:
            temp = histos['2M'].GetBinContent(2)
            histos['2M'].SetBinContent(2, temp + d4c2_1.projs_cuts[run].GetEntries())
        else:
            temp = histos['5M'].GetBinContent(2)
            histos['5M'].SetBinContent(2, temp + d4c2_1.projs_cuts[run].GetEntries())
    for run in d4c3_1.runs:
        if run in d4c3_1.run2k:
            temp = histos['2k'].GetBinContent(3)
            histos['2k'].SetBinContent(3, temp + d4c3_1.projs_cuts[run].GetEntries())
        elif run in d4c3_1.run20k:
            temp = histos['20k'].GetBinContent(3)
            histos['20k'].SetBinContent(3, temp + d4c3_1.projs_cuts[run].GetEntries())
        elif run in d4c3_1.run80k:
            temp = histos['80k'].GetBinContent(3)
            histos['80k'].SetBinContent(3, temp + d4c3_1.projs_cuts[run].GetEntries())
        elif run in d4c3_1.run200k:
            temp = histos['200k'].GetBinContent(3)
            histos['200k'].SetBinContent(3, temp + d4c3_1.projs_cuts[run].GetEntries())
        elif run in d4c3_1.run600k:
            temp = histos['600k'].GetBinContent(3)
            histos['600k'].SetBinContent(3, temp + d4c3_1.projs_cuts[run].GetEntries())
        elif run in d4c3_1.run20k:
            temp = histos['2M'].GetBinContent(3)
            histos['2M'].SetBinContent(3, temp + d4c3_1.projs_cuts[run].GetEntries())
        else:
            temp = histos['5M'].GetBinContent(3)
            histos['5M'].SetBinContent(3, temp + d4c3_1.projs_cuts[run].GetEntries())
    for run in d4c4_1.runs:
        if run in d4c4_1.run2k:
            temp = histos['2k'].GetBinContent(4)
            histos['2k'].SetBinContent(4, temp + d4c4_1.projs_cuts[run].GetEntries())
        elif run in d4c4_1.run20k:
            temp = histos['20k'].GetBinContent(4)
            histos['20k'].SetBinContent(4, temp + d4c4_1.projs_cuts[run].GetEntries())
        elif run in d4c4_1.run80k:
            temp = histos['80k'].GetBinContent(4)
            histos['80k'].SetBinContent(4, temp + d4c4_1.projs_cuts[run].GetEntries())
        elif run in d4c4_1.run200k:
            temp = histos['200k'].GetBinContent(4)
            histos['200k'].SetBinContent(4, temp + d4c4_1.projs_cuts[run].GetEntries())
        elif run in d4c4_1.run600k:
            temp = histos['600k'].GetBinContent(4)
            histos['600k'].SetBinContent(4, temp + d4c4_1.projs_cuts[run].GetEntries())
        elif run in d4c4_1.run20k:
            temp = histos['2M'].GetBinContent(4)
            histos['2M'].SetBinContent(4, temp + d4c4_1.projs_cuts[run].GetEntries())
        else:
            temp = histos['5M'].GetBinContent(4)
            histos['5M'].SetBinContent(4, temp + d4c4_1.projs_cuts[run].GetEntries())

    graphClus4_1.SetPoint(0, d4c1_1.mean_fluxes['2k'], histos['2k'].GetMean())
    graphClus4_1.SetPointError(0, 0.1*d4c1_1.mean_fluxes['2k'], histos['2k'].GetStdDev())
    graphClus4_1.SetPoint(1, d4c1_1.mean_fluxes['20k'], histos['20k'].GetMean())
    graphClus4_1.SetPointError(1, 0.1*d4c1_1.mean_fluxes['20k'], histos['20k'].GetStdDev())
    graphClus4_1.SetPoint(2, d4c1_1.mean_fluxes['80k'], histos['80k'].GetMean())
    graphClus4_1.SetPointError(2, d4c1_1.mean_fluxes['80k']*0.1, histos['80k'].GetStdDev())
    if type is 2:
        graphClus4_1.SetPoint(3, d4c1_1.mean_fluxes['200k'], histos['200k'].GetMean())
        graphClus4_1.SetPointError(3, 0.1*d4c1_1.mean_fluxes['200k'], histos['200k'].GetStdDev())
        graphClus4_1.SetPoint(4, d4c1_1.mean_fluxes['600k'], histos['600k'].GetMean())
        graphClus4_1.SetPointError(4, 0.1*d4c1_1.mean_fluxes['600k'], histos['600k'].GetStdDev())
        graphClus4_1.SetPoint(5, d4c1_1.mean_fluxes['2M'], histos['2M'].GetMean())
        graphClus4_1.SetPointError(5, d4c1_1.mean_fluxes['2M']*0.1, histos['2M'].GetStdDev())
        graphClus4_1.SetPoint(6, d4c1_1.mean_fluxes['5M'], histos['5M'].GetMean())
        graphClus4_1.SetPointError(6, 0.1*d4c1_1.mean_fluxes['5M'], histos['5M'].GetStdDev())
    else:
        graphClus4_1.SetPoint(3, d4c1_1.mean_fluxes['600k'], histos['600k'].GetMean())
        graphClus4_1.SetPointError(3, 0.1*d4c1_1.mean_fluxes['600k'], histos['600k'].GetStdDev())
        graphClus4_1.SetPoint(4, d4c1_1.mean_fluxes['2M'], histos['2M'].GetMean())
        graphClus4_1.SetPointError(4, d4c1_1.mean_fluxes['2M']*0.1, histos['2M'].GetStdDev())
        graphClus4_1.SetPoint(5, d4c1_1.mean_fluxes['5M'], histos['5M'].GetMean())
        graphClus4_1.SetPointError(5, 0.1*d4c1_1.mean_fluxes['5M'], histos['5M'].GetStdDev())


class RateAnalysisExtraction:
    def __init__(self, ini, fin, diamond=4, campaign='1510', exclude='', size=0):
        self.ini = ini
        self.fin = fin
        self.dia = diamond
        self.campaign = campaign
        self.excl = exclude.replace(',', '_')
        self.exclude = exclude.split(',')
        self.cl_size = size;
        self.dirs = self.get_list_dirs()
        self.runs = self.get_runs()
        self.runs = sorted(self.runs)
        self.runs = self.runs[self.runs.index(self.ini):self.runs.index(self.fin)+1]
        self.dir_roots = 'Ini_{i}_Fin_{f}_Exc_{e}/Root'.format(i=self.ini, f=self.fin, e=self.excl)
        self.dir_plots = 'Ini_{i}_Fin_{f}_Exc_{e}/Plots'.format(i=self.ini, f=self.fin, e=self.excl)
        if not os.path.isdir(self.dir_roots):
            os.makedirs(self.dir_roots)
        if not os.path.isdir(self.dir_plots):
            os.makedirs(self.dir_plots)
        if self.exclude[0] is not '':
            for excl in self.exclude:
                if int(excl) in self.runs:
                    self.runs.remove(int(excl))
                    self.dirs.remove('{c}_{e}'.format(c=self.campaign, e=excl))
        self.dirs = self.get_dirs_dict()
        if self.dia is 4:
            self.projs = {run: self.get_histo(run, 'c_phROC4_all') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phROC4_1cl') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phROC4_2cl') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phROC4_3cl') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phROC4_M4cl') for run in self.runs}
            self.projs_cuts = {run: self.get_histo(run, 'c_phROC4_all_cuts') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phROC4_1cl_cuts') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phROC4_2cl_cuts') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phROC4_3cl_cuts') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phROC4_M4cl_cuts') for run in self.runs}
            self.projs_syst = {run: self.get_histo(run, 'c_phROC6_all') for run in self.runs}
            self.projs_syst_cuts = {run: self.get_histo(run, 'c_phROC6_all_cuts') for run in self.runs}
        elif self.dia is 5:
            self.projs = {run: self.get_histo(run, 'c_phROC5_all') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phROC5_1cl') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phROC5_2cl') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phROC5_3cl') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phROC5_M4cl') for run in self.runs}
            self.projs_cuts = {run: self.get_histo(run, 'c_phROC5_all_cuts') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phROC5_1cl_cuts') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phROC5_2cl_cuts') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phROC5_3cl_cuts') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phROC5_M4cl_cuts') for run in self.runs}
            self.projs_syst = {run: self.get_histo(run, 'c_phROC6_all') for run in self.runs}
            self.projs_syst_cuts = {run: self.get_histo(run, 'c_phROC6_all_cuts') for run in self.runs}
        elif self.dia is 6:
            self.projs = {run: self.get_histo(run, 'c_phROC6_all') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phROC6_1cl') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phROC6_2cl') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phROC6_3cl') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phROC6_M4cl') for run in self.runs}
            self.projs_cuts = {run: self.get_histo(run, 'c_phROC6_all_cuts') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phROC6_1cl_cuts') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phROC6_2cl_cuts') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phROC6_3cl_cuts') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phROC6_M4cl_cuts') for run in self.runs}

        # self.projs = {run: self.do_projection(run, False) for run in self.runs}
        # self.projs_cuts = {run: self.do_projection(run, True) for run in self.runs}
        self.means = {run: self.projs[run].GetMean() for run in self.runs}
        self.means_cuts = {run: self.projs_cuts[run].GetMean() for run in self.runs}
        self.sigmas = {run: self.projs[run].GetMeanError() for run in self.runs}  # Micha said GetMeanError(). Trying GetStdDev()
        self.sigmas_cuts = {run: self.projs_cuts[run].GetMeanError() for run in self.runs}  # Micha said GetMeanError(). Trying GetStdDev()
        self.mpvs = {run: self.projs[run].GetBinCenter(self.projs[run].GetMaximumBin()) for run in self.runs}
        self.mpvs_cuts = {run: self.projs_cuts[run].GetBinCenter(self.projs_cuts[run].GetMaximumBin()) for run in self.runs}
        self.runinfo = self.get_runinfo()
        self.run2k, self.run20k, self.run80k, self.run200k, self.run600k, self.run2M, self.run5M = [], [], [], [], [], [], []
        self.fluxes = {run: self.calc_flux(run) for run in self.runs}
        self.systematics = {'2k': [], '20k': [], '80k': [], '200k': [], '600k': [], '2M': [], '5M': []}
        self.systematics_cuts = {'2k': [], '20k': [], '80k': [], '200k': [], '600k': [], '2M': [], '5M': []}
        if self.dia is not 6:
            self.get_systematics()
            for run in self.runs:
                self.sigmas[run] = sqrt(self.sigmas[run]**2 + std(self.systematics['2k'])**2) if run in self.run2k else sqrt(self.sigmas[run]**2 + std(self.systematics['20k'])**2) if run in self.run20k else sqrt(self.sigmas[run]**2 + std(self.systematics['80k'])**2) if run in self.run80k else sqrt(self.sigmas[run]**2 + std(self.systematics['200k'])**2) if run in self.run200k else sqrt(self.sigmas[run]**2 + std(self.systematics['600k'])**2) if run in self.run600k else sqrt(self.sigmas[run]**2 + std(self.systematics['2M'])**2) if run in self.run2M else sqrt(self.sigmas[run]**2 + std(self.systematics['5M'])**2) 
                self.sigmas_cuts[run] = sqrt(self.sigmas_cuts[run]**2 + std(self.systematics_cuts['2k'])**2) if run in self.run2k else sqrt(self.sigmas_cuts[run]**2 + std(self.systematics_cuts['20k'])**2) if run in self.run20k else sqrt(self.sigmas_cuts[run]**2 + std(self.systematics_cuts['80k'])**2) if run in self.run80k else sqrt(self.sigmas_cuts[run]**2 + std(self.systematics_cuts['200k'])**2) if run in self.run200k else sqrt(self.sigmas_cuts[run]**2 + std(self.systematics_cuts['600k'])**2) if run in self.run600k else sqrt(self.sigmas_cuts[run]**2 + std(self.systematics_cuts['2M'])**2) if run in self.run2M else sqrt(self.sigmas_cuts[run]**2 + std(self.systematics_cuts['5M'])**2)
        self.mean_fluxes = {'2k': self.calculate_mean_fluxees(self.run2k), '20k': self.calculate_mean_fluxees(self.run20k), '80k': self.calculate_mean_fluxees(self.run80k), '200k': self.calculate_mean_fluxees(self.run200k), '600k': self.calculate_mean_fluxees(self.run600k), '2M': self.calculate_mean_fluxees(self.run2M), '5M': self.calculate_mean_fluxees(self.run5M)}

        # self.clusters = {run:  0 for run in self.runs}
        # self.clusters_cuts = {run: 0 for run in self.runs}
        # if self.dia is not 6:
        #     self.get_cluster_analysis()
    def calculate_mean_fluxees(self, runrate):
        temp = 0
        for run in runrate:
            temp += self.fluxes[run]/float(len(runrate))
        return temp

    def get_list_dirs(self):
        return glob('{c}*'.format(c=self.campaign))

    def get_runs(self):
        return [int(dir.split('_')[-1]) for dir in self.dirs]

    def get_dirs_dict(self):
        return {run: '{camp}_{r}'.format(camp=self.campaign, r=run) for run in self.runs}

    def get_histo(self, run, name):
        histo = (TFile('{c}_{r}/Root/{name}.root'.format(c= self.campaign, r=run, name=name), 'READ').Get(name)).GetPrimitive(name.replace('c_', ''))
        name_old = histo.GetName()
        histo.SetName('{name}_{r}'.format(name=name_old, r=run))
        return histo

    def do_projection(self, run, with_cuts=True):
        name = self.histos_cuts[run].GetName() if with_cuts else self.histos[run].GetName()
        if with_cuts:
            return self.histos_cuts[run].ProjectionY('p_{n}'.format(n=name))
        else:
            return self.histos[run].ProjectionY('p_{n}'.format(n=name))
        
    # def get_cluster_analysis(self):
    #     for run in self.runs:
    #         self.clusters[run] += self.projs[run].GetEntries()
    #         self.clusters_cuts[run] += self.projs_cuts[run].GetEntries()

    def get_systematics(self):
        for run in self.runs:
            if run in self.run2k:
                self.systematics['2k'].append(self.projs_syst[run].GetMean())
                self.systematics_cuts['2k'].append(self.projs_syst_cuts[run].GetMean())
            elif run in self.run20k:
                self.systematics['20k'].append(self.projs_syst[run].GetMean())
                self.systematics_cuts['20k'].append(self.projs_syst_cuts[run].GetMean())
            elif run in self.run80k:
                self.systematics['80k'].append(self.projs_syst[run].GetMean())
                self.systematics_cuts['80k'].append(self.projs_syst_cuts[run].GetMean())
            elif run in self.run200k:
                self.systematics['200k'].append(self.projs_syst[run].GetMean())
                self.systematics_cuts['200k'].append(self.projs_syst_cuts[run].GetMean())
            elif run in self.run600k:
                self.systematics['600k'].append(self.projs_syst[run].GetMean())
                self.systematics_cuts['600k'].append(self.projs_syst_cuts[run].GetMean())
            elif run in self.run2M:
                self.systematics['2M'].append(self.projs_syst[run].GetMean())
                self.systematics_cuts['2M'].append(self.projs_syst_cuts[run].GetMean())
            else:
                self.systematics['5M'].append(self.projs_syst[run].GetMean())
                self.systematics_cuts['5M'].append(self.projs_syst_cuts[run].GetMean())

    def plot_rate_results(self, name='PH_Flux_Scan_with_cuts', color=kBlue, paint=True, scale=2):
        self.graph = TGraphErrors(len(self.means))
        self.graph.SetNameTitle('fluxes', 'fluxes')
        self.graph_cuts = TGraphErrors(len(self.means_cuts))
        self.graph_cuts.SetNameTitle(name, name)
        for i in xrange(len(self.means)):
            self.graph.SetPoint(i, self.fluxes[self.runs[i]], self.means[self.runs[i]])
            self.graph.SetPointError(i, self.fluxes[self.runs[i]] * 0.1, self.sigmas[self.runs[i]])
        for i in xrange(len(self.means_cuts)):
            self.graph_cuts.SetPoint(i, self.fluxes[self.runs[i]], self.means_cuts[self.runs[i]])
            self.graph_cuts.SetPointError(i, self.fluxes[self.runs[i]] * 0.1, self.sigmas_cuts[self.runs[i]])
        self.set_1D_options(self.graph, 'PH_Flux_Scan', 'Flux (kHz/cm^{2})', 'mean Charge (a.u.)',color, 20, scale, 0, self.graph.GetMaximum())
        self.set_1D_options(self.graph_cuts, name, 'Flux (kHz/cm^{2})', 'mean Charge (a.u.)',color, 20, scale, 0, self.graph_cuts.GetMaximum())
        self.CreateCanvas('Charge_Vs_Flux', scale)
        self.CreateCanvas('Charge_Vs_Flux_with_cuts', scale)
        self.startPoint = TGraph(1)
        self.startPoint_cuts = TGraph(1)
        self.startPoint.SetName('Run Start w/o cuts')
        self.startPoint_cuts.SetName('Run Start')
        self.startPoint.SetTitle('Run Start w/o cuts')
        self.startPoint_cuts.SetTitle('Run Start')
        self.finishPoint = TGraph(1)
        self.finishPoint_cuts = TGraph(1)
        self.finishPoint.SetName('Run Finish w/o cuts')
        self.finishPoint_cuts.SetName('Run Finish')
        self.finishPoint.SetTitle('Run Finish w/o cuts')
        self.finishPoint_cuts.SetTitle('Run Finish')
        self.startPoint.SetPoint(0, self.fluxes[self.runs[0]], self.means[self.runs[0]])
        self.startPoint_cuts.SetPoint(0, self.fluxes[self.runs[0]], self.means_cuts[self.runs[0]])
        self.finishPoint.SetPoint(0, self.fluxes[self.runs[-1]], self.means[self.runs[-1]])
        self.finishPoint_cuts.SetPoint(0, self.fluxes[self.runs[-1]], self.means_cuts[self.runs[-1]])
        self.startPoint.SetMarkerStyle(22)
        self.startPoint_cuts.SetMarkerStyle(22)
        self.startPoint.SetMarkerColor(kGreen)
        self.startPoint_cuts.SetMarkerColor(kGreen)
        self.startPoint.SetMarkerSize(gStyle.GetMarkerSize()*1.7)
        self.startPoint_cuts.SetMarkerSize(gStyle.GetMarkerSize()*1.7)
        self.finishPoint.SetMarkerStyle(23)
        self.finishPoint_cuts.SetMarkerStyle(23)
        self.finishPoint.SetMarkerColor(kRed)
        self.finishPoint_cuts.SetMarkerColor(kRed)
        self.finishPoint.SetMarkerSize(gStyle.GetMarkerSize()*1.7)
        self.finishPoint_cuts.SetMarkerSize(gStyle.GetMarkerSize()*1.7)
        self.graph.SetFillColor(kWhite)
        self.graph_cuts.SetFillColor(kWhite)
        self.startPoint.SetFillColor(kWhite)
        self.startPoint_cuts.SetFillColor(kWhite)
        self.startPoint.SetLineColor(kWhite)
        self.startPoint_cuts.SetLineColor(kWhite)
        self.startPoint.SetLineColorAlpha(kWhite, 1)
        self.startPoint_cuts.SetLineColorAlpha(kWhite, 1)
        self.finishPoint.SetFillColor(kWhite)
        self.finishPoint_cuts.SetFillColor(kWhite)
        self.finishPoint.SetLineColor(kWhite)
        self.finishPoint_cuts.SetLineColor(kWhite)
        self.finishPoint.SetLineColorAlpha(kWhite, 1)
        self.finishPoint_cuts.SetLineColorAlpha(kWhite, 1)
        if paint:
            self.cCharge_Vs_Flux.cd()
            self.graph.Draw('APL')
            self.startPoint.Draw('P')
            self.finishPoint.Draw('P')
            self.cCharge_Vs_Flux.BuildLegend(0.65, 0.7, 0.9, 0.9)
            if self.cl_size is 0:
                self.cCharge_Vs_Flux.SaveAs('{dir}/Charge_Vs_Flux_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(dir=self.dir_roots, d=self.dia, i=self.ini, f=self.fin, e=self.excl))
                self.cCharge_Vs_Flux.SaveAs('{dir}/Charge_Vs_Flux_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(dir=self.dir_plots, d=self.dia, i=self.ini, f=self.fin, e=self.excl))
            elif self.cl_size < 4:
                self.cCharge_Vs_Flux.SaveAs('{dir}/Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(dir=self.dir_roots, d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
                self.cCharge_Vs_Flux.SaveAs('{dir}/Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(dir=self.dir_plots, d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
            else:
                self.cCharge_Vs_Flux.SaveAs('{dir}/Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(dir=self.dir_roots, d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))
                self.cCharge_Vs_Flux.SaveAs('{dir}/Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(dir=self.dir_plots, d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))
            self.cCharge_Vs_Flux_with_cuts.cd()
            self.graph_cuts.Draw('APL')
            self.startPoint_cuts.Draw('P')
            self.finishPoint_cuts.Draw('P')
            self.cCharge_Vs_Flux_with_cuts.BuildLegend(0.65, 0.7, 0.9, 0.9)
            if self.cl_size is 0:
                self.cCharge_Vs_Flux_with_cuts.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(dir=self.dir_roots, d=self.dia, i=self.ini, f=self.fin, e=self.excl))
                self.cCharge_Vs_Flux_with_cuts.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(dir=self.dir_plots, d=self.dia, i=self.ini, f=self.fin, e=self.excl))
            elif self.cl_size < 4:
                self.cCharge_Vs_Flux_with_cuts.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(dir=self.dir_roots, d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
                self.cCharge_Vs_Flux_with_cuts.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(dir=self.dir_plots, d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
            else:
                self.cCharge_Vs_Flux_with_cuts.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(dir=self.dir_roots, d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))
                self.cCharge_Vs_Flux_with_cuts.SaveAs('{dir}/Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(dir=self.dir_plots, d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))

    def get_runinfo(self):
        self.config_file = ConfigParser()
        self.config_file.read('../Configuration/RunConfig_20{tc}.cfg'.format(tc=self.campaign))
        self.mask_dir = self.config_file.get('BASIC', 'maskfilepath')
        runinfo_file = self.config_file.get('BASIC', 'runinfofile')
        f = open(runinfo_file, 'r')
        runinfo = json.load(f)
        f.close()
        return runinfo

    def set_1D_options(self, histo, name, nameX, nameY, color=kBlack, markerStyle=3, scale=2, min=0, max=-999):
        # self.CreateCanvas(name, scale)
        histo.SetTitle(name)
        histo.SetLineColor(color)
        histo.GetXaxis().SetTitle(nameX)
        histo.GetYaxis().SetTitle(nameY)
        histo.SetLineWidth(scale*gStyle.GetLineWidth())
        histo.SetMarkerStyle(markerStyle)
        histo.SetMarkerColor(color)
        histo.SetMarkerSize(float(scale)/3*gStyle.GetMarkerSize())
        if max is not -999: histo.SetMaximum(max)
        histo.SetMinimum(min)
        histo.GetYaxis().SetTitleOffset(1.35)

    def CreateCanvas(self, name, scale=2):
        width, height = 700*scale, 500*scale
        exec('self.c{n} = TCanvas("c1_{n}", "c1_{n}", {w}, {h})'.format(n=name, w=width, h=height))
        exec('self.c{n}.cd()'.format(n=name))
        exec('self.c{n}.SetLeftMargin(0.1)'.format(n=name))
        exec('self.c{n}.SetRightMargin(0.1)'.format(n=name))
        exec('self.c{n}.SetLogx(1)'.format(n=name))

    def calc_flux(self, run):
        if 'for1' not in self.runinfo[str(run)] or self.runinfo[str(run)]['for1'] == 0:
            if 'measuredflux' in self.runinfo[str(run)]:
                return self.runinfo[str(run)]['measuredflux']
        self.maskfile = '{dir}/{m}'.format(dir=self.mask_dir, m=self.runinfo[str(run)]['maskfile'])
        if os.path.isfile(self.maskfile):
            f = open(self.maskfile, 'r')
        else:
            print 'Error: Could not read maskfile!'
            return
        data = []
        for line in f:
            if len(line) > 3:
                line = line.split()
                data.append([int(line[2])] + [int(line[3])])
        f.close()
        pixel_size = 0.01 * 0.015
        area = [(data[1][0] - data[0][0]) * (data[1][1] - data[0][1]) * pixel_size, (data[3][0] - data[2][0]) * (data[3][1] - data[2][1]) * pixel_size]
        flux = [self.runinfo[str(run)]['for{0}'.format(i + 1)] / area[i] / 1000. for i in xrange(2)]
        if mean(flux) <= 15:
            self.run2k.append(run)
        elif mean(flux) <= 50:
            self.run20k.append(run)
        elif mean(flux) <= 130:
            self.run80k.append(run)
        elif mean(flux) <= 500:
            self.run200k.append(run)
        elif mean(flux) <= 1100:
            self.run600k.append(run)
        elif mean(flux) <= 3000:
            self.run2M.append(run)
        else:
            self.run5M.append(run)
        return mean(flux)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', '--initial', dest='ini', type='int', default=303, help='Number of initial run: e.g. 303')
    parser.add_option('-f', '--finish', dest='fin', type='int', default=313, help='Number of last runs: e.g. 313')
    parser.add_option('-d', '--diamond', dest='dia', type='int', default=4, help='ROC number corresponding to the device under test: e.g. 4')
    parser.add_option('-c', '--campaign', dest='camp', type='string', default='1510', help='String of testcampaign: e.g. \'1510\'')
    parser.add_option('-e', '--exclude', dest='e', type='string', default='', help='String with the list of excluded runs in between: e.g. \'304,305,312\'')
    parser.add_option('-s', '--size', dest='s', type='int', default=0, help='Selects the cluster size to analyse. i.e. 0-> all, 1-> 1cl, ...')
    parser.add_option('-a', '--action', dest='a', type='int', default=0, help='Action: 0 do all (default), 1 or other, select diamond and size')
    (options, args) = parser.parse_args()
    ini = int(options.ini)
    fin = int(options.fin)
    dia = int(options.dia)
    camp = str(options.camp)
    e = str(options.e)
    s = int(options.s)
    a = int(options.a)
    if a is 0:
        gROOT.SetBatch(True)
        gROOT.ProcessLine("gErrorIgnoreLevel = {f};".format(f=kError))
        for di in xrange(4,7):
            for cs in xrange(0,5):
                exec('d{di}_s{cs} = RateAnalysisExtraction(ini, fin, {di}, camp, e, {cs})'.format(di=di, cs=cs))
                exec('d{di}_s{cs}.plot_rate_results()'.format(di=di, cs=cs))
        gROOT.SetBatch(False)
    elif a is 10:
        gROOT.SetBatch(True)
        gROOT.ProcessLine("gErrorIgnoreLevel = {f};".format(f=kError))
        for di in xrange(4,6):
            exec('d{di}_s0_1 = RateAnalysisExtraction(297, 313, {di}, camp, "286,288", 0)'.format(di=di))
            exec('d{di}_s0_2 = RateAnalysisExtraction(314, 356, {di}, camp, "317,325,340", 0)'.format(di=di))
        plotAll(d4_s0_1, d4_s0_2, d5_s0_1, d5_s0_2)
        gROOT.SetBatch(False)

    elif a is 20:
        gROOT.SetBatch(True)
        gROOT.ProcessLine("gErrorIgnoreLevel = {f};".format(f=kError))
        for di in xrange(4,6):
            for cs in xrange(1,5):
                exec('d{di}_s{cs}_1 = RateAnalysisExtraction(297, 313, {di}, camp, "286,288", {cs})'.format(di=di, cs=cs))
                exec('d{di}_s{cs}_2 = RateAnalysisExtraction(314, 356, {di}, camp, "317,325,340", {cs})'.format(di=di, cs=cs))
        doClusterAll(d4_s1_1, d4_s2_1, d4_s3_1, d4_s4_1, d5_s1_1, d5_s2_1, d5_s3_1, d5_s4_1, d4_s1_2, d4_s2_2, d4_s3_2, d4_s4_2, d5_s1_2, d5_s2_2, d5_s3_2, d5_s4_2)
        gROOT.SetBatch(False)

    else:
        z = RateAnalysisExtraction(ini, fin, dia, camp, e, s)
        z.plot_rate_results()


