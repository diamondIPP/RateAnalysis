from ROOT import TFile, TCanvas, TProfile, TH2D, kBlue, kRed, TH1D, gStyle, TF1, kWhite, gROOT, kBlack, TGraphErrors, TGraph, kGreen
from optparse import OptionParser
from ConfigParser import ConfigParser
from array import array
from numpy import mean
import os
import json
from glob import glob

__author__ = "Diego Alejandro"

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
        if self.exclude[0] is not '':
            for excl in self.exclude:
                if int(excl) in self.runs:
                    self.runs.remove(int(excl))
                    self.dirs.remove('{c}_{e}'.format(c=self.campaign, e=excl))
        self.dirs = self.get_dirs_dict()
        if self.dia is 4:
            self.histos = {run: self.get_histo(run, 'c_phAllVsEventROC4') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phCl1VsEventROC4') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phCl2VsEventROC4') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phCl3VsEventROC4') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phClM4VsEventROC4') for run in self.runs}
            self.histos_cuts = {run: self.get_histo(run, 'c_phAllVsEventROC4_cuts') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phCl1VsEventROC4_cuts') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phCl2VsEventROC4_cuts') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phCl3VsEventROC4_cuts') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phClM4VsEventROC4_cuts') for run in self.runs}
        elif self.dia is 5:
            self.histos = {run: self.get_histo(run, 'c_phAllVsEventROC5') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phCl1VsEventROC5') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phCl2VsEventROC5') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phCl3VsEventROC5') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phClM4VsEventROC5') for run in self.runs}
            self.histos_cuts = {run: self.get_histo(run, 'c_phAllVsEventROC5_cuts') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phCl1VsEventROC5_cuts') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phCl2VsEventROC5_cuts') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phCl3VsEventROC5_cuts') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phClM4VsEventROC5_cuts') for run in self.runs}
        elif self.dia is 6:
            self.histos = {run: self.get_histo(run, 'c_phAllVsEventROC6') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phCl1VsEventROC6') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phCl2VsEventROC6') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phCl3VsEventROC6') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phClM4VsEventROC6') for run in self.runs}
            self.histos_cuts = {run: self.get_histo(run, 'c_phAllVsEventROC6_cuts') for run in self.runs} if size is 0 else {run: self.get_histo(run, 'c_phCl1VsEventROC6_cuts') for run in self.runs} if size is 1 else {run: self.get_histo(run, 'c_phCl2VsEventROC6_cuts') for run in self.runs} if size is 2 else {run: self.get_histo(run, 'c_phCl3VsEventROC6_cuts') for run in self.runs} if size is 3 else {run: self.get_histo(run, 'c_phClM4VsEventROC6_cuts') for run in self.runs}
        self.projs = {run: self.do_projection(run) for run in self.runs}
        self.projs_cuts = {run: self.do_projection(run, True) for run in self.runs}
        self.means = {run: self.projs[run].GetMean() for run in self.runs}
        self.means_cuts = {run: self.projs_cuts[run].GetMean() for run in self.runs}
        self.sigmas = {run: self.projs[run].GetMeanError() for run in self.runs}
        self.sigmas_cuts = {run: self.projs_cuts[run].GetMeanError() for run in self.runs}
        self.mpvs = {run: self.projs[run].GetBinCenter(self.projs[run].GetMaximumBin()) for run in self.runs}
        self.mpvs_cuts = {run: self.projs_cuts[run].GetBinCenter(self.projs_cuts[run].GetMaximumBin()) for run in self.runs}
        self.runinfo = self.get_runinfo()
        self.fluxes = {run: self.calc_flux(run) for run in self.runs}

    def get_list_dirs(self):
        return glob('{c}*'.format(c=self.campaign))

    def get_runs(self):
        return [int(dir.split('_')[-1]) for dir in self.dirs]

    def get_dirs_dict(self):
        return {run: '{camp}_{r}'.format(camp=self.campaign, r=run) for run in self.runs}

    def get_histo(self, run, name):
        histo = (TFile('{c}_{r}/{name}.root'.format(c= self.campaign, r=run, name=name), 'READ').Get(name)).GetPrimitive(name.replace('c_', ''))
        name_old = histo.GetName()
        histo.SetName('{name}_{r}'.format(name=name_old, r=run))
        return histo

    def do_projection(self, run, with_cuts=True):
        name = self.histos_cuts[run].GetName() if with_cuts else self.histos[run].GetName()
        if with_cuts:
            return self.histos_cuts[run].ProjectionY('p_{n}'.format(n=name))
        else:
            return self.histos[run].ProjectionY('p_{n}'.format(n=name))

    def plot_results(self):
        self.graph = TGraphErrors(len(self.means))
        self.graph.SetNameTitle('fluxes', 'fluxes')
        self.graph_cuts = TGraphErrors(len(self.means_cuts))
        self.graph_cuts.SetNameTitle('fluxes_with_cuts', 'fluxes_with_cuts')
        for i in xrange(len(self.means)):
            self.graph.SetPoint(i, self.fluxes[self.runs[i]], self.means[self.runs[i]])
            self.graph.SetPointError(i, self.fluxes[self.runs[i]] * 0.1, self.sigmas[self.runs[i]])
        for i in xrange(len(self.means_cuts)):
            self.graph_cuts.SetPoint(i, self.fluxes[self.runs[i]], self.means_cuts[self.runs[i]])
            self.graph_cuts.SetPointError(i, self.fluxes[self.runs[i]] * 0.1, self.sigmas_cuts[self.runs[i]])
        self.set_1D_options(self.graph, 'PH_Flux_Scan', 'Flux (kHz/cm^{2}/s)', 'mean Charge (e)',kBlue, 20, 2, 0, self.graph.GetMaximum())
        self.set_1D_options(self.graph_cuts, 'PH_Flux_Scan', 'Flux (kHz/cm^{2}/s)', 'mean Charge (e)',kBlue, 20, 2, 0, self.graph_cuts.GetMaximum())
        self.CreateCanvas('Charge_Vs_Flux', 2)
        self.CreateCanvas('Charge_Vs_Flux_with_cuts', 2)
        self.startPoint = TGraph(1)
        self.startPoint_cuts = TGraph(1)
        self.startPoint.SetName('Run Start')
        self.startPoint_cuts.SetName('Run Start wc')
        self.startPoint.SetTitle('Run Start')
        self.startPoint_cuts.SetTitle('Run Start wc')
        self.finishPoint = TGraph(1)
        self.finishPoint_cuts = TGraph(1)
        self.finishPoint.SetName('Run Finish')
        self.finishPoint_cuts.SetName('Run Finish wc')
        self.finishPoint.SetTitle('Run Finish')
        self.finishPoint_cuts.SetTitle('Run Finish wc')
        self.startPoint.SetPoint(0, self.fluxes[self.runs[0]], self.means[self.runs[0]])
        self.startPoint_cuts.SetPoint(0, self.fluxes[self.runs[0]], self.means_cuts[self.runs[0]])
        self.finishPoint.SetPoint(0, self.fluxes[self.runs[-1]], self.means[self.runs[-1]])
        self.finishPoint_cuts.SetPoint(0, self.fluxes[self.runs[-1]], self.means_cuts[self.runs[-1]])
        self.startPoint.SetMarkerStyle(22)
        self.startPoint_cuts.SetMarkerStyle(22)
        self.startPoint.SetMarkerColor(kGreen)
        self.startPoint_cuts.SetMarkerColor(kGreen)
        self.startPoint.SetMarkerSize(gStyle.GetMarkerSize()*1.5)
        self.startPoint_cuts.SetMarkerSize(gStyle.GetMarkerSize()*1.5)
        self.finishPoint.SetMarkerStyle(23)
        self.finishPoint_cuts.SetMarkerStyle(23)
        self.finishPoint.SetMarkerColor(kRed)
        self.finishPoint_cuts.SetMarkerColor(kRed)
        self.finishPoint.SetMarkerSize(gStyle.GetMarkerSize()*1.5)
        self.finishPoint_cuts.SetMarkerSize(gStyle.GetMarkerSize()*1.5)
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
        self.cCharge_Vs_Flux.cd()
        self.graph.Draw('APL')
        self.startPoint.Draw('P')
        self.finishPoint.Draw('P')
        self.cCharge_Vs_Flux.BuildLegend(0.65, 0.2, 0.9, 0.4)
        if self.cl_size is 0:
            self.cCharge_Vs_Flux.SaveAs('Charge_Vs_Flux_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(d=self.dia, i=self.ini, f=self.fin, e=self.excl))
            self.cCharge_Vs_Flux.SaveAs('Charge_Vs_Flux_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(d=self.dia, i=self.ini, f=self.fin, e=self.excl))
        elif self.cl_size < 4:
            self.cCharge_Vs_Flux.SaveAs('Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
            self.cCharge_Vs_Flux.SaveAs('Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
        else:
            self.cCharge_Vs_Flux.SaveAs('Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))
            self.cCharge_Vs_Flux.SaveAs('Charge_Vs_Flux_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))
        self.cCharge_Vs_Flux_with_cuts.cd()
        self.graph_cuts.Draw('APL')
        self.startPoint_cuts.Draw('P')
        self.finishPoint_cuts.Draw('P')
        self.cCharge_Vs_Flux_with_cuts.BuildLegend(0.65, 0.2, 0.9, 0.4)
        if self.cl_size is 0:
            self.cCharge_Vs_Flux_with_cuts.SaveAs('Charge_Vs_Flux_with_cuts_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(d=self.dia, i=self.ini, f=self.fin, e=self.excl))
            self.cCharge_Vs_Flux_with_cuts.SaveAs('Charge_Vs_Flux_with_cuts_ROC{d}_All_Cl_sizes_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(d=self.dia, i=self.ini, f=self.fin, e=self.excl))
        elif self.cl_size < 4:
            self.cCharge_Vs_Flux_with_cuts.SaveAs('Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
            self.cCharge_Vs_Flux_with_cuts.SaveAs('Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(d=self.dia, s=self.cl_size, i=self.ini, f=self.fin, e=self.excl))
        else:
            self.cCharge_Vs_Flux_with_cuts.SaveAs('Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.root'.format(d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))
            self.cCharge_Vs_Flux_with_cuts.SaveAs('Charge_Vs_Flux_with_cuts_ROC{d}_{s}_pix_Cl_size_Ini_{i}_Fin_{f}_Exc_{e}.png'.format(d=self.dia, s='4More', i=self.ini, f=self.fin, e=self.excl))

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
        histo.GetYaxis().SetTitleOffset(1.2)

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
        return mean(flux)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-i', '--initial', dest='ini', type='int', default=303, help='Number of initial run: e.g. 303')
    parser.add_option('-f', '--finish', dest='fin', type='int', default=313, help='Number of last runs: e.g. 313')
    parser.add_option('-d', '--diamond', dest='dia', type='int', default=4, help='ROC number corresponding to the device under test: e.g. 4')
    parser.add_option('-c', '--campaign', dest='camp', type='string', default='1510', help='String of testcampaign: e.g. \'1510\'')
    parser.add_option('-e', '--exclude', dest='e', type='string', default='', help='String with the list of excluded runs in between: e.g. \'304,305,312\'')
    parser.add_option('-s', '--size', dest='s', type='int', default=0, help='Selects the cluster size to analyse. i.e. 0-> all, 1-1cl, ...')
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
        for di in xrange(4,7):
            for cs in xrange(0,5):
                exec('d{di}_s{cs} = RateAnalysisExtraction(ini, fin, {di}, camp, e, {cs})'.format(di=di, cs=cs))
                exec('d{di}_s{cs}.plot_results()'.format(di=di, cs=cs))
    else:
        z = RateAnalysisExtraction(ini, fin, dia, camp, e, s)
        z.plot_results()
