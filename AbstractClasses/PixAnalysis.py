#!/usr/bin/env python
# --------------------------------------------------------
#       Main class for Rate Pixel Analysis
# created some time in 2016 by D. Sanz (sandiego@phys.ethz.ch), maintained by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TH2D, TH1D, gROOT, TFormula, TCut, TH1I, TProfile, THStack, TProfile2D, TF1, TGraph, TPie, gRandom, TH3D, TMultiGraph, TCutG, TH1F
from argparse import ArgumentParser
from collections import OrderedDict, Counter
from copy import deepcopy
from math import ceil
from numpy import array, mean, corrcoef, arange
from os.path import join as joinpath
from time import sleep

from CurrentInfo import Currents
from CutPix import CutPix
from Elementary import Elementary
from TelescopeAnalysis import Analysis
from Utils import *

__author__ = 'DA & Micha'


# ==============================================
# MAIN CLASS
# ==============================================
class PixAnalysis(Analysis):
    def __init__(self, run, dut=1, load_tree=True, verbose=False, binning=10000):

        Analysis.__init__(self, run, verbose=verbose, binning=binning, load_tree=load_tree)

        # main
        self.RunNumber = run
        self.DiamondName = self.load_diamond_name(dut)
        self.DiamondNumber = dut
        self.Bias = self.run.Bias[dut - 1]
        self.Dut = dut + 3
        self.save_dir = '{dia}/{run}/'.format(run=str(self.RunNumber).zfill(3), dia=self.DiamondName)
        self.NRocs = self.load_n_rocs()

        if load_tree:
            # stuff
            self.plots.save_dir = self.save_dir

            # cuts
            self.Settings = self.plots.Settings
            self.Cut = CutPix(self, dut)

            # alignment
            # self.IsAligned = self.check_alignment()

            # pulse height calibrations
            self.Fit = None
            self.Parameters = None
            self.Vcals = None
            self.Points = None
            self.get_calibration_data()

        # currents
        self.Currents = Currents(self)

        self.stuff = []

    def __del__(self):
        for c in gROOT.GetListOfCanvases():
            c.Close()
        for lst in self.histos.itervalues():
            if not type(lst) is list:
                lst = [lst]
            for obj in lst:
                self.del_rootobj(obj)

    def draw_current(self, relative_time=True, volt_range=None):
        self.Currents.draw_indep_graphs(rel_time=relative_time, v_range=volt_range)

    # ==========================================================================
    # region INIT
    def load_diamond_name(self, dut):
        assert dut in [1, 2, 3], 'You have to choose either dut 1, 2 or 3'
        return self.run.DiamondNames[dut - 1]

    def load_n_rocs(self):
        if self.ana_config_parser.has_option('BASIC', 'n_rocs'):
            return self.ana_config_parser.getint('BASIC', 'n_rocs')
        else:
            return 7
    # endregion

    def do_analysis(self, do_cut_dist=False, do_res_ana=False, do_cut_ana=False, do_occupancy=True, do_pulse_height=True):
        """ Does automatic analysis with the selected options """
        TFormula.SetMaxima(1000000, 10000, 10000000)  # (1000,1000,1000)

        if do_cut_dist:
            self.Cut.do_cuts_distributions()
        if do_res_ana:
            self.Cut.do_res_analysis()
        if do_cut_ana:
            self.do_cuts_analysis(do_occupancy, do_pulse_height, True)

    def do_cuts_analysis(self, do_occupancy, do_pulse_height, normalize_ph_plots=False):
        """ Calls do_cut_analysis in the Cut method """
        self.Cut.do_cuts_analysis(do_occupancy, do_pulse_height, normalize_ph_plots)

    def draw_occupancy(self, cut=None, show=True, fid=False, prnt=True, adc=None, roc=None, tel_coods=False, res=sqrt(12)):
        """ Does the occupancy of a roc with the specified cut and it is saved on the given histogram. If none is given, it will create a histogram and return a deepcopy of it """
        roc = self.Dut if roc is None else roc
        cut_string = self.Cut.generate_special_cut(excluded='fiducial' if not fid else [], cluster=False) if cut is None else TCut(cut)
        cut_string += '{p}=={d}'.format(d=roc, p='cluster_plane' if tel_coods else 'plane')
        cut_string += self.Cut.add_adc_cut(adc)
        self.set_root_output(False)
        h = TH2D('h_oc', 'Occupancy {d}'.format(d=self.DiamondName), *(self.Settings['2DBins'] if not tel_coods else self.plots.get_global_bins(res)))
        draw_var = 'row:col' if not tel_coods else 'cluster_ypos_tel:cluster_xpos_tel'
        self.tree.Draw('{d} >> {n}'.format(n='h_oc', d=draw_var), cut_string, 'goff')
        save_name = 'Occupancy{c}'.format(c=make_cut_string(cut, self.Cut.NCuts))
        set_statbox(x=.81, entries=8, opt=1000000010)
        xtit, ytit = ('col', 'row') if not tel_coods else ('x [cm]', 'y [cm]')
        self.format_histo(h, x_tit=xtit, y_tit=ytit, z_tit='Number of Entries', y_off=1.3, z_off=1.5)
        self.save_histo(h, save_name, show, rm=.17, lm=.13, draw_opt='colz', prnt=prnt)
        return h

    def draw_time_occupancy(self, cut=None, roc=None, fid=False, binning=10000):
        self.set_bin_size(binning)
        roc = self.Dut if roc is None else roc
        cut_string = self.Cut.generate_special_cut(excluded='fiducial' if not fid else [], cluster=False) if cut is None else TCut(cut)
        cut_string += 'plane=={r}'.format(r=roc)
        h = TH3D('h_to', 'to', len(self.time_binning) - 1, array([t / 1000. for t in self.time_binning], 'd'), *self.plots.get_arrays(self.Settings['2DBins']))
        self.tree.Draw('row:col:time/1000.>>h_to', cut_string, 'goff')
        self.format_histo(h, y_tit='col', z_tit='row')
        gStyle.SetNumberContours(20)
        titles = ['Mean X', 'Sigma X', 'Mean Y', 'Sigma Y']
        graphs = [self.make_tgrapherrors('g_to{i}'.format(i=i), titles[i]) for i in xrange(4)]
        for ibin in xrange(h.GetNbinsX() - 1):
            h.GetXaxis().SetRange(ibin, ibin + 1)
            p = h.Project3D('zy')
            self.draw_histo(p, draw_opt='colz') if not ibin else p.Draw('samecolz')
            c = gROOT.GetListOfCanvases()[-1]
            c.Update()
            c.Modified()
            sleep(.5)
        for ibin in xrange(h.GetNbinsX() - 1):
            t = h.GetXaxis().GetBinCenter(ibin)
            h.GetXaxis().SetRange(ibin, ibin + 1)
            p = h.Project3D('zy')
            px = p.ProjectionX()
            py = p.ProjectionY()
            fits = [px.Fit('gaus', 'qs', '', .5, 25), py.Fit('gaus', 'qs', '', 15, 60)]
            for i in xrange(4):
                graphs[i].SetPoint(ibin, t, fits[i % 2].Parameter(1 + i % 2))
                graphs[i].SetPointError(ibin, 0, fits[i % 2].ParError(1 + i % 2))
        for i in xrange(4):
            set_time_axis(graphs[i], off=self.run.startTime / 1000 + 3600)
            self.draw_histo(graphs[i], draw_opt='alp')
        return h

    def draw_pulse_height_vs_event(self, cut=None, show=True, adc=False):
        """ Pulse height analysis vs event for a given cut. If no cut is provided it will take all. """
        if cut is None:
            cut_string = self.Cut.all_cut if not adc else self.Cut.HitMapCut
        else:
            cut_string = TCut(cut)
        cut_string += 'charge_all_ROC{d} > 0'.format(d=self.Dut) if not adc else 'plane == {n}'.format(n=self.Dut)
        ybins = [self.Settings['ph1DbinsD{n}'.format(n=self.Dut)], self.Settings['ph1DminD{n}'.format(n=self.Dut)], self.Settings['ph1DmaxD{n}'.format(n=self.Dut)]]
        ybins = [255, 0, 255] if adc else ybins
        self.set_root_output(False)
        h = TH2D('h_ph', 'Pulse Height {d}'.format(d=self.DiamondName), self.n_bins - 1, array([t / 1000 for t in self.time_binning]), *ybins)
        cut_var = 'charge_all_ROC{d}'.format(d=self.Dut) if not adc else 'adc'
        self.tree.Draw('{v}:time / 1000. >> h_ph'.format(v=cut_var), cut_string, 'goff')
        self.format_histo(h, x_tit='Time [hh:mm]', y_tit='Cluster Charge [e]' if not adc else 'adc', z_tit='Number of Entries', y_off=2.05, z_off=1.3, stats=0)
        set_time_axis(h, off=self.run.startTime / 1000 + 3600)
        h.SetNdivisions(520)
        self.save_histo(h, 'PulseHeightVsEvent', show, rm=.15, lm=.16, draw_opt='colz', save=show)
        return h

    def draw_adc_vs_event(self, cut=None, show=True):
        h = self.draw_pulse_height_vs_event(cut, show=False, adc=True)
        self.format_histo(h, title='ADC vs Time - {d}'.format(d=self.DiamondName))
        self.save_histo(h, 'ADCvsEvent', show, rm=.15, lm=.16, draw_opt='colz', logz=True)
        return h

    def draw_adc_map(self, show=True, cut=None, adc=None):
        cut_string = deepcopy(self.Cut.HitMapCut) if cut is None else TCut(cut)
        cut_string += 'plane == {n}'.format(n=self.Dut)
        cut_string += '' if adc is None else 'adc>0'
        self.set_root_output(False)
        h = TProfile2D('p_am', 'ADC Map', *self.Settings['2DBins'])
        self.tree.Draw('adc:row:col >> p_am', cut_string, 'goff')
        set_statbox(x=.81, entries=8, opt=1000000010)
        self.format_histo(h, x_tit='col', y_tit='row', z_tit='ADC', y_off=1.3, z_off=1.5)
        self.save_histo(h, 'ADCMap', show, rm=.17, lm=.13, draw_opt='colz')

    def draw_zero_contribution(self, show=True, cut=None, entries=False, one_hit=False):
        cut_string = deepcopy(self.Cut.HitMapCut) if cut is None else TCut(cut)
        cut_string += 'plane == {n}'.format(n=self.Dut)
        cut_string += 'adc>0' if entries else ''
        cut_string += 'n_hits[4]==1' if one_hit else ''
        h = TProfile2D('p_zc', 'Good Events', *self.Settings['2DBins']) if not entries else TH2D('h_zc', 'Good Events', *self.Settings['2DBins'])
        self.tree.Draw('((adc>0)*100):row:col >> p_zc' if not entries else 'row:col >> h_zc', cut_string, 'goff')
        set_statbox(x=.81, entries=8, opt=1000000010)
        self.format_histo(h, x_tit='col', y_tit='row', z_tit='Good Events{p}'.format(p='' if entries else ' [%]'), y_off=1.3, z_off=1.5)
        self.save_histo(h, 'ZeroContribution{e}'.format(e='Entries' if entries else ''), show, rm=.17, lm=.13, draw_opt='colz')

    def get_calibration_data(self):
        f = open(joinpath(self.run.converter.TrackingDir, 'calibration_lists', 'GKCalibrationList_Telescope{n}.txt'.format(n=self.run.converter.TelescopeID)))
        lines = f.readlines()
        f.close()
        # calibration fit
        file_names = [joinpath(self.run.converter.TrackingDir, lines[0].strip('./\n'), line.strip('\n')) for i, line in enumerate(lines) if i]
        fit = None
        params = [[[0 for _ in xrange(self.Settings['nRows'])] for _ in xrange(self.Settings['nCols'])] for _ in xrange(self.NRocs)]
        for roc, file_name in enumerate(file_names):
            f = open(file_name)
            f.readline()
            fit_string = f.readline().replace('par', '').strip('\n')
            fit = TF1('ErFit', fit_string, 0, 255 * 7)
            f.readline()
            for line in f.readlines():
                line = line.split()
                params[roc][int(line[-2])][int(line[-1])] = [float(line[i]) for i in xrange(4)]
            f.close()
        # calibration points
        calib_files = [joinpath(self.run.converter.TrackingDir, lines[0].strip('./\n'), 'phCalibration_C{n}.dat'.format(n=n)) for n in xrange(self.NRocs)]
        vcals = None
        points = [[[0 for _ in xrange(self.Settings['nRows'])] for _ in xrange(self.Settings['nCols'])] for _ in xrange(self.NRocs)]
        for roc, file_name in enumerate(calib_files, 0):
            try:
                f = open(file_name)
                f.readline()
                vcals = [int(i) for i in f.readline().split()[2:]] + [int(i) * 7 for i in f.readline().split()[2:]]
                f.readline()
                for line in f.readlines():
                    line = line.split()
                    points[roc][int(line[-2])][int(line[-1])] = [int(line[i]) for i in xrange(len(vcals))]
            except IOError as e:
                log_warning(e)

        if self.Fit is None:
            self.Fit = fit
            self.Parameters = params
            self.Vcals = vcals
            self.Points = points

    def draw_calibration_fit(self, col, row, show=True, roc=None):
        roc = self.Dut if roc is None else roc
        self.Fit.SetParameters(*self.Parameters[roc][col][row])
        self.format_histo(self.Fit, title='Calibration Fit for Pix {c} {r}'.format(c=col, r=row), x_tit='vcal', y_tit='adc', y_off=1.4, color=632, lw=2)
        gr = TGraph(len(self.Vcals), array(self.Vcals, 'd'), array(self.Points[roc][col][row], 'd'))
        self.format_histo(gr, marker=20, name='gr_fp', title='Calibration Fit for Pix {c} {r}'.format(c=col, r=row), x_tit='vcal', y_tit='adc', y_off=1.4)
        self.draw_histo(gr, lm=.12, draw_opt='ap')
        self.Fit.Draw('same')
        self.save_plots('CalFit{c}{r}'.format(c=col, r=row), show=show)

    def draw_threshold_map(self, show=True, vcal=True):
        h = TProfile2D('p_tm', 'Artificial Threshold Map', *self.Settings['2DBins'])
        for (col, row), thresh in self.get_thresholds(vcal=vcal).iteritems():
            h.Fill(col, row, thresh if vcal else thresh / 1000.)
        self.format_histo(h, x_tit='col', y_tit='row', z_tit='Artificial Treshold [{u}]'.format(u='vcal' if vcal else 'ke'), y_off=1.3, z_off=1.5, stats=0)
        self.save_histo(h, 'ThresholdMap', show, rm=.17, lm=.13, draw_opt='colz')
        return h

    def get_thresholds(self, cols=None, pix=None, vcal=True):
        columns, rows = self.Cut.CutConfig['FidRegionLocal'][:2], self.Cut.CutConfig['FidRegionLocal'][2:]
        columns = cols if cols is not None else columns
        columns = [columns, columns] if type(columns) == int else columns
        columns, rows = ([pix[0]] * 2, [pix[1]] * 2) if pix is not None else (columns, rows)
        dic = {}
        for col in xrange(columns[0], columns[1] + 1):
            for row in xrange(rows[0], rows[1] + 1):
                self.Fit.SetParameters(*self.Parameters[self.Dut][col][row])
                dic[(col, row)] = self.Fit.GetX(0) * (47 if not vcal else 1)
        return dic

    def fit_erf(self, col, row, roc=None, show=True):
        roc = self.Dut if roc is None else roc
        # fit.SetParameters(*self.Parameters[roc][col][row])
        points = OrderedDict(sorted({vcal: value for vcal, value in zip(self.Vcals, self.Points[roc][col][row])}.iteritems()))
        good_points = OrderedDict(sorted({vcal: value for vcal, value in zip(self.Vcals, self.Points[roc][col][row]) if value}.iteritems()))
        start = 0
        for vcal, value in points.iteritems():
            if value:
                start = vcal
                break
        if len(good_points) > 3:
            fit = deepcopy(self.Fit)
            fit.SetParameters(309.2062, 112.8961, 1.022439, 35.89524)
        else:
            fit = TF1('fit', 'pol1', 0, 3000)
        if len(good_points) > 1:
            gr = TGraph(len(points), array(points.keys(), 'd'), array(points.values(), 'd'))
            gr.Fit(fit, 'q', '', start, 3000 if len(good_points) > 3 else 1000)
            self.format_histo(gr, marker=20, x_tit='vcal', y_tit='adc', y_off=1.3, title='Calibration Fit for Pix {c} {r}'.format(c=col, r=row))
            self.draw_histo(gr, draw_opt='ap', show=show, lm=.12)
            return fit

    def draw_adc_fixed_vcal_map(self, show=True, vcal=200, roc=None):
        roc = self.Dut if roc is None else roc
        h = TProfile2D('p_pm', 'Pulse Height Map for Vcal {v}'.format(v=vcal), *self.Settings['2DBins'])
        cols, rows = self.Cut.CutConfig['MaskCols'], self.Cut.CutConfig['MaskRows']
        for col in xrange(cols[0][1] + 1, cols[1][0]):
            for row in xrange(rows[0][1], self.Settings['nRows']):
                self.Fit.SetParameters(*self.Parameters[roc][col][row])
                h.Fill(col, row, self.Fit(vcal))
        self.format_histo(h, x_tit='col', y_tit='row', z_tit='Pulse Height [adc]', y_off=1.3, z_off=1.5, stats=0)
        self.save_histo(h, 'PulseHeightMap{v}'.format(v=vcal), show, rm=.17, lm=.13, draw_opt='colz')

    def draw_adc_disto(self, cut=None, show=True, col=None, pix=None):
        h = TH1I('h_adc', 'ADC Distribution {d}'.format(d=self.DiamondName), 255, 0, 255)
        cut_string = deepcopy(self.Cut.HitMapCut) if cut is None else TCut(cut)
        cut_string += 'plane == {n}'.format(n=self.Dut)
        cut_string += 'col=={c}'.format(c=col) if col is not None else ''
        cut_string += 'col=={c}&&row=={r}'.format(c=pix[0], r=pix[1]) if pix is not None else ''
        self.set_root_output(False)
        self.tree.Draw('adc>>h_adc', cut_string, 'goff')
        set_statbox(entries=8, opt=1000000010)
        self.format_histo(h, x_tit='adc', y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor)
        self.save_histo(h, 'ADCDisto', show, lm=0.13, logy=True)
        return h

    def draw_vcal_disto(self, cut=None, show=True, col=None, pix=None, rnge=None, electrons=False):
        h = TH1F('h_vcal', 'vcal Distribution {d}'.format(d=self.DiamondName), 200, 0, 1000 if not electrons else 100000)
        cut_string = deepcopy(self.Cut.HitMapCut) if cut is None else TCut(cut)
        cut_string += 'plane == {n}'.format(n=self.Dut)
        cut_string += 'n_hits[{n}] == 1'.format(n=self.Dut)
        cut_string += 'col=={c}'.format(c=col) if col is not None else ''
        cut_string += 'col=={c}&&row=={r}'.format(c=pix[0], r=pix[1]) if pix is not None else ''
        cut_string += 'col>={c1}&&col<={c2}&&row>={r1}&&row<={r2}'.format(c1=rnge[0], c2=rnge[1], r1=rnge[2], r2=rnge[3]) if rnge is not None else ''
        self.set_root_output(False)
        n = self.tree.Draw('adc:col:row', cut_string, 'goff')
        for i in xrange(n):
            row = int(self.tree.GetV3()[i])
            col = int(self.tree.GetV2()[i])
            adc = self.tree.GetV1()[i]
            self.Fit.SetParameters(*self.Parameters[self.Dut][col][row])
            h.Fill(self.Fit.GetX(adc) * (47 if electrons else 1))
        set_statbox(entries=8, opt=1000000010)
        self.format_histo(h, x_tit='Pulse Height [{u}]'.format(u='vcal' if not electrons else 'e'), y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor)
        self.save_histo(h, '{p}Disto'.format(p='Ph' if electrons else 'Vcal'), show, lm=0.13)
        return h

    def check_adc(self):
        for i in xrange(10000, 12000):
            self.tree.GetEntry(i)
            if 4 in self.tree.plane:
                ind = list(self.tree.plane).index(4)
                if self.tree.adc[ind]:
                    print i, self.tree.adc[ind], list(self.tree.charge_all_ROC4)

    def draw_pulse_height_disto(self, cut=None, show=True, prnt=True, sup_zero=False, pix=None, roc=None, vcal=False, redo=False):
        area_string = '_'.join(str(pi) for pi in pix) if pix is not None else 'all'
        roc = self.Dut if roc is None else roc
        pickle_path = self.make_pickle_path('PulseHeight', run=self.RunNumber, suf='{r}_{a}{z}'.format(r=roc, a=area_string, z='_0' if sup_zero else ''))
        cut_string = deepcopy(self.Cut.all_cut) if cut is None else TCut(cut)
        cut_string += 'cluster_charge>0'.format(d=self.Dut) if sup_zero else ''
        cut_string += 'cluster_plane=={r}'.format(r=roc)
        if pix is not None:
            pix = [pix, pix] if not type(pix[0]) == list else pix
            cut_string += 'cluster_col>={c1}&&cluster_col<={c2}&&cluster_row>={r1}&&cluster_row<={r2}'.format(c1=pix[0][0], c2=pix[1][0], r1=pix[0][1], r2=pix[1][1], d=self.Dut)

        def func():
            self.set_root_output(False)
            h1 = TH1D('h_phd', 'Pulse Height Distribution - {d}'.format(d=self.DiamondName), *self.Settings['phBins' if not vcal else 'vcalBins'])
            self.tree.Draw('cluster_charge{v}>>h_phd'.format(d=self.Dut, v='/47.5 + 427.4/47.5' if vcal else ''), cut_string, 'goff')
            set_statbox(entries=8, opt=1000000010, x=.92)
            self.format_histo(h1, x_tit='Pulse Height [{u}]'.format(u='vcal' if vcal else 'e'), y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor)
            self.save_histo(h1, 'PulseHeightDisto{c}'.format(c=make_cut_string(cut, self.Cut.NCuts)), show, lm=.13, prnt=prnt, rm=.06)
            return h1

        h = func() if redo else None
        h = self.do_pickle(pickle_path, func, h)
        self.draw_histo(h, show=show, lm=.13, prnt=prnt, rm=.06)
        return h

    def draw_cluster_disto(self, n=1, cut=None):
        cut_string = deepcopy(self.Cut.all_cut) if cut is None else TCut(cut) + TCut('clusters_per_plane[{r}] == {n}'.format(r=self.Dut, n=n))
        self.draw_pulse_height_disto(cut=cut_string + TCut('clusters_per_plane[{r}] == {n}'.format(r=self.Dut, n=n)), redo=True)

    def draw_cluster_ph_distos(self, cut=None, show=True):
        cut_string = deepcopy(self.Cut.all_cut) if cut is None else TCut(cut)
        hs = OrderedDict([(str(n), self.draw_pulse_height_disto(cut=cut_string + TCut('clusters_per_plane[{r}] == {n}'.format(r=self.Dut, n=n)), show=False, redo=True)) for n in xrange(1, 4)])
        hs['>3'] = self.draw_pulse_height_disto(cut=cut_string + TCut('clusters_per_plane[{r}] > 3'.format(r=self.Dut)), show=False, redo=True)
        stack = THStack('h_cph', 'Pulser Height per Cluster Size')
        l = self.make_legend(y2=.5, nentries=5, x1=.59)
        for name, h in hs.iteritems():
            self.format_histo(h, color=self.get_color())
            stack.Add(h, name)
            l.AddEntry(h, name, 'h')
        self.format_histo(stack, x_tit='Pulse Height [e]', y_tit='Number of Entries', y_off=1.4, draw_first=True)
        self.save_histo(stack, 'ClusterPulseHeight', show, lm=.13, l=l, draw_opt='nostack', gridy=True)
        self.reset_colors()

    def draw_pulse_height_map(self, show=True, cut=None, roc=None, sup_zero=True):
        roc = self.Dut if roc is None else roc
        cut_string = self.Cut.all_cut if cut is None else TCut(cut)
        cut_string += 'cluster_plane=={r}'.format(r=roc)
        cut_string += 'cluster_charge>0'.format(d=self.Dut) if sup_zero else ''
        self.set_root_output(False)
        h = TProfile2D('p_phm', 'Pulse Height Map', *self.Settings['2DBins'])
        self.tree.Draw('cluster_charge:cluster_row:cluster_col>>p_phm'.format(d=self.Dut), cut_string, 'goff')
        set_statbox(entries=8, opt=1000000010)
        self.format_histo(h, x_tit='col', y_tit='row', z_tit='Pulse Height [e]', z_off=1.5, y_off=1.4)
        self.save_histo(h, 'PulseHeightMap', show, lm=.13, rm=.15, draw_opt='colz')

    def draw_hit_efficiency(self, roc=None, save=True, cut='all', vs_time=True, binning=5000, n=1e9, start=0, show=True):
        roc = self.Dut if roc is None else roc
        self.set_root_output(False)
        suffix = 'ROC {n}'.format(n=roc) if roc < 4 else self.load_diamond_name(roc - 3)
        h = TProfile('h_he', 'Hit Efficiency {s}'.format(s=suffix), *(self.get_time_bins(binning) if vs_time else self.get_bins(binning)))
        cut_string = self.Cut.generate_special_cut(excluded=['masks', 'rhit']) if cut == 'all' else TCut(cut)
        x_var = 'time / 1000' if vs_time else 'event_number'
        self.tree.Draw('(n_hits[{r}]>0)*100:{x} >> h_he'.format(r=roc, x=x_var), cut_string, 'goff', int(n), start)
        set_time_axis(h, off=self.run.startTime / 1000 + 3600) if vs_time else do_nothing()
        self.format_histo(h, x_tit='Time [hh:mm]' if vs_time else 'Event Number', y_tit='Efficiency [%]', y_off=1.4, ndiv=505, y_range=[-5, 105], stats=0)
        self.save_histo(h, 'HitEfficiencyROC{n}'.format(n=roc), show, lm=.13, save=save, gridy=True)
        return h

    def fit_hit_efficiency(self, roc=None, show=True, save=True, cut='all', vs_time=True, n=1e9, start=0, binning=5000):
        pickle_path = self.make_pickle_path('Efficiency', run=self.RunNumber, suf='{r}{c}_{t}'.format(r=roc, c='_Cuts' if cut else '', t='Time' if vs_time else 'EvntNr'))

        def func():
            set_statbox(y=.37, only_fit=True, entries=1.5, form='4.1f')
            h = self.draw_hit_efficiency(roc, show=False, save=False, cut=cut, vs_time=vs_time, n=n, start=start, binning=binning)
            if h.GetEntries() < 100:
                return FitRes()
            self.format_histo(h, stats=1, name='Fit Result')
            fit = h.Fit('pol0', 'qs')
            self.save_histo(h, 'HitEfficiencyROC{n}Fit'.format(n=roc), show, lm=.13, save=save, gridy=True, prnt=show)
            return FitRes(fit)

        fit_res = func() if show or save else None
        return self.do_pickle(pickle_path, func, fit_res)

    def draw_all_efficiencies(self, show=True):
        stack = THStack('s_he', 'Raw Hit Efficiencies')
        l = self.make_legend(y2=.5, nentries=5, x1=.59)
        for roc in xrange(self.NRocs):
            h = self.draw_hit_efficiency(roc, show=False)
            fit = self.fit_hit_efficiency(roc, show=False)
            self.format_histo(h, color=self.get_color())
            stack.Add(h, 'ROC{n}'.format(n=roc))
            leg_string = 'ROC{n}'.format(n=roc) if roc < 4 else self.load_diamond_name(roc - 3)
            leg_string += ' ({v:5.2f}%)'.format(v=fit.Parameter(0))
            l.AddEntry(h, leg_string, 'pl')
        self.format_histo(stack, x_tit='Time [hh:mm]', y_tit='Efficiency [%]', y_off=1.4, ndiv=505, y_range=[-5, 105], stats=0, draw_first=True)
        set_time_axis(stack, off=self.run.startTime / 1000 + 3600)
        self.save_histo(stack, 'HitEfficiencies', show, lm=.13, l=l, draw_opt='nostack', gridy=True)
        self.reset_colors()
        return stack

    def draw_hits_dut(self, cut=None, show=True):
        cut_string = z.Cut.all_cut if cut is None else TCut(cut)
        h = TH2D('h_hd', 'Hits in Dia vs Hits in Silicon', 40, 0, 40, 40, 0, 40)
        self.tree.Draw('n_hits[5]:n_hits[4]>>h_hd', cut_string, 'goff')
        set_statbox(entries=4, opt=1000000010, x=.81)
        self.format_histo(h, x_tit='Hits in Diamond', y_tit='Hits in Silicon', y_off=1.3, z_tit='Number of Entries', z_off=1.4)
        self.save_histo(h, 'HitsDiaSil', show, draw_opt='colz', rm=0.17, lm=.13, logz=True)

    def draw_hit_pie(self):
        zero = z.tree.GetEntries('n_hits[4]==0&&n_hits[5]==0')
        dia = z.tree.GetEntries('n_hits[4]>0&&n_hits[5]==0')
        sil = z.tree.GetEntries('n_hits[4]==0&&n_hits[5]>0')
        both = z.tree.GetEntries('n_hits[4]>0&&n_hits[5]>0')
        names = ['No Hits', 'Diamond Hit', 'Silicon Hit', 'Both Hits']
        values = [zero, dia, sil, both]
        colors = [self.get_color() for _ in xrange(1, len(values) + 1)]
        self.reset_colors()
        pie = TPie('pie', 'Hit Contributions', len(values), array(values, 'f'), array(colors, 'i'))
        for i, label in enumerate(names):
            pie.SetEntryRadiusOffset(i, .05)
            pie.SetEntryLabel(i, label.title())
        pie.SetHeight(.04)
        pie.SetRadius(.2)
        pie.SetTextSize(.025)
        pie.SetAngle3D(70)
        pie.SetAngularOffset(250)
        self.draw_histo(pie, draw_opt='3drsc')

    @staticmethod
    def draw_fid_cut():
        cut = gROOT.FindObject('fid')
        cut.Draw()

    def draw_efficiency_map(self, res=5, cut='', show=True):
        cut_string = TCut(cut) + self.Cut.CutStrings['tracks']
        cut_string = self.Cut.generate_special_cut(excluded=['masks', 'fiducial', 'rhit']) if cut == 'all' else cut_string
        p = TProfile2D('p_em', 'Efficiency Map {d}'.format(d=self.DiamondName), *self.plots.get_global_bins(res=res))
        self.tree.Draw('(n_hits[{r}]>0)*100:diam{r1}_track_y:diam{r1}_track_x>>p_em'.format(r=self.Dut, r1=self.Dut - 3), cut_string, 'goff')
        set_statbox(entries=4, opt=1000000010, x=.81)
        self.format_histo(p, x_tit='Track x [cm]', y_tit='Track y [cm]', z_tit='Efficiency [%]', y_off=1.4, z_off=1.5)
        self.save_histo(p, 'Efficiency Map', show, lm=.13, rm=.17, draw_opt='colz')

    def draw_track_occupancy(self, cut='', show=True, res=2):
        cut_string = self.Cut.all_cut if cut is None else TCut(cut) + self.Cut.CutStrings['tracks']
        h = TH2D('h_to', 'Track Occupancy {d}'.format(d=self.DiamondName), *self.plots.get_global_bins(res=res))
        self.tree.Draw('diam{nr}_track_y:diam{nr}_track_x>>h_to'.format(nr=self.Dut - 3), cut_string, 'goff')
        set_statbox(entries=4, opt=1000000010, x=.81)
        self.format_histo(h, x_tit='Track x [cm]', y_tit='Track y [cm]', z_tit='Number of Entries', y_off=1.4, z_off=1.5)
        self.save_histo(h, 'TrackOccupancy', show, lm=.12, rm=.17, draw_opt='colz')

    def draw_trigphase_offset(self, cut=None, show=True):
        cut_string = deepcopy(self.Cut.all_cut) if cut is None else TCut(cut)
        h = TH1I('h_tp', 'Trigger Phase Offset', 19, -9, 10)
        self.tree.Draw('trigger_phase[1] - trigger_phase[0]>>h_tp', cut_string, 'goff')
        set_statbox(entries=4, opt=1000000010, y=0.88)
        self.format_histo(h, x_tit='Trigger Phase', y_tit='Number of Entries', y_off=1.8, fill_color=self.FillColor, ndiv=20)
        self.save_histo(h, 'TriggerPhase', show, lm=.16)

    def draw_hit_eff_vs_trigphase(self, roc=None, show=True, n=1e9, start=0):
        roc = self.Dut if roc is None else roc
        x = range(10)
        cut_string = self.Cut.generate_special_cut(excluded=['masks', 'rhit', 'trigger_phase'])
        y = [self.fit_hit_efficiency(roc=roc, show=False, cut=cut_string + TCut('trigger_phase[1]=={v}'.format(v=i)), n=n, start=start).Parameter(0) for i in xrange(10)]
        y = [0 if i is None else i for i in y]
        gr = self.make_tgrapherrors('gr_etp', 'Efficiency per Trigger Phase', x=x, y=y)
        gr.GetXaxis().SetLimits(-1, 10)
        self.format_histo(gr, fill_color=self.FillColor, x_tit='Trigger Phase', y_tit='Efficiency [%]', y_off=1.4)
        self.save_histo(gr, 'EffVsTrigPhase', show, draw_opt='ba', lm=.13)

    def induce_threshold(self):
        # todo
        pass

    def find_landau(self, aver=10, m1=2500, m2=5000, s1=500, s2=1600, pix=([18, 51], [18, 53])):
        seed = self.draw_pulse_height_disto(show=False, sup_zero=False, pix=pix)
        h = deepcopy(seed)
        m_range = range(m1, m2 + 1, 100)
        s_range = range(s1, s2 + 1, 50)
        p = TProfile2D('g_fl', 'Find Landau', len(m_range) - 1, m_range[0], m_range[-1], len(s_range) - 1, s_range[0], s_range[-1])
        self.start_pbar(len(m_range) * len(s_range) * aver)
        i = 0
        r_min, r_max = .5, 1.5
        for _ in xrange(aver):

            for m in m_range:
                for s in s_range:
                    i += 1
                    self.ProgressBar.update(i)
                    if r_min < m / 4. / s < r_max:
                        diff = self.model_landau(seed, h, m, s, show=False, thresh=True)
                        p.Fill(m, s, diff)
        self.ProgressBar.finish()
        set_statbox(entries=8, opt=1000000010, x=.82)
        self.format_histo(p, x_tit='MPV [e]', y_tit='Sigma [e]', z_tit='#chi^{2} to Seed Function', y_off=1.7, z_off=1.3)
        self.draw_histo(p, draw_opt='colz', lm=.13, rm=0.16)
        self.draw_ms_ratios(r_min, r_max, m1, m2, s1, s2)
        self.save_plots('FindLandau')
        self.find_working_point(p)

    def draw_ms_ratios(self, r_min, r_max, m1, m2, s1, s2, step=.1):
        ratios = arange(r_min, r_max + step, step)
        off = .01
        for i, ratio in enumerate(ratios):
            cut = TCutG('ms{n}'.format(n=i), 2, array([0, 20000], 'd'), array([0, 20000 / ratio / 4]))
            x_pos = m1 if m1 / ratio / 4 > s1 else 4 * ratio * s1
            y_pos = m1 / ratio / 4 if m1 / ratio / 4 > s1 else s1
            self.draw_tlatex(x_pos + off * (m2 - m1), y_pos + off * (s2 - s1), text='{0:3.1f}'.format(ratio), size=.02, align=11)
            cut.Draw('same')
            self.RootObjects.append(cut)

    @staticmethod
    def find_working_point(h):
        ps = [h.ProfileY(), h.ProfileX()]
        fits = [TF1('f{n}'.format(n=i), 'pol2', 0, 10000) for i in xrange(2)]
        for fit, p in zip(fits, ps):
            p.Fit(fit, 'qs0')
        mins = [fit.GetMinimumX() for fit in fits]
        print mins, mins[1] / mins[0] / 4

    def model_landau(self, seed=None, h=None, m=10000, s=1000, show=True, thresh=False, pix=([18, 51], [18, 53])):
        # seed = self.draw_pulse_height_disto(show=False, sup_zero=False, col=col) if seed is None else seed
        seed = self.draw_pulse_height_disto(show=False, sup_zero=False,  pix=pix) if seed is None else seed
        h = deepcopy(seed) if h is None else h
        h.SetName('h_ml')
        n = seed.GetEntries()
        h.Reset()
        # thresholds = self.get_thresholds(cols=col, vcal=False).values()
        thresholds = [47 * i for i in [150, 160, 170]]
        for _ in xrange(int(n)):
            v = gRandom.Landau(m, s)
            threshold = thresholds[int(gRandom.Rndm() * len(thresholds))]
            h.Fill(v if v > threshold else 0 if thresh else v)
            # h.Fill(v if v > 11421 else 0 if thresh else v)
        diff = mean([(h.GetBinContent(i) - seed.GetBinContent(i)) ** 2 for i in xrange(h.GetNbinsX()) if i is not h.FindBin(0)])
        seed.SetFillColor(2)
        h.SetFillColor(self.FillColor)
        seed.SetFillStyle(3017)
        h.SetFillStyle(3018)
        l = self.make_legend(y2=.76)
        l.AddEntry(h, 'Simulation', 'f')
        l.AddEntry(seed, 'Original', 'f')
        self.draw_histo(h, show=show, l=l)
        self.draw_histo(seed, show=show, draw_opt='same', canvas=gROOT.GetListOfCanvases()[-1])
        seed.Draw('same')
        return diff

    def landau_vid(self, save=False, mpv=5000, sigma=820):
        h = self.draw_pulse_height_disto(sup_zero=False)
        h.GetYaxis().SetRangeUser(0, 2500)
        zero_bin = h.FindBin(0)
        zeros = int(h.GetBinContent(zero_bin))
        entries = int(h.GetEntries())
        print entries
        c = gROOT.GetListOfCanvases()[-1]
        thresholds = self.get_thresholds(vcal=False)
        for i in xrange(entries):
            h.SetBinContent(zero_bin, zeros)
            v = gRandom.Landau(mpv, sigma)
            threshold = thresholds[int(gRandom.Rndm() * len(thresholds))]
            if v < threshold:
                h.Fill(v)
                zeros -= 1
            if i % 100 == 0:
                c.Update()
                c.Modified()
            if i % 100 == 0 and save:
                self.save_canvas(c, name='l{i:04d}'.format(i=i), show=False, print_names=False)

    def draw_correlation(self, plane1=2, plane2=None, mode='y', chi2=1, show=True, start=0, evts=1000000000):
        plane2 = self.Dut if plane2 is None else plane2
        h = TH2D('h_pc', 'Plane Correlation', *self.plots.get_global_bins(mode=mode, res=sqrt(12)))
        draw_var = 'cluster_{m}pos_tel'.format(m=mode)
        cut_string = 'clusters_per_plane[{p1}]==1&&clusters_per_plane[{p2}]==1&&cluster_plane=={{p}}'.format(p1=plane1, p2=plane2)
        n = self.tree.Draw(draw_var, TCut(cut_string.format(p=plane1)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff', evts, start)
        x1 = [self.tree.GetV1()[i] for i in xrange(n)]
        n = self.tree.Draw(draw_var, TCut(cut_string.format(p=plane2)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff', evts, start)
        x2 = [self.tree.GetV1()[i] for i in xrange(n)]
        for i, j in zip(x1, x2):
            h.Fill(i, j)
        self.log_info('Correlation Factor: {f:4.3f}'.format(f=h.GetCorrelationFactor()))
        self.format_histo(h, x_tit='{m} Plane {p}'.format(p=plane1, m=mode), y_tit='{m} Plane {p}'.format(p=plane2, m=mode), y_off=1.5, stats=0, z_tit='Number of Entries', z_off=1.5)
        self.save_histo(h, 'PlaneCorrelation{m}{p1}{p2}'.format(m=mode.title(), p1=plane1, p2=plane2), show,  lm=.13, draw_opt='colz', rm=.17)

    def draw_alignment(self, plane1=2, plane2=None, mode='y', binning=5000, chi2=1, show=True, vs_time=True):
        plane2 = self.Dut if plane2 is None else plane2
        picklepath = self.make_pickle_path('Alignment', run=self.RunNumber, suf='{m}_{p1}{p2}_{b}_{t}'.format(m=mode, p1=plane1, p2=plane2, b=binning, t='Time' if vs_time else 'EvtNr'))

        def func():
            start = self.log_info('Checking for alignment between plane {p1} and {p2} ... '.format(p1=plane1, p2=plane2), next_line=False)
            self.set_bin_size(binning)
            h = TH3D('h_pa', 'pa', len(self.time_binning) - 1, array([t / 1000. for t in self.time_binning], 'd'), *self.plots.get_global_bins(res=sqrt(12), mode=mode, arrays=True))
            if not vs_time:
                h = TH3D('h_pae', 'pa', len(self.binning) - 1, array(self.binning, 'd'), *self.plots.get_global_bins(res=sqrt(12), mode=mode, arrays=True))
            draw_var = 'cluster_{m}pos_tel'.format(m=mode)
            cut_string = 'clusters_per_plane[{p1}]==1&&clusters_per_plane[{p2}]==1&&cluster_plane=={{p}}'.format(p1=plane1, p2=plane2)
            n = self.tree.Draw(draw_var, TCut(cut_string.format(p=plane1)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff')
            y = [self.tree.GetV1()[i] for i in xrange(n)]
            n = self.tree.Draw('{v}:{t}'.format(v=draw_var, t='time' if vs_time else 'event_number'), TCut(cut_string.format(p=plane2)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff')
            t = [self.tree.GetV1()[i] for i in xrange(n)]
            x = [self.tree.GetV2()[i] / 1000. if vs_time else self.tree.GetV2()[i] for i in xrange(n)]
            for i, j, k in zip(x, y, t):
                h.Fill(i, j, k)
            g = self.make_tgrapherrors('g_pa', 'Plane Correlation {p1} {p2}'.format(p1=plane1, p2=plane2), marker_size=.5)
            for ibin in xrange(h.GetNbinsX() - 1):
                h.GetXaxis().SetRange(ibin, ibin + 1)
                p = h.Project3D('yz')
                g.SetPoint(ibin, h.GetXaxis().GetBinCenter(ibin), p.GetCorrelationFactor())
            if vs_time:
                set_time_axis(g, off=self.run.startTime / 1000 + 3600)
            self.format_histo(g, x_tit='Time [hh::mm]' if vs_time else 'Event Number', y_tit='Correlation Factor', y_off=1.5, y_range=[0, 1])
            self.add_info(start)
            return g

        gr = self.do_pickle(picklepath, func)
        self.save_histo(gr, 'PixelAligment', show, draw_opt='alp', lm=.13, prnt=show)
        return gr

    def check_alignment(self):
        gr = self.draw_alignment(show=False)
        values = [gr.GetY()[i_ev] for i_ev in xrange(gr.GetN())]
        mean_, sigma = calc_mean(values)
        if mean_ < .4:
            log_warning('Planes are not correlated!')
        elif sigma > .05:
            log_warning('Large fluctuations in correlation!')
        return mean_ > .3

    def draw_event_offsets(self, evnts=1000, start=0, pnts=None, rnge=1, show=True):
        pnts = int(ceil(self.run.n_entries / evnts)) if pnts is None else pnts
        n_graphs = 2 * rnge + 1
        graphs = [self.make_tgrapherrors('g_to{i}'.format(i=i), '', color=self.get_color(), marker_size=.5) for i in xrange(n_graphs)]
        l = self.make_legend(x1=.8, y2=.5, nentries=n_graphs - 1)
        self.start_pbar(pnts * evnts - start)
        p1 = self.Dut
        p2 = 2
        for k in xrange(pnts):
            x, y = OrderedDict(), OrderedDict()
            n = self.tree.Draw('plane:row:event_number', 'n_hits[2]==1 || n_hits[{p1}]==1'.format(p1=p1), 'goff', evnts, start + k * evnts)
            planes = [int(self.tree.GetV1()[i]) for i in xrange(n)]
            rows = [int(self.tree.GetV2()[i]) for i in xrange(n)]
            nrs = Counter([int(self.tree.GetV3()[i]) for i in xrange(n)])
            for ev, size in sorted(nrs.iteritems()):
                self.ProgressBar.update(ev + 1 - start)
                plane = [planes.pop(0) for _ in xrange(size)]
                row = [rows.pop(0) for _ in xrange(size)]
                if plane.count(p1) == 1:
                    x[ev] = row[plane.index(p1)]
                if plane.count(p2) == 1:
                    y[ev] = row[plane.index(p2)]
            xts = [[] for _ in xrange(n_graphs)]
            yts = [[] for _ in xrange(n_graphs)]
            for i in xrange(-rnge, rnge + 1):
                j = i + rnge
                for ev, row in x.iteritems():
                    if ev + i in y:
                        xts[j].append(row)
                        yts[j].append(y[ev + i])
                corr = corrcoef(xts[j], yts[j])[0][1]
                graphs[j].SetPoint(k, k * evnts + start, corr)
        self.ProgressBar.finish()
        mg = TMultiGraph('m_eo', 'Correlations')
        for i, gr in enumerate(graphs):
            l.AddEntry(gr, 'offset: {i}'.format(i=i - rnge), 'pl')
            mg.Add(gr, 'pl')
        self.format_histo(mg, x_tit='Event Number', y_tit='Correlation Factor', y_off=1.5, y_range=[-.2, 1], draw_first=True)
        self.save_histo(mg, 'EventOffsets', show, draw_opt='ap', l=l, lm=.13)
        self.reset_colors()

    def draw_cluster_size(self, show=True, cut=''):
        h = TH1I('h_cs', 'Cluster Size {d}'.format(d=self.DiamondName), 10, 0, 10)
        self.tree.Draw('clusters_per_plane[{d}]>>h_cs'.format(d=self.Dut), TCut(cut), 'goff')
        set_statbox(only_entries=True)
        self.format_histo(h, x_tit='Cluster Size', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.save_histo(h, 'ClusterSize', logy=True)

    # ==========================================================================

    def __placeholder(self):
        pass


if __name__ == '__main__':
    st = time()
    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=489, type=int, help='Run to be analysed {e.g.334}')
    parser.add_argument('dut', nargs='?', default=1, type=int, help='Number of the DUT to analyse (either 1, 2 or 3)')
    # parser.add_argument('-t', '--doTelescope', action='store_true', dest='doTelscp', default=False, help='set with -t or with --doTelescope to do telescope analysis')
    parser.add_argument('-d', '--doDUTs', action='store_true', dest='doDUTs', default=False, help='set with -d or with --doDUTs to do DUTs analysis')
    parser.add_argument('-u', '--doCutDist', action='store_true', dest='doCutDist', default=False,
                        help='set with -u or with --doCutDist to do Cuts distributions on selected devices (DUTs and/or telescope)')
    parser.add_argument('-e', '--doResolution', action='store_true', dest='doResolution', default=False,
                        help='set with -e or with --doResolution to do resolution analysis on selected devices (DUTs and/or telescope)')
    parser.add_argument('-c', '--doCutAna', action='store_true', dest='doCutAna', default=False, help='set with -c or with --doCutAna to do Cuts analysis on selected devices (DUTs and/or telescope)')
    parser.add_argument('-o', '--doOccupancy', action='store_true', dest='doOccupancy', default=False,
                        help='set with -o or with --doOccupancy to do occupancies (hit maps) on selected devices (DUTs and/or telescope)')
    parser.add_argument('-x', '--doCorrel', action='store_true', dest='doCorrel', default=False, help='set with -x or with --doCorrel to do correlations between important planes')
    parser.add_argument('-g', '--doCharge', action='store_true', dest='doCharge', default=False,
                        help='set with -g or with --doCharge to do pulse height analysis on selected devices (DUTs and/or telescope)')
    parser.add_argument('-p', '--progBar', action='store_true', dest='progBar', default=False, help='show progress bar')
    parser.add_argument('-v', '--verb', action='store_true', dest='verb', default=True, help='show verbose')
    parser.add_argument('-a', '--analyse', action='store_true', dest='doAna', default=False, help='run the whole analysis with the options entered')
    parser.add_argument('-tc', '--testcampaign', nargs='?', default='')
    parser.add_argument('-t', '--tree', default=True, action='store_false')

    args = parser.parse_args()

    # doTelscp = bool(args.doTelscp)
    doDUTs = bool(args.doDUTs)
    doCutDist = bool(args.doCutDist)
    doResolution = bool(args.doResolution)
    doCutAna = bool(args.doCutAna)
    doHitMap = bool(args.doOccupancy)
    doCorrel = bool(args.doCorrel)
    doPH = bool(args.doCharge)
    pbar = bool(args.progBar)
    verb = bool(args.verb)
    doAna = bool(args.doAna)

    command = 'Analysing run ' + str(args.run) + ' with:'
    # command += ' telescope,' if doTelscp else ' no telescope,'
    command += ' DUTs,' if doDUTs else ' no DUTs,'
    command += ' cuts distributions,' if doCutDist else ' no cuts distributions,'
    command += ' resolution analysis,' if doResolution else ' no resolution analysis,'
    command += ' cuts analysis,' if doCutAna else ' no cuts analysis,'
    command += ' hitmaps,' if doHitMap else ' no hitmpas,'
    command += ' correlations,' if doCorrel else ' no correlations,'
    command += ' pulse heights,' if doPH else ' no heights,'
    command += ' progress bar,' if pbar else ' no progress bar,'
    command += ' verbose' if pbar else ' no verbose,'
    command += ' with automatic analysis' if doAna else '. Start the Analysis by typing "z.do_analysis(doTelscp, doDUTs, doCutDist, doCutAna, doHitMap, doCorrel, doPH, pbar, verb)" when ready'
    command += '\n'

    tc = args.testcampaign if args.testcampaign.startswith('201') else None
    el = Elementary(tc)
    el.print_testcampaign()
    el.print_banner('STARTING PIXEL-ANALYSIS OF RUN {0}'.format(args.run))
    # print command
    print
    z = PixAnalysis(args.run, args.dut, args.tree, args.verb)
    z.print_elapsed_time(st, 'Instantiation')

    if doAna:
        print 'Starting automatic analysis...'
        z.do_analysis(doCutDist, doResolution, doCutAna, doHitMap, doPH)
        print 'Finished automatic analysis :)'
