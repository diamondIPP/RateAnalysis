# ==============================================
# IMPORTS
# ==============================================
from ROOT import TH2D, TH1D, gROOT, TFormula, TCut, TH1I, TProfile, THStack, TProfile2D, TF1, TGraph, TPie, gRandom, TH3D
from TelescopeAnalysis import Analysis
# from CurrentInfo import Currents
from argparse import ArgumentParser
from time import time
from copy import deepcopy
from CutPix import CutPix
from Elementary import Elementary
from numpy import array
from Utils import *
from os.path import join as joinpath
from collections import OrderedDict

__author__ = 'DA'


# ==============================================
# MAIN CLASS
# ==============================================
class PixAnalysis(Analysis):
    def __init__(self, run, dut=1, verbose=False, binning=10000):

        Analysis.__init__(self, run, verbose=verbose, binning=binning)

        # main
        self.RunNumber = run
        self.DiamondName = self.load_diamond_name(dut)
        self.Bias = self.run.bias[dut - 1]
        self.Dut = dut + 3
        self.save_dir = '{dia}/{run}/'.format(run=str(self.RunNumber).zfill(3), dia=self.DiamondName)
        self.NRocs = self.load_n_rocs()

        # stuff
        self.plots.save_dir = self.save_dir

        # cuts
        self.Settings = self.plots.Settings
        self.Cut = CutPix(self, dut)

        # alignment
        self.IsAligned = self.check_alignment(show=False)

        # pulse height calibrations
        self.Fit = None
        self.Parameters = None
        self.Vcals = None
        self.Points = None
        self.get_calibration_data()

        self.stuff = []

    def __del__(self):
        for c in gROOT.GetListOfCanvases():
            c.Close()
        for lst in self.histos.itervalues():
            if not type(lst) is list:
                lst = [lst]
            for obj in lst:
                self.del_rootobj(obj)

    def load_diamond_name(self, dut):
        assert dut in [1, 2, 3], 'You have to choose either dut 1, 2 or 3'
        return self.run.diamond_names[dut - 1]

    def load_n_rocs(self):
        if self.ana_config_parser.has_option('BASIC', 'n_rocs'):
            return self.ana_config_parser.getint('BASIC', 'n_rocs')
        else:
            return 7

    def print_info(self, event):
        # self.log_info('Event {e}:'.format(e=event))
        self.tree.GetEntry(event)
        adc = [self.tree.adc[i] for i, j in enumerate(self.tree.plane) if (j == 4 and ord(self.tree.n_tracks) == 1)]
        cluster_charge = [self.tree.cluster_charge[i] for i, j in enumerate(self.tree.cluster_plane) if (j == 4 and ord(self.tree.n_tracks) == 1)]
        if adc:
            print event, adc, cluster_charge

    def print_infos(self, i, j):
        for n in xrange(i, j):
            self.print_info(n)

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

    def draw_occupancy(self, cut=None, show=True, fid=False, prnt=True, adc=None, roc=None, tel_coods=False):
        """ Does the occupancy of a roc with the specified cut and it is saved on the given histogram. If none is given, it will create a histogram and return a deepcopy of it """
        roc = self.Dut if roc is None else roc
        cut_string = self.Cut.generate_special_cut(excluded='fiducial' if not fid else [], cluster=False) if cut is None else TCut(cut)
        cut_string += '{p}=={d}'.format(d=roc, p='cluster_plane' if tel_coods else 'plane')
        cut_string += self.Cut.add_adc_cut(adc)
        self.set_root_output(False)
        h = TH2D('h_oc', 'Occupancy {d}'.format(d=self.DiamondName), *(self.Settings['2DBins'] if not tel_coods else self.plots.get_global_bins(sqrt(12))))
        draw_var = 'row:col' if not tel_coods else 'cluster_ypos_tel:cluster_xpos_tel'
        self.tree.Draw('{d} >> {n}'.format(n='h_oc', d=draw_var), cut_string, 'goff')
        save_name = 'Occupancy{c}'.format(c=make_cut_string(cut, self.Cut.NCuts))
        set_statbox(x=.81, entries=8, opt=1000000010)
        xtit, ytit = ('col', 'row') if not tel_coods else ('x [cm]', 'y [cm]')
        self.format_histo(h, x_tit=xtit, y_tit=ytit, z_tit='Number of Entries', y_off=1.3, z_off=1.5)
        self.save_histo(h, save_name, show, rm=.17, lm=.13, draw_opt='colz', prnt=prnt)
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
        params = [[[0 for _ in xrange(self.Settings['nRows'])] for _ in xrange(self.Settings['nRows'])] for _ in xrange(self.NRocs)]
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
        calib_files = [joinpath(self.run.converter.TrackingDir, lines[0].strip('./\n'), 'phCalibration_{e}'.format(e=line.strip('\n')[-6:])) for line in lines[5:]]
        vcals = None
        points = [[[0 for _ in xrange(self.Settings['nRows'])] for _ in xrange(self.Settings['nRows'])] for _ in xrange(self.NRocs)]
        for roc, file_name in enumerate(calib_files, 4):
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
            h.Fill(col, row, thresh)
        self.format_histo(h, x_tit='col', y_tit='row', z_tit='Artificial Treshold [vcal]', y_off=1.3, z_off=1.5, stats=0)
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

    def check_adc(self):
        for i in xrange(10000, 12000):
            self.tree.GetEntry(i)
            if 4 in self.tree.plane:
                ind = list(self.tree.plane).index(4)
                if self.tree.adc[ind]:
                    print i, self.tree.adc[ind], list(self.tree.charge_all_ROC4)

    def draw_pulse_height_disto(self, cut=None, show=True, prnt=True, sup_zero=True, col=None, pix=None, roc=None, vcal=False):
        roc = self.Dut if roc is None else roc
        cut_string = deepcopy(self.Cut.generate_special_cut(excluded=['fiducial', 'trigger_phase'])) if cut is None else TCut(cut)
        cut_string += 'cluster_charge>0'.format(d=self.Dut) if sup_zero else ''
        cut_string += 'cluster_col=={c}'.format(c=col, d=self.Dut) if col is not None else ''
        cut_string += 'cluster_col=={c}&&cluster_row=={r}'.format(c=pix[0], r=pix[1], d=self.Dut) if pix is not None else ''
        cut_string += 'cluster_plane=={r}'.format(r=roc)
        self.set_root_output(False)
        h = TH1D('h_phd', 'Pulse Height Distribution - {d}'.format(d=self.DiamondName), *self.Settings['phBins' if not vcal else 'vcalBins'])
        self.tree.Draw('cluster_charge{v}>>h_phd'.format(d=self.Dut, v='/47.5 + 427.4/47.5' if vcal else ''), cut_string, 'goff')
        set_statbox(entries=8, opt=1000000010, x=.92)
        self.format_histo(h, x_tit='Pulse Height [{u}]'.format(u='vcal' if vcal else 'e'), y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor)
        self.save_histo(h, 'PulseHeightDisto{c}'.format(c=make_cut_string(cut, self.Cut.NCuts)), show, lm=.13, prnt=prnt, rm=.06)
        return h

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

    def draw_hit_efficiency(self, roc=None, show=True, save=True, cut=''):
        roc = self.Dut if roc is None else roc
        self.set_root_output(False)
        suffix = 'ROC {n}'.format(n=roc) if roc < 4 else self.load_diamond_name(roc - 3)
        h = TProfile('h_he', 'Hit Efficiency {s}'.format(s=suffix), int(self.run.n_entries / 5000), z.run.startTime / 1000, self.run.endTime / 1000.)
        cut_string = self.Cut.generate_special_cut(excluded=['masks', 'rhit']) if cut == 'all' else TCut(cut)
        self.tree.Draw('(n_hits[{r}]>0)*100:time / 1000 >> h_he'.format(r=roc), cut_string, 'goff')
        set_time_axis(h, off=self.run.startTime / 1000 + 3600)
        self.format_histo(h, x_tit='Time [hh:mm]', y_tit='Efficiency [%]', y_off=1.4, ndiv=505, y_range=[-5, 105], stats=0)
        self.save_histo(h, 'HitEfficiencyROC{n}'.format(n=roc), show, lm=.13, save=save, gridy=True)
        return h

    def fit_hit_efficiency(self, roc=None, show=True, save=True, cut=''):
        pickle_path = self.make_pickle_path('Efficiency', run=self.RunNumber, suf='{r}{c}'.format(r=roc, c='_Cuts' if cut else ''))

        def func():
            set_statbox(y=.37, only_fit=True, entries=1.5)
            h = self.draw_hit_efficiency(roc, show=False, save=False, cut=cut)
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

    def draw_hit_eff_vs_trigphase(self, roc=None, show=True):
        roc = self.Dut if roc is None else roc
        x = range(10)
        cut_string = self.Cut.generate_special_cut(excluded=['masks', 'rhit', 'trigger_phase'])
        y = [self.fit_hit_efficiency(roc=roc, show=False, cut=cut_string + TCut('trigger_phase[1]=={v}'.format(v=i))).Parameter(0) for i in xrange(10)]
        y = [0 if i is None else i for i in y]
        gr = self.make_tgrapherrors('gr_etp', 'Efficiency per Trigger Phase', x=x, y=y)
        gr.GetXaxis().SetLimits(-1, 10)
        self.format_histo(gr, fill_color=self.FillColor, x_tit='Trigger Phase', y_tit='Efficiency [%]', y_off=1.4)
        self.save_histo(gr, 'EffVsTrigPhase', show, draw_opt='ba', lm=.13)


    def draw_correlation(self, plane1=1, plane2=None, mode='y', chi2=1, show=True):
        plane2 = self.Dut if plane2 is None else plane2
        h = TH2D('h_pc', 'Plane Correlation', *self.plots.get_global_bins(mode=mode, res=sqrt(12)))
        draw_var = 'cluster_{m}pos_tel'.format(m=mode)
        cut_string = 'clusters_per_plane[{p1}]==1&&clusters_per_plane[{p2}]==1&&cluster_plane=={{p}}'.format(p1=plane1, p2=plane2)
        n = self.tree.Draw(draw_var, TCut(cut_string.format(p=plane1)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff')
        x1 = [self.tree.GetV1()[i] for i in xrange(n)]
        n = self.tree.Draw(draw_var, TCut(cut_string.format(p=plane2)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff')
        x2 = [self.tree.GetV1()[i] for i in xrange(n)]
        for i, j in zip(x1, x2):
            h.Fill(i, j)
        self.log_info('Correlation Factor: {f:4.3f}'.format(f=h.GetCorrelationFactor()))
        self.format_histo(h, x_tit='{m} Plane {p}'.format(p=plane1, m=mode), y_tit='{m} Plane {p}'.format(p=plane2, m=mode), y_off=1.5, stats=0, z_tit='Number of Entries', z_off=1.5)
        self.save_histo(h, 'PlaneCorrelation{m}{p1}{p2}'.format(m=mode.title(), p1=plane1, p2=plane2), show,  lm=.13, draw_opt='colz', rm=.17)

    def check_alignment(self, plane1=1, plane2=None, mode='y', binning=5000, chi2=1, show=True):
        plane2 = self.Dut if plane2 is None else plane2
        picklepath = self.make_pickle_path('Alignment', run=self.RunNumber, suf='{m}_{p1}{p2}_{b}'.format(m=mode, p1=plane1, p2=plane2, b=binning))

        def func():
            start = self.log_info('Checking for alignment between plane {p1} and {p2} ... '.format(p1=plane1, p2=plane2), next_line=False)
            self.set_bin_size(binning)
            h = TH3D('h_pa', 'pa', len(self.time_binning) - 1, array([t / 1000. for t in self.time_binning], 'd'), *self.plots.get_global_bins(res=sqrt(12), mode=mode, arrays=True))
            draw_var = 'cluster_{m}pos_tel'.format(m=mode)
            cut_string = 'clusters_per_plane[{p1}]==1&&clusters_per_plane[{p2}]==1&&cluster_plane=={{p}}'.format(p1=plane1, p2=plane2)
            n = self.tree.Draw(draw_var, TCut(cut_string.format(p=plane1)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff')
            x1 = [self.tree.GetV1()[i] for i in xrange(n)]
            n = self.tree.Draw('{v}:time'.format(v=draw_var), TCut(cut_string.format(p=plane2)) + TCut(self.Cut.generate_chi2(mode, chi2)), 'goff')
            x2 = [self.tree.GetV1()[i] for i in xrange(n)]
            t = [self.tree.GetV2()[i] / 1000. for i in xrange(n)]
            for i, j, k in zip(x1, x2, t):
                h.Fill(k, i, j)
            g = self.make_tgrapherrors('g_pa', 'Plane Correlation {p1} {p2}'.format(p1=plane1, p2=plane2), marker_size=.5)
            for ibin in xrange(h.GetNbinsX() - 1):
                h.GetXaxis().SetRange(ibin, ibin + 1)
                p = h.Project3D('yz')
                g.SetPoint(ibin, h.GetXaxis().GetBinCenter(ibin), p.GetCorrelationFactor())
            set_time_axis(g, off=self.run.startTime / 1000 + 3600)
            self.format_histo(g, x_tit='Time [hh::mm]', y_tit='Correlation Factor', y_off=1.5, y_range=[0, 1])
            self.add_info('Done', start)
            return g

        gr = self.do_pickle(picklepath, func)
        self.save_histo(gr, 'PixelAligment', show, draw_opt='alp', lm=.13, prnt=show)
        values = [gr.GetY()[i_ev] for i_ev in xrange(gr.GetN())]
        mean_, sigma = calc_mean(values)
        if sigma > .05:
            log_warning('Large fluctuations in correlation!')
        if mean_ < .4:
            log_warning('Planes are not correlated!')
        return mean_ > .3


    def do_pulse_height_roc(self, roc=4, num_clust=1, cut='', histoevent=None, histo=None):
        """
        Does the pulse height extraction for the specified roc with the specified cluster size and applying the given
        cummulative cuts. The results are saved in the given histograms (histoevent and histo). If no histogram is given,
        it will create its own. Returns a dictionary with deepcopies of the 2D histogram ph vs event and the 1D histogram
        of ph
        :param roc: roc number of the device (e.g. 4, 5)
        :param num_clust: number of pixels in the cluster to take into adcount
        :param cut: name of the cut up to which the cuts should be applied (e.g. 'fiducial' would apply all the cuts up to fiducial cut)
        :param histoevent: histogram to save the pulse height vs event. If None, it will create its own
        :param histo: histogram to save the 1D pulse height. If None, it will create its own
        :return: a dictionary with deepcopies of the 2D histogram ph vs event, and the 1D histogram of the ph
        """
        phbins = {self.roc_diam1: self.plot_settings['ph1DbinsD4'], self.roc_diam2: self.plot_settings['ph1DbinsD5'], self.roc_si: self.plot_settings['ph1DbinsSi']}
        phmin = {self.roc_diam1: self.plot_settings['ph1DminD4'], self.roc_diam2: self.plot_settings['ph1DminD5'], self.roc_si: self.plot_settings['ph1DminSi']}
        phmax = {self.roc_diam1: self.plot_settings['ph1DmaxD4'], self.roc_diam2: self.plot_settings['ph1DmaxD5'], self.roc_si: self.plot_settings['ph1DmaxSi']}
        phdelta = {self.roc_diam1: phmax[self.roc_diam1] - phmin[self.roc_diam1], self.roc_diam2: phmax[self.roc_diam2] - phmin[self.roc_diam2],
                   self.roc_si: phmax[self.roc_si] - phmin[self.roc_si]}
        event_histo = self.do_pulse_height_event_roc(roc, num_clust, cut, histoevent)
        if cut != '' and num_clust != 0:
            if histo is None: name = 'ph{n}_roc{r}_{c}'.format(n=num_clust, r=roc, c=cut)
            if self.verbose: print 'Doing pulse height', num_clust, 'pix cluster histogram for ROC', roc, 'upto', cut, 'cut'
        elif cut != '':
            if histo is None: name = 'ph_roc{r}_{c}'.format(r=roc, c=cut)
            if self.verbose: print 'Doing pulse height histogram for ROC', roc, 'upto', cut, 'cut'
        elif num_clust != 0:
            if histo is None: name = 'ph{n}_roc{r}'.format(n=num_clust, r=roc)
            if self.verbose: print 'Doing pulse height', num_clust, 'pix cluster histogram for ROC', roc, 'without cuts'
        else:
            if histo is None: name = 'ph_roc{r}'.format(r=roc)
            if self.verbose: print 'Doing pulse height histogram for ROC', roc, 'without cuts'
        if histo is None:
            histo = TH1D(name, name, phbins[roc] + 1, phmin[roc] - phdelta[roc] / (2 * float(phbins[roc])), phmax[roc] + phdelta[roc] / float(2 * phbins[roc]))
        gROOT.SetBatch(True)
        name = histo.GetName()
        histoevent.ProjectionY(name, 0, -1, 'e')
        gROOT.SetBatch(False)
        return {'event_histo': deepcopy(histoevent), 'histo': deepcopy(histo)}

    def do_pulse_height_roc_map(self, roc=4, num_clust=1, cut='', histo=None):
        """
        Does the ph map of the specified roc with the given cluster size after applying the specified cuts. It is saved on
        the histogram given. If none is given, it will create its own. a deepcopy of the 2D average map is returned
        :param roc: roc number of the device (e.g. 4, 5)
        :param num_clust: number of pixels in the cluster to take into adcount
        :param cut: name of the cut up to which the cuts should be applied (e.g. 'fiducial' would apply all the cuts up to fiducial cut)
        :param histo: histogram to save the 2D average pulse height. If None, it will create its own
        :return: a deepcopy of the 2D average histogram of the ph
        """
        if cut != '' and num_clust != 0:
            if histo is None:
                name = 'ph{n}_map_roc{r}_{c}'.format(n=num_clust, r=roc, c=cut)
                ZTitle = 'ph {n} pix clusters(e)'.format(n=num_clust)
            cut_string = 'cluster_size_ROC{r}=={n}&&{c}'.format(r=roc, n=num_clust, c=self.Cut.cuts_pixelated_roc_incr[roc][self.Cut.dict_cuts[cut]])
            if self.verbose: print 'Doing pulse height', num_clust, 'pix cluster map for ROC', roc, 'upto', cut, 'cut'
        elif cut != '':
            if histo is None:
                name = 'ph_map_roc{r}_{c}'.format(r=roc, c=cut)
                ZTitle = 'ph all pix clusters(e)'
            cut_string = '{c}'.format(c=self.Cut.cuts_pixelated_roc_incr[roc][self.Cut.dict_cuts[cut]])
            if self.verbose: print 'Doing pulse height map for ROC', roc, 'upto', cut, 'cut'
        elif num_clust != 0:
            if histo is None:
                name = 'ph{n}_map_roc{r}'.format(n=num_clust, r=roc)
                ZTitle = 'ph {n} pix clusters(e)'.format(n=num_clust)
            cut_string = 'cluster_size_ROC{r}=={n}'.format(r=roc, n=num_clust)
            if self.verbose: print 'Doing pulse height', num_clust, 'pix cluster map for ROC', roc, 'without cuts'
        else:
            if histo is None:
                name = 'ph_map_roc{r}'.format(r=roc)
                ZTitle = 'ph all pix clusters(e)'
            cut_string = ''
            if self.verbose: print 'Doing pulse height map for ROC', roc, 'without cuts'
        if histo is None:
            histo = self.plots.create_2D_profile('spatial', name, name, 'x(um)', 'y(um)', ZTitle, 'auto', -1)
        gROOT.SetBatch(True)
        name = histo.GetName()
        self.tree.Draw('charge_all_ROC{r}:10000*(residual_ROC{r}_Local_Y+cluster_pos_ROC{r}_Local_Y):10000*(residual_ROC{r}_Local_X+cluster_pos_ROC{r}_Local_X)>>{n}'.format(r=roc, n=name),
                       cut_string, 'goff prof')
        gROOT.SetBatch(False)
        return deepcopy(histo)

    def fill_occupancy(self, show_progressBar=False, do_tlscp=False, verbosity=False):
        # for i in xrange(len(self.plane)):
        #     if 0 <= self.plane[i] < self.num_devices:
        #         self.plots.hitMap[self.plane[i]].Fill(self.col[i], self.row[i])
        self.print_banner('Starting Occupancy Analysis...', '%')
        if show_progressBar:
            widgets = [
                progressbar.Percentage(),
                ' ', progressbar.Bar(marker='>'),
                ' ', progressbar.Timer(),
                ' ', progressbar.ETA()  # , DA: this two work great!
                # ' ', progressbar.AdaptiveETA(),
                # ' ', progressbar.AdaptiveTransferSpeed(),
            ]
            bar = progressbar.ProgressBar(widgets=widgets, max_value=self.num_devices) if do_tlscp else progressbar.ProgressBar(widgets=widgets, max_value=self.num_devices - len(self.roc_tel))
            bar.start()
        devini = 0 if do_tlscp else self.num_devices - len(self.roc_tel) + 1
        for iROC in xrange(devini, self.num_devices):
            if iROC not in self.roc_tel:
                if verbosity: self.print_banner('Analysing ROC {r}...'.format(r=iROC))
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}'.format(n=iROC)):
                    self.tree.Draw('row:col >> hitMapROC{n}'.format(n=iROC), 'plane == {n} && {mask}'.format(n=iROC, mask=self.Cut.mask_hitmap_roc[iROC].GetTitle()), 'goff')
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}_cuts'.format(n=iROC)):
                    self.tree.Draw('row:col >> hitMapROC{n}_cuts'.format(n=iROC),
                                   'plane == {n} && {mask} && fidcut_hitmap_roc{n}'.format(n=iROC, mask=self.Cut.cuts_hitmap_roc_incr[iROC][self.Cut.num_cuts - 1].GetTitle()), 'goff')
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}'.format(n=iROC)):
                    self.plots.save_individual_plots(self.plots.hitMap[iROC], self.plots.hitMap[iROC].GetName(), self.plots.hitMap[iROC].GetTitle(), self.Cut.fid_cut_hitmap_roc[iROC], 'colz', 0,
                                                     self.save_dir, verbosity)
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}_cuts'.format(n=iROC)):
                    self.plots.save_individual_plots(self.plots.hitMap_cuts[iROC], self.plots.hitMap_cuts[iROC].GetName(), self.plots.hitMap_cuts[iROC].GetTitle(), self.Cut.fid_cut_hitmap_roc[iROC],
                                                     'colz', 0, self.save_dir, verbosity)
            elif do_tlscp:
                self.tree.Draw('row:col >> hitMapROC{n}'.format(n=iROC), 'plane == {n}'.format(n=iROC), 'goff')
                # self.tree.Draw('row:col >> hitMapROC{n}_cuts'.format(n=iROC), 'plane == {n}'.format(n=iROC), 'goff')
            if show_progressBar: bar.update(iROC + 1 - devini)
            if verbosity: self.print_banner('ROC {r} analysis -> Done'.format(r=iROC))
        if show_progressBar: bar.finish()
        self.print_banner('Occupancy Analysis -> Done', '%')

    def do_correlations(self, planeTel1=1, planeTel2=2, DUT1=4, DUT2=6, pbar=False, verbosity=False):
        self.print_banner('Creating main correlations ...', '%')
        if pbar:
            widgets = [
                progressbar.Percentage(),
                ' ', progressbar.Bar(marker='>'),
                ' ', progressbar.Timer(),
                ' ', progressbar.ETA()  # , DA: this two work great!
                # ' ', progressbar.AdaptiveETA(),
                # ' ', progressbar.AdaptiveTransferSpeed(),
            ]
            bar = progressbar.ProgressBar(widgets=widgets, max_value=4)
            bar.start()
        created_Tel1Dut1_col = self.correlate_planes(planeTel1, DUT1, 'col', verbosity)
        created_Tel1Dut1_x = self.correlate_planes(planeTel1, DUT1, 'x', verbosity)
        if pbar: bar.update(1)
        created_Tel1Dut1_row = self.correlate_planes(planeTel1, DUT1, 'row', verbosity)
        created_Tel1Dut1_y = self.correlate_planes(planeTel1, DUT1, 'y', verbosity)
        if pbar: bar.update(2)
        created_Tel2Dut2_col = self.correlate_planes(planeTel2, DUT2, 'col', verbosity)
        created_Tel2Dut2_x = self.correlate_planes(planeTel2, DUT2, 'x', verbosity)
        if pbar: bar.update(3)
        created_Tel2Dut2_row = self.correlate_planes(planeTel2, DUT2, 'row', verbosity)
        created_Tel2Dut2_y = self.correlate_planes(planeTel2, DUT2, 'y', verbosity)
        if pbar: bar.update(4)
        if pbar: bar.finish()
        self.print_banner('Creating transposed correlations ...')
        if created_Tel1Dut1_col and created_Tel1Dut1_row and created_Tel1Dut1_x and created_Tel1Dut1_y:
            exec (
            "self.plots.clone_correlation_histograms(planeTel1, DUT1, verbosity, 1001, self.corrf_col_rocs_{rx}_{ry}, self.corrf_row_rocs_{rx}_{ry}, self.corrf_x_rocs_{rx}_{ry}, self.corrf_y_rocs_{rx}_{ry})".format(
                rx=planeTel1, ry=DUT1))
        if created_Tel2Dut2_col and created_Tel2Dut2_row and created_Tel2Dut2_x and created_Tel2Dut2_y:
            exec (
            "self.plots.clone_correlation_histograms(planeTel2, DUT2, verbosity, 1001, self.corrf_col_rocs_{rx}_{ry}, self.corrf_row_rocs_{rx}_{ry}, self.corrf_x_rocs_{rx}_{ry}, self.corrf_y_rocs_{rx}_{ry})".format(
                rx=planeTel2, ry=DUT2))
        self.print_banner('Transposed correlation creation -> Done')
        self.print_banner('Main correlations creation -> Done', '%')

    def correlate_planes(self, rocx=1, rocy=4, var='col', verbosity=False):
        if not self.plots.check_plot_existence(self.save_dir, 'c_corr_{rx}_{ry}_{v}'.format(v=var, rx=rocx, ry=rocy)):
            if var is 'col' or var is 'row':
                self.tree.Draw('cluster_{v}_ROC{rx}:cluster_{v}_ROC{ry} >> corr_{rx}_{ry}_{v}'.format(v=var, rx=rocx, ry=rocy), self.Cut.cuts_pixelated_roc_incr[rocy][self.Cut.num_cuts - 1], 'goff')
            elif var is 'x' or var is 'y':
                self.tree.Draw('cluster_pos_ROC{rx}_Telescope_{vv}*10:cluster_pos_ROC{ry}_Telescope_{vv}*10 >> corr_{rx}_{ry}_{v}'.format(v=var, rx=rocx, ry=rocy, vv=var.title()),
                               self.Cut.cuts_pixelated_roc_incr[rocy][self.Cut.num_cuts - 1], 'goff')
            gROOT.SetBatch(True)
            exec ("self.fit_{x}_rocs_{rx}_{ry} = self.plots.correl_{x}[rocx][rocy].Fit('pol1', 'qsfm')".format(x=var, rx=rocx, ry=rocy))
            gROOT.SetBatch(False)
            exec ("self.fit_{x}_rocs_{rx}_{ry}_angle = self.fit_{x}_rocs_{rx}_{ry}.Parameter(1)".format(x=var, rx=rocx, ry=rocy))
            exec ("self.corrf_{x}_rocs_{rx}_{ry} = self.plots.correl_{x}[rocx][rocy].GetCorrelationFactor()".format(x=var, rx=rocx, ry=rocy))
            exec (
            "self.plots.save_individual_plots(self.plots.correl_{v}[rocx][rocy], self.plots.correl_{v}[rocx][rocy].GetName(), self.plots.correl_{v}[rocx][rocy].GetTitle(), None, 'colz', 1000000011, self.save_dir, verbosity, 1001, self.corrf_{v}_rocs_{rx}_{ry})".format(
                v=var, rx=rocx, ry=rocy))
            return True
        else:
            return False

    def show_current(self, relative_time=True):
        self.Currents.draw_graphs(relative_time=relative_time)

    # ==========================================================================

    def __placeholder(self):
        pass


if __name__ == '__main__':
    st = time()
    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=489, type=int, help='Run to be analysed {e.g.334}')
    parser.add_argument('dut', nargs='?', default=1, type=int, help='Number of the DUT to analyse (either 1, 2 or 3)')
    parser.add_argument('-t', '--doTelescope', action='store_true', dest='doTelscp', default=False, help='set with -t or with --doTelescope to do telescope analysis')
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

    args = parser.parse_args()

    doTelscp = bool(args.doTelscp)
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
    command += ' telescope,' if doTelscp else ' no telescope,'
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
    z = PixAnalysis(args.run, args.dut, args.verb)
    z.print_elapsed_time(st, 'Instantiation')

    if doAna:
        print 'Starting automatic analysis...'
        z.do_analysis(doCutDist, doResolution, doCutAna, doHitMap, doPH)
        print 'Finished automatic analysis :)'
