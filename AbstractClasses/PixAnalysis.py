# ==============================================
# IMPORTS
# ==============================================
from ROOT import TH2D, TH1D, gROOT, TFormula, TCut, TH1I, TProfile, THStack, TProfile2D, TF1
from TelescopeAnalysis import Analysis
# from CurrentInfo import Currents
from argparse import ArgumentParser
from time import time
from copy import deepcopy
import progressbar
from CutPix import CutPix
from numpy import array
from Utils import *
from os.path import join as joinpath

__author__ = 'DA'


# ==============================================
# MAIN CLASS
# ==============================================
class SignalPixAnalysis(Analysis):
    def __init__(self, run, dut=1, verbose=False, binning=10000):

        Analysis.__init__(self, run, verbose=verbose, binning=binning)

        # main
        self.RunNumber = run
        self.DiamondName = self.load_diamond_name(dut)
        self.Bias = self.run.bias[dut - 1]
        self.Dut = dut + 3
        self.save_dir = '{dia}/{run}/'.format(run=str(self.run_number).zfill(3), dia=self.DiamondName)
        self.NRocs = self.load_n_rocs()

        # stuff
        self.plots.save_dir = self.save_dir

        # cuts
        self.Cut = CutPix(self, dut)
        self.Settings = self.plots.plot_settings

        # pulse height calibrations
        self.Fit = None
        self.Parameters = None
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

    def draw_occupancy(self, cut=None, show=True, fid=False, prnt=True, adc=None):
        """ Does the occupancy of a roc with the specified cut and it is saved on the given histogram. If none is given, it will create a histogram and return a deepcopy of it """
        cut_string = self.Cut.generate_special_cut(excluded='fiducial' if not fid else [], cluster=False) if cut is None else TCut(cut)
        cut_string += 'plane=={d}'.format(d=self.Dut)
        cut_string += self.Cut.add_adc_cut(adc)
        self.set_root_output(False)
        h = TH2D('h_oc', 'Occupancy {d}'.format(d=self.DiamondName), *self.Settings['2DBins'])
        self.tree.Draw('row:col >> {n}'.format(n='h_oc'), cut_string, 'goff')
        save_name = 'Occupancy{c}'.format(c=make_cut_string(cut, self.Cut.NCuts))
        set_statbox(x=.81, entries=8, opt=1000000010)
        self.format_histo(h, x_tit='col', y_tit='row', z_tit='Number of Entries', y_off=1.3, z_off=1.5)
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

    def draw_adc_disto(self, cut=None, show=True):
        h = TH1I('h_adc', 'ADC Distribution {d}'.format(d=self.DiamondName), 255, 0, 255)
        cut_string = self.Cut.HitMapCut if cut is None else TCut(cut)
        cut_string += 'plane == {n}'.format(n=self.Dut)

    def get_calibration_data(self):
        f = open(joinpath(self.run.converter.TrackingDir, 'calibration_lists', 'GKCalibrationList_Telescope{n}.txt'.format(n=self.run.converter.TelescopeID)))
        lines = f.readlines()
        f.close()
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
        if self.Fit is None:
            self.Fit = fit
            self.Parameters = params

    def draw_calibration_fit(self, col, row, show=True, roc=None):
        roc = self.Dut if roc is None else roc
        self.Fit.SetParameters(*self.Parameters[roc][col][row])
        self.format_histo(self.Fit, title='Calibration Fitf for Pix {c} {r}'.format(c=col, r=row), x_tit='vcal', y_tit='adc', y_off=1.4, color=632, lw=2)
        self.save_histo(self.Fit, 'CalFit{c}{r}'.format(c=col, r=row), show, lm=0.12)

        self.set_root_output(False)
        self.tree.Draw('adc>>h_adc', cut_string, 'goff')
        self.format_histo(h, x_tit='adc', y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor, stats=0)
        self.save_histo(h, 'ADCDisto', show, lm=0.13, logy=True)
        return h

    def draw_pulse_height_disto(self, cut=None, show=True, prnt=True, sup_zero=True):
        cut_string = self.Cut.all_cut if cut is None else TCut(cut)
        cut_string += 'charge_all_ROC{d}!=0'.format(d=self.Dut) if sup_zero else ''
        self.set_root_output(False)
        h = TH1D('h_phd', 'Pulse Height Distribution - {d}'.format(d=self.DiamondName), *self.Settings['phBinsD{n}'.format(n=self.Dut)])
        self.tree.Draw('charge_all_ROC{d}>>h_phd'.format(d=self.Dut), cut_string, 'goff')
        self.format_histo(h, x_tit='Pulse Height [e]', y_tit='Number of Entries', y_off=1.4, stats=0, fill_color=self.FillColor)
        self.save_histo(h, 'PulseHeightDisto{c}'.format(c=make_cut_string(cut, self.Cut.NCuts)), show, lm=.13, prnt=prnt)
        return h

    def draw_hit_efficiency(self, roc, show=True, save=True, cut=''):
        self.set_root_output(False)
        suffix = 'ROC {n}'.format(n=roc) if roc < 4 else self.load_diamond_name(roc - 3)
        h = TProfile('h_he', 'Hit Efficiency {s}'.format(s=suffix), int(self.run.n_entries / 5000), z.run.startTime / 1000, self.run.endTime / 1000.)
        cut_string = self.Cut.generate_special_cut(excluded=['fiducial', 'masks', 'rhit']) if cut == 'all' else TCut(cut)
        self.tree.Draw('(clusters_per_plane[{r}]>0)*100:time / 1000 >> h_he'.format(r=roc), cut_string, 'goff')
        set_time_axis(h, off=self.run.startTime / 1000 + 3600)
        self.format_histo(h, x_tit='Time [hh:mm]', y_tit='Efficiency [%]', y_off=1.4, ndiv=505, y_range=[-5, 105], stats=0)
        self.save_histo(h, 'HitEfficiencyROC{n}'.format(n=roc), show, lm=.13, save=save, gridy=True)
        return h

    def fit_hit_efficiency(self, roc, show=True, save=True, cut=''):
        pickle_path = self.make_pickle_path('Efficiency', run=self.RunNumber, suf='{r}{c}'.format(r=roc, c='_Cuts' if cut else ''))

        def func():
            set_statbox(y=.37, only_fit=True)
            h = self.draw_hit_efficiency(roc, show=False, save=False, cut=cut)
            self.format_histo(h, stats=1, name='Fit Results')
            fit = h.Fit('pol0', 'qs')
            self.save_histo(h, 'HitEfficiencyROC{n}Fit'.format(n=roc), show, lm=.13, save=save, gridy=True)
            return fit

        fit_res = func() if show else None
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

    command = '\nAnalysing run ' + str(args.run) + ' with:'
    command = command + ' telescope,' if doTelscp else command + ' no telescope,'
    command = command + ' DUTs,' if doDUTs else command + ' no DUTs,'
    command = command + ' cuts distributions,' if doCutDist else command + ' no cuts distributions,'
    command = command + ' resolution analysis,' if doResolution else command + ' no resolution analysis,'
    command = command + ' cuts analysis,' if doCutAna else command + ' no cuts analysis,'
    command = command + ' hitmaps,' if doHitMap else command + ' no hitmpas,'
    command = command + ' correlations,' if doCorrel else command + ' no correlations,'
    command = command + ' pulse heights,' if doPH else command + ' no heights,'
    command = command + ' progress bar,' if pbar else command + ' no progress bar,'
    command = command + ' verbose' if pbar else command + ' no verbose,'
    command = command + ' with automatic analysis' if doAna else command + '. Start the Analysis by typing "z.do_analysis(doTelscp, doDUTs, doCutDist, doCutAna, doHitMap, doCorrel, doPH, pbar, verb)" when ready'
    command = command + '\n'

    print command

    z = SignalPixAnalysis(args.run, args.dut, args.verb)
    z.print_elapsed_time(st, 'Instantiation')

    if doAna:
        print 'Starting automatic analysis...'
        z.do_analysis(doTelscp, doDUTs, doCutDist, doResolution, doCutAna, doHitMap, doCorrel, doPH, pbar, verb)
        print 'Finished automatic analysis :)'
