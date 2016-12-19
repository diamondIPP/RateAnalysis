# ==============================================
# IMPORTS
# ==============================================
from ROOT import TH2D, TH1D, gROOT, TFormula
from TelescopeAnalysis import Analysis
# from CurrentInfo import Currents
from argparse import ArgumentParser
from time import time
from copy import deepcopy
import progressbar
# from CutPix import CutPix

__author__ = 'DA'


# ==============================================
# MAIN CLASS
# ==============================================
class SignalPixAnalysis(Analysis):
    def __init__(self, run, dut=1, verbose=False):

        Analysis.__init__(self, run, verbose=verbose)

        # main
        self.RunNumber = run
        self.DiamondName = self.load_diamond_name(dut)
        self.Bias = self.run.bias[dut - 1]
        self.save_dir = '{dia}/{run}/'.format(run=str(self.run_number).zfill(3), dia=self.DiamondName)

        # stuff
        self.roc_diam1 = 4
        self.roc_diam2 = 5 if self.TESTCAMPAIGN != '201610' else -1
        self.roc_si = 6 if self.TESTCAMPAIGN != '201610' else 5
        self.roc_tel = [0, 1, 2, 3]
        self.plots.roc_tel = self.roc_tel
        self.plots.roc_d1 = self.roc_diam1
        self.plots.roc_d2 = self.roc_diam2
        self.plots.roc_si = self.roc_si
        self.plots.save_dir = self.save_dir
        self.devices = {'tel': [], 'dut': []}
        self.Cut.reset_cut_dicts()
        self.Cut.do_cuts()
        self.plot_settings = self.plots.plot_settings
        self.stuff = []

    def __del__(self):
        for c in gROOT.GetListOfCanvases():
            c.Close()
        for lst in self.histos.itervalues():
            if not type(lst) is list:
                lst = [lst]
            for obj in lst:
                self.del_rootobj(obj)

    def change_roc_ids(self, tel0, tel1, tel2, tel3, dia1, dia2, si):
        """
        Changes the id of the ROCs. By default telescope are [0, 1, 2, 3], dia1 = 4, dia2 = 5, si = 6
        :param tel0:
        :param tel1:
        :param tel2:
        :param tel3:
        :param dia1:
        :param dia2:
        :param si:
        :return: No return
        """
        self.roc_tel = [tel0, tel1, tel2, tel3]
        self.roc_diam1 = dia1
        self.roc_diam2 = dia2
        self.roc_si = si
        self.dut_names = {self.roc_diam1: self.run.diamond_names[0], self.roc_diam2: self.run.diamond_names[3], self.roc_si: 'Si'}
        self.set_cuts_rocs()

    def set_cuts_rocs(self):
        """
        Changes the roc id's in the Cut instance of CutPix
        :return:
        """
        self.Cut.roc_tel = self.roc_tel
        self.Cut.roc_diam1 = self.roc_diam1
        self.Cut.roc_diam2 = self.roc_diam2
        self.Cut.roc_si = self.roc_si
        self.Cut.dut_names = self.dut_names
        self.Cut.reset_cuts_dicts()

    def add_telescope_device(self):
        self.devices['tel'] = self.roc_tel

    def add_duts_device(self):
        self.devices['dut'] = [self.roc_diam1, self.roc_diam2, self.roc_si] if self.TESTCAMPAIGN != '201610' else [self.roc_diam1, self.roc_si]

    def do_analysis(self, do_tlscp=False, do_duts=True, do_cut_dist=False, do_res_ana=False, do_cut_ana=False, do_occupancy=True, do_correlations=False, do_pulse_height=True, show_progressBar=False, verbosity=False):
        """
        Does automatic analysis with the selected options
        :param do_tlscp: if true, analyses the telescope planes too (not implemented yet)
        :param do_duts: if true, analyses the DUTs
        :param do_cut_dist: if true, calculates and saves the plot and root file of the cut distributions (done in the CutPix method)
        :param do_res_ana: if true, calculates and saves the plots and root files of the resolution analysis
        :param do_cut_ana: if true, does the cut analysis on the selected devices (currently only on DUTs)
        :param do_occupancy: if true, does occupancy analysis on selected devices
        :param do_correlations: if true, does correlation analysis
        :param do_pulse_height: if true, does pulse height analysis
        :param show_progressBar: if true, shows progress bar (not implmented)
        :param verbosity: if true, prints extra verbosity
        :return: smiley face when finished
        """
        gROOT.SetBatch(True)
        gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
        gROOT.SetBatch(False)
        TFormula.SetMaxima(1000000,10000,10000000)  # (1000,1000,1000)
        self.kmax = int(self.plots.plot_settings['num_diff_cluster_sizes'] + 1)
        self.deltaX = self.plots.plot_settings['deltaX']
        self.deltaY = self.plots.plot_settings['deltaY']

        if do_tlscp:
            self.add_telescope_device()
        if do_duts:
            self.add_duts_device()
        if do_cut_dist:
            self.Cut.do_cuts_distributions()
        if do_res_ana:
            self.Cut.do_res_analysis()
        if do_cut_ana:
            self.do_cuts_analysis(do_occupancy, do_pulse_height, True)

    def do_cuts_analysis(self, do_occupancy, do_pulse_height, normalize_ph_plots=False):
        """
        Calls do_cut_analysis in the Cut method
        :param do_occupancy:
        :param do_pulse_height:
        :param normalize_ph_plots:
        :return:
        """
        self.Cut.do_cuts_analysis(do_occupancy, do_pulse_height, normalize_ph_plots)

    def do_occupancy_roc(self, roc=4, cut='', histo=None):
        """
        Does the occupancy of a roc with the specified cut and it is saved on the given histogram. If none is given, it will
        create a histogram and return a deepcopy of it
        :param roc: roc number of the device (e.g. 4, 5)
        :param cut: name of the cut up to which the cuts should be applied (e.g. 'fiducial' would apply all the cuts up to fiducial cut)
        :param histo: histogram to save the occupancy map. If None, it will create its own
        :return: returns a deepcopy of the filled histogram after applying the specified cuts
        """
        if cut != '':
            if histo is None: name = 'hitmap_roc_{r}_{c}'.format(r=roc, c=cut)
            cut_string = 'plane=={r}&&{c}'.format(r=roc, c=self.Cut.cuts_hitmap_roc_incr[roc][self.Cut.dict_cuts[cut]])
            if self.verbose: print 'Doing occupancy for ROC', roc, 'upto', cut, 'cut'
        else:
            if histo is None: name = 'hitmap_roc_{r}'.format(r=roc)
            cut_string = 'plane=={r}'.format(r=roc)
            if self.verbose: print 'Doing occupancy for ROC', roc, 'without cuts'
        if histo is None:
            histo = TH2D(name, name, self.plot_settings['nBinCol']+1,
                         self.plot_settings['minCol']-(self.plot_settings['maxCol']-self.plot_settings['minCol'])/(
                             2*float(self.plot_settings['nBinCol'])), self.plot_settings['maxCol']+(
                             self.plot_settings['maxCol']-self.plot_settings['minCol'])/(
                             2*float(self.plot_settings['nBinCol'])), self.plot_settings['nBinRow']+1,
                         self.plot_settings['minRow']-(self.plot_settings['maxRow']-self.plot_settings['minRow'])/(
                             2*float(self.plot_settings['nBinRow'])), self.plot_settings['maxRow']+(
                             self.plot_settings['maxRow']-self.plot_settings['minRow'])/(
                             2*float(self.plot_settings['nBinRow'])))
        gROOT.SetBatch(True)
        name = histo.GetName()
        self.tree.Draw('row:col >> {n}'.format(n=name), cut_string, 'goff')
        gROOT.SetBatch(False)
        return deepcopy(histo)

    def do_pulse_height_event_roc(self, roc=4, num_clust=1, cut='', histo=None):
        """
        Does the pulse height analysis vs event for the specified roc for a specified number of cluster with a specified cut,
        on the given histogram. If none is given, it will create one. A deepcopy of the histogram is returned
        :param roc: roc number of the device (e.g. 4, 5)
        :param num_clust: number of pixels in the cluster to take into adcount
        :param cut: name of the cut up to which the cuts should be applied (e.g. 'fiducial' would apply all the cuts up to fiducial cut)
        :param histo: histogram to save the pulse height vs event. If None, it will create its own
        :return: a deepcopy of the filled histogram after applying the specified cuts
        """
        phbins = {self.roc_diam1: self.plot_settings['ph1DbinsD4'], self.roc_diam2: self.plot_settings['ph1DbinsD5'], self.roc_si: self.plot_settings['ph1DbinsSi']}
        phmin = {self.roc_diam1: self.plot_settings['ph1DminD4'], self.roc_diam2: self.plot_settings['ph1DminD5'], self.roc_si: self.plot_settings['ph1DminSi']}
        phmax = {self.roc_diam1: self.plot_settings['ph1DmaxD4'], self.roc_diam2: self.plot_settings['ph1DmaxD5'], self.roc_si: self.plot_settings['ph1DmaxSi']}
        phdelta = {self.roc_diam1: phmax[self.roc_diam1] - phmin[self.roc_diam1], self.roc_diam2: phmax[self.roc_diam2] - phmin[self.roc_diam2],
                   self.roc_si: phmax[self.roc_si] - phmin[self.roc_si]}
        if cut != '' and num_clust != 0:
            if histo is None: name = 'ph{n}_evt_roc{r}_{c}'.format(n=num_clust, r=roc, c=cut)
            cut_string = 'cluster_size_ROC{r}=={n}&&{cu}'.format(r=roc,n=num_clust,cu=self.Cut.cuts_pixelated_roc_incr[roc][self.Cut.dict_cuts[cut]])
            if self.verbose: print 'Doing pulse height', num_clust, 'pix cluster Vs. event for ROC', roc, 'upto', cut, 'cut'
        elif cut != '':
            if histo is None: name = 'ph_evt_roc{r}_{c}'.format(r=roc, c=cut)
            cut_string = '{cu}'.format(cu=self.Cut.cuts_pixelated_roc_incr[roc][self.Cut.dict_cuts[cut]])
            if self.verbose: print 'Doing pulse height Vs. event for ROC', roc, 'upto', cut, 'cut'
        elif num_clust != 0:
            if histo is None: name = 'ph{n}_evt_roc{r}'.format(n=num_clust, r=roc)
            cut_string = 'cluster_size_ROC{r}=={n}'.format(r=roc,n=num_clust)
            if self.verbose: print 'Doing pulse height', num_clust, 'pix cluster Vs. event for ROC', roc, 'without cuts'
        else:
            if histo is None: name = 'ph_evt_roc{r}'.format(r=roc)
            cut_string = ''
            if self.verbose: print 'Doing pulse height Vs. event for ROC', roc, 'without cuts'
        if histo is None:
            histo = TH2D(name, name, self.plot_settings['event_bins']+1, self.plot_settings['event_min']-(
                             self.plot_settings['event_max']-self.plot_settings['event_min'])/(
                             2*float(self.plot_settings['event_bins'])), self.plot_settings['event_max']+(
                             self.plot_settings['event_max']-self.plot_settings['event_min'])/(
                             2*float(self.plot_settings['event_bins'])), phbins[roc]+1, phmin[roc]-phdelta[roc]/(
                             2*float(phbins[roc])), phmax[roc]+phdelta[roc]/float(2*phbins[roc]))
        gROOT.SetBatch(True)
        name = histo.GetName()
        self.tree.Draw('charge_all_ROC{r}:event_number >> {n}'.format(r=roc, n=name), cut_string, 'goff')
        gROOT.SetBatch(False)
        return deepcopy(histo)

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
            if self.verbose: print 'Doing pulse height', num_clust,'pix cluster histogram for ROC', roc, 'without cuts'
        else:
            if histo is None: name = 'ph_roc{r}'.format(r=roc)
            if self.verbose: print 'Doing pulse height histogram for ROC', roc, 'without cuts'
        if histo is None:
            histo = TH1D(name, name, phbins[roc]+1, phmin[roc]-phdelta[roc]/(2*float(phbins[roc])), phmax[roc]+phdelta[roc]/float(2*phbins[roc]))
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
            if self.verbose: print 'Doing pulse height', num_clust,'pix cluster map for ROC', roc, 'without cuts'
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
        self.tree.Draw('charge_all_ROC{r}:10000*(residual_ROC{r}_Local_Y+cluster_pos_ROC{r}_Local_Y):10000*(residual_ROC{r}_Local_X+cluster_pos_ROC{r}_Local_X)>>{n}'.format(r=roc,n=name),
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
                ' ', progressbar.ETA()  #, DA: this two work great!
                # ' ', progressbar.AdaptiveETA(),
                # ' ', progressbar.AdaptiveTransferSpeed(),
                ]
            bar = progressbar.ProgressBar(widgets=widgets, max_value=self.num_devices) if do_tlscp else progressbar.ProgressBar(widgets=widgets, max_value=self.num_devices-len(self.roc_tel))
            bar.start()
        devini = 0 if do_tlscp else self.num_devices - len(self.roc_tel) + 1
        for iROC in xrange(devini, self.num_devices):
            if iROC not in self.roc_tel:
                if verbosity: self.print_banner('Analysing ROC {r}...'.format(r=iROC))
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}'.format(n=iROC)):
                    self.tree.Draw('row:col >> hitMapROC{n}'.format(n=iROC), 'plane == {n} && {mask}'.format(n=iROC, mask=self.Cut.mask_hitmap_roc[iROC].GetTitle()), 'goff')
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}_cuts'.format(n=iROC)):
                    self.tree.Draw('row:col >> hitMapROC{n}_cuts'.format(n=iROC), 'plane == {n} && {mask} && fidcut_hitmap_roc{n}'.format(n=iROC, mask=self.Cut.cuts_hitmap_roc_incr[iROC][self.Cut.num_cuts-1].GetTitle()), 'goff')
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}'.format(n=iROC)):
                    self.plots.save_individual_plots(self.plots.hitMap[iROC], self.plots.hitMap[iROC].GetName(), self.plots.hitMap[iROC].GetTitle(), self.Cut.fid_cut_hitmap_roc[iROC], 'colz', 0, self.save_dir, verbosity)
                if not self.plots.check_plot_existence(self.save_dir, 'c_hitMapROC{n}_cuts'.format(n=iROC)):
                    self.plots.save_individual_plots(self.plots.hitMap_cuts[iROC], self.plots.hitMap_cuts[iROC].GetName(), self.plots.hitMap_cuts[iROC].GetTitle(), self.Cut.fid_cut_hitmap_roc[iROC], 'colz', 0, self.save_dir, verbosity)
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
                ' ', progressbar.ETA()  #, DA: this two work great!
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
            exec("self.plots.clone_correlation_histograms(planeTel1, DUT1, verbosity, 1001, self.corrf_col_rocs_{rx}_{ry}, self.corrf_row_rocs_{rx}_{ry}, self.corrf_x_rocs_{rx}_{ry}, self.corrf_y_rocs_{rx}_{ry})".format(rx=planeTel1, ry=DUT1))
        if created_Tel2Dut2_col and created_Tel2Dut2_row and created_Tel2Dut2_x and created_Tel2Dut2_y:
            exec("self.plots.clone_correlation_histograms(planeTel2, DUT2, verbosity, 1001, self.corrf_col_rocs_{rx}_{ry}, self.corrf_row_rocs_{rx}_{ry}, self.corrf_x_rocs_{rx}_{ry}, self.corrf_y_rocs_{rx}_{ry})".format(rx=planeTel2, ry=DUT2))
        self.print_banner('Transposed correlation creation -> Done')
        self.print_banner('Main correlations creation -> Done', '%')

    def correlate_planes(self, rocx=1, rocy=4, var='col', verbosity=False):
        if not self.plots.check_plot_existence(self.save_dir, 'c_corr_{rx}_{ry}_{v}'.format(v=var, rx=rocx, ry=rocy)):
            if var is 'col' or var is 'row':
                self.tree.Draw('cluster_{v}_ROC{rx}:cluster_{v}_ROC{ry} >> corr_{rx}_{ry}_{v}'.format(v=var, rx=rocx, ry=rocy), self.Cut.cuts_pixelated_roc_incr[rocy][self.Cut.num_cuts-1], 'goff')
            elif var is 'x' or var is 'y':
                self.tree.Draw('cluster_pos_ROC{rx}_Telescope_{vv}*10:cluster_pos_ROC{ry}_Telescope_{vv}*10 >> corr_{rx}_{ry}_{v}'.format(v=var, rx=rocx, ry=rocy, vv=var.title()), self.Cut.cuts_pixelated_roc_incr[rocy][self.Cut.num_cuts-1], 'goff')
            gROOT.SetBatch(True)
            exec("self.fit_{x}_rocs_{rx}_{ry} = self.plots.correl_{x}[rocx][rocy].Fit('pol1', 'qsfm')".format(x=var, rx=rocx, ry=rocy))
            gROOT.SetBatch(False)
            exec("self.fit_{x}_rocs_{rx}_{ry}_angle = self.fit_{x}_rocs_{rx}_{ry}.Parameter(1)".format(x=var, rx=rocx, ry=rocy))
            exec("self.corrf_{x}_rocs_{rx}_{ry} = self.plots.correl_{x}[rocx][rocy].GetCorrelationFactor()".format(x=var, rx=rocx, ry=rocy))
            exec("self.plots.save_individual_plots(self.plots.correl_{v}[rocx][rocy], self.plots.correl_{v}[rocx][rocy].GetName(), self.plots.correl_{v}[rocx][rocy].GetTitle(), None, 'colz', 1000000011, self.save_dir, verbosity, 1001, self.corrf_{v}_rocs_{rx}_{ry})".format(v=var, rx=rocx, ry=rocy))
            return True
        else: return False

    def show_current(self, relative_time=True):
        self.Currents.draw_graphs(relative_time=relative_time)

    # ==========================================================================

    def __placeholder(self):
        pass

if __name__ == "__main__":
    st = time()
    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=392, type=int, help='Run to be analysed {e.g.334}')
    parser.add_argument('-t', '--doTelescope', action='store_true', dest='doTelscp', default=False, help='set with -t or with --doTelescope to do telescope analysis')
    parser.add_argument('-d', '--doDUTs', action='store_true', dest='doDUTs', default=False, help='set with -d or with --doDUTs to do DUTs analysis')
    parser.add_argument('-u', '--doCutDist', action='store_true', dest='doCutDist', default=False, help='set with -u or with --doCutDist to do Cuts distributions on selected devices (DUTs and/or telescope)')
    parser.add_argument('-e', '--doResolution', action='store_true', dest='doResolution', default=False, help='set with -e or with --doResolution to do resolution analysis on selected devices (DUTs and/or telescope)')
    parser.add_argument('-c', '--doCutAna', action='store_true', dest='doCutAna', default=False, help='set with -c or with --doCutAna to do Cuts analysis on selected devices (DUTs and/or telescope)')
    parser.add_argument('-o', '--doOccupancy', action='store_true', dest='doOccupancy', default=False, help='set with -o or with --doOccupancy to do occupancies (hit maps) on selected devices (DUTs and/or telescope)')
    parser.add_argument('-x', '--doCorrel', action='store_true', dest='doCorrel', default=False, help='set with -x or with --doCorrel to do correlations between important planes')
    parser.add_argument('-g', '--doCharge', action='store_true', dest='doCharge', default=False, help='set with -g or with --doCharge to do pulse height analysis on selected devices (DUTs and/or telescope)')
    parser.add_argument('-p', '--progBar', action='store_true', dest='progBar', default=False, help='show progress bar')
    parser.add_argument('-v', '--verb', action='store_true', dest='verb', default=False, help='show verbose')
    parser.add_argument('-a', '--analyse', action='store_true', dest='doAna', default=False, help='run the whole analysis with the options entered')

    args = parser.parse_args()
    run = int(args.run)
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

    command = '\nAnalysing run ' + str(run) + ' with:'
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

    z = SignalPixAnalysis(run)
    z.print_elapsed_time(st, 'Instantiation')

    if doAna:
        print 'Starting automatic analysis...'
        z.do_analysis(doTelscp, doDUTs, doCutDist, doResolution, doCutAna, doHitMap, doCorrel, doPH, pbar, verb)
        print 'Finished automatic analysis :)'
