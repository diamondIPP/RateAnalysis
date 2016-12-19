import os

import sys
from numpy import array
from ROOT import TCut, gROOT, TH1F, kRed, TCutG, gDirectory, kBlue, TH2D, TH1D, kGreen
from math import ceil
from Cut import Cut
from json import loads

class CutPix(Cut):
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, analysis, dut=0):
        Cut.__init__(self, analysis, skip=True)
        self.__dict__.update(analysis.Cut.__dict__)

        self.Dut = dut + 4

        self.plot_settings = self.analysis.plots.plot_settings

        self.cuts_hitmap_roc = {}  # each cut separately
        self.cuts_pixelated_roc = {}  # each cut separately
        self.cuts_hitmap_roc_incr = {}  # each cut incremental. last position used for all cuts
        self.cuts_pixelated_roc_incr = {}  # each cut incremental. last position used for all cuts

        self.events_after_cuts = {}
        self.events_after_cuts_roc = {}
        self.events_after_cuts_incr = {}
        self.events_after_cuts_roc_incr = {}

        self.load_pixel_config()

        self.plots = self.analysis.plots

    def load_run_config(self):
        return self.load_run_configs(self.RunNumber)

    def calculate_initial_events(self):
        gROOT.SetBatch(1)
        self.events_after_cuts['None'] = self.analysis.tree.GetEntries()
        for iROC in self.duts_list:
            self.analysis.tree.Draw('time>>temp0', 'plane[{r}]'.format(r=iROC), 'goff')
            bla = gDirectory.Get('temp0')
            self.events_after_cuts_roc['None'] = bla.GetEntries()

        self.analysis.tree.Draw('time>>temp0', 'plane[{')

    def do_cuts(self):
        """
        Calculates or gets the cut strings to apply in the analysis for each of the cuts. Each of these cuts generation,
        fills the self.cuts_hitmap_roc, self.cuts_hitmap_roc_incr, self.cuts_pixelated_roc, self.cuts_pixelated_roc_incr
        :return:
        """
        # generate cut strings
        self.print_banner('Generating Cut strings...')
        print 'The following cuts will be implemented:', self.cut_names
        if 'ini_fin' in self.cut_names:
            self.generate_ini_fin_cuts()
        if 'beam' in self.cut_names:
            self.generate_beam_interruption_cut()
        if 'tracks' in self.cut_names:
            self.generate_tracks_cut()
        if 'hit' in self.cut_names:
            self.generate_hit_cut()
        if 'masks' in self.cut_names:
            self.generate_masks()
        if 'fiducial' in self.cut_names:
            self.generate_fid_cuts()
        if 'chi2x' in self.cut_names or 'chi2y' in self.cut_names:
            self.generate_chi2_cuts()
        if 'anglex' in self.cut_names or 'angley' in self.cut_names:
            self.generate_angle_cuts()
        if 'rhit' in self.cut_names:
            self.generate_rhit_cuts()
        # self.gen_incr_vect_cuts()
        self.cuts_done = True
        self.print_banner('Finished generating Cut stringss')

    def list_duts(self):
        """
        creates a list with the roc numbers of the DUTs
        :return:
        """
        self.duts_list = [self.roc_diam1, self.roc_diam2, self.roc_si] if self.TESTCAMPAIGN != '201610' else [self.roc_diam1, self.roc_si]

    def reset_cuts_dicts(self):
        """
        Resets the lists, dictionaries and the numbering of the DUTs when they are changed in analysis. This method should be called by PixAnalysis
        :return:
        """
        self.list_duts()
        self.cuts_hitmap_roc = {iROC: {} for iROC in self.duts_list}  # each cut separately
        self.cuts_pixelated_roc = {iROC: {} for iROC in self.duts_list}  # each cut separately
        self.cuts_hitmap_roc_incr = {iROC: {} for iROC in self.duts_list}  # each cut incremental. last position used for all cuts
        self.cuts_pixelated_roc_incr = {iROC: {} for iROC in self.duts_list}  # each cut incremental. last position used for all cuts
        self.num_cuts = 0
        self.events_after_cuts = {}
        self.events_after_cuts_roc = {iROC: {} for iROC in self.duts_list}
        self.events_after_cuts_incr = {}
        self.events_after_cuts_roc_incr = {iROC: {} for iROC in self.duts_list}

    def do_cuts_distributions(self):
        """
        Does the cuts distribution for beam interruption, chi2x, chi2y, anglex, angley and rhit for each of the ROCs to be analysed
        :return:
        """
        self.print_banner('Doing cuts distributions...')
        self.do_beam_distribution()
        self.do_chi2_distributions()
        self.do_angle_distributions()
        self.do_rhit_distribution()
        self.print_banner('Finished with distribution cuts')

    def do_res_analysis(self):
        """
        Calculates and saves the plots of res Y vs res X, rhit vs res x, rhit vs res y, chi2 vs resx, chi2 vs resy,
        chi2x vs resx, chi2y vs resx, resx vs Y predicted hit position, resy vs X predicted hit position
        :return:
        """
        self.h_resy_resx = {}
        self.h_rhit_resx = {}
        self.h_rhit_resy = {}
        self.h_chi2_resx = {}
        self.h_chi2_resy = {}
        self.h_chi2x_resx = {}
        self.h_chi2y_resx = {}
        self.h_resx_hitposy = {}
        self.h_resy_hitposx = {}
        self.print_banner('Doing resolution plots...')
        if self.verbose: print 'Res_Y Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resy_resx[iroc] = TH2D('h_resy_resx_roc{r}'.format(r=iroc), 'h_resy_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 21, -1050, 1050)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*residual_ROC{r}_Local_X>>h_resy_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*residual_ROC{r}_Local_X>>h_resy_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_resy_resx[iroc], 'Res_X(um)', 'Res_y(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resy_resx[iroc], 'h_resy_resx_roc{r}'.format(r=iroc), 'Res_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Rhit Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_resx[iroc] = TH2D('h_rhit_resx_roc{r}'.format(r=iroc), 'h_rhit_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 101, -0.5, 100.5)
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_X>>h_rhit_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_X>>h_rhit_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_rhit_resx[iroc], 'Res_X(um)', 'R_Hit(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_rhit_resx[iroc], 'h_rhit_resx_roc{r}'.format(r=iroc), 'R_Hit Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Rhit Vs Res_Y...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_resy[iroc] = TH2D('h_rhit_resy_roc{r}'.format(r=iroc), 'h_rhit_resy_roc{r}'.format(r=iroc), 21, -1050, 1050, 101, -0.5, 100.5)
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_Y>>h_rhit_resy_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_Y>>h_rhit_resy_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_rhit_resy[iroc], 'Res_Y(um)', 'R_Hit(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_rhit_resy[iroc], 'h_rhit_resy_roc{r}'.format(r=iroc), 'R_Hit Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2 Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2_resx[iroc] = TH2D('h_chi2_resx_roc{r}'.format(r=iroc), 'h_chi2_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_X>>h_chi2_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_X>>h_chi2_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2_resx[iroc], 'Res_X(um)', 'Chi2', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2_resx[iroc], 'h_chi2_resx_roc{r}'.format(r=iroc), 'Chi2 Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2 Vs Res_Y...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2_resy[iroc] = TH2D('h_chi2_resy_roc{r}'.format(r=iroc), 'h_chi2_resy_roc{r}'.format(r=iroc), 21, -1050, 1050, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_Y>>h_chi2_resy_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_Y>>h_chi2_resy_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2_resy[iroc], 'Res_Y(um)', 'Chi2', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2_resy[iroc], 'h_chi2_resy_roc{r}'.format(r=iroc), 'Chi2 Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2_X Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2x_resx[iroc] = TH2D('h_chi2x_resx_roc{r}'.format(r=iroc), 'h_chi2x_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_x:10000*residual_ROC{r}_Local_X>>h_chi2x_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_x:10000*residual_ROC{r}_Local_X>>h_chi2x_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2x_resx[iroc], 'Res_X(um)', 'Chi2_X', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2x_resx[iroc], 'h_chi2x_resx_roc{r}'.format(r=iroc), 'Chi2_X Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2_Y Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2y_resx[iroc] = TH2D('h_chi2y_resx_roc{r}'.format(r=iroc), 'h_chi2y_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_y:10000*residual_ROC{r}_Local_X>>h_chi2y_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_y:10000*residual_ROC{r}_Local_X>>h_chi2y_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2y_resx[iroc], 'Res_X(um)', 'Chi2_Y', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2y_resx[iroc], 'h_chi2y_resx_roc{r}'.format(r=iroc), 'Chi2_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Res_X Vs Hit_Y...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resx_hitposy[iroc] = TH2D('h_resx_hitposy_roc{r}'.format(r=iroc), 'h_resx_hitposy_roc{r}'.format(r=iroc), 161, -4025, 4025, 21, -1575, 1575)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_X:10000*(residual_ROC{r}_Local_Y+cluster_pos_ROC{r}_Local_Y)>>h_resx_hitposy_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('10000*residual_ROC{r}_Local_X:10000*(residual_ROC{r}_Local_Y+cluster_pos_ROC{r}_Local_Y)>>h_resx_hitposy_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_resx_hitposy[iroc], 'Hit_Y(um)', 'Res_X(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resx_hitposy[iroc], 'h_resx_hitposy_roc{r}'.format(r=iroc), 'Hit_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Res_Y Vs Hit_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resy_hitposx[iroc] = TH2D('h_resy_hitposx_roc{r}'.format(r=iroc), 'h_resy_hitposx_roc{r}'.format(r=iroc), 105, -3937.5, 3937.5, 21, -1050, 1050)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*(residual_ROC{r}_Local_X+cluster_pos_ROC{r}_Local_X)>>h_resy_hitposx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*(residual_ROC{r}_Local_X+cluster_pos_ROC{r}_Local_X)>>h_resy_hitposx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_resy_hitposx[iroc], 'Hit_X(um)', 'Res_Y(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resy_hitposx[iroc], 'h_resy_hitposx_roc{r}'.format(r=iroc), 'Hit_X Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        self.print_banner('Finished with resolution plots')

    def do_occupancy_analysis(self):
        """
        Does the occupancy analysis for each roc applying each of the cummulative cuts and saving them on the corresponding
        dictionaries
        :return:
        """
        self.print_banner('Starting occupancy cut analysis...')
        self.h_hitmaps_cuts = {}
        maxz_hitmap = {iroc: -10000000 for iroc in self.duts_list}
        for iroc in self.duts_list:
            self.h_hitmaps_cuts[iroc] = {}
            for cut in self.cut_names:
                if self.verbose: print 'Analysing ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut),; sys.stdout.flush()
                self.h_hitmaps_cuts[iroc][cut] = TH2D('hitmap_roc{r}_{c}'.format(r=iroc,c=cut), 'hitmap_roc{r}_{c}'.
                                                      format(r=iroc,c=cut), self.plot_settings['nBinCol']+1,
                                                      self.plot_settings['minCol']-(
                                                          self.plot_settings['maxCol']-self.plot_settings['minCol'])/(
                                                          2*float(self.plot_settings['nBinCol'])),
                                                      self.plot_settings['maxCol']+(
                                                          self.plot_settings['maxCol']-self.plot_settings['minCol'])/(
                                                          2*float(self.plot_settings['nBinCol'])),
                                                      self.plot_settings['nBinRow']+1, self.plot_settings['minRow']-(
                                                          self.plot_settings['maxRow']-self.plot_settings['minRow'])/(
                                                          2*float(self.plot_settings['nBinRow'])),
                                                      self.plot_settings['maxRow']+(
                                                          self.plot_settings['maxRow']-self.plot_settings['minRow'])/(
                                                          2*float(self.plot_settings['nBinRow'])))
                self.h_hitmaps_cuts[iroc][cut] = self.analysis.do_occupancy_roc(iroc, cut, self.h_hitmaps_cuts[iroc][cut])
                if maxz_hitmap[iroc] < self.h_hitmaps_cuts[iroc][cut].GetBinContent(self.h_hitmaps_cuts[iroc][cut].GetMaximumBin()):
                    maxz_hitmap[iroc] = self.h_hitmaps_cuts[iroc][cut].GetBinContent(self.h_hitmaps_cuts[iroc][cut].GetMaximumBin())
                if self.verbose: print 'Done'
                if self.verbose: print 'Saving for ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut), ; sys.stdout.flush()
                self.plots.set_2D_options(self.h_hitmaps_cuts[iroc][cut], 'col', 'row', 'entries', max_val=maxz_hitmap[iroc])
                self.plots.save_individual_plots(self.h_hitmaps_cuts[iroc][cut], 'hitmap_roc{r}_{c}'.format(r=iroc,c=cut),
                                                 'hitmap_roc{r}_{c}'.format(r=iroc,c=cut), self.fid_cut_hitmap_roc[iroc],
                                                 'colz', 1, self.plots.save_dir+'/cuts', doLogZ=True)
                if self.verbose: print 'Done'
        self.print_banner('Finished occupancy cut analysis')

    def do_pulse_height_analysis(self, normalize_ph_plots=True):
        """
        does the pulse height analysis to all the DUTs for each cummulative cut. The histograms ph vs event, ph 2D average map,
        and ph 1D are saved on the respective dictionaries
        :param normalize_ph_plots: if true, the pulse height maps have the same limits in Z to compare them between each other
        :return:
        """
        self.print_banner('Starting pulse height cut analysis...')
        self.h_ph1_evt_cuts = {}
        self.h_ph1_cuts = {}
        self.h_ph2_evt_cuts = {}
        self.h_ph2_cuts = {}
        self.h_ph2_cuts = {}
        self.h_ph1_map_cuts = {}
        self.h_ph2_map_cuts = {}
        maxz_ph1 = -10000000
        maxz_ph2 = -10000000
        minz_ph1 = 10000000
        minz_ph2 = 10000000
        max_ph1_map = {iroc: -10000000 for iroc in self.duts_list}
        max_ph2_map = {iroc: -10000000 for iroc in self.duts_list}
        min_ph1_map = {iroc: 10000000 for iroc in self.duts_list}
        min_ph2_map = {iroc: 10000000 for iroc in self.duts_list}
        phbins = {self.roc_diam1: self.plot_settings['ph1DbinsD4'], self.roc_diam2: self.plot_settings['ph1DbinsD5'],
                  self.roc_si: self.plot_settings['ph1DbinsSi']}
        phmin = {self.roc_diam1: self.plot_settings['ph1DminD4'], self.roc_diam2: self.plot_settings['ph1DminD5'],
                 self.roc_si: self.plot_settings['ph1DminSi']}
        phmax = {self.roc_diam1: self.plot_settings['ph1DmaxD4'], self.roc_diam2: self.plot_settings['ph1DmaxD5'],
                 self.roc_si: self.plot_settings['ph1DmaxSi']}
        phdelta = {self.roc_diam1: phmax[self.roc_diam1] - phmin[self.roc_diam1],
                   self.roc_diam2: phmax[self.roc_diam2] - phmin[self.roc_diam2],
                   self.roc_si: phmax[self.roc_si] - phmin[self.roc_si]}
        for iroc in self.duts_list:
            self.h_ph1_evt_cuts[iroc] = {}
            self.h_ph1_cuts[iroc] = {}
            self.h_ph2_evt_cuts[iroc] = {}
            self.h_ph2_cuts[iroc] = {}
            self.h_ph1_map_cuts[iroc] = {}
            self.h_ph2_map_cuts[iroc] = {}
            for cut in self.cut_names:
                if self.verbose: print 'Analysing ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut),; sys.stdout.flush()
                self.h_ph1_map_cuts[iroc][cut] = self.plots.create_2D_profile('spatial',
                                                                              'ph1_map_roc{r}_{c}'.format(r=iroc,c=cut),
                                                                              'ph1_map_roc{r}_{c}'.format(r=iroc,c=cut),
                                                                              'x(um)', 'y(um)', 'ph 1 pix cluster(e)',
                                                                              'auto', -1)
                self.h_ph2_map_cuts[iroc][cut] = self.plots.create_2D_profile('spatial',
                                                                              'ph2_map_roc{r}_{c}'.format(r=iroc,c=cut),
                                                                              'ph2_map_roc{r}_{c}'.format(r=iroc,c=cut),
                                                                              'x(um)', 'y(um)', 'ph 2 pix cluster(e)',
                                                                              'auto', -1)
                self.h_ph1_evt_cuts[iroc][cut] = TH2D('ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut),
                                                      'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), self.plot_settings[
                                                          'event_bins']+1, self.plot_settings['event_min']-(
                                                          self.plot_settings['event_max']-self.plot_settings[
                                                              'event_min'])/(2*float(self.plot_settings['event_bins'])),
                                                      self.plot_settings['event_max']+(
                                                          self.plot_settings['event_max']-self.plot_settings[
                                                              'event_min'])/(2*float(self.plot_settings['event_bins'])),
                                                      phbins[iroc]+1, phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])),
                                                      phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.h_ph1_cuts[iroc][cut] = TH1D('ph1_roc{r}_{c}'.format(r=iroc,c=cut),
                                                  'ph1_roc{r}_{c}'.format(r=iroc,c=cut), phbins[iroc]+1,
                                                  phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])),
                                                  phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.h_ph2_evt_cuts[iroc][cut] = TH2D('ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut),
                                                      'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), self.plot_settings[
                                                          'event_bins']+1, self.plot_settings['event_min']-(
                                                          self.plot_settings['event_max']-self.plot_settings[
                                                              'event_min'])/(2*float(self.plot_settings['event_bins'])),
                                                      self.plot_settings['event_max']+(
                                                          self.plot_settings['event_max']-self.plot_settings[
                                                              'event_min'])/(2*float(self.plot_settings['event_bins'])),
                                                      phbins[iroc]+1, phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])),
                                                      phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.h_ph2_cuts[iroc][cut] = TH1D('ph2_roc{r}_{c}'.format(r=iroc,c=cut),
                                                  'ph2_roc{r}_{c}'.format(r=iroc,c=cut), phbins[iroc]+1,
                                                  phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])),
                                                  phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.h_ph1_map_cuts[iroc][cut] = self.analysis.do_pulse_height_roc_map(iroc, 1, cut, self.h_ph1_map_cuts[iroc][cut])
                self.h_ph2_map_cuts[iroc][cut] = self.analysis.do_pulse_height_roc_map(iroc, 2, cut, self.h_ph2_map_cuts[iroc][cut])
                tempPh1 = self.analysis.do_pulse_height_roc(iroc, 1, cut, self.h_ph1_evt_cuts[iroc][cut], self.h_ph1_cuts[iroc][cut])
                self.h_ph1_evt_cuts[iroc][cut] = tempPh1['event_histo']
                self.h_ph1_cuts[iroc][cut] = tempPh1['histo']
                tempPh2 = self.analysis.do_pulse_height_roc(iroc, 2, cut, self.h_ph2_evt_cuts[iroc][cut], self.h_ph2_cuts[iroc][cut])
                self.h_ph2_evt_cuts[iroc][cut] = tempPh2['event_histo']
                self.h_ph2_cuts[iroc][cut] = tempPh2['histo']
                if maxz_ph1 < self.h_ph1_evt_cuts[iroc][cut].GetBinContent(self.h_ph1_evt_cuts[iroc][cut].GetMaximumBin()):
                    maxz_ph1 = self.h_ph1_evt_cuts[iroc][cut].GetBinContent(self.h_ph1_evt_cuts[iroc][cut].GetMaximumBin())
                if (self.dict_cuts[cut] >= self.dict_cuts[self.cut_near_fiducial()]) and (
                            minz_ph1 > self.h_ph1_evt_cuts[iroc][cut].GetBinContent(
                            self.h_ph1_evt_cuts[iroc][cut].GetMinimumBin())):
                    minz_ph1 = self.h_ph1_evt_cuts[iroc][cut].GetBinContent(self.h_ph1_evt_cuts[iroc][cut].GetMinimumBin())
                if maxz_ph2 < self.h_ph2_evt_cuts[iroc][cut].GetBinContent(self.h_ph2_evt_cuts[iroc][cut].GetMaximumBin()):
                    maxz_ph2 = self.h_ph2_evt_cuts[iroc][cut].GetBinContent(self.h_ph2_evt_cuts[iroc][cut].GetMaximumBin())
                if (self.dict_cuts[cut] >= self.dict_cuts[self.cut_near_fiducial()]) and (
                            minz_ph2 > self.h_ph2_evt_cuts[iroc][cut].GetBinContent(
                            self.h_ph2_evt_cuts[iroc][cut].GetMinimumBin())):
                    minz_ph2 = self.h_ph2_evt_cuts[iroc][cut].GetBinContent(self.h_ph2_evt_cuts[iroc][cut].GetMinimumBin())
                if max_ph1_map[iroc] < self.h_ph1_map_cuts[iroc][cut].GetBinContent(
                        self.h_ph1_map_cuts[iroc][cut].GetMaximumBin()) and self.dict_cuts[cut] > 3:
                    max_ph1_map[iroc] = self.h_ph1_map_cuts[iroc][cut].GetBinContent(self.h_ph1_map_cuts[iroc][cut].GetMaximumBin())
                if (self.dict_cuts[cut] >= self.dict_cuts[self.cut_near_fiducial()]) and (
                            min_ph1_map[iroc] > self.h_ph1_map_cuts[iroc][cut].GetBinContent(
                            self.h_ph1_map_cuts[iroc][cut].GetMinimumBin())) and self.dict_cuts[cut] > 3:
                    min_ph1_map[iroc] = self.h_ph1_map_cuts[iroc][cut].GetBinContent(self.h_ph1_map_cuts[iroc][cut].GetMinimumBin())
                if max_ph2_map[iroc] < self.h_ph2_map_cuts[iroc][cut].GetBinContent(
                        self.h_ph2_map_cuts[iroc][cut].GetMaximumBin()) and self.dict_cuts[cut] > 3:
                    max_ph2_map[iroc] = self.h_ph2_map_cuts[iroc][cut].GetBinContent(self.h_ph2_map_cuts[iroc][cut].GetMaximumBin())
                if (self.dict_cuts[cut] >= self.dict_cuts[self.cut_near_fiducial()]) and (
                            min_ph2_map[iroc] > self.h_ph2_map_cuts[iroc][cut].GetBinContent(
                            self.h_ph2_map_cuts[iroc][cut].GetMinimumBin())) and self.dict_cuts[cut] > 3:
                    min_ph2_map[iroc] = self.h_ph2_map_cuts[iroc][cut].GetBinContent(self.h_ph2_map_cuts[iroc][cut].GetMinimumBin())
                if self.verbose: print 'Done'
                if not normalize_ph_plots:
                    if self.verbose: print 'Saving for ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut), ; sys.stdout.flush()
                    self.plots.set_2D_options(self.h_ph1_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries')
                    self.plots.set_1D_options('ph',self.h_ph1_cuts[iroc][cut],'ph 1 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph2_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries')
                    self.plots.set_1D_options('ph',self.h_ph2_cuts[iroc][cut],'ph 2 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph1_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 1 pix cluster(e)')
                    self.plots.set_2D_options(self.h_ph2_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 2 pix cluster(e)')

                    self.plots.save_individual_plots(self.h_ph1_evt_cuts[iroc][cut], 'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph2_evt_cuts[iroc][cut], 'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph1_cuts[iroc][cut], 'ph1_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_roc{r}_{c}'.format(r=iroc,c=cut), None, '', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph2_cuts[iroc][cut], 'ph2_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_roc{r}_{c}'.format(r=iroc,c=cut), None, '', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph1_map_cuts[iroc][cut], 'ph1_map_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_map_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph2_map_cuts[iroc][cut], 'ph2_map_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_map_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    if self.verbose: print 'Done'

        if normalize_ph_plots:
            min_ph1_map[iroc] = min(min_ph1_map[iroc], 0)
            min_ph2_map[iroc] = min(min_ph2_map[iroc], 0)
            minz_ph1 = min(minz_ph1, 0)
            minz_ph2 = min(minz_ph2, 0)
            for iroc in self.duts_list:
                for cut in self.cut_names:
                    if self.verbose: print 'Saving for ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut), ; sys.stdout.flush()
                    self.plots.set_2D_options(self.h_ph1_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries', min_val=minz_ph1, max_val=maxz_ph1)
                    self.plots.set_1D_options('ph',self.h_ph1_cuts[iroc][cut],'ph 1 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph2_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries', min_val=minz_ph2, max_val=maxz_ph2)
                    self.plots.set_1D_options('ph',self.h_ph2_cuts[iroc][cut],'ph 2 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph1_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 1 pix cluster(e)', min_val=min_ph1_map[iroc], max_val=max_ph1_map[iroc])
                    self.plots.set_2D_options(self.h_ph2_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 2 pix cluster(e)', min_val=min_ph2_map[iroc], max_val=max_ph2_map[iroc])

                    self.plots.save_individual_plots(self.h_ph1_evt_cuts[iroc][cut], 'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph2_evt_cuts[iroc][cut], 'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph1_cuts[iroc][cut], 'ph1_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_roc{r}_{c}'.format(r=iroc,c=cut), None, '', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph2_cuts[iroc][cut], 'ph2_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_roc{r}_{c}'.format(r=iroc,c=cut), None, '', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph1_map_cuts[iroc][cut], 'ph1_map_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_map_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    self.plots.save_individual_plots(self.h_ph2_map_cuts[iroc][cut], 'ph2_map_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_map_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                    if self.verbose: print 'Done'


    def do_cuts_analysis(self, do_occupancy=True, do_pulse_height=False, normalize_ph_plots=True):
        """
        calls the occupancy analysis and the ph analysis
        :param do_occupancy: if true, the method will call the occupancy analysis
        :param do_pulse_height: if true, the method will call the ph analysis
        :param normalize_ph_plots: If true, the ph analysis plots will have the same limits in Z
        :return:
        """
        self.print_banner('Starting Cuts Analysis...')
        self.print_banner('Creating histograms with cuts...')
        if do_occupancy: self.do_occupancy_analysis()
        if do_pulse_height: self.do_pulse_height_analysis(normalize_ph_plots)
        self.print_banner('Finished Cut Analysis', ':)')

    def cut_near_fiducial(self):
        """
        finds and returns the nearest cut that is enabled near the fiducial cut ('fiducial')
        :return: the nearest cut to 'fiducial'
        """
        for cut in ['fiducial', 'chi2x', 'anglex', 'rhit', 'masks', 'hit', 'tracks', 'beam']:
            if cut in self.cut_names:
                return cut
        return 'ini_fin'

    # ==============================================
    # region GET CONFIG

    def load_dut_type(self):
        """
        loads wether it is a Pad analysis or a Pixel analysis
        :return: dut type
        """
        dut_type = self.run_config_parser.get("BASIC", "type")
        assert dut_type.lower() in ["pixel", "pad"], "The DUT type {0} should be 'pixel' or 'pad'".format(dut_type)
        return dut_type

    def load_config(self):
        """
        Loads the configuration parameters from the config file
        :return:
        """
        # self.CutConfig['IndividualChCut'] = ''
        self.CutConfig['ExcludeFirst'] = self.ana_config_parser.getint('CUT', 'excludefirst') if self.ana_config_parser. \
            has_option('CUT', 'excludefirst') else 0
        self.CutConfig['ExcludeBeforeJump'] = self.ana_config_parser.getint('CUT', 'excludeBeforeJump') if self.ana_config_parser. \
            has_option('CUT', 'excludeBeforeJump') else 0
        self.CutConfig['ExcludeAfterJump'] = self.ana_config_parser.getint('CUT', 'excludeAfterJump') if self.ana_config_parser. \
            has_option('CUT', 'excludeAfterJump') else 0
        # self.CutConfig['EventRange'] = self.load_event_range(json.loads(self.ana_config_parser.get('CUT', 'EventRange')))
        self.CutConfig['chi2X'] = self.ana_config_parser.getint('CUT', 'chi2X') if self.ana_config_parser. \
            has_option('CUT', 'chi2X') else ''
        self.CutConfig['chi2Y'] = self.ana_config_parser.getint('CUT', 'chi2Y') if self.ana_config_parser. \
            has_option('CUT', 'chi2Y') else ''
        self.CutConfig['rhit'] = self.ana_config_parser.getint('CUT', 'r_hit') if self.ana_config_parser. \
            has_option('CUT', 'r_hit') else ''
        self.CutConfig['track_angle'] = self.ana_config_parser.getfloat('CUT', 'track_angle') if self.ana_config_parser. \
            has_option('CUT', 'track_angle') else ''
        self.CutConfig['MaskRowsROC4'] = self.ana_config_parser.get('CUT', 'MaskRowsROC4') if self.ana_config_parser. \
            has_option('CUT', 'MaskRowsROC4') else ''
        self.CutConfig['MaskRowsROC5'] = self.ana_config_parser.get('CUT', 'MaskRowsROC5') if self.ana_config_parser. \
            has_option('CUT', 'MaskRowsROC5') else ''
        self.CutConfig['MaskRowsROC6'] = self.ana_config_parser.get('CUT', 'MaskRowsROC6') if self.ana_config_parser. \
            has_option('CUT', 'MaskRowsROC6') else ''
        self.CutConfig['MaskColsROC4'] = self.ana_config_parser.get('CUT', 'MaskColsROC4') if self.ana_config_parser. \
            has_option('CUT', 'MaskColsROC4') else ''
        self.CutConfig['MaskColsROC5'] = self.ana_config_parser.get('CUT', 'MaskColsROC5') if self.ana_config_parser. \
            has_option('CUT', 'MaskColsROC5') else ''
        self.CutConfig['MaskColsROC6'] = self.ana_config_parser.get('CUT', 'MaskColsROC6') if self.ana_config_parser. \
            has_option('CUT', 'MaskColsROC6') else ''
        self.CutConfig['MaskPixelsROC4'] = self.ana_config_parser.get('CUT', 'MaskPixelsROC4') if self.ana_config_parser. \
            has_option('CUT', 'MaskPixelsROC4') else ''
        self.CutConfig['MaskPixelsROC5'] = self.ana_config_parser.get('CUT', 'MaskPixelsROC5') if self.ana_config_parser. \
            has_option('CUT', 'MaskPixelsROC5') else ''
        self.CutConfig['MaskPixelsROC6'] = self.ana_config_parser.get('CUT', 'MaskPixelsROC6') if self.ana_config_parser. \
            has_option('CUT', 'MaskPixelsROC6') else ''
        self.CutConfig['FidRegionROC4'] = self.ana_config_parser.get('CUT', 'FidRegionROC4') if self.ana_config_parser. \
            has_option('CUT', 'FidRegionROC4') else ''
        self.CutConfig['FidRegionROC5'] = self.ana_config_parser.get('CUT', 'FidRegionROC5') if self.ana_config_parser. \
            has_option('CUT', 'FidRegionROC5') else ''
        self.CutConfig['FidRegionROC6'] = self.ana_config_parser.get('CUT', 'FidRegionROC6') if self.ana_config_parser. \
            has_option('CUT', 'FidRegionROC6') else ''
        self.CutConfig['IniFin'] = self.ana_config_parser.getboolean('CUT', 'IniFin') if self.ana_config_parser.has_option('CUT', 'IniFin') else False
        self.CutConfig['Beam'] = self.ana_config_parser.getboolean('CUT', 'Beam') if self.ana_config_parser.has_option('CUT', 'Beam') else False
        self.CutConfig['Tracks'] = self.ana_config_parser.getboolean('CUT', 'Tracks') if self.ana_config_parser.has_option('CUT', 'Tracks') else False
        self.CutConfig['Hit'] = self.ana_config_parser.getboolean('CUT', 'Hit') if self.ana_config_parser.has_option('CUT', 'Hit') else False
        self.CutConfig['Mask'] = self.ana_config_parser.getboolean('CUT', 'Mask') if self.ana_config_parser.has_option('CUT', 'Mask') else False
        self.CutConfig['Fiducial'] = self.ana_config_parser.getboolean('CUT', 'Fiducial') if self.ana_config_parser.has_option('CUT', 'Fiducial') else False
        self.CutConfig['Chi2'] = self.ana_config_parser.getboolean('CUT', 'Chi2') if self.ana_config_parser.has_option('CUT', 'Chi2') else False
        self.CutConfig['Angle'] = self.ana_config_parser.getboolean('CUT', 'Angle') if self.ana_config_parser.has_option('CUT', 'Angle') else False
        self.CutConfig['RHit'] = self.ana_config_parser.getboolean('CUT', 'RHit') if self.ana_config_parser.has_option('CUT', 'RHit') else False

    def is_first_cut(self):
        """
        tells if it is the first cut to be applied
        :return: returns True, if it is the first cut to be applied
        """
        return self.num_cuts == 0

    def generate_ini_fin_cuts(self):
        """
        generates the ini_fin cut and accumulates its result on the accumulated dictionaries cut_{pixelated/hitmap}_roc(_incr).
        :return:
        """
        if self.verbose: print 'Creating cut for initial and final', abs(self.CutConfig['ExcludeFirst']), 'seconds...', ; sys.stdout.flush()
        picklepath = 'Configuration/Individual_Configs/IniFin/{tc}_{r}.pickle'.format(tc=self.TESTCAMPAIGN, r=self.run_number)
        if not os.path.isdir('Configuration/Individual_Configs/IniFin'): os.makedirs('Configuration/Individual_Configs/IniFin')
        def func0():
            nentries = self.analysis.tree.GetEntries()
            self.analysis.tree.GetEntry(0)
            first_t = self.analysis.tree.time
            self.analysis.tree.GetEntry(nentries-1)
            last_t = self.analysis.tree.time
            ini_fin_cut = 'time>{ini}&&time<{fin}'.format(ini=first_t + abs(self.CutConfig['ExcludeFirst'])*1000, fin=last_t - abs(self.CutConfig['ExcludeFirst'])*1000)
            return deepcopy(str(ini_fin_cut))
        nentries = self.analysis.tree.GetEntries()
        self.analysis.tree.GetEntry(0)
        first_t = self.analysis.tree.time
        self.analysis.tree.GetEntry(nentries-1)
        last_t = self.analysis.tree.time
        self.ini_fin_cut = self.do_pickle(picklepath, func0)
        for iroc in self.duts_list:
            self.gen_vect_cuts(self.ini_fin_cut, self.ini_fin_cut, iroc)
        self.num_cuts += 1
        if self.verbose: print 'Done'

    def generate_hit_cut(self): # TODO implement cut! must change tree structure from tracking telescope
        """
        Needs to be implemented. Have to change trackingTelescope for this
        :return:
        """
        if self.verbose: print 'Creating cut to require at least one hit for each DUT plane...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            self.gen_vect_cuts('1==1', '1==1', iroc)
        self.num_cuts += 1
        if self.verbose: print 'Bla'

    def generate_beam_interruption_cut(self):
        """
        locates the beam interruptions and cuts some seconds before and some seconds after which are specified in the config file.
        overlapping segments are fused to one to reduce the size of the string. An average of the number of event per time bin is calculated
        and every bin below 90% or above 20% is excluded
        :return:
        """
        # time is in ms. good results found with bin size of 5 seconds
        picklepath = 'Configuration/Individual_Configs/Beam/{tc}_{r}.pickle'.format(tc=self.TESTCAMPAIGN, r=self.run_number)
        if not os.path.isdir('Configuration/Individual_Configs/Beam'): os.makedirs('Configuration/Individual_Configs/Beam')
        print 'Generating beam interruption and overshoots cut for run', self.run_number, '...', ; sys.stdout.flush()

        def func0():
            nentries = self.analysis.tree.GetEntries()
            self.analysis.tree.GetEntry(0)
            first_t = self.analysis.tree.time
            self.analysis.tree.GetEntry(nentries-1)
            last_t = self.analysis.tree.time
            bins = int((last_t-first_t)/float(5000))
            vector_interr = {}
            max_interr = 0
            gROOT.SetBatch(True)
            h1 = TH1F('h_beam_time_', 'h_beam_time_', bins+1, first_t-(last_t-first_t)/float(2*bins), last_t+(last_t-first_t)/float(2*bins))
            self.analysis.tree.Draw('time>>h_beam_time_', self.cuts_pixelated_roc_incr[self.duts_list[0]][self.num_cuts-1], 'goff') if not self.is_first_cut() else self.analysis.tree.Draw('time>>h_beam_time_', '', 'goff')
            gROOT.SetBatch(False)
            mean = h1.Integral()/float(h1.GetNbinsX())
            beam_interr_cut = TCut('beam_interruptions_cut', '')
            vector_interr = {}
            max_interr = 0
            for t in xrange(1, h1.GetNbinsX()+1):
                if h1.GetBinContent(t) < mean*0.9 or h1.GetBinContent(t) > mean*1.2:
                    if max_interr != 0:
                        if h1.GetBinLowEdge(t) - abs(self.CutConfig['ExcludeBeforeJump'])*1000 < vector_interr[max_interr-1]['f']:
                            vector_interr[max_interr-1]['f'] = h1.GetBinLowEdge(t)+h1.GetBinWidth(t)+abs(self.CutConfig['ExcludeAfterJump'])*1000
                        else:
                            vector_interr[max_interr] = {'i': h1.GetBinLowEdge(t)-abs(self.CutConfig['ExcludeBeforeJump'])*1000, 'f': h1.GetBinLowEdge(t)+h1.GetBinWidth(t)+abs(self.CutConfig['ExcludeAfterJump'])*1000}
                            max_interr += 1
                    else:
                        vector_interr[max_interr] = {'i': h1.GetBinLowEdge(t)-abs(self.CutConfig['ExcludeBeforeJump'])*1000, 'f': h1.GetBinLowEdge(t)+h1.GetBinWidth(t)+abs(self.CutConfig['ExcludeAfterJump'])*1000}
                        max_interr += 1
            for interr in xrange(max_interr):
                beam_interr_cut = beam_interr_cut + TCut('bi{i}'.format(i=interr), 'time<{low}||time>{high}'.format(low=vector_interr[interr]['i'], high=vector_interr[interr]['f']))
            return deepcopy(str(beam_interr_cut.GetTitle()))

        self.beam_interr_cut = self.do_pickle(picklepath, func0)
        for iroc in self.duts_list:
            self.gen_vect_cuts(self.beam_interr_cut, self.beam_interr_cut, iroc)
        self.num_cuts += 1; nentries = self.analysis.tree.GetEntries(); self.analysis.tree.GetEntry(0); first_t = self.analysis.tree.time; self.analysis.tree.GetEntry(nentries-1); last_t = self.analysis.tree.time
        if self.verbose: print 'Done'

    def generate_rhit_cuts(self):
        if self.verbose: print 'Generatin R-hit cut for distances greater than', self.CutConfig['rhit'], 'um between predicted position from the track and the seed cluster position...', ; sys.stdout.flush()
        self.rhit_cut = {}
        for iroc in self.duts_list:
            self.generate_rhit_cuts_DUT(iroc)
            self.gen_vect_cuts(self.rhit_cut[iroc], self.rhit_cut[iroc], iroc)
        self.num_cuts += 1
        if self.verbose: print 'Done'

    def generate_rhit_cuts_DUT(self, dut):
        picklepath = 'Configuration/Individual_Configs/RHitRoc{ir}/{tc}_{r}.pickle'.format(ir=dut, tc=self.TESTCAMPAIGN, r=self.run_number)
        if not os.path.isdir('Configuration/Individual_Configs/RHitRoc{ir}'.format(ir=dut)): os.makedirs('Configuration/Individual_Configs/RHitRoc{ir}'.format(ir=dut))
        def func0():
            value = self.CutConfig['rhit']
            string = '((10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))<{val}&&(sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>=0)'.format(n=dut, val=value)
            return string
        self.rhit_cut[dut] = self.do_pickle(picklepath, func0)

    def generate_tracks_cut(self):
        if self.verbose: print 'Generating tracks cut...',
        self.cut_tracks = TCut('cut_tracks', 'n_tracks==1')
        for iroc in self.duts_list:
            self.gen_vect_cuts(self.cut_tracks.GetTitle(), self.cut_tracks.GetTitle(), iroc)
        self.num_cuts += 1
        if self.verbose: print 'Done'

    def generate_chi2_cuts(self):
        if self.verbose: print 'Generating cut for a chi2 with a percentile of', self.CutConfig['chi2X'], 'in X and', self.CutConfig['chi2Y'], 'in Y...', ; sys.stdout.flush()
        self.chi2_cut = {i:{} for i in self.duts_list}
        self.h_chi2 = {i:{} for i in self.duts_list}
        self.h_chi2_cut = {i:{} for i in self.duts_list}
        num_prev_cut = self.num_cuts - 1
        for iroc in self.duts_list:
            self.generate_chi2('x', num_prev_cut, iroc)
            self.generate_chi2('y', num_prev_cut, iroc)
            self.gen_vect_cuts(self.chi2_cut[iroc]['x'], self.chi2_cut[iroc]['x'], iroc)
        self.num_cuts += 1
        for iroc in self.duts_list:
            self.gen_vect_cuts(self.chi2_cut[iroc]['y'], self.chi2_cut[iroc]['y'], iroc)
        self.num_cuts += 1
        if self.verbose: print 'Done'

    def generate_chi2(self, mode='x', num_prev_cut=4, iroc=4):
        picklepath = 'Configuration/Individual_Configs/Chi2Roc{ir}{m}/{tc}_{r}.pickle'.format(ir=iroc, m=mode, tc=self.TESTCAMPAIGN, r=self.run_number)
        if not os.path.isdir('Configuration/Individual_Configs/Chi2Roc{ir}{m}'.format(ir=iroc, m=mode)): os.makedirs('Configuration/Individual_Configs/Chi2Roc{ir}{m}'.format(ir=iroc, m=mode))
        def func0():
            gROOT.SetBatch(1)
            h_chi2 = TH1F('h_chi2_roc{r}_{m}_'.format(r=iroc,m=mode), 'h_chi2_roc{r}_{m}_'.format(r=iroc,m=mode), 201, -0.1, 40.1)
            nq = 100
            chi2s = zeros(nq)
            xq = array([(i + 1) / float(nq) for i in xrange(nq)])
            self.analysis.tree.Draw('chi2_{mod}>>h_chi2_roc{r}_{mod}_'.format(r=iroc,mod=mode), self.cuts_pixelated_roc_incr[iroc][num_prev_cut], 'goff') if not self.is_first_cut() else self.analysis.tree.Draw('chi2_{mod}>>h_chi2_roc{r}_{mod}_'.format(r=iroc,mod=mode), '', 'goff')
            h_chi2.GetQuantiles(nq, chi2s, xq)
            gROOT.SetBatch(0)
            quantile = self.CutConfig['chi2{mod}'.format(mod=mode.title())]
            assert type(quantile) is int and 0 < quantile <= 100, 'chi2 quantile has to be an integer between (0, 100]'
            string = 'chi2_{mod}<{val}&&chi2_{mod}>=0'.format(val=chi2s[quantile], mod=mode)
            return string
        cut = self.do_pickle(picklepath, func0)
        self.chi2_cut[iroc][mode] = cut

    def generate_fid_cuts(self):
        if self.verbose: print 'Generating fiducial cuts...', ; sys.stdout.flush()
        self.fid_cut_hitmap_roc = {}
        self.fid_cut_pixelated_roc = {}
        # self.fid_cut_local_roc = {}
        # self.fid_cut_telescope_roc = {}
        for iroc in self.duts_list:
            self.generate_fid_cuts_DUT(iroc)
            self.gen_vect_cuts('fidcut_hitmap_roc{r}'.format(r=iroc), 'fidcut_pixelated_roc{r}'.format(r=iroc), iroc)
        self.num_cuts += 1
        if self.verbose: print 'Done'

    def generate_fid_cuts_DUT(self, roc):
        fidcutstring = self.CutConfig['FidRegionROC{d}'.format(d=roc)]
        varx = 'col'
        vary = 'row'
        xmin = self.plot_settings['minCol']
        xmax = self.plot_settings['maxCol']
        ymin = self.plot_settings['minRow']
        ymax = self.plot_settings['maxRow']
        deltax = 1
        deltay = 1
        if fidcutstring is not '':
            fidcutstring = fidcutstring.replace(']', '').replace('[', '')
            if fidcutstring is not '':
                fidcutstring = fidcutstring.split(';')
                roc = int(fidcutstring[0])
                xmin = int(fidcutstring[1].split(',')[0])
                xmax = int(fidcutstring[1].split(',')[1])
                ymin = int(fidcutstring[2].split(',')[0])
                ymax = int(fidcutstring[2].split(',')[1])
        self.fid_cut_hitmap_roc[roc] = self.make_fid_cut('fidcut_hitmap_roc{r}'.format(r=roc), varx, vary, xmin, xmax,
                                                         deltax, ymin, ymax, deltay)
        varx2 = 'cluster_col_ROC{n}'.format(n=roc)
        vary2 = 'cluster_row_ROC{n}'.format(n=roc)
        self.fid_cut_pixelated_roc[roc] = self.make_fid_cut('fidcut_pixelated_roc{r}'.format(r=roc), varx2, vary2, xmin,
                                                            xmax, deltax, ymin, ymax, deltay)
        # self.fid_cut_local_roc[roc] = self.make_fid_cut('fidcut_local_roc{r}'.format(r=roc), varx, vary, xmin,
        #                                                     xmax, deltax, ymin, ymax, deltay)
        # self.fid_cut_telescope_roc[roc] = self.make_fid_cut('fidcut_telescope_roc{r}'.format(r=roc), varx, vary, xmin,
        #                                                     xmax, deltax, ymin, ymax, deltay)

    def make_fid_cut(self, name='fidcut', varx='col', vary='row', xmin=0, xmax=51, deltax=1, ymin=0, ymax=79, deltay=1):
        cut = TCutG(name, 5)
        cut.SetVarX(varx)
        cut.SetVarY(vary)
        cut.SetPoint(0, float(xmin) - float(deltax)/2, float(ymin) - float(deltay)/2)
        cut.SetPoint(1, float(xmin) - float(deltax)/2, float(ymax) + float(deltay)/2)
        cut.SetPoint(2, float(xmax) + float(deltax)/2, float(ymax) + float(deltay)/2)
        cut.SetPoint(3, float(xmax) + float(deltax)/2, float(ymin) - float(deltay)/2)
        cut.SetPoint(4, float(xmin) - float(deltax)/2, float(ymin) - float(deltay)/2)
        cut.SetLineColor(kRed)
        cut.SetLineWidth(3*3)
        return cut

    def gen_vect_cuts(self, cut_hitmap, cut_pixelated, roc=4):
        self.cuts_hitmap_roc[roc][self.num_cuts] = cut_hitmap
        self.cuts_pixelated_roc[roc][self.num_cuts] = cut_pixelated
        self.accum_incr_vect_cuts(roc)

    def accum_incr_vect_cuts(self, roc=4):
        self.cuts_hitmap_roc_incr[roc][self.num_cuts] = self.cuts_hitmap_roc_incr[roc][self.num_cuts-1] + '&&(' + self.cuts_hitmap_roc[roc][self.num_cuts] + ')' if self.num_cuts != 0 else '(' + self.cuts_hitmap_roc[roc][self.num_cuts] + ')'
        self.cuts_pixelated_roc_incr[roc][self.num_cuts] = self.cuts_pixelated_roc_incr[roc][self.num_cuts-1] + '&&(' + self.cuts_pixelated_roc[roc][self.num_cuts] + ')' if self.num_cuts != 0 else '(' + self.cuts_pixelated_roc[roc][self.num_cuts] + ')'

    def generate_masks(self):
        if self.verbose: print 'Generating masks cuts...', ; sys.stdout.flush()
        self.mask_hitmap_roc = {}
        self.mask_pixelated_roc = {}
        self.generate_col_masks()
        self.generate_row_masks()
        self.generate_pixel_masks()
        for roc in self.duts_list:
            self.mask_hitmap_roc[roc] = self.mask_hitmap_roc[roc].GetTitle()
            self.mask_pixelated_roc[roc] = self.mask_pixelated_roc[roc].GetTitle()
            self.gen_vect_cuts(self.mask_hitmap_roc[roc], self.mask_pixelated_roc[roc], roc)
        self.num_cuts += 1
        if self.verbose: print 'Done'

    def generate_col_masks(self):
        self.col_mask_hitmap_roc = {}
        self.col_mask_pixelate_roc = {}
        for iroc in self.duts_list:
            self.generate_col_masks_DUT(iroc)

    def generate_col_masks_DUT(self, roc):
        maskcolstring = self.CutConfig['MaskColsROC{d}'.format(d=roc)]
        title = ''
        name = 'mask_col_hitmap_roc{r}'.format(r=roc)
        mask_col_hitmap_temp = TCut('temp0', title)
        if maskcolstring is not '':
            maskcolstring = maskcolstring.replace('[', '').replace(']', '')
            if maskcolstring is not '':
                maskcolstring = maskcolstring.split(';')
                roc = int(maskcolstring[0])
                for i in xrange(1, len(maskcolstring)):
                    if ':' in maskcolstring[i]:
                        tempstring = maskcolstring[i].split(':')
                        tempmask = TCut('temp', '(col<{inf}||col>{sup})'.format(inf=tempstring[0], sup=tempstring[1]))
                    else:
                        tempmask = TCut('temp', '(col!={val})'.format(val=maskcolstring[i]))
                    mask_col_hitmap_temp = mask_col_hitmap_temp + tempmask
        self.col_mask_hitmap_roc[roc] = TCut(name, '')
        self.col_mask_hitmap_roc[roc] = self.col_mask_hitmap_roc[roc] + mask_col_hitmap_temp
        name2 = 'mask_col_pixelated_roc{r}'.format(r=roc)
        self.col_mask_pixelate_roc[roc] = TCut(name2, self.col_mask_hitmap_roc[roc].GetTitle().replace('col', 'cluster_col_ROC{n}'.format(n=roc)))
        name3 = 'mask_hitmap_roc{r}'.format(r=roc)
        self.mask_hitmap_roc[roc] = TCut(name3, '')
        self.mask_hitmap_roc[roc] = self.mask_hitmap_roc[roc] + self.col_mask_hitmap_roc[roc]
        name4 = 'mask_pixelated_roc{r}'.format(r=roc)
        self.mask_pixelated_roc[roc] = TCut(name4, '')
        self.mask_pixelated_roc[roc] = self.mask_pixelated_roc[roc] + self.col_mask_pixelate_roc[roc]

    def generate_row_masks(self):
        self.row_mask_hitmap_roc = {}
        self.row_mask_pixelated_roc = {}
        for iroc in self.duts_list:
            self.generate_row_masks_DUT(iroc)

    def generate_row_masks_DUT(self, roc):
        maskrowstring = self.CutConfig['MaskRowsROC{d}'.format(d=roc)]
        title = ''
        name = 'mask_row_hitmap_roc{r}'.format(r=roc)
        mask_row_hitmap_temp = TCut('temp0', title)
        if maskrowstring is not '':
            maskrowstring = maskrowstring.replace('[', '').replace(']', '')
            if maskrowstring is not '':
                maskrowstring = maskrowstring.split(';')
                roc = int(maskrowstring[0])
                for i in xrange(1,len(maskrowstring)):
                    if ':' in maskrowstring[i]:
                        tempstring = maskrowstring[i].split(':')
                        tempmask = TCut('temp', '(row<{inf}||row>{sup})'.format(inf=tempstring[0], sup=tempstring[1]))
                    else:
                        tempmask = TCut('temp', '(row!={val})'.format(val=maskrowstring[i]))
                    mask_row_hitmap_temp = mask_row_hitmap_temp + tempmask
        self.row_mask_hitmap_roc[roc] = TCut(name, '')
        self.row_mask_hitmap_roc[roc] = self.row_mask_hitmap_roc[roc] + mask_row_hitmap_temp
        name2 = 'mask_row_pixelated_roc{r}'.format(r=roc)
        self.row_mask_pixelated_roc[roc] = TCut(name2, self.row_mask_hitmap_roc[roc].GetTitle().replace('row', 'cluster_row_ROC{n}'.format(n=roc)))
        self.mask_hitmap_roc[roc] = self.mask_hitmap_roc[roc] + self.row_mask_hitmap_roc[roc]
        self.mask_pixelated_roc[roc] = self.mask_pixelated_roc[roc] + self.row_mask_pixelated_roc[roc]

    def generate_pixel_masks(self):
        self.pixel_mask_hitmap_roc = {}
        self.pixel_mask_pixelated_roc = {}
        for iroc in self.duts_list:
            self.generate_pixel_masks_DUT(iroc)

    def generate_pixel_masks_DUT(self, roc):
        maskpixelstring = self.CutConfig['MaskPixelsROC{d}'.format(d=roc)]
        title = ''
        name = 'mask_pixel_hitmap_roc{r}'.format(r=roc)
        mask_pixel_hitmap_temp = TCut('temp0', title)
        if maskpixelstring is not '':
            maskpixelstring = maskpixelstring.replace('[', '').replace(']', '')
            if maskpixelstring is not '':
                maskpixelstring = maskpixelstring.split(';')
                roc = int(maskpixelstring[0])
                for i in xrange(1, len(maskpixelstring)):
                    if ',' in maskpixelstring[i]:
                        tempstring = maskpixelstring[i].split(',')
                        tempmask = TCut('temp', '(col!={x}||row!={y})'.format(x=tempstring[0], y=tempstring[1]))
                        mask_pixel_hitmap_temp = mask_pixel_hitmap_temp + tempmask
        self.pixel_mask_hitmap_roc[roc] = TCut(name, '')
        self.pixel_mask_hitmap_roc[roc] = self.pixel_mask_hitmap_roc[roc] + mask_pixel_hitmap_temp
        name2 = 'mask_pixel_pixelated_roc{r}'.format(r=roc)
        self.pixel_mask_pixelated_roc[roc] = TCut(name2, self.pixel_mask_hitmap_roc[roc].GetTitle().replace('row', 'cluster_row_ROC{n}'.format(n=roc)).replace('col', 'cluster_col_ROC{n}'.format(n=roc)))
        self.mask_hitmap_roc[roc] = self.mask_hitmap_roc[roc] + self.pixel_mask_hitmap_roc[roc]
        self.mask_pixelated_roc[roc] = self.mask_pixelated_roc[roc] + self.pixel_mask_pixelated_roc[roc]

    def generate_angle_cuts(self):
        if self.verbose: print 'Generating angle cuts of', self.CutConfig['track_angle'], 'deg in X and Y...', ; sys.stdout.flush()
        self.angle_cut = {i: {} for i in self.duts_list}
        self.h_angle = {i: {} for i in self.duts_list}
        self.h_angle_cut = {i: {} for i in self.duts_list}
        prev_num_cut = self.num_cuts - 1
        for iroc in self.duts_list:
            self.generate_angle('x', prev_num_cut, iroc)
            self.generate_angle('y', prev_num_cut, iroc)
            self.gen_vect_cuts(self.angle_cut[iroc]['x'], self.angle_cut[iroc]['x'], iroc)
        self.num_cuts += 1
        for iroc in self.duts_list:
            self.gen_vect_cuts(self.angle_cut[iroc]['y'], self.angle_cut[iroc]['y'], iroc)
        self.num_cuts += 1
        if self.verbose: print 'Done'
        # for iroc in self.duts_list:
        #     gROOT.SetBatch(1)
        #     self.analysis.tree.Draw('angle_x>>h_angle_roc{r}_x_cut'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
        #     self.analysis.tree.Draw('angle_y>>h_angle_roc{r}_y_cut'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
        #     self.plots.set_1D_options('angle', self.h_angle[iroc]['x'], 'angleX', 'entries', kBlue)
        #     self.plots.set_1D_options('angle', self.h_angle[iroc]['y'], 'angleY', 'entries', kBlue)
        #     self.plots.set_1D_options('angle', self.h_angle_cut[iroc]['x'], 'angleX', 'entries', kRed)
        #     self.plots.set_1D_options('angle', self.h_angle_cut[iroc]['y'], 'angleY', 'entries', kRed)
        #     gROOT.SetBatch(0)
        #     self.plots.save_cuts_distributions(self.h_angle[iroc]['x'], self.h_angle_cut[iroc]['x'], 'angle_roc{r}_x_cut_overlay'.format(r=iroc), 'Angle ROC{r} x Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
        #     self.plots.save_cuts_distributions(self.h_angle[iroc]['y'], self.h_angle_cut[iroc]['y'], 'angle_roc{r}_y_cut_overlay'.format(r=iroc), 'Angle ROC{r} y Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)

    def generate_angle(self, mode='x', prev_num_cut=6, iroc=4):
        picklepath = 'Configuration/Individual_Configs/AngleRoc{ir}{m}/{tc}_{run}.pickle'.format(ir=iroc, m=mode, tc=self.TESTCAMPAIGN, run=self.analysis.lowest_rate_run)
        if not os.path.isdir('Configuration/Individual_Configs/AngleRoc{ir}{m}'.format(ir=iroc, m=mode)): os.makedirs('Configuration/Individual_Configs/AngleRoc{ir}{m}'.format(ir=iroc, m=mode))
        angle = self.CutConfig['track_angle']

        # def func():
        #     print 'generating slope cut for run {run}...'.format(run=self.analysis.run_number)
        # fit the slope to get the mean
        # gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        def func0():
            gROOT.SetBatch(1)
            gROOT.ProcessLine('gErrorIgnoreLevel = 6000;')
            h = TH1F('h_angle_', 'h_angle_', 61, -3.05, 3.05)
            self.analysis.tree.Draw('angle_{x}>>h_angle_'.format(x=mode), self.cuts_pixelated_roc_incr[iroc][prev_num_cut], 'goff') if not self.is_first_cut() else self.analysis.tree.Draw('angle_{x}>>h_angle_'.format(x=mode), '', 'goff')
            h_angle = TH1F('h_angle_roc_', 'h_angle_roc_', 51, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) - 3*h.GetRMS() - 3*h.GetRMS()/50, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) + 3*h.GetRMS() + 3*h.GetRMS()/50)
            # self.h_angle_cut[iroc][mode] = TH1F('h_angle_roc{r}_{m}_cut'.format(r=iroc,m=mode), 'h_angle_roc{r}_{m}_cut'.format(r=iroc,m=mode), 51, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) - 2*h.GetRMS() - 2*h.GetRMS()/50, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) + 2*h.GetRMS() + 2*h.GetRMS()/50)
            self.analysis.tree.Draw('angle_{m}>>h_angle_roc_'.format(m=mode), self.cuts_pixelated_roc_incr[iroc][prev_num_cut], 'goff') if not self.is_first_cut() else self.analysis.tree.Draw('angle_{m}>>h_angle_roc_'.format(m=mode), '', 'goff')
            fit_result = h_angle.Fit('gaus', 'qs')
            # fit_result = h.Fit('gaus', 'qs')# , '', xmin, xmax)
            x_mean = fit_result.Parameter(1)
            angles = [x_mean - angle, x_mean + angle]
            # c = gROOT.FindObject('c1')
            # c.Close()
            gROOT.ProcessLine('gErrorIgnoreLevel = -1;')
            gROOT.SetBatch(0)
            return angles

        angles = self.do_pickle(picklepath, func0)
        # create the cut string
        string = 'angle_{x}>={minx}&&angle_{x}<={maxx}'.format(x=mode, minx=angles[0], maxx=angles[1])
        self.angle_cut[iroc][mode] = string if angle > 0 else ''

    def get_nearest_prev_existing_cut_key(self, cutname):
        allcuts = {0: 'ini_fin', 1: 'beam', 2: 'tracks', 3: 'hit', 4: 'masks', 5: 'fiducial', 6: 'chi2x', 7: 'chi2y', 8: 'anglex', 9: 'angley', 10: 'rhit'}
        for key in allcuts.keys():
            if cutname == allcuts[key]:
                cutkey = key
                break
        for key in xrange(cutkey-1, allcuts.keys()[0], -1):
            if allcuts[key] in self.cut_names:
                return key
        return allcuts.keys()[0]

    def do_beam_distribution(self):
        if self.verbose:
            print 'Beam interruption...', ; sys.stdout.flush()
        nentries = self.analysis.tree.GetEntries()
        self.analysis.tree.GetEntry(0)
        first_t = self.analysis.tree.time
        self.analysis.tree.GetEntry(nentries-1)
        last_t = self.analysis.tree.time
        bins = int((last_t-first_t)/float(5000))
        gROOT.SetBatch(True)
        self.h_beam_time = TH1F('h_beam_time', 'h_beam_time', bins+1, first_t-(last_t-first_t)/float(2*bins), last_t+(last_t-first_t)/float(2*bins))
        self.h_beam_time_cut = TH1F('h_beam_time_cut', 'h_beam_time_cut', bins+1, first_t-(last_t-first_t)/float(2*bins), last_t+(last_t-first_t)/float(2*bins))
        self.h_beam_mean_cut = TH1F('h_beam_mean_cut', 'h_beam_mean_cut', bins+1, first_t-(last_t-first_t)/float(2*bins), last_t+(last_t-first_t)/float(2*bins))
        if 'beam' in self.cut_names:
            self.analysis.tree.Draw('time>>h_beam_time', self.cuts_pixelated_roc_incr[self.duts_list[0]][self.dict_cuts['beam']-1],'goff') if self.dict_cuts['beam'] != 0 else self.analysis.tree.Draw('time>>h_beam_time', '','goff')
            self.analysis.tree.Draw('time>>h_beam_time_cut', self.cuts_pixelated_roc_incr[self.duts_list[0]][self.dict_cuts['beam']],'goff')
        else:
            key = self.get_nearest_prev_existing_cut_key('beam')
            self.analysis.tree.Draw('time>>h_beam_time', self.cuts_pixelated_roc_incr[self.duts_list[0]][key],'goff') if key != -1 else self.analysis.tree.Draw('time>>h_beam_time', '','goff')
            self.analysis.tree.Draw('time>>h_beam_time_cut', self.cuts_pixelated_roc_incr[self.duts_list[0]][key],'goff') if key != -1 else self.analysis.tree.Draw('time>>h_beam_time_cut', '','goff')
        self.mean_events_5sec = self.h_beam_time.Integral()/float(self.h_beam_time.GetNbinsX())
        binsEvents = int(ceil(nentries/float(self.mean_events_5sec)))
        self.plot_settings['event_bins'] = binsEvents
        for bin in xrange(1, bins+2):
            self.h_beam_mean_cut.SetBinContent(bin, self.mean_events_5sec)
        gROOT.SetBatch(False)
        self.plots.set_1D_options('time', self.h_beam_time, 'time(ms)', 'entries', kBlue)
        self.plots.set_1D_options('time', self.h_beam_time_cut, 'time(ms)', 'entries', kRed)
        self.plots.set_1D_options('time', self.h_beam_mean_cut, 'time(ms)', 'entries', color=kGreen)
        self.plots.save_cuts_distributions(self.h_beam_time, self.h_beam_time_cut, 'beam_time_cut_overlay', 'Beam cut overlay', '', 1000000011, self.plots.save_dir+'/cuts', False, self.h_beam_mean_cut)
        if self.verbose: print 'Done'

    def do_chi2_distributions(self):
        self.h_chi2x_dist = {}
        self.h_chi2y_dist = {}
        self.h_chi2x_cut_dist = {}
        self.h_chi2y_cut_dist = {}
        if self.verbose: print 'Chi2...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2x_dist[iroc] = TH1F('h_chi2x_roc{r}'.format(r=iroc), 'h_chi2x_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            self.h_chi2y_dist[iroc] = TH1F('h_chi2y_roc{r}'.format(r=iroc), 'h_chi2y_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            self.h_chi2x_cut_dist[iroc] = TH1F('h_chi2x_cut_roc{r}'.format(r=iroc), 'h_chi2x_cut_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            self.h_chi2y_cut_dist[iroc] = TH1F('h_chi2y_cut_roc{r}'.format(r=iroc), 'h_chi2y_cut_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            if 'chi2x' in self.cut_names and 'chi2y' in self.cut_names:
                self.analysis.tree.Draw('chi2_x>>h_chi2x_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2x'] - 1], 'goff') if self.dict_cuts['chi2x'] != 0 else self.analysis.tree.Draw('chi2_x>>h_chi2x_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('chi2_y>>h_chi2y_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2x'] - 1], 'goff') if self.dict_cuts['chi2x'] != 0 else self.analysis.tree.Draw('chi2_y>>h_chi2y_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('chi2_x>>h_chi2x_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2y']], 'goff')
                self.analysis.tree.Draw('chi2_y>>h_chi2y_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2y']], 'goff')
            else:
                key = self.get_nearest_prev_existing_cut_key('chi2x')
                self.analysis.tree.Draw('chi2_x>>h_chi2x_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('chi2_x>>h_chi2x_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('chi2_y>>h_chi2y_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('chi2_y>>h_chi2y_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('chi2_x>>h_chi2x_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('chi2_x>>h_chi2x_cut_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('chi2_y>>h_chi2y_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('chi2_y>>h_chi2y_cut_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_1D_options('chi2', self.h_chi2x_dist[iroc], 'chi2X', 'entries', kBlue)
            self.plots.set_1D_options('chi2', self.h_chi2y_dist[iroc], 'chi2Y', 'entries', kBlue)
            self.plots.set_1D_options('chi2', self.h_chi2x_cut_dist[iroc], 'chi2X', 'entries', kRed)
            self.plots.set_1D_options('chi2', self.h_chi2y_cut_dist[iroc], 'chi2Y', 'entries', kRed)
            gROOT.SetBatch(False)
            self.plots.save_cuts_distributions(self.h_chi2x_dist[iroc], self.h_chi2x_cut_dist[iroc], 'chi2_roc{r}_x_cut_overlay'.format(r=iroc), 'Chi2 roc{r} x Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
            self.plots.save_cuts_distributions(self.h_chi2y_dist[iroc], self.h_chi2y_cut_dist[iroc], 'chi2_roc{r}_y_cut_overlay'.format(r=iroc), 'Chi2 roc{r} Y Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
        if self.verbose: print 'Done'

    def do_angle_distributions(self):
        if self.verbose: print 'Angle...', ; sys.stdout.flush()
        self.h_anglex_dist = {}
        self.h_angley_dist = {}
        self.h_anglex_cut_dist = {}
        self.h_angley_cut_dist = {}
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_anglex_dist[iroc] = TH1F('h_anglex_roc{r}'.format(r=iroc), 'h_anglex_roc{r}'.format(r=iroc), 121, -3.025, 3.025)
            self.h_angley_dist[iroc] = TH1F('h_angley_roc{r}'.format(r=iroc), 'h_angley_roc{r}'.format(r=iroc), 121, -3.025, 3.025)
            self.h_anglex_cut_dist[iroc] = TH1F('h_anglex_cut_roc{r}'.format(r=iroc), 'h_anglex_cut_roc{r}'.format(r=iroc), 121, -3.025, 3.025)
            self.h_angley_cut_dist[iroc] = TH1F('h_angley_cut_roc{r}'.format(r=iroc), 'h_angley_cut_roc{r}'.format(r=iroc), 121, -3.025, 3.025)
            if 'anglex' in self.cut_names and 'angley' in self.cut_names:
                self.analysis.tree.Draw('angle_x>>h_anglex_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['anglex']-1], 'goff') if self.dict_cuts['anglex'] != 0 else self.analysis.tree.Draw('angle_x>>h_anglex_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('angle_y>>h_angley_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['anglex']-1], 'goff') if self.dict_cuts['anglex'] != 0 else self.analysis.tree.Draw('angle_y>>h_angley_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('angle_x>>h_anglex_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['angley']], 'goff')
                self.analysis.tree.Draw('angle_y>>h_angley_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['angley']], 'goff')
            else:
                key = self.get_nearest_prev_existing_cut_key('anglex')
                self.analysis.tree.Draw('angle_x>>h_anglex_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('angle_x>>h_anglex_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('angle_y>>h_angley_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('angle_y>>h_angley_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('angle_x>>h_anglex_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('angle_x>>h_anglex_cut_roc{r}'.format(r=iroc), '', 'goff')
                self.analysis.tree.Draw('angle_y>>h_angley_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('angle_y>>h_angley_cut_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_1D_options('angle', self.h_anglex_dist[iroc], 'angleX(deg)', 'entries', kBlue)
            self.plots.set_1D_options('angle', self.h_angley_dist[iroc], 'angleY(deg)', 'entries', kBlue)
            self.plots.set_1D_options('angle', self.h_anglex_cut_dist[iroc], 'angleX(deg)', 'entries', kRed)
            self.plots.set_1D_options('angle', self.h_angley_cut_dist[iroc], 'angleY(deg)', 'entries', kRed)
            gROOT.SetBatch(False)
            self.plots.save_cuts_distributions(self.h_anglex_dist[iroc], self.h_anglex_cut_dist[iroc], 'angle_roc{r}_x_cut_overlay'.format(r=iroc), 'Angle roc{r} X Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
            self.plots.save_cuts_distributions(self.h_angley_dist[iroc], self.h_angley_cut_dist[iroc], 'angle_roc{r}_y_cut_overlay'.format(r=iroc), 'Angle roc{r} Y Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
        if self.verbose: print 'Done'

    def do_rhit_distribution(self):
        if self.verbose: print 'R_hit...', ; sys.stdout.flush()
        self.h_rhit_dist = {}
        self.h_rhit_cut_dist = {}
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_dist[iroc] = TH1F('h_rhit_roc{r}'.format(r=iroc), 'h_rhit_roc{r}'.format(r=iroc), 101, -0.5, 100.5)
            self.h_rhit_cut_dist[iroc] = TH1F('h_rhit_cut_roc{r}'.format(r=iroc), 'h_rhit_cut_roc{r}'.format(r=iroc), 101, -0.5, 100.5)
            if 'rhit' in self.cut_names:
                self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_roc{d}'.format(n=iroc, d=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['rhit'] - 1], 'goff') if self.dict_cuts['rhit'] != 0 else self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_roc{d}'.format(n=iroc, d=iroc), '', 'goff')
                self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_cut_roc{d}'.format(n=iroc, d=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['rhit']], 'goff')
            else:
                key = self.get_nearest_prev_existing_cut_key('rhit')
                self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_roc{d}'.format(n=iroc, d=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_roc{d}'.format(n=iroc, d=iroc), '', 'goff')
                self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_cut_roc{d}'.format(n=iroc, d=iroc), self.cuts_pixelated_roc_incr[iroc][key], 'goff') if key != -1 else self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_roc{d}'.format(n=iroc, d=iroc), '', 'goff')
            self.plots.set_1D_options('rhit', self.h_rhit_dist[iroc], 'R_Hit(um)', 'entries', kBlue, 0.1)
            self.plots.set_1D_options('rhit', self.h_rhit_cut_dist[iroc], 'R_Hit(um)', 'entries', kRed, 0.1)
            gROOT.SetBatch(False)
            self.plots.save_cuts_distributions(self.h_rhit_dist[iroc], self.h_rhit_cut_dist[iroc], 'rhit_roc{r}_x_cut_overlay'.format(r=iroc), 'Rhit roc{r} x Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False, '', True)
        if self.verbose: print 'Done'
