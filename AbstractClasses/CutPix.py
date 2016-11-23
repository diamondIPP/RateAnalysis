import os
import pickle
import json

import sys
from numpy import array, zeros, arange, delete
from Elementary import Elementary
from ROOT import TCut, gROOT, TH1F, kRed, TCutG, gDirectory, kBlue, TH2D, TH2F, TH1D, kGreen, kMagenta
from collections import OrderedDict
from copy import deepcopy
from math import ceil

class CutPix(Elementary):
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, parent_analysis, verbose=True, skip=False):

        if not skip:
            Elementary.__init__(self, verbose=verbose)
            self.analysis = parent_analysis
            self.run_number = self.analysis.run_number
            self.plot_settings = self.analysis.plots.plot_settings
            # saving stuff
            self.picklepath = 'Configuration/Individual_Configs/cuts/{tc}_{r}.pickle'.format(tc=self.TESTCAMPAIGN, r=self.run_number)
            self.histos = {}
            self.roc_diam1 = 4
            self.roc_diam2 = 5 if self.TESTCAMPAIGN != '201610' else -1
            self.roc_si = 6 if self.TESTCAMPAIGN != '201610' else 5
            self.duts_list = []
            self.list_duts()
            self.roc_tel = [0, 1, 2, 3]
            self.dut_names = {self.roc_diam1: 'ROC4', self.roc_diam2: 'ROC5', self.roc_si: 'ROC6'} if self.TESTCAMPAIGN != '201610' else {self.roc_diam1: 'ROC4', self.roc_si: 'ROC5'}
            # config
            self.DUTType = self.load_dut_type()
            # self.beaminterruptions_folder = self.ana_config_parser.get('CUT', 'beaminterruptions_folder')
            # self.exclude_before_jump = self.ana_config_parser.getint('CUT', 'excludeBeforeJump')
            # self.exclude_after_jump = self.ana_config_parser.getint('CUT', 'excludeAfterJump')
            self.CutConfig = {}
            self.cuts_hitmap_roc = {iROC: {} for iROC in self.duts_list}  # each cut separately
            self.cuts_pixelated_roc = {iROC: {} for iROC in self.duts_list}  # each cut separately
            self.cuts_hitmap_roc_incr = {iROC: {} for iROC in self.duts_list}  # each cut incremental. last position used for all cuts
            self.cuts_pixelated_roc_incr = {iROC: {} for iROC in self.duts_list}  # each cut incremental. last position used for all cuts
            self.num_cuts = 0

            self.events_after_cuts = {}
            self.events_after_cuts_roc = {iROC: {} for iROC in self.duts_list}
            self.events_after_cuts_incr = {}
            self.events_after_cuts_roc_incr = {iROC: {} for iROC in self.duts_list}

            self.cut_names = ['ini_fin', 'beam', 'tracks', 'masks', 'fiducial', 'chi2x', 'chi2y', 'anglex', 'angley', 'rhit']
            self.dict_cuts = {'ini_fin': 0, 'beam': 1, 'tracks': 2, 'masks': 3, 'fiducial': 4, 'chi2x': 5, 'chi2y': 6, 'anglex': 7, 'angley':8, 'rhit': 9}

            self.load_config()
            self.cuts_done = False

            self.plots = self.analysis.plots

            # define cut strings
            # self.EasyCutStrings = self.init_easy_cutstrings()
            # self.CutStrings = self.define_cutstrings()

            # self.region_cut = TCut('region_cut', '')
            # self.JumpCut = TCut('JumpCut', '')

            # beam interrupts
            # self.jumps = None
            # self.jump_ranges = None

    def calculate_initial_events(self):
        gROOT.SetBatch(1)
        self.events_after_cuts['None'] = self.analysis.tree.GetEntries()
        for iROC in self.duts_list:
            self.analysis.tree.Draw('time>>temp0', 'plane[{r}]'.format(r=iROC), 'goff')
            bla = gDirectory.Get('temp0')
            self.events_after_cuts_roc['None'] = bla.GetEntries()

        self.analysis.tree.Draw('time>>temp0', 'plane[{')

    def do_cuts(self):
        # generate cut strings
        self.print_banner('Generating Cut strings...')
        self.generate_ini_fin_cuts()
        self.generate_beam_interruption_cut()
        self.generate_tracks_cut()
        self.generate_masks()
        self.generate_fid_cuts()
        self.generate_chi2_cuts()
        self.generate_angle_cuts()
        self.generate_rhit_cuts()
        # self.gen_incr_vect_cuts()
        self.cuts_done = True
        self.print_banner('Finished generating Cut stringss')
        # self.add_cuts()

        # self.generate_cut_string()  # DA TODO
        # self.all_cut = self.generate_all_cut()  # DA TODO
    def list_duts(self):
        self.duts_list = [self.roc_diam1, self.roc_diam2, self.roc_si] if self.TESTCAMPAIGN != '201610' else [self.roc_diam1, self.roc_si]

    def reset_cuts_dicts(self):
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
        self.print_banner('Doing cuts distributions...')
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
        self.analysis.tree.Draw('time>>h_beam_time', self.cuts_pixelated_roc_incr[self.duts_list[0]][self.dict_cuts['beam']-1],'goff')
        self.analysis.tree.Draw('time>>h_beam_time_cut', self.cuts_pixelated_roc_incr[self.duts_list[0]][self.dict_cuts['beam']],'goff')
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
        self.h_chi2x_dist = {}
        self.h_chi2y_dist = {}
        self.h_chi2x_cut_dist = {}
        self.h_chi2y_cut_dist = {}
        self.h_resy_resx = {}
        self.h_rhit_resx = {}
        self.h_rhit_resy = {}
        self.h_chi2_resx = {}
        self.h_chi2_resy = {}
        self.h_chi2x_resx = {}
        self.h_chi2y_resx = {}
        self.h_resx_hitposy = {}
        self.h_resy_hitposx = {}

        if self.verbose: print 'Chi2...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2x_dist[iroc] = TH1F('h_chi2x_roc{r}'.format(r=iroc), 'h_chi2x_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            self.h_chi2y_dist[iroc] = TH1F('h_chi2y_roc{r}'.format(r=iroc), 'h_chi2y_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            self.h_chi2x_cut_dist[iroc] = TH1F('h_chi2x_cut_roc{r}'.format(r=iroc), 'h_chi2x_cut_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            self.h_chi2y_cut_dist[iroc] = TH1F('h_chi2y_cut_roc{r}'.format(r=iroc), 'h_chi2y_cut_roc{r}'.format(r=iroc), 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_x>>h_chi2x_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2x'] - 1], 'goff')
            self.analysis.tree.Draw('chi2_y>>h_chi2y_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2x'] - 1], 'goff')
            self.analysis.tree.Draw('chi2_x>>h_chi2x_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2y']], 'goff')
            self.analysis.tree.Draw('chi2_y>>h_chi2y_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['chi2y']], 'goff')
            self.plots.set_1D_options('chi2', self.h_chi2x_dist[iroc], 'chi2X', 'entries', kBlue)
            self.plots.set_1D_options('chi2', self.h_chi2y_dist[iroc], 'chi2Y', 'entries', kBlue)
            self.plots.set_1D_options('chi2', self.h_chi2x_cut_dist[iroc], 'chi2X', 'entries', kRed)
            self.plots.set_1D_options('chi2', self.h_chi2y_cut_dist[iroc], 'chi2Y', 'entries', kRed)
            gROOT.SetBatch(False)
            self.plots.save_cuts_distributions(self.h_chi2x_dist[iroc], self.h_chi2x_cut_dist[iroc], 'chi2_roc{r}_x_cut_overlay'.format(r=iroc), 'Chi2 roc{r} x Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
            self.plots.save_cuts_distributions(self.h_chi2y_dist[iroc], self.h_chi2y_cut_dist[iroc], 'chi2_roc{r}_y_cut_overlay'.format(r=iroc), 'Chi2 roc{r} Y Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
        if self.verbose: print 'Done'
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
            self.analysis.tree.Draw('angle_x>>h_anglex_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['anglex']-1], 'goff')
            self.analysis.tree.Draw('angle_y>>h_angley_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['anglex']-1], 'goff')
            self.analysis.tree.Draw('angle_x>>h_anglex_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['angley']], 'goff')
            self.analysis.tree.Draw('angle_y>>h_angley_cut_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['angley']], 'goff')
            self.plots.set_1D_options('angle', self.h_anglex_dist[iroc], 'angleX(deg)', 'entries', kBlue)
            self.plots.set_1D_options('angle', self.h_angley_dist[iroc], 'angleY(deg)', 'entries', kBlue)
            self.plots.set_1D_options('angle', self.h_anglex_cut_dist[iroc], 'angleX(deg)', 'entries', kRed)
            self.plots.set_1D_options('angle', self.h_angley_cut_dist[iroc], 'angleY(deg)', 'entries', kRed)
            gROOT.SetBatch(False)
            self.plots.save_cuts_distributions(self.h_anglex_dist[iroc], self.h_anglex_cut_dist[iroc], 'angle_roc{r}_x_cut_overlay'.format(r=iroc), 'Angle roc{r} X Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
            self.plots.save_cuts_distributions(self.h_angley_dist[iroc], self.h_angley_cut_dist[iroc], 'angle_roc{r}_y_cut_overlay'.format(r=iroc), 'Angle roc{r} Y Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
        if self.verbose: print 'Done'
        if self.verbose: print 'R_hit...', ; sys.stdout.flush()
        self.h_rhit_dist = {}
        self.h_rhit_cut_dist = {}
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_dist[iroc] = TH1F('h_rhit_roc{r}'.format(r=iroc), 'h_rhit_roc{r}'.format(r=iroc), 101, -0.25, 50.25)
            self.h_rhit_cut_dist[iroc] = TH1F('h_rhit_cut_roc{r}'.format(r=iroc), 'h_rhit_cut_roc{r}'.format(r=iroc), 101, -0.25, 50.25)
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_roc{d}'.format(n=iroc, d=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['rhit'] - 1], 'goff')
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>>h_rhit_cut_roc{d}'.format(n=iroc, d=iroc), self.cuts_pixelated_roc_incr[iroc][self.dict_cuts['rhit']], 'goff')
            self.plots.set_1D_options('rhit', self.h_rhit_dist[iroc], 'R_Hit(um)', 'entries', kBlue)
            self.plots.set_1D_options('rhit', self.h_rhit_cut_dist[iroc], 'R_Hit(um)', 'entries', kRed)
            gROOT.SetBatch(False)
            self.plots.save_cuts_distributions(self.h_rhit_dist[iroc], self.h_rhit_cut_dist[iroc], 'rhit_roc{r}_x_cut_overlay'.format(r=iroc), 'Rhit roc{r} x Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
        if self.verbose: print 'Done'
        self.print_banner('Finished with distribution cuts')

        self.print_banner('Doing resolution plots...')
        if self.verbose: print 'Res_Y Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resy_resx[iroc] = TH2D('h_resy_resx_roc{r}'.format(r=iroc), 'h_resy_resx_roc{r}'.format(r=iroc), 41, -1537.5, 1537.5, 41, -1025, 1025)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*residual_ROC{r}_Local_X>>h_resy_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_resy_resx[iroc], 'Res_X(um)', 'Res_y(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resy_resx[iroc], 'h_resy_resx_roc{r}'.format(r=iroc), 'Res_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Rhit Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_resx[iroc] = TH2D('h_rhit_resx_roc{r}'.format(r=iroc), 'h_rhit_resx_roc{r}'.format(r=iroc), 41, -1537.5, 1537.5, 101, -0.25, 50.25)
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_X>>h_rhit_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_rhit_resx[iroc], 'Res_X(um)', 'R_Hit(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_rhit_resx[iroc], 'h_rhit_resx_roc{r}'.format(r=iroc), 'R_Hit Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Rhit Vs Res_Y...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_resy[iroc] = TH2D('h_rhit_resy_roc{r}'.format(r=iroc), 'h_rhit_resy_roc{r}'.format(r=iroc), 41, -1025, 1025, 101, -0.25, 50.25)
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_Y>>h_rhit_resy_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_rhit_resy[iroc], 'Res_Y(um)', 'R_Hit(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_rhit_resy[iroc], 'h_rhit_resy_roc{r}'.format(r=iroc), 'R_Hit Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2 Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2_resx[iroc] = TH2D('h_chi2_resx_roc{r}'.format(r=iroc), 'h_chi2_resx_roc{r}'.format(r=iroc), 41, -1537.5, 1537.5, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_X>>h_chi2_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_chi2_resx[iroc], 'Res_X(um)', 'Chi2', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2_resx[iroc], 'h_chi2_resx_roc{r}'.format(r=iroc), 'Chi2 Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2 Vs Res_Y...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2_resy[iroc] = TH2D('h_chi2_resy_roc{r}'.format(r=iroc), 'h_chi2_resy_roc{r}'.format(r=iroc), 41, -1025, 1025, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_Y>>h_chi2_resy_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_chi2_resy[iroc], 'Res_Y(um)', 'Chi2', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2_resy[iroc], 'h_chi2_resy_roc{r}'.format(r=iroc), 'Chi2 Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2_X Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2x_resx[iroc] = TH2D('h_chi2x_resx_roc{r}'.format(r=iroc), 'h_chi2x_resx_roc{r}'.format(r=iroc), 41, -1537.5, 1537.5, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_x:10000*residual_ROC{r}_Local_X>>h_chi2x_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_chi2x_resx[iroc], 'Res_X(um)', 'Chi2_X', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2x_resx[iroc], 'h_chi2x_resx_roc{r}'.format(r=iroc), 'Chi2_X Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2_Y Vs Res_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2y_resx[iroc] = TH2D('h_chi2y_resx_roc{r}'.format(r=iroc), 'h_chi2y_resx_roc{r}'.format(r=iroc), 41, -1537.5, 1537.5, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_y:10000*residual_ROC{r}_Local_X>>h_chi2y_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_chi2y_resx[iroc], 'Res_X(um)', 'Chi2_Y', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2y_resx[iroc], 'h_chi2y_resx_roc{r}'.format(r=iroc), 'Chi2_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Res_X Vs Hit_Y...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resx_hitposy[iroc] = TH2D('h_resx_hitposy_roc{r}'.format(r=iroc), 'h_resx_hitposy_roc{r}'.format(r=iroc), 161, -4025, 4025, 41, -1537.5, 1537.5)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_X:10000*(residual_ROC{r}_Local_Y+cluster_pos_ROC{r}_Local_Y)>>h_resx_hitposy_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_resx_hitposy[iroc], 'Hit_Y(um)', 'Res_X(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resx_hitposy[iroc], 'h_resx_hitposy_roc{r}'.format(r=iroc), 'Hit_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Res_Y Vs Hit_X...', ; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resy_hitposx[iroc] = TH2D('h_resy_hitposx_roc{r}'.format(r=iroc), 'h_resy_hitposx_roc{r}'.format(r=iroc), 105, -3937.5, 3937.5, 41, -1025, 1025)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*(residual_ROC{r}_Local_X+cluster_pos_ROC{r}_Local_X)>>h_resy_hitposx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts-1], 'goff')
            self.plots.set_2D_options(self.h_resy_hitposx[iroc], 'Hit_X(um)', 'Res_Y(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resy_hitposx[iroc], 'h_resy_hitposx_roc{r}'.format(r=iroc), 'Hit_X Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir+'/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        self.print_banner('Finished with resolution plots')

    def do_cuts_analysis(self):
        self.print_banner('Starting Cuts Analysis...')
        self.do_cuts_distributions()
        self.print_banner('Creating histograms with cuts...')
        self.h_hitmaps_cuts = {}
        self.h_ph1_evt_cuts = {}
        self.h_ph1_cuts = {}
        self.h_ph2_evt_cuts = {}
        self.h_ph2_cuts = {}
        self.h_ph2_cuts = {}
        maxz_hitmap = 0
        maxz_ph1 = 0
        maxz_ph2 = 0

        for iroc in self.duts_list:
            self.h_hitmaps_cuts[iroc] = {}
            self.h_ph1_evt_cuts[iroc] = {}
            self.h_ph1_cuts[iroc] = {}
            self.h_ph2_evt_cuts[iroc] = {}
            self.h_ph2_cuts[iroc] = {}
            phbins = {self.roc_diam1: self.plot_settings['ph1DbinsD4'], self.roc_diam2: self.plot_settings['ph1DbinsD5'], self.roc_si: self.plot_settings['ph1DbinsSi']}
            phmin = {self.roc_diam1: self.plot_settings['ph1DminD4'], self.roc_diam2: self.plot_settings['ph1DminD5'], self.roc_si: self.plot_settings['ph1DminSi']}
            phmax = {self.roc_diam1: self.plot_settings['ph1DmaxD4'], self.roc_diam2: self.plot_settings['ph1DmaxD5'], self.roc_si: self.plot_settings['ph1DmaxSi']}
            phdelta = {self.roc_diam1: phmax[self.roc_diam1] - phmin[self.roc_diam1], self.roc_diam2: phmax[self.roc_diam2] - phmin[self.roc_diam2], self.roc_si: phmax[self.roc_si] - phmin[self.roc_si]}
            for cut in self.cut_names:
                if self.verbose: print 'Analysing ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut), ; sys.stdout.flush()
                gROOT.SetBatch(1)
                self.h_hitmaps_cuts[iroc][cut] = TH2D('hitmap_roc{r}_{c}'.format(r=iroc,c=cut), 'hitmap_roc{r}_{c}'.format(r=iroc,c=cut), self.plot_settings['nBinCol']+1, self.plot_settings['minCol']-(self.plot_settings['maxCol']-self.plot_settings['minCol'])/(2*float(self.plot_settings['nBinCol'])), self.plot_settings['maxCol']+(self.plot_settings['maxCol']-self.plot_settings['minCol'])/(2*float(self.plot_settings['nBinCol'])), self.plot_settings['nBinRow']+1, self.plot_settings['minRow']-(self.plot_settings['maxRow']-self.plot_settings['minRow'])/(2*float(self.plot_settings['nBinRow'])), self.plot_settings['maxRow']+(self.plot_settings['maxRow']-self.plot_settings['minRow'])/(2*float(self.plot_settings['nBinRow'])))
                self.h_ph1_evt_cuts[iroc][cut] = TH2D('ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), self.plot_settings['event_bins']+1, self.plot_settings['event_min']-(self.plot_settings['event_max']-self.plot_settings['event_min'])/(2*float(self.plot_settings['event_bins'])), self.plot_settings['event_max']+(self.plot_settings['event_max']-self.plot_settings['event_min'])/(2*float(self.plot_settings['event_bins'])), phbins[iroc]+1, phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])), phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.h_ph1_cuts[iroc][cut] = TH1D('ph1_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_roc{r}_{c}'.format(r=iroc,c=cut), phbins[iroc]+1, phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])), phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.h_ph2_evt_cuts[iroc][cut] = TH2D('ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), self.plot_settings['event_bins']+1, self.plot_settings['event_min']-(self.plot_settings['event_max']-self.plot_settings['event_min'])/(2*float(self.plot_settings['event_bins'])), self.plot_settings['event_max']+(self.plot_settings['event_max']-self.plot_settings['event_min'])/(2*float(self.plot_settings['event_bins'])), phbins[iroc]+1, phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])), phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.h_ph2_cuts[iroc][cut] = TH1D('ph2_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_roc{r}_{c}'.format(r=iroc,c=cut), phbins[iroc]+1, phmin[iroc]-phdelta[iroc]/(2*float(phbins[iroc])), phmax[iroc]+phdelta[iroc]/float(2*phbins[iroc]))
                self.analysis.tree.Draw('row:col >> hitmap_roc{r}_{c}'.format(r=iroc,c=cut), 'plane=={r}&&{cu}'.format(r=iroc,cu=self.cuts_hitmap_roc_incr[iroc][self.dict_cuts[cut]]), 'goff')
                if maxz_hitmap < self.h_hitmaps_cuts[iroc][cut].GetBinContent(self.h_hitmaps_cuts[iroc][cut].GetMaximumBin()): maxz_hitmap = self.h_hitmaps_cuts[iroc][cut].GetBinContent(self.h_hitmaps_cuts[iroc][cut].GetMaximumBin())
                self.analysis.tree.Draw('charge_all_ROC{r}:event_number >> ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'cluster_size_ROC{r}==1&&{cu}'.format(r=iroc,cu=self.cuts_pixelated_roc_incr[iroc][self.dict_cuts[cut]]), 'goff')
                if maxz_ph1 < self.h_ph1_evt_cuts[iroc][cut].GetBinContent(self.h_ph1_evt_cuts[iroc][cut].GetMaximumBin()): maxz_ph1 = self.h_ph1_evt_cuts[iroc][cut].GetBinContent(self.h_ph1_evt_cuts[iroc][cut].GetMaximumBin())
                self.analysis.tree.Draw('charge_all_ROC{r}:event_number >> ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'cluster_size_ROC{r}==2&&{cu}'.format(r=iroc,cu=self.cuts_pixelated_roc_incr[iroc][self.dict_cuts[cut]]), 'goff')
                if maxz_ph2 < self.h_ph2_evt_cuts[iroc][cut].GetBinContent(self.h_ph2_evt_cuts[iroc][cut].GetMaximumBin()): maxz_ph2 = self.h_ph2_evt_cuts[iroc][cut].GetBinContent(self.h_ph2_evt_cuts[iroc][cut].GetMaximumBin())
                self.h_ph1_evt_cuts[iroc][cut].ProjectionY('ph1_roc{r}_{c}'.format(r=iroc,c=cut),0,-1,'e')
                self.h_ph2_evt_cuts[iroc][cut].ProjectionY('ph2_roc{r}_{c}'.format(r=iroc,c=cut),0,-1,'e')
                gROOT.SetBatch(0)
                if self.verbose: print 'Done'
        for iroc in self.duts_list:
            for cut in self.cut_names:
                if self.verbose: print 'Saving for ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut), ; sys.stdout.flush()
                self.plots.set_2D_options(self.h_hitmaps_cuts[iroc][cut], 'col', 'row', 'entries', max_val=maxz_hitmap)
                self.plots.set_2D_options(self.h_ph1_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries', max_val=maxz_ph1)
                self.plots.set_1D_options('ph',self.h_ph1_cuts[iroc][cut],'ph 1 pix cl (e)', 'entries')
                self.plots.set_2D_options(self.h_ph2_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries', max_val=maxz_ph2)
                self.plots.set_1D_options('ph',self.h_ph2_cuts[iroc][cut],'ph 2 pix cl (e)', 'entries')
                self.plots.save_individual_plots(self.h_hitmaps_cuts[iroc][cut], 'hitmap_roc{r}_{c}'.format(r=iroc,c=cut), 'hitmap_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts', doLogZ=True)
                self.plots.save_individual_plots(self.h_ph1_evt_cuts[iroc][cut], 'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_evt_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                self.plots.save_individual_plots(self.h_ph2_evt_cuts[iroc][cut], 'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_evt_roc{r}_{c}'.format(r=iroc,c=cut), None, 'colz', 1, self.plots.save_dir+'/cuts')
                self.plots.save_individual_plots(self.h_ph1_cuts[iroc][cut], 'ph1_roc{r}_{c}'.format(r=iroc,c=cut), 'ph1_roc{r}_{c}'.format(r=iroc,c=cut), None, '', 1, self.plots.save_dir+'/cuts')
                self.plots.save_individual_plots(self.h_ph2_cuts[iroc][cut], 'ph2_roc{r}_{c}'.format(r=iroc,c=cut), 'ph2_roc{r}_{c}'.format(r=iroc,c=cut), None, '', 1, self.plots.save_dir+'/cuts')
                if self.verbose: print 'Done'
        self.print_banner('Finished Cut Analysis', ':)')

    def generate_all_cut(self):
        cut = TCut('all_cuts', '')
        for key, value in self.CutStrings.iteritems():
            if not key.startswith('old') and not key.startswith('all_cut'):
                cut += value
        return cut

    def get_included_events(self, maxevent=None):
        """
        :param maxevent:
        :return: list of included event numbers not excluded by: excludeFirst, EventRange or BeamInterruptions
        """
        minevent = self.get_min_event()
        maxevent = self.get_max_event() if maxevent is None else maxevent

        excluded = [i for i in arange(0, minevent)]  # first events
        for start, stop in zip(self.jump_ranges['start'], self.jump_ranges['stop']):
            excluded += [i for i in xrange(start, stop + 1)]  # events around jumps
        excluded.sort()
        all_events = arange(0, maxevent)
        included = delete(all_events, excluded)
        return included

    @staticmethod
    def init_easy_cutstrings():
        dic = OrderedDict()
        dic['IndividualChCut'] = ''
        dic['EventRange'] = ''
        dic['noPulser'] = ''
        dic['ExcludeFirst'] = ''
        dic['notSaturated'] = ''
        dic['noBeamInter'] = ''
        dic['Tracks'] = ''
        dic['peakPos_high'] = ''
        dic['spread_low'] = ''
        dic['absMedian_high'] = ''
        dic['pedestalsigma'] = ''
        return dic

    @staticmethod
    def define_cutstrings():
        dic = OrderedDict()
        dic['raw'] = TCut('raw', '')
        dic['pulser'] = TCut('pulser', '')
        dic['event_range'] = TCut('event_range', '')
        # waveform
        dic['beam_interruptions'] = TCut('beam_interruptions', '')
        dic['ped_sigma'] = TCut('ped_sigma', '')
        dic['spread_low'] = TCut('spread_low', '')
        dic['median'] = TCut('median', '')
        # tracks
        dic['tracks'] = TCut('tracks', '')
        dic['chi2X'] = TCut('chi2X', '')
        dic['chi2Y'] = TCut('chi2Y', '')
        dic['track_angle'] = TCut('track_angle', '')
        # waveform
        dic['saturated'] = TCut('saturated', '')
        dic['signal_peak_pos'] = TCut('signal_peak_pos', '')
        dic['trigger_cell'] = TCut('trigger_cell', '')
        dic['old_bucket'] = TCut('old_bucket', '')
        dic['bucket'] = TCut('bucket', '')
        dic['all_cuts'] = TCut('all_cuts', '')
        return dic

    # ==============================================
    # region GET CONFIG

    def load_dut_type(self):
        dut_type = self.run_config_parser.get("BASIC", "type")
        assert dut_type.lower() in ["pixel", "pad"], "The DUT type {0} should be 'pixel' or 'pad'".format(dut_type)
        return dut_type

    def load_config(self):
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
        self.CutConfig['rhit'] = self.ana_config_parser.getint('CUT', 'rhit') if self.ana_config_parser. \
            has_option('CUT', 'rhit') else ''
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

    # def add_cuts(self):
    #     for iROC in xrange(4, 7):
    #         self.cuts_hitmap_roc[iROC] = self.mask_hitmap_roc[iROC] + self.chi2x_cut + self.chi2y_cut \
    #                                      + self.cut_tracks + self.angle_x_cut + self.angle_y_cut + self.ini_fin_cut \
    #                                      + beam_interr_cut
    #         self.cuts_pixelated_roc[iROC] = self.mask_pixelated_roc[iROC] + self.chi2x_cut + self.chi2y_cut \
    #                                         + self.cut_tracks + self.angle_x_cut + self.angle_y_cut + self.ini_fin_cut \
    #                                         + beam_interr_cut + self.rhit_cut[iROC]

    def generate_ini_fin_cuts(self):
        if self.verbose: print 'Creating cut for initial and final', abs(self.CutConfig['ExcludeFirst']), 'seconds...', ; sys.stdout.flush()
        picklepath = 'Configuration/Individual_Configs/IniFin/{tc}_{r}.pickle'.format(tc=self.TESTCAMPAIGN, r=self.run_number)
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

    def generate_beam_interruption_cut(self):
        # time is in ms. good results found with bin size of 5 seconds
        picklepath = 'Configuration/Individual_Configs/Beam/{tc}_{r}.pickle'.format(tc=self.TESTCAMPAIGN, r=self.run_number)
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
            self.analysis.tree.Draw('time>>h_beam_time_', self.cuts_pixelated_roc_incr[self.duts_list[0]][self.num_cuts-1],'goff')
            gROOT.SetBatch(False)
            mean = h1.Integral()/float(h1.GetNbinsX())
            beam_interr_cut = TCut('beam_interruptions_cut', '')
            vector_interr = {}
            max_interr = 0
            for t in xrange(1, h1.GetNbinsX()+1):
                if h1.GetBinContent(t) < mean*0.9 or h1.GetBinContent(t) > mean*1.2:
                    if t != 1:
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
        # self.h_beam_time_cut = TH1F('h_beam_time_cut', 'h_beam_time_cut', self.h_beam_time.GetXaxis().GetNbins(), self.h_beam_time.GetXaxis().GetXmin(), self.h_beam_time.GetXaxis().GetXmax())
        # self.analysis.tree.Draw('time>>h_beam_time_cut', self.cuts_pixelated_roc_incr[self.duts_list[0]][self.num_cuts-1],'goff')
        # self.plots.set_1D_options('time', self.h_beam_time, 'time (ms)', 'events', kBlue)
        # self.plots.set_1D_options('time', self.h_beam_time_cut, 'time (ms)', 'events', kRed)
        # gROOT.SetBatch(False)
        # self.plots.save_cuts_distributions(self.h_beam_time, self.h_beam_time_cut, 'Beam_Interruptions_cut', 'Beam_Interruptions_cut', '', 1000000011, self.plots.save_dir+'/cuts', False)


    def generate_rhit_cuts(self):
        if self.verbose: print 'Generatin R-hit cut for distances greater than', self.CutConfig['rhit'], 'um between predicted position from the track and the seed cluster position...', ; sys.stdout.flush()
        self.rhit_cut = {}
        # self.h_rhit = {}
        # self.h_rhit_cut = {}
        for iroc in self.duts_list:
            self.generate_rhit_cuts_DUT(iroc)
            self.gen_vect_cuts(self.rhit_cut[iroc], self.rhit_cut[iroc], iroc)
        self.num_cuts += 1
        if self.verbose: print 'Done'
        # for iroc in self.duts_list:
        #     gROOT.SetBatch(1)
        #     self.analysis.tree.Draw('(10000*sqrt((track_x_ROC{n}-cluster_pos_ROC{n}_Telescope_X)**2+(track_y_ROC{n}-cluster_pos_ROC{n}_Telescope_Y)**2))>>h_rhit_ROC{d}_cut'.format(n=iroc, d=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1], 'goff')
        #     gROOT.SetBatch(0)
        #     self.plots.save_cuts_distributions(self.h_rhit[iroc], self.h_rhit_cut[iroc], 'rhit_ROC{r}cut_overlay'.format(r=iroc), 'R_Hit ROC{r} cuts Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)


    def generate_rhit_cuts_DUT(self, dut):
        picklepath = 'Configuration/Individual_Configs/RHitRoc{ir}/{tc}_{r}.pickle'.format(ir=dut, tc=self.TESTCAMPAIGN, r=self.run_number)
        def func0():
            # gROOT.SetBatch(1)
            # h_rhit = TH1F('h_rhit_', 'h_rhit_', 201, -5, 2005)
            # self.h_rhit_cut[dut] = TH1F('h_rhit_ROC{d}_cut'.format(d=dut), 'h_rhit_ROC{d}_cut'.format(d=dut), 201, -5, 2005)
            # self.plots.set_1D_options('rhit', self.h_rhit[dut], 'R_hit(um)', 'entries', kBlue)
            # self.plots.set_1D_options('rhit', self.h_rhit_cut[dut], 'R_hit(um)', 'entries', kRed)
            # self.analysis.tree.Draw('(10000*sqrt((track_x_ROC{n}-cluster_pos_ROC{n}_Telescope_X)**2+(track_y_ROC{n}-cluster_pos_ROC{n}_Telescope_Y)**2))>>h_rhit_ROC{d}'.format(n=dut, d=dut), self.cuts_pixelated_roc_incr[dut][self.num_cuts - 1], 'goff')
            # self.analysis.tree.Draw('sqrt((10*(track_x_ROC{n}-cluster_pos_ROC{n}_Telescope_X))**2+(10*(track_x_ROC{n}-cluster_pos_ROC{n}_Telescope_Y))**2)>>h'.format(n=dut),'','goff')
            # h.GetQuantiles(nq, rhits, xq)
            # gROOT.SetBatch(0)
            value = self.CutConfig['rhit']
            # string=''
            string = '((10000*sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))<{val}&&(sqrt((residual_ROC{n}_Local_X)**2+(residual_ROC{n}_Local_Y)**2))>=0)'.format(n=dut, val=value)
            return string
        self.rhit_cut[dut] = self.do_pickle(picklepath, func0)

    def generate_tracks_cut(self):
        if self.verbose: print 'Creating tracks cut...',
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
        # for iroc in self.duts_list:
        #     gROOT.SetBatch(1)
        #     self.analysis.tree.Draw('chi2_x>>h_chi2_roc{r}_x_cut'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1], 'goff')
        #     self.analysis.tree.Draw('chi2_y>>h_chi2_roc{r}_y_cut'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1], 'goff')
        #     self.plots.set_1D_options('chi2', self.h_chi2[iroc]['x'], 'chi2X', 'entries', kBlue)
        #     self.plots.set_1D_options('chi2', self.h_chi2[iroc]['y'], 'chi2Y', 'entries', kBlue)
        #     self.plots.set_1D_options('chi2', self.h_chi2_cut[iroc]['x'], 'chi2X', 'entries', kRed)
        #     self.plots.set_1D_options('chi2', self.h_chi2_cut[iroc]['y'], 'chi2Y', 'entries', kRed)
        #     gROOT.SetBatch(0)
        #     self.plots.save_cuts_distributions(self.h_chi2[iroc]['x'], self.h_chi2_cut[iroc]['x'], 'chi2_roc{r}_x_cut_overlay'.format(r=iroc), 'Chi2 roc{r} x Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)
        #     self.plots.save_cuts_distributions(self.h_chi2[iroc]['y'], self.h_chi2_cut[iroc]['y'], 'chi2_roc{r}_y_cut_overlay'.format(r=iroc), 'Chi2 roc{r} y Cut Overlay'.format(r=iroc), '', 1000000011, self.plots.save_dir+'/cuts', False)

    def generate_chi2(self, mode='x', num_prev_cut=4, iroc=4):
        picklepath = 'Configuration/Individual_Configs/Chi2Roc{ir}{m}/{tc}_{r}.pickle'.format(ir=iroc, m=mode, tc=self.TESTCAMPAIGN, r=self.run_number)
        def func0():
            gROOT.SetBatch(1)
            h_chi2 = TH1F('h_chi2_roc{r}_{m}_'.format(r=iroc,m=mode), 'h_chi2_roc{r}_{m}_'.format(r=iroc,m=mode), 201, -0.1, 40.1)
            # h_chi2_cut[iroc][mode] = TH1F('h_chi2_roc{r}_{m}_cut'.format(r=iroc,m=mode), 'h_chi2_roc{r}_{m}_cut'.format(r=iroc,m=mode), 201, -0.1, 40.1)
            nq = 100
            chi2s = zeros(nq)
            xq = array([(i + 1) / float(nq) for i in xrange(nq)])
            self.analysis.tree.Draw('chi2_{mod}>>h_chi2_roc{r}_{mod}_'.format(r=iroc,mod=mode), self.cuts_pixelated_roc_incr[iroc][num_prev_cut], 'goff')
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
        # if cut_hitmap != 0:
        self.cuts_hitmap_roc[roc][self.num_cuts] = cut_hitmap
        self.cuts_pixelated_roc[roc][self.num_cuts] = cut_pixelated
        self.accum_incr_vect_cuts(roc)

    def accum_incr_vect_cuts(self, roc=4):
        self.cuts_hitmap_roc_incr[roc][self.num_cuts] = self.cuts_hitmap_roc_incr[roc][self.num_cuts-1] + '&&(' + self.cuts_hitmap_roc[roc][self.num_cuts] + ')' if self.num_cuts != 0 else '(' + self.cuts_hitmap_roc[roc][self.num_cuts] + ')'
        self.cuts_pixelated_roc_incr[roc][self.num_cuts] = self.cuts_pixelated_roc_incr[roc][self.num_cuts-1] + '&&(' + self.cuts_pixelated_roc[roc][self.num_cuts] + ')' if self.num_cuts != 0 else '(' + self.cuts_pixelated_roc[roc][self.num_cuts] + ')'

    # def gen_incr_vect_cuts(self):  # creates incremental cuts in a vector including the fiducial cut. For this reason it is given as a string and not as a TCut
    #     ifiducial = self.num_cuts
    #     for roc in self.duts_list:
    #         for i in range(self.num_cuts):
    #             if self.dict_cuts['fiducial'] != i:
    #                 self.cuts_hitmap_roc_incr[roc][i] = TCut('cut_hit_incr_roc_{r}_pos_{ii}'.format(r=roc, ii=i), '')
    #                 self.cuts_pixelated_roc_incr[roc][i] = TCut('cut_pix_incr_roc_{r}_pos_{ii}'.format(r=roc, ii=i), '')
    #                 for j in range(i+1):
    #                     # if i != self.num_cuts - 1:
    #                     if self.dict_cuts['fiducial'] != j:
    #                         self.cuts_hitmap_roc_incr[roc][i] = self.cuts_hitmap_roc_incr[roc][i] + self.cuts_hitmap_roc[roc][j]
    #                         self.cuts_pixelated_roc_incr[roc][i] = self.cuts_pixelated_roc_incr[roc][i] + self.cuts_pixelated_roc[roc][j]
    #             else:
    #                 ifiducial = i
    #     for roc in self.duts_list:
    #         if ifiducial != self.num_cuts:
    #             for i in range(ifiducial):
    #                 self.cuts_hitmap_roc_incr[roc][i] = self.cuts_hitmap_roc_incr[roc][i].GetTitle()
    #                 self.cuts_pixelated_roc_incr[roc][i] = self.cuts_pixelated_roc_incr[roc][i].GetTitle()
    #             self.cuts_hitmap_roc_incr[roc][ifiducial] = self.cuts_hitmap_roc_incr[roc][ifiducial - 1].GetTitle()
    #             self.cuts_hitmap_roc_incr[roc][ifiducial] = self.cuts_hitmap_roc_incr[roc][ifiducial] + '&&(fidcut_hitmap_roc{n})'.format(n=roc)
    #             self.cuts_pixelated_roc_incr[roc][ifiducial] = self.cuts_pixelated_roc_incr[roc][ifiducial - 1].GetTitle()
    #             self.cuts_pixelated_roc_incr[roc][ifiducial] = self.cuts_pixelated_roc_incr[roc][ifiducial] + '&&(fidcut_pixelated_roc{n})'.format(n=roc)
    #             for i in range(ifiducial + 1, self.num_cuts):
    #                 self.cuts_hitmap_roc_incr[roc][i] = self.cuts_hitmap_roc_incr[roc][i].GetTitle() + '&&(fidcut_hitmap_roc{n})'.format(n=roc)
    #                 self.cuts_pixelated_roc_incr[roc][i] = self.cuts_pixelated_roc_incr[roc][i].GetTitle() + '&&(fidcut_pixelated_roc{n})'.format(n=roc)
    #         else:
    #             for i in range(self.num_cuts):
    #                 self.cuts_hitmap_roc_incr[roc][i] = self.cuts_hitmap_roc_incr[roc][i].GetTitle()
    #                 self.cuts_pixelated_roc_incr[roc][i] = self.cuts_pixelated_roc_incr[roc][i].GetTitle()

    def generate_masks(self):
        if self.verbose: print 'generating masks cuts...', ; sys.stdout.flush()
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
        angle = self.CutConfig['track_angle']

        # def func():
        #     print 'generating slope cut for run {run}...'.format(run=self.analysis.run_number)
        # fit the slope to get the mean
        # gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        def func0():
            gROOT.SetBatch(1)
            h = TH1F('h_angle_', 'h_angle_', 61, -3.05, 3.05)
            self.analysis.tree.Draw('angle_{x}>>h_angle_'.format(x=mode), self.cuts_pixelated_roc_incr[iroc][prev_num_cut], 'goff')
            h_angle = TH1F('h_angle_roc_', 'h_angle_roc_', 51, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) - 3*h.GetRMS() - 3*h.GetRMS()/50, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) + 3*h.GetRMS() + 3*h.GetRMS()/50)
            # self.h_angle_cut[iroc][mode] = TH1F('h_angle_roc{r}_{m}_cut'.format(r=iroc,m=mode), 'h_angle_roc{r}_{m}_cut'.format(r=iroc,m=mode), 51, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) - 2*h.GetRMS() - 2*h.GetRMS()/50, h.GetXaxis().GetBinCenter(h.GetMaximumBin()) + 2*h.GetRMS() + 2*h.GetRMS()/50)
            self.analysis.tree.Draw('angle_{m}>>h_angle_roc_'.format(m=mode), self.cuts_pixelated_roc_incr[iroc][prev_num_cut], 'goff')
            fit_result = h_angle.Fit('gaus', 'qs', '')
            # fit_result = h.Fit('gaus', 'qs')# , '', xmin, xmax)
            x_mean = fit_result.Parameter(1)
            angles = [x_mean - angle, x_mean + angle]
            # c = gROOT.FindObject('c1')
            # c.Close()
            gROOT.SetBatch(0)
            # gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
            return angles

        angles = self.do_pickle(picklepath, func0)
        # create the cut string
        string = 'angle_{x}>={minx}&&angle_{x}<={maxx}'.format(x=mode, minx=angles[0], maxx=angles[1])
        self.angle_cut[iroc][mode] = string if angle > 0 else ''

    def load_event_range(self, event_range=None):
        """
        Gets the event range cut. If the arguments are negative, they are interpreted as time in minutes. Therefore, e.g.
        load_event_range(-10, 700000) means that only events are considered, which fulfill: >10 minutes after run start event number < 700000
        :param event_range:
        :return: event range
        """
        if event_range is None:
            event_range = [0, 0]
        for i, value in enumerate(event_range):
            if value < 0:
                event_range[i] = self.analysis.get_event_at_time(time_sec=-1 * value * 60)
        if not event_range[1]:
            event_range[1] = self.analysis.get_event_at_time(-1)
        if not event_range[0]:
            event_range[0] = self.CutConfig['ExcludeFirst']
        return event_range

    def set_event_range(self, event_range):
        self.CutConfig['EventRange'] = self.load_event_range(event_range)

    def load_exclude_first(self, value):
        """
        Sets how many events at the very beginning of the run should be excluded. if the argument is negative, it will be interpreted as time in minutes. For a positive argument it is interpreted as
        maximum event number.
        :param value: events or time in minutes
        :return:
        """
        if value > 0:
            self.EasyCutStrings['ExcludeFirst'] = str(int(value) / 1000) + 'k+'
            return value
        elif value == 0:
            self.EasyCutStrings['ExcludeFirst'] = ''
            return 0
        else:
            self.EasyCutStrings['ExcludeFirst'] = str(-1 * value) + 'min+'
            seconds = -1 * value * 60
            event = self.analysis.get_event_at_time(seconds)
            return event

    def set_exclude_first(self, value):
        self.CutConfig['ExcludeFirst'] = self.load_exclude_first(value)

    def load_peakpos_high(self, high):
        if high > 0:
            self.EasyCutStrings['peakPos_high'] = 'peakPos<{high}'.format(high=high)
            return high
        else:
            return -1

    def set_peakpos_high(self, value):
        self.CutConfig['peakPos_high'] = self.load_peakpos_high(value)

    def load_spread_low(self, value):
        if value > 0:
            self.EasyCutStrings['spread_low'] = 'spread>{low}'.format(low=value)
            return value
        else:
            return -1

    # endregion

    def get_event_range(self):
        """
        Returns a the lowest and highest event numbers to consider in the analysis.
        :return: cut eventrange as list, empty if no cut applied
        """
        return self.CutConfig["EventRange"]

    def get_min_event(self):
        """ :return: the smallest event number satisfying the cut conditions. """
        return self.CutConfig["EventRange"][0]

    def get_n_events(self):
        """ :return: number of events in EventRange """
        total_events = self.analysis.get_event_at_time(-1)
        return total_events if not self.CutConfig["EventRange"] else self.CutConfig["EventRange"][1] - self.CutConfig["EventRange"][0]

    def get_max_event(self):
        """ :return: maximum event number """
        return self.CutConfig["EventRange"][1]

    # ==============================================
    # region GENERATE CUT STRINGS
    def generate_event_range(self):
        if self.CutConfig['EventRange']:
            self.CutStrings['event_range'] += '(event_number<={max}&&event_number>={min})'.format(min=self.CutConfig['EventRange'][0], max=self.CutConfig['EventRange'][1])
        elif self.CutConfig['ExcludeFirst']:
            self.CutStrings['event_range'] += 'event_number>={min}'.format(min=self.CutConfig['ExcludeFirst'])

    # def generate_chi2(self, mode='x'):
    #     picklepath = 'Configuration/Individual_Configs/Chi2/{tc}_{run}_{mod}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.run.run_number, mod=mode.title())
    #
    #     def func():
    #         print 'generating chi2 cut in {mod} for run {run}...'.format(run=self.analysis.run_number, mod=mode)
    #         gROOT.SetBatch(1)
    #         h = TH1F('h', '', 200, 0, 100)
    #         nq = 100
    #         chi2s = zeros(nq)
    #         xq = array([(i + 1) / float(nq) for i in range(nq)])
    #         self.analysis.tree.Draw('chi2_{mod}>>h'.format(mod=mode), '', 'goff')
    #         h.GetQuantiles(nq, chi2s, xq)
    #         gROOT.SetBatch(0)
    #         return chi2s
    #
    #     chi2 = self.do_pickle(picklepath, func)
    #     quantile = self.CutConfig['chi2{mod}'.format(mod=mode.title())]
    #     assert type(quantile) is int and 0 < quantile <= 100, 'chi2 quantile has to be and integer between 0 and 100'
    #     string = 'chi2_{mod}<{val}&&chi2_{mod}>=0'.format(val=chi2[quantile], mod=mode)
    #     return string if quantile > 0 else ''



    def generate_cut_string(self):
        """ Creates the cut string. """
        gROOT.SetBatch(1)

        # --TRACKS --
        # self.CutStrings['chi2X'] += self.generate_chi2('x')
        # self.CutStrings['chi2Y'] += self.generate_chi2('y')
        # self.CutStrings['track_angle'] += self.generate_slope()
        # self.CutStrings['tracks'] += 'n_tracks'

        # -- EVENT RANGE CUT --
        # self.generate_event_range()
        # if self.CutConfig['EventRange']:
        #     self.EasyCutStrings['EventRange'] = 'Evts.{min}k-{max}k'.format(min=int(self.CutConfig['EventRange'][0]) / 1000, max=int(self.CutConfig['EventRange'][1]) / 1000)
        #     self.EasyCutStrings['ExcludeFirst'] = 'Evts.{min}k+'.format(min=int(self.CutConfig['ExcludeFirst']) / 1000) if self.CutConfig['ExcludeFirst'] > 0 else ''

        # -- PULSER CUT --
        # self.CutStrings['pulser'] += '!pulser'

        # -- BEAM INTERRUPTION CUT --
        # self.__generate_beam_interruptions()
        # self.EasyCutStrings['noBeamInter'] = 'BeamOn'
        # self.generate_jump_cut()

        # -- MASK PIXELS --


        # -- FIDUCIAL REGION --

        gROOT.SetBatch(0)

    def __generate_beam_interruptions(self):
        """
        This adds the restrictions to the cut string such that beam interruptions are excluded each time the cut is applied.
        """
        self.get_beam_interruptions()

        njumps = len(self.jump_ranges["start"])
        cut_string = ''
        start_event = self.CutConfig['EventRange'][0]
        for i in xrange(njumps):
            upper = self.jump_ranges["stop"][i]
            lower = self.jump_ranges["start"][i]
            if upper > start_event:
                lower = start_event if lower < start_event else lower
                string = "!(event_number<={up}&&event_number>={low})".format(up=upper, low=lower)
                # new separate strings
                if cut_string != '':
                    cut_string += '&&'
                cut_string += string
        self.CutStrings['beam_interruptions'] += cut_string
    # endregion

    # ==============================================
    # region BEAM INTERRUPTS
    def generate_jump_cut(self):
        cut_string = ''
        start_event = self.CutConfig['EventRange'][0]
        for tup in self.jumps:
            if tup[1] > start_event:
                low = start_event if tup[0] < start_event else tup[0]
                cut_string += '&&' if cut_string else ''
                cut_string += '!(event_number<={up}&&event_number>={low})'.format(up=tup[1], low=low)
        self.JumpCut += cut_string

    def find_beam_interruptions(self):
        return self.find_pad_beam_interruptions() if self.DUTType == 'pad' else self.find_pixel_beam_interruptions()

    def find_pad_beam_interruptions(self):
        """
        Looking for the beam interruptions by investigating the pulser rate.
        :return: interrupt list
        """
        print 'Searching for beam interruptions...'
        binning = 200
        nbins = int(self.analysis.run.tree.GetEntries()) / binning
        rate = []
        for i in xrange(nbins):
            pulserevents = self.analysis.run.tree.Draw('1', 'pulser', 'goff', binning, i * binning)
            rate.append(100 * pulserevents / binning)
        interrupts = []
        last_rate = 0
        tup = [0, 0]
        cut = 30  # if rate goes higher than n %
        for i, value in enumerate(rate):
            if value > cut > last_rate:
                tup[0] = i * binning
            elif value < cut < last_rate:
                tup[1] = i * binning
                interrupts.append(tup)
                tup = [0, 0]
            last_rate = value
        return interrupts

    def find_pixel_beam_interruptions(self):
        # todo DA, just return the same format as find_pad_beam_interruptions does
        pass

    def __save_beaminterrupts(self):
        # check if directories exist
        if not os.path.exists(self.beaminterruptions_folder):
            os.mkdir(self.beaminterruptions_folder)
        if not os.path.exists(self.beaminterruptions_folder + '/data'):
            os.mkdir(self.beaminterruptions_folder + '/data')

        # save jump list to file
        jumpfile = open(self.beaminterruptions_folder + '/data/{testcampaign}Run_{run}.pickle'.format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number), 'wb')
        pickle.dump(self.jumps, jumpfile)
        jumpfile.close()

    def __create_jump_ranges(self):
        if self.jump_ranges is None and len(self.jumps) > 0:
            print 'generating jump ranges...'
            start = []
            stop = []
            time_offset = self.analysis.run.get_time_at_event(0)
            t_max = (self.analysis.run.get_time_at_event(-1) - time_offset) / 1000.
            last_stop = 0
            for tup in self.jumps:
                t_start = (self.analysis.run.get_time_at_event(tup[0]) - time_offset) / 1000.
                t_stop = (self.analysis.run.get_time_at_event(tup[1]) - time_offset) / 1000.
                # add offsets from config file
                t_start -= -1 * self.exclude_before_jump if t_start >= -1 * self.exclude_before_jump else 0
                t_stop = t_stop + -1 * self.exclude_after_jump if t_stop + -1 * self.exclude_after_jump <= t_max else t_max
                if t_start < last_stop:
                    stop[-1] = self.analysis.get_event_at_time(t_stop)
                    last_stop = t_stop
                    continue
                start.append(self.analysis.get_event_at_time(t_start))
                stop.append(self.analysis.get_event_at_time(t_stop))
                last_stop = t_stop

            self.jump_ranges = {"start": start,
                                "stop": stop}

        return [self.exclude_before_jump, self.exclude_after_jump, self.jump_ranges]

    def get_beam_interruptions(self):
        """
        If beam interruption data exist in beaminterruptions/data/, it will load it in order to account for beam interruptions. The data is stored as a list of jumps, dumped into a pickle file.
        If no pickle file exists, it will perform a beam interruption analysis in order to identify the beam interruptions. The found interruptions are stored in a list at .jumps and dumped into
        a pickle file.
        :return: list of events where beam interruptions occures
        """
        if self.jump_ranges is None:
            jumps_pickle = self.beaminterruptions_folder + "/data/{testcampaign}Run_{run}.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
            range_pickle = self.beaminterruptions_folder + "/data/{testcampaign}_{run}_Jump_Ranges.pickle".format(testcampaign=self.TESTCAMPAIGN, run=self.analysis.run.run_number)
            self.jumps = self.do_pickle(jumps_pickle, self.find_beam_interruptions)
            ranges = self.do_pickle(range_pickle, self.__create_jump_ranges)
            # redo range pickle if config parameters have changed
            if ranges[0] != self.exclude_before_jump or ranges[1] != self.exclude_after_jump:
                os.remove(range_pickle)
                ranges = self.do_pickle(range_pickle, self.__create_jump_ranges)
            self.jump_ranges = ranges[2]
        return self.jumps
    # endregion

    def get_easy_cutstring(self):
        """
        Returns a short, more user-friendly cut string, which can be used to display the cut configuration as terminal prompt or inside a canvas.
        :return:
        """
        string_ = ""
        for type_ in self.EasyCutStrings.keys():
            if self.EasyCutStrings[type_] != "":
                string_ += self.EasyCutStrings[type_] + ", "
        if string_ != "":
            string_ = string_[:-2]
        return string_

    def reset_cut(self, name):
        if name in self.CutStrings:
            self.CutStrings[name].SetTitle('')
        else:
            print 'There is no cut with the name "{name}"!'.format(name=name)
        self.all_cut = self.generate_all_cut()

    def show_cuts(self, easy=True):
        cuts = self.EasyCutStrings if easy else self.CutStrings
        max_len = max(len(key) for key, value in cuts.iteritems() if str(value))
        for key, value in cuts.iteritems():
            if not key == 'all_cuts' and str(value):
                print '{key}:'.format(key=key.rjust(max_len)), value
        return
