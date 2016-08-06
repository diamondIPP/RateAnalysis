import os
import pickle
import json
from numpy import array, zeros, arange, delete
from Elementary import Elementary
from ROOT import TCut, gROOT, TH1F, kRed, TCutG
from collections import OrderedDict


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
            self.plot_settings = self.analysis.plots.plot_settings
            # saving stuff
            self.histos = {}

            # config
            self.DUTType = self.load_dut_type()
            # self.beaminterruptions_folder = self.ana_config_parser.get('CUT', 'beaminterruptions_folder')
            # self.exclude_before_jump = self.ana_config_parser.getint('CUT', 'excludeBeforeJump')
            # self.exclude_after_jump = self.ana_config_parser.getint('CUT', 'excludeAfterJump')
            self.CutConfig = {}

            # define cut strings
            # self.EasyCutStrings = self.init_easy_cutstrings()
            # self.CutStrings = self.define_cutstrings()

            # self.region_cut = TCut('region_cut', '')
            # self.JumpCut = TCut('JumpCut', '')

            # beam interrupts
            # self.jumps = None
            # self.jump_ranges = None

            self.load_config()
            # generate cut strings
            self.generate_fid_cuts()
            self.generate_masks()
            self.generate_chi2_cuts()
            self.generate_tracks_cut()
            self.generate_slope_cuts()
            self.generate_ini_fin_cuts()
            self.generate_beam_interruption_cut()
            self.add_cuts()

            # self.generate_cut_string()  # DA TODO
            # self.all_cut = self.generate_all_cut()  # DA TODO

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
        self.CutConfig['ExcludeFirst'] = self.ana_config_parser.getint('CUT', 'excludefirst') if self.ana_config_parser.\
            has_option('CUT', 'excludefirst') else 0
        self.CutConfig['ExcludeBeforeJump'] = self.ana_config_parser.getint('CUT', 'excludeBeforeJump') if self.ana_config_parser.\
            has_option('CUT', 'excludeBeforeJump') else 0
        self.CutConfig['ExcludeAfterJump'] = self.ana_config_parser.getint('CUT', 'excludeAfterJump') if self.ana_config_parser.\
            has_option('CUT', 'excludeAfterJump') else 0
        # self.CutConfig['EventRange'] = self.load_event_range(json.loads(self.ana_config_parser.get('CUT', 'EventRange')))
        self.CutConfig['chi2X'] = self.ana_config_parser.getint('CUT', 'chi2X') if self.ana_config_parser.\
            has_option('CUT', 'chi2X') else ''
        self.CutConfig['chi2Y'] = self.ana_config_parser.getint('CUT', 'chi2Y') if self.ana_config_parser.\
            has_option('CUT', 'chi2Y') else ''
        self.CutConfig['track_angle'] = self.ana_config_parser.getfloat('CUT', 'track_angle') if self.ana_config_parser.\
            has_option('CUT', 'track_angle') else ''
        self.CutConfig['MaskRowsDUT1'] = self.ana_config_parser.get('CUT', 'MaskRowsDUT1') if self.ana_config_parser.\
            has_option('CUT', 'MaskRowsDUT1') else ''
        self.CutConfig['MaskRowsDUT2'] = self.ana_config_parser.get('CUT', 'MaskRowsDUT2') if self.ana_config_parser.\
            has_option('CUT', 'MaskRowsDUT2') else ''
        self.CutConfig['MaskRowsDUT3'] = self.ana_config_parser.get('CUT', 'MaskRowsDUT3') if self.ana_config_parser.\
            has_option('CUT', 'MaskRowsDUT3') else ''
        self.CutConfig['MaskColsDUT1'] = self.ana_config_parser.get('CUT', 'MaskColsDUT1') if self.ana_config_parser.\
            has_option('CUT', 'MaskColsDUT1') else ''
        self.CutConfig['MaskColsDUT2'] = self.ana_config_parser.get('CUT', 'MaskColsDUT2') if self.ana_config_parser.\
            has_option('CUT', 'MaskColsDUT2') else ''
        self.CutConfig['MaskColsDUT3'] = self.ana_config_parser.get('CUT', 'MaskColsDUT3') if self.ana_config_parser.\
            has_option('CUT', 'MaskColsDUT3') else ''
        self.CutConfig['MaskPixelsDUT1'] = self.ana_config_parser.get('CUT', 'MaskPixelsDUT1') if self.ana_config_parser.\
            has_option('CUT', 'MaskPixelsDUT1') else ''
        self.CutConfig['MaskPixelsDUT2'] = self.ana_config_parser.get('CUT', 'MaskPixelsDUT2') if self.ana_config_parser.\
            has_option('CUT', 'MaskPixelsDUT2') else ''
        self.CutConfig['MaskPixelsDUT3'] = self.ana_config_parser.get('CUT', 'MaskPixelsDUT3') if self.ana_config_parser.\
            has_option('CUT', 'MaskPixelsDUT3') else ''
        self.CutConfig['FidRegionDUT1'] = self.ana_config_parser.get('CUT', 'FidRegionDUT1') if self.ana_config_parser.\
            has_option('CUT', 'FidRegionDUT1') else ''
        self.CutConfig['FidRegionDUT2'] = self.ana_config_parser.get('CUT', 'FidRegionDUT2') if self.ana_config_parser.\
            has_option('CUT', 'FidRegionDUT2') else ''
        self.CutConfig['FidRegionDUT3'] = self.ana_config_parser.get('CUT', 'FidRegionDUT3') if self.ana_config_parser.\
            has_option('CUT', 'FidRegionDUT3') else ''

    def add_cuts(self):
        self.cuts_hitmap_roc = {}
        self.cuts_pixelated_roc = {}
        for iROC in xrange(4, 7):
            self.cuts_hitmap_roc[iROC] = self.mask_hitmap_roc[iROC] + self.chi2x_cut + self.chi2y_cut \
                                         + self.cut_tracks + self.angle_x_cut + self.angle_y_cut + self.ini_fin_cut \
                                         + self.beam_interr_cut
            self.cuts_pixelated_roc[iROC] = self.mask_pixelated_roc[iROC] + self.chi2x_cut + self.chi2y_cut \
                                            + self.cut_tracks + self.angle_x_cut + self.angle_y_cut + self.ini_fin_cut \
                                            + self.beam_interr_cut

    def generate_ini_fin_cuts(self):
        nentries = self.analysis.tree.GetEntries()
        self.analysis.tree.GetEntry(0)
        first_t = self.analysis.tree.time
        self.analysis.tree.GetEntry(nentries-1)
        last_t = self.analysis.tree.time
        self.ini_fin_cut = TCut('initial_final_cuts', 'time>{ini}&&time<{fin}'.format(ini=first_t + abs(self.CutConfig['ExcludeFirst'])*1000, fin=last_t - abs(self.CutConfig['ExcludeFirst'])*1000))

    def generate_beam_interruption_cut(self):
        nentries = self.analysis.tree.GetEntries()
        self.analysis.tree.GetEntry(0)
        first_t = self.analysis.tree.time
        self.analysis.tree.GetEntry(nentries-1)
        last_t = self.analysis.tree.time
        bins = int((last_t-first_t)/float(5000))
        gROOT.SetBatch(True)
        h = TH1F('h', 'h', bins+1, first_t-(last_t-first_t)/float(2*bins), last_t+(last_t-first_t)/float(2*bins))
        mean = h.Integral()/float(h.GetNbinsX())
        self.beam_interr_cut = TCut('beam_interruptions_cut', '')
        for t in xrange(1, h.GetNbinsX()+1):
            if h.GetBinContent(t) < mean*0.9:
                self.beam_interr_cut=self.beam_interr_cut+TCut('bi{i}'.format(i=t), 'time<{low}&&time>{high}'.format(low=h.GetBinLowEdge(t)-abs(self.CutConfig['ExcludeBeforeJump'])*1000, high=h.GetBinLowEdge(t)+h.GetBinWidth(t)+abs(self.CutConfig['ExcludeAfterJump'])*1000))
        gROOT.SetBatch(False)

    def generate_tracks_cut(self):
        self.cut_tracks = TCut('cut_tracks', 'n_tracks')

    def generate_chi2_cuts(self):
        self.generate_chi2('x')
        self.generate_chi2('y')

    def generate_chi2(self, mode='x'):
        gROOT.SetBatch(1)
        h = TH1F('h', 'h', 2001, -0.05, 200.05)
        nq = 100
        chi2s = zeros(nq)
        xq = array([(i + 1) / float(nq) for i in xrange(nq)])
        self.analysis.tree.Draw('chi2_{mod}>>h'.format(mod=mode), '', 'goff')
        h.GetQuantiles(nq, chi2s, xq)
        gROOT.SetBatch(0)
        quantile = self.CutConfig['chi2{mod}'.format(mod=mode.title())]
        assert type(quantile) is int and 0 < quantile <= 100, 'chi2 quantile has to be an integer between (0, 100]'
        string = 'chi2_{mod}<{val}&&chi2_{mod}>=0'.format(val=chi2s[quantile], mod=mode)
        if quantile > 0:
            exec('self.chi2{mode}_cut = TCut("chi2_cut_{mode}", "{cut}")'.format(mode=mode, cut=string))
        else:
            exec('self.chi2{mode}_cut = TCut("chi2_cut_{mode}", {cut})'.format(mode=mode, cut=''))

    def generate_fid_cuts(self):
        self.fid_cut_hitmap_roc = {}
        self.fid_cut_pixelated_roc = {}
        # self.fid_cut_local_roc = {}
        # self.fid_cut_telescope_roc = {}
        self.generate_fid_cuts_DUT(1)
        self.generate_fid_cuts_DUT(2)
        self.generate_fid_cuts_DUT(3)

    def generate_fid_cuts_DUT(self, dut):
        fidcutstring = self.CutConfig['FidRegionDUT{d}'.format(d=dut)]
        roc = 4 if dut is 1 else 5 if dut is 2 else 6
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

    def generate_masks(self):
        self.mask_hitmap_roc = {}
        self.mask_pixelated_roc = {}
        self.generate_col_masks()
        self.generate_row_masks()
        self.generate_pixel_masks()

    def generate_col_masks(self):
        self.col_mask_hitmap_roc = {}
        self.col_mask_pixelate_roc = {}
        self.generate_col_masks_DUT(1)
        self.generate_col_masks_DUT(2)
        self.generate_col_masks_DUT(3)

    def generate_col_masks_DUT(self, dut):
        maskcolstring = self.CutConfig['MaskColsDUT{d}'.format(d=dut)]
        roc = 4 if dut is 1 else 5 if dut is 2 else 6
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
        self.generate_row_masks_DUT(1)
        self.generate_row_masks_DUT(2)
        self.generate_row_masks_DUT(3)

    def generate_row_masks_DUT(self, dut):
        maskrowstring = self.CutConfig['MaskRowsDUT{d}'.format(d=dut)]
        roc = 4 if dut is 1 else 5 if dut is 2 else 6
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
        self.generate_pixel_masks_DUT(1)
        self.generate_pixel_masks_DUT(2)
        self.generate_pixel_masks_DUT(3)

    def generate_pixel_masks_DUT(self, dut):
        maskpixelstring = self.CutConfig['MaskPixelsDUT{d}'.format(d=dut)]
        roc = 4 if dut is 1 else 5 if dut is 2 else 6
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

    def generate_slope_cuts(self):
        self.generate_slope('x')
        self.generate_slope('y')

    def generate_slope(self, mode='x'):
        # picklepath = 'Configuration/Individual_Configs/Slope/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.analysis.lowest_rate_run)
        angle = self.CutConfig['track_angle']

        # def func():
        #     print 'generating slope cut for run {run}...'.format(run=self.analysis.run_number)
            # fit the slope to get the mean
            # gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        h = TH1F('h', 'h', 61, -3.05, 3.05)
        self.analysis.tree.Draw('slope_{x}>>h'.format(x=mode), '', 'goff')
        fit_result = h.Fit('gaus', 'qs')
        x_mean = fit_result.Parameter(1)
        x_sigm = fit_result.Parameter(2)
        xmin = x_mean - x_sigm
        xmax = x_mean + x_sigm
        h1 = TH1F('h1', 'h1', 51, xmin - 2*x_sigm/100, xmax + 2*x_sigm/100)
        self.analysis.tree.Draw('slope_{x}>>h1'.format(x=mode), '', 'goff')
        fit_result1 = h1.Fit('gaus', 'qs')
        x_mean = fit_result1.Parameter(1)
        slopes = [x_mean - angle, x_mean + angle]
            # c = gROOT.FindObject('c1')
            # c.Close()
        gROOT.SetBatch(0)
            # gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
            # return slopes

        # slope = self.do_pickle(picklepath, func)
        # create the cut string
        string = 'slope_{x}>={minx}&&slope_{x}<={maxx}'.format(x=mode, minx=slopes[0], maxx=slopes[1])
        if angle > 0:
            exec('self.angle_{x}_cut = TCut("angle_cut_{x}", "{cut}")'.format(x=mode, cut=string))
        else:
            exec('self.angle_{x}_cut = TCut("angle_cut_{x}", {cut})'.format(x=mode, cut=''))

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
