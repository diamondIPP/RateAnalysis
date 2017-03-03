import sys
from numpy import array
from ROOT import TCut, gROOT, TH1F, kRed, TCutG, TH2D, TH1D, THStack
from Cut import Cut
from json import loads
from collections import OrderedDict


class CutPix(Cut):
    """
    A cut contains all cut settings which corresponds to a single diamond in a single run. Thus, an Analysis object holds two Cut instances, one for each diamond. The default configuration
    is loaded from the Analysis config file, whereas the individual cut settings are loaded from a JSON file located at Configuration/Individual_Configs. The JSON files are generated
    by the Analysis method SetIndividualCuts().
    """
    def __init__(self, analysis, dut=1):
        Cut.__init__(self, analysis, skip=True)
        self.__dict__.update(analysis.Cut.__dict__)

        self.Dut = dut + 3

        self.Settings = self.analysis.Settings
        self.plots = self.analysis.plots

        self.load_pixel_config()

        self.generate_pixel_cutstrings()
        self.all_cut = self.generate_all_cut()
        self.CutStrings['all_cuts'] = self.all_cut
        self.HitMapCut = self.generate_hitmap_cutstrings()

        self.ConsecutiveCuts = self.generate_consecutive_cuts()
        self.ConsecutiveHitMapCuts = self.generate_consecutive_cuts(cluster=False)

    def load_run_config(self):
        return self.load_run_configs(self.RunNumber)

    def generate_pixel_cutstrings(self):
        """ Generates the cut strings to apply in the analysis for each of the cuts. """
        self.CutStrings['hit'] += self.generate_hit()
        self.CutStrings['masks'] += self.generate_masks()
        self.CutStrings['fiducial'] += self.generate_fiducial()
        self.CutStrings['rhit'] += self.generate_rhit()
        self.CutStrings['trigger_phase'] += self.generate_trigger_phase()
        self.CutStrings['alignment'] += self.generate_alignment()

    def generate_hitmap_cutstrings(self):
        self.set_hitmap_cuts()
        cut = self.generate_all_cut()
        self.set_hitmap_cuts(False)
        return cut

    def reload_cuts(self):
        self.ana_config_parser = self.load_ana_config()
        self.load_pixel_config()
        self.generate_pixel_cutstrings()

    def set_hitmap_cuts(self, on=True):
        self.set_cut('masks', self.generate_masks(cluster=not on))

    def generate_consecutive_cuts(self, cluster=True):
        self.set_hitmap_cuts(not cluster)
        cuts = OrderedDict({'raw': TCut('0', '')})
        n = 1
        for key, value in self.CutStrings.iteritems():
            if str(value) and key != 'all_cuts' and not key.startswith('old'):
                new_cut = cuts.values()[n - 1] + value
                key = 'beam stop' if key.startswith('beam') else key
                key = key.replace('track_', '')
                cuts['+ {n}'.format(n=key)] = TCut('{n}'.format(n=n), new_cut.GetTitle())
                n += 1
        self.set_hitmap_cuts(False)
        return cuts

    def do_cuts_distributions(self):
        """ Does the cuts distribution for all other cuts """
        # todo implement this...
        self.print_banner('Doing cuts distributions...')
        self.print_banner('Finished with distribution cuts')

    def load_pixel_config(self):
        """ Loads the pixel configuration parameters from the config file. """
        self.CutConfig['rhit'] = self.get_config('r_hit')
        self.CutConfig['trigPhase'] = loads(self.get_config('trigger_phase'))[str(self.Dut)]
        self.CutConfig['MaskRows'] = self.load_mask('MaskRows')
        self.CutConfig['MaskCols'] = self.load_mask('MaskCols')
        self.CutConfig['MaskPixels'] = self.load_mask('MaskPixels')
        self.CutConfig['FidRegion'] = self.load_fiducial()
        self.CutConfig['FidRegionLocal'] = self.load_fiducial('FidPixLocal')

    def get_config(self, option):
        return self.ana_config_parser.get('CUT', option) if self.ana_config_parser.has_option('CUT', option) else None

    def load_mask(self, name):
        string = self.get_config('{n}ROC{d}'.format(n=name, d=self.Dut))
        if string == '[]' or string == 'None':
            return []
        lst = []
        if string is not None:
            string = string.replace('[', '').replace(']', '')
            tuples = string.split(';')
            for tup in tuples:
                lst.append([int(i) for i in tup.split('{s}'.format(s=',' if 'Pix' in name else ':'))])
        return lst if lst else None

    def load_fiducial(self, name='FidPix'):
        if self.ana_config_parser.has_option('CUT', name):
            dic = loads(self.ana_config_parser.get('CUT', name))
            dut = str(self.Dut)
            if dut in dic:
                return dic[dut]

    def generate_special_cut(self, excluded=None, included=None, name='special_cut', cluster=True):
        cut = TCut(name, '')
        self.NCuts = 0
        excluded = [excluded] if type(excluded) is not list else excluded
        self.set_hitmap_cuts()
        for key, value in self.CutStrings.iteritems():
            if excluded and key in excluded:
                continue
            if included and key not in included:
                continue
            if key.startswith('old') or key.startswith('all_cut'):
                continue
            if value.GetTitle() == '':
                continue
            cut += value
            self.NCuts += 1
        self.log_info('generated {name} cut with {num} cuts'.format(name=name, num=self.NCuts))
        self.set_hitmap_cuts(False)
        return cut

    @staticmethod
    def generate_hit():  # TODO implement cut! must change tree structure from tracking telescope
        """ Needs to be implemented. Have to change trackingTelescope for this """
        cut_string = ''
        return cut_string

    def generate_rhit(self):
        value = self.CutConfig['rhit']
        string = '(10000*sqrt((residuals_x[{n}])**2+(residuals_y[{n}])**2))<{val}'.format(n=self.Dut, val=value)
        return string

    def generate_trigger_phase(self):
        mi, ma = self.CutConfig['trigPhase']
        return 'trigger_phase[1]>={min}&&trigger_phase[1]<={max}'.format(min=mi, max=ma)

    def generate_fiducial(self, name='fid'):
        xy = self.CutConfig['FidRegion']
        if xy is not None:
            d = 0
            x = array([xy[0] - d, xy[0] - d, xy[1] + d, xy[1] + d, xy[0] - d], 'd')
            y = array([xy[2] - d, xy[3] + d, xy[3] + d, xy[2] - d, xy[2] - d], 'd')
            cut = TCutG(name, 5, x, y)
            cut.SetVarX('diam{n}_track_x'.format(n=self.Dut - 3))
            cut.SetVarY('diam{n}_track_y'.format(n=self.Dut - 3))
            self.ROOTObjects.append(cut)
            cut.SetLineColor(kRed)
            cut.SetLineWidth(3)
        return name

    def generate_alignment(self):
        """ This adds the restrictions to the cut string such that misalignments are excluded each time the cut is applied. """
        interruptions = self.find_misaligments()
        cut_string = TCut('')
        for beg, end in zip(interruptions['begin'], interruptions['end']):
            cut_string += TCut('event_number<{low}||event_number>{high}'.format(low=beg, high=end))
        self.EasyCutStrings['alignment'] = 'AlignOn'
        return cut_string.GetTitle()

    def find_misaligments(self):
        picklepath = self.make_pickle_path('Cuts', 'Alignment', run=self.RunNumber, suf=self.analysis.Dut)

        def func():
            n = 1000
            start = self.log_info('Generating aligment cut ... ', next_line=False)
            g = self.analysis.draw_alignment(show=False, binning=n, vs_time=False)
            h = TH1F('h_p', 'Pull', g.GetN() / 3, 0, 1)
            for i in xrange(g.GetN()):
                h.Fill(g.GetY()[i])
            max_x = h.GetBinCenter(h.GetMaximumBin())
            fit = h.Fit('gaus', 'qs', '', max_x - .2, max_x + .2)
            mean_, sigma = fit.Parameter(1), fit.Parameter(2)
            mis_events = {'begin': [], 'end': []}
            cut = mean_ - 5 * sigma
            min_events = OrderedDict(sorted({i: int(g.GetX()[i]) for i in xrange(g.GetN()) if g.GetY()[i] < cut}.iteritems()))
            if not min_events:
                return mis_events
            mis_events['begin'].append(min_events.values()[0] - n)
            n *= 2
            for i in xrange(len(min_events) - 1):
                if not min_events.keys()[i] == min_events.keys()[i + 1] - 1:
                    mis_events['end'].append(min_events.values()[i] + n)
                    mis_events['begin'].append(min_events.values()[i + 1] - n)
            mis_events['end'].append(min_events.values()[-1] + n)
            self.draw_histo(h)
            self.add_info(start)
            return mis_events

        events = self.do_pickle(picklepath, func)
        return events

    def generate_masks(self, cluster=True):
        # t = self.log_info('Generating mask cuts ...', False)
        cut_string = TCut('')
        cut_string += self.generate_line_mask('col', cluster)
        cut_string += self.generate_line_mask('row', cluster)
        cut_string += self.generate_pixel_mask(cluster)
        # self.add_info('Done', t)
        return cut_string.GetTitle()

    def generate_line_mask(self, line, cluster=True):
        cut_string = ''
        cut_var = 'cluster_{l}'.format(n=self.Dut, l=line) if cluster else line
        for tup in self.CutConfig['Mask{l}s'.format(l=line.title())]:
            cut_string += '||' if cut_string else ''
            if len(tup) == 2:
                cut_string += '{v}>={i}&&{v}<={f}'.format(i=tup[0], f=tup[1], v=cut_var)
            elif len(tup) == 1:
                cut_string += '{v}=={i}'.format(i=tup[0], v=cut_var)
        if not cut_string:
            return ''
        plane_var = 'plane' if not cluster else 'cluster_plane'
        cut_string = TCut(cut_string) + TCut('{p}=={r}'.format(r=self.Dut, p=plane_var))
        return '!({c})'.format(c=cut_string.GetTitle())

    def generate_pixel_mask(self, cluster=True):
        cut_string = ''
        cut_var1 = 'cluster_col'.format(n=self.Dut) if cluster else 'col'
        cut_var2 = 'cluster_row'.format(n=self.Dut) if cluster else 'row'
        for tup in self.CutConfig['MaskPixels']:
            cut_string += '||' if cut_string else ''
            cut_string += '{v1}=={x}&&{v2}=={y}'.format(x=tup[0], y=tup[1], v1=cut_var1, v2=cut_var2)
        if not cut_string:
            return ''
        plane_var = 'plane' if not cluster else 'cluster_plane'
        cut_string = TCut(cut_string) + TCut('{p}=={r}'.format(r=self.Dut, p=plane_var))
        return '!({c})'.format(c=cut_string.GetTitle())

    @staticmethod
    def add_adc_cut(adc):
        if adc is not None:
            return 'adc>0' if adc else 'adc==0'
        else:
            return ''

    # region ANALYSIS

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
        if self.verbose: print 'Res_Y Vs Res_X...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resy_resx[iroc] = TH2D('h_resy_resx_roc{r}'.format(r=iroc), 'h_resy_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 21, -1050, 1050)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*residual_ROC{r}_Local_X>>h_resy_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1],
                                    'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*residual_ROC{r}_Local_X>>h_resy_resx_roc{r}'.format(r=iroc), '',
                                                                                               'goff')
            self.plots.set_2D_options(self.h_resy_resx[iroc], 'Res_X(um)', 'Res_y(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resy_resx[iroc], 'h_resy_resx_roc{r}'.format(r=iroc), 'Res_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir + '/cuts',
                                             False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Rhit Vs Res_X...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_resx[iroc] = TH2D('h_rhit_resx_roc{r}'.format(r=iroc), 'h_rhit_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 101, -0.5, 100.5)
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_X>>h_rhit_resx_roc{r}'.format(r=iroc),
                                    self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw(
                '(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_X>>h_rhit_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_rhit_resx[iroc], 'Res_X(um)', 'R_Hit(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_rhit_resx[iroc], 'h_rhit_resx_roc{r}'.format(r=iroc), 'R_Hit Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir + '/cuts',
                                             False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Rhit Vs Res_Y...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_rhit_resy[iroc] = TH2D('h_rhit_resy_roc{r}'.format(r=iroc), 'h_rhit_resy_roc{r}'.format(r=iroc), 21, -1050, 1050, 101, -0.5, 100.5)
            self.analysis.tree.Draw('(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_Y>>h_rhit_resy_roc{r}'.format(r=iroc),
                                    self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw(
                '(10000*sqrt((residual_ROC{r}_Local_X)**2+(residual_ROC{r}_Local_Y)**2)):10000*residual_ROC{r}_Local_Y>>h_rhit_resy_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_rhit_resy[iroc], 'Res_Y(um)', 'R_Hit(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_rhit_resy[iroc], 'h_rhit_resy_roc{r}'.format(r=iroc), 'R_Hit Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir + '/cuts',
                                             False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2 Vs Res_X...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2_resx[iroc] = TH2D('h_chi2_resx_roc{r}'.format(r=iroc), 'h_chi2_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_X>>h_chi2_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1],
                                    'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_X>>h_chi2_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2_resx[iroc], 'Res_X(um)', 'Chi2', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2_resx[iroc], 'h_chi2_resx_roc{r}'.format(r=iroc), 'Chi2 Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir + '/cuts',
                                             False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2 Vs Res_Y...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2_resy[iroc] = TH2D('h_chi2_resy_roc{r}'.format(r=iroc), 'h_chi2_resy_roc{r}'.format(r=iroc), 21, -1050, 1050, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_Y>>h_chi2_resy_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1],
                                    'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_tracks:10000*residual_ROC{r}_Local_Y>>h_chi2_resy_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2_resy[iroc], 'Res_Y(um)', 'Chi2', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2_resy[iroc], 'h_chi2_resy_roc{r}'.format(r=iroc), 'Chi2 Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0, self.plots.save_dir + '/cuts',
                                             False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2_X Vs Res_X...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2x_resx[iroc] = TH2D('h_chi2x_resx_roc{r}'.format(r=iroc), 'h_chi2x_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_x:10000*residual_ROC{r}_Local_X>>h_chi2x_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1],
                                    'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_x:10000*residual_ROC{r}_Local_X>>h_chi2x_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2x_resx[iroc], 'Res_X(um)', 'Chi2_X', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2x_resx[iroc], 'h_chi2x_resx_roc{r}'.format(r=iroc), 'Chi2_X Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0,
                                             self.plots.save_dir + '/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Chi2_Y Vs Res_X...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_chi2y_resx[iroc] = TH2D('h_chi2y_resx_roc{r}'.format(r=iroc), 'h_chi2y_resx_roc{r}'.format(r=iroc), 21, -1575, 1575, 51, -0.1, 10.1)
            self.analysis.tree.Draw('chi2_y:10000*residual_ROC{r}_Local_X>>h_chi2y_resx_roc{r}'.format(r=iroc), self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1],
                                    'goff') if self.num_cuts != 0 else self.analysis.tree.Draw('chi2_y:10000*residual_ROC{r}_Local_X>>h_chi2y_resx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_chi2y_resx[iroc], 'Res_X(um)', 'Chi2_Y', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_chi2y_resx[iroc], 'h_chi2y_resx_roc{r}'.format(r=iroc), 'Chi2_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0,
                                             self.plots.save_dir + '/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Res_X Vs Hit_Y...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resx_hitposy[iroc] = TH2D('h_resx_hitposy_roc{r}'.format(r=iroc), 'h_resx_hitposy_roc{r}'.format(r=iroc), 161, -4025, 4025, 21, -1575, 1575)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_X:10000*(residual_ROC{r}_Local_Y+cluster_pos_ROC{r}_Local_Y)>>h_resx_hitposy_roc{r}'.format(r=iroc),
                                    self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw(
                '10000*residual_ROC{r}_Local_X:10000*(residual_ROC{r}_Local_Y+cluster_pos_ROC{r}_Local_Y)>>h_resx_hitposy_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_resx_hitposy[iroc], 'Hit_Y(um)', 'Res_X(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resx_hitposy[iroc], 'h_resx_hitposy_roc{r}'.format(r=iroc), 'Hit_Y Vs. Res_X roc{r}'.format(r=iroc), None, 'colz', 0,
                                             self.plots.save_dir + '/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        if self.verbose: print 'Res_Y Vs Hit_X...',; sys.stdout.flush()
        for iroc in self.duts_list:
            gROOT.SetBatch(True)
            self.h_resy_hitposx[iroc] = TH2D('h_resy_hitposx_roc{r}'.format(r=iroc), 'h_resy_hitposx_roc{r}'.format(r=iroc), 105, -3937.5, 3937.5, 21, -1050, 1050)
            self.analysis.tree.Draw('10000*residual_ROC{r}_Local_Y:10000*(residual_ROC{r}_Local_X+cluster_pos_ROC{r}_Local_X)>>h_resy_hitposx_roc{r}'.format(r=iroc),
                                    self.cuts_pixelated_roc_incr[iroc][self.num_cuts - 1], 'goff') if self.num_cuts != 0 else self.analysis.tree.Draw(
                '10000*residual_ROC{r}_Local_Y:10000*(residual_ROC{r}_Local_X+cluster_pos_ROC{r}_Local_X)>>h_resy_hitposx_roc{r}'.format(r=iroc), '', 'goff')
            self.plots.set_2D_options(self.h_resy_hitposx[iroc], 'Hit_X(um)', 'Res_Y(um)', '# entries', 0, -1)
            gROOT.SetBatch(False)
            self.plots.save_individual_plots(self.h_resy_hitposx[iroc], 'h_resy_hitposx_roc{r}'.format(r=iroc), 'Hit_X Vs. Res_Y roc{r}'.format(r=iroc), None, 'colz', 0,
                                             self.plots.save_dir + '/cuts', False, 0, doLogZ=True)
        if self.verbose: print 'Done'

        self.print_banner('Finished with resolution plots')

    def do_occupancy_analysis(self):
        """ Does the occupancy analysis for each roc applying each of the cummulative cuts """
        self.print_banner('Starting occupancy cut analysis...')
        self.start_pbar(len(self.ConsecutiveHitMapCuts))
        for i, cut in enumerate(self.ConsecutiveHitMapCuts):
            self.NCuts = i
            self.analysis.draw_occupancy(cut, show=False, fid=True, prnt=False)
            self.ProgressBar.update(i + 1)

    def do_pulse_height_analysis(self, show=True):
        self.print_banner('Starting pulse height cut analysis...')
        self.start_pbar(len(self.ConsecutiveHitMapCuts))
        stack = THStack('s_ph', 'Pulse Height Distribution with Consecutive Cuts')
        legend = self.make_legend(.71, .88, nentries=len(self.ConsecutiveCuts))
        for i, (name, cut) in enumerate(self.ConsecutiveCuts.iteritems()):
            self.NCuts = i
            h = self.analysis.draw_pulse_height_disto(cut, show=False, prnt=False)
            color = self.get_color()
            self.format_histo(h, color=color, fill_color=color)
            stack.Add(h)
            legend.AddEntry(h, name, 'f')
            self.ProgressBar.update(i + 1)
        self.reset_colors()
        self.format_histo(stack, draw_first=True, x_tit='Pulse Height [e]', y_tit='Number of Entries', y_off=1.5, stats=0)
        self.save_histo(stack, 'ConsecutivePulseHeights', show, draw_opt='nostack', l=legend, lm=.14)

    def do_pulse_heights_analysis(self, normalize_ph_plots=True):
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
        phbins = {self.roc_diam1: self.Settings['ph1DbinsD4'], self.roc_diam2: self.Settings['ph1DbinsD5'],
                  self.roc_si: self.Settings['ph1DbinsSi']}
        phmin = {self.roc_diam1: self.Settings['ph1DminD4'], self.roc_diam2: self.Settings['ph1DminD5'],
                 self.roc_si: self.Settings['ph1DminSi']}
        phmax = {self.roc_diam1: self.Settings['ph1DmaxD4'], self.roc_diam2: self.Settings['ph1DmaxD5'],
                 self.roc_si: self.Settings['ph1DmaxSi']}
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
                                                                              'ph1_map_roc{r}_{c}'.format(r=iroc, c=cut),
                                                                              'ph1_map_roc{r}_{c}'.format(r=iroc, c=cut),
                                                                              'x(um)', 'y(um)', 'ph 1 pix cluster(e)',
                                                                              'auto', -1)
                self.h_ph2_map_cuts[iroc][cut] = self.plots.create_2D_profile('spatial',
                                                                              'ph2_map_roc{r}_{c}'.format(r=iroc, c=cut),
                                                                              'ph2_map_roc{r}_{c}'.format(r=iroc, c=cut),
                                                                              'x(um)', 'y(um)', 'ph 2 pix cluster(e)',
                                                                              'auto', -1)
                self.h_ph1_evt_cuts[iroc][cut] = TH2D('ph1_evt_roc{r}_{c}'.format(r=iroc, c=cut),
                                                      'ph1_evt_roc{r}_{c}'.format(r=iroc, c=cut), self.Settings[
                                                          'event_bins'] + 1, self.Settings['event_min'] - (
                                                          self.Settings['event_max'] - self.Settings[
                                                              'event_min']) / (2 * float(self.Settings['event_bins'])),
                                                      self.Settings['event_max'] + (
                                                          self.Settings['event_max'] - self.Settings[
                                                              'event_min']) / (2 * float(self.Settings['event_bins'])),
                                                      phbins[iroc] + 1, phmin[iroc] - phdelta[iroc] / (2 * float(phbins[iroc])),
                                                      phmax[iroc] + phdelta[iroc] / float(2 * phbins[iroc]))
                self.h_ph1_cuts[iroc][cut] = TH1D('ph1_roc{r}_{c}'.format(r=iroc, c=cut),
                                                  'ph1_roc{r}_{c}'.format(r=iroc, c=cut), phbins[iroc] + 1,
                                                  phmin[iroc] - phdelta[iroc] / (2 * float(phbins[iroc])),
                                                  phmax[iroc] + phdelta[iroc] / float(2 * phbins[iroc]))
                self.h_ph2_evt_cuts[iroc][cut] = TH2D('ph2_evt_roc{r}_{c}'.format(r=iroc, c=cut),
                                                      'ph2_evt_roc{r}_{c}'.format(r=iroc, c=cut), self.Settings[
                                                          'event_bins'] + 1, self.Settings['event_min'] - (
                                                          self.Settings['event_max'] - self.Settings[
                                                              'event_min']) / (2 * float(self.Settings['event_bins'])),
                                                      self.Settings['event_max'] + (
                                                          self.Settings['event_max'] - self.Settings[
                                                              'event_min']) / (2 * float(self.Settings['event_bins'])),
                                                      phbins[iroc] + 1, phmin[iroc] - phdelta[iroc] / (2 * float(phbins[iroc])),
                                                      phmax[iroc] + phdelta[iroc] / float(2 * phbins[iroc]))
                self.h_ph2_cuts[iroc][cut] = TH1D('ph2_roc{r}_{c}'.format(r=iroc, c=cut),
                                                  'ph2_roc{r}_{c}'.format(r=iroc, c=cut), phbins[iroc] + 1,
                                                  phmin[iroc] - phdelta[iroc] / (2 * float(phbins[iroc])),
                                                  phmax[iroc] + phdelta[iroc] / float(2 * phbins[iroc]))
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
                    if self.verbose: print 'Saving for ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut),; sys.stdout.flush()
                    self.plots.set_2D_options(self.h_ph1_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries')
                    self.plots.set_1D_options('ph', self.h_ph1_cuts[iroc][cut], 'ph 1 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph2_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries')
                    self.plots.set_1D_options('ph', self.h_ph2_cuts[iroc][cut], 'ph 2 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph1_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 1 pix cluster(e)')
                    self.plots.set_2D_options(self.h_ph2_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 2 pix cluster(e)')

                    self.plots.save_individual_plots(self.h_ph1_evt_cuts[iroc][cut], 'ph1_evt_roc{r}_{c}'.format(r=iroc, c=cut), 'ph1_evt_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph2_evt_cuts[iroc][cut], 'ph2_evt_roc{r}_{c}'.format(r=iroc, c=cut), 'ph2_evt_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph1_cuts[iroc][cut], 'ph1_roc{r}_{c}'.format(r=iroc, c=cut), 'ph1_roc{r}_{c}'.format(r=iroc, c=cut), None, '', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph2_cuts[iroc][cut], 'ph2_roc{r}_{c}'.format(r=iroc, c=cut), 'ph2_roc{r}_{c}'.format(r=iroc, c=cut), None, '', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph1_map_cuts[iroc][cut], 'ph1_map_roc{r}_{c}'.format(r=iroc, c=cut), 'ph1_map_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph2_map_cuts[iroc][cut], 'ph2_map_roc{r}_{c}'.format(r=iroc, c=cut), 'ph2_map_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
                    if self.verbose: print 'Done'

        if normalize_ph_plots:
            min_ph1_map[iroc] = min(min_ph1_map[iroc], 0)
            min_ph2_map[iroc] = min(min_ph2_map[iroc], 0)
            minz_ph1 = min(minz_ph1, 0)
            minz_ph2 = min(minz_ph2, 0)
            for iroc in self.duts_list:
                for cut in self.cut_names:
                    if self.verbose: print 'Saving for ROC {r} with cummulative cut {c}...'.format(r=iroc, c=cut),; sys.stdout.flush()
                    self.plots.set_2D_options(self.h_ph1_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries', min_val=minz_ph1, max_val=maxz_ph1)
                    self.plots.set_1D_options('ph', self.h_ph1_cuts[iroc][cut], 'ph 1 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph2_evt_cuts[iroc][cut], 'event', 'ph(e)', 'entries', min_val=minz_ph2, max_val=maxz_ph2)
                    self.plots.set_1D_options('ph', self.h_ph2_cuts[iroc][cut], 'ph 2 pix cl (e)', 'entries')
                    self.plots.set_2D_options(self.h_ph1_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 1 pix cluster(e)', min_val=min_ph1_map[iroc], max_val=max_ph1_map[iroc])
                    self.plots.set_2D_options(self.h_ph2_map_cuts[iroc][cut], 'x(um)', 'y(um)', 'ph 2 pix cluster(e)', min_val=min_ph2_map[iroc], max_val=max_ph2_map[iroc])

                    self.plots.save_individual_plots(self.h_ph1_evt_cuts[iroc][cut], 'ph1_evt_roc{r}_{c}'.format(r=iroc, c=cut), 'ph1_evt_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph2_evt_cuts[iroc][cut], 'ph2_evt_roc{r}_{c}'.format(r=iroc, c=cut), 'ph2_evt_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph1_cuts[iroc][cut], 'ph1_roc{r}_{c}'.format(r=iroc, c=cut), 'ph1_roc{r}_{c}'.format(r=iroc, c=cut), None, '', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph2_cuts[iroc][cut], 'ph2_roc{r}_{c}'.format(r=iroc, c=cut), 'ph2_roc{r}_{c}'.format(r=iroc, c=cut), None, '', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph1_map_cuts[iroc][cut], 'ph1_map_roc{r}_{c}'.format(r=iroc, c=cut), 'ph1_map_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
                    self.plots.save_individual_plots(self.h_ph2_map_cuts[iroc][cut], 'ph2_map_roc{r}_{c}'.format(r=iroc, c=cut), 'ph2_map_roc{r}_{c}'.format(r=iroc, c=cut), None, 'colz', 1,
                                                     self.plots.save_dir + '/cuts')
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
        # endregion
