from argparse import ArgumentParser
from copy import deepcopy
from time import sleep

from ROOT import TCanvas, TH2F, gROOT, TProfile, TH1F, TLegend, gStyle, kGreen, TArrow, kOrange, kViolet, kCyan, TText, TCut
from numpy import array, zeros

from Elementary import Elementary
from RunClass import Run
from Cut import Cut


class Analysis(Elementary):
    """ Class for the analysis of the non-channel specific stuff of a single run. """

    def __init__(self, run, diamonds=3, verbose=False, high_low_rate=None):
        """
        Parent class for all analyses, which contains all the basic stuff about the Telescope.
        :param run:             run object of type "Run" or integer run number
        :param diamonds:        an integer number defining the diamonds activated for analysis: 0x1=ch0 (diamond 1) 0x2=ch3 (diamond 2)
        :param verbose:         if True, verbose printing is activated
        :param high_low_rate:   list of highest and lowest rate runs for an analysis collection
        """
        Elementary.__init__(self, verbose=verbose)
        self.histos = []
        self.RootObjects = []

        # basics
        self.diamonds = diamonds
        self.run = self.init_run(run)
        self.run.analysis = self
        self.run_number = self.run.run_number
        self.RunInfo = deepcopy(self.run.RunInfo)
        self.lowest_rate_run = high_low_rate['min'] if high_low_rate is not None else self.run.run_number
        self.highest_rate_run = high_low_rate['max'] if high_low_rate is not None else self.run.run_number
        self.PickleDir = self.get_program_dir() + self.ana_config_parser.get('SAVE', 'pickle_dir')
        # self.saveMCData = self.ana_config_parser.getboolean("SAVE", "SaveMCData")
        self.ana_save_dir = '{run}'.format(run=self.run.run_number)
        self.set_root_titles()

        # DUT
        self.DUTType = self.run.DUTType
        
        # tree
        self.tree = self.run.tree

        # miscellaneous
        self.channel = self.channel if hasattr(self, 'channel') else None

        # general for pads and pixels
        self.Cut = Cut(self)
        self.StartEvent = self.Cut.CutConfig['EventRange'][0]
        self.EndEvent = self.Cut.CutConfig['EventRange'][1]

        # save histograms // canvases
        self.signal_canvas = None

        # alignment
        self.IsAligned = self.check_alignment(draw=False, save_plot=False)

    # ============================================================================================
    # region INIT

    def init_run(self, run):
        if not isinstance(run, Run):
            assert type(run) is int, 'run has to be either a Run instance or an integer run number'
            return Run(run, self.diamonds)
        else:
            assert run.run_number is not None, 'No run selected, choose run.SetRun(run_nr) before you pass the run object'
            return run

    def set_root_titles(self):
        if self.MainConfigParser.has_option('SAVE', 'activate_title'):
            gStyle.SetOptTitle(self.MainConfigParser.getboolean('SAVE', 'activate_title'))
    # endregion

    # ============================================================================================
    # region REGIONS AND PEAK INTEGRAL

    def __draw_single_wf(self, event=None, show=True):
        start = self.StartEvent if event is None else event
        if hasattr(self, 'draw_waveforms') and self.run.wf_exists(self.channel):
            h = self.draw_waveforms(n=1, show=show, start_event=start)[0]
        else:
            h = TH2F('regions', '', 1024, 0, 511, 1000, -200, 50)
            if self.run.wf_exists(0):
                self.tree.Draw('wf0:Iteration$/2>>regions', self.Cut.all_cut, 'goff', 1, start)
        h.GetXaxis().SetNdivisions(26)
        self.format_histo(h, markersize=0.3, x_tit='Time [ns]', y_tit='Signal [au]', stats=0)
        self.RootObjects.append(self.save_histo(h, 'Regions', show, self.ana_save_dir, lm=.075, rm=.045, x=2000, y=1000))
        return h

    def draw_regions(self, ped=True, event=None):
        h = self.__draw_single_wf(event=event, show=False)
        c = TCanvas('c1', 'Regions', 1000, 500)
        c.SetMargin(.075, .045, .2, .1)
        c.SetGrid()
        h.Draw()
        tit = 'Pedestal Regions' if ped else 'Signal Regions'
        h.SetTitle(tit)
        lines = []
        starts = []
        titles = []
        regions = self.run.pedestal_regions if ped else self.run.signal_regions
        gr = self.make_tgrapherrors('gr', '', color=2, marker_size=0, width=3)
        i = 0
        gStyle.SetEndErrorSize(4)
        sleep(.5)
        for reg, lst in regions.iteritems():
            if len(reg) < 3:
                offset = 40 if not lst[0] in starts else 20
                if lst[1] - lst[0] > 1:
                    gr.SetPoint(i, (lst[1] + lst[0]) / 4., c.GetUymax() - offset)
                    gr.SetPointError(i, (lst[1] - lst[0]) / 4., 0)
                    l = self.make_tlatex(gr.GetX()[i], gr.GetY()[i] + 3, '{sig}{reg}'.format(reg=reg, sig='p' if ped else 's'), color=2, size=.04)
                    gr.GetListOfFunctions().Add(l)
                    i += 1
                l1 = self.make_tgaxis(lst[0] / 2, c.GetUymin(), c.GetUymax() - offset, '', 2)
                l2 = self.make_tgaxis(lst[1] / 2, c.GetUymin(), c.GetUymax() - offset, '', 2) if lst[1] - lst[0] > 1 else 0
                if not lst[1] - lst[0] > 1:
                    l1.SetLineColor(4)
                    l1.SetLineWidth(2)
                    l1.SetTitleColor(4)
                    l1.SetY2(c.GetUymax() - 100)
                    tit = self.make_tlatex(lst[0] / 2, c.GetUymax() - 97, '{sig}{reg}'.format(reg=reg, sig='p' if ped else 's'), size=.04, color=4)
                    tit.Draw()
                    titles.append(tit)
                l1.Draw()
                l2.Draw() if l2 else self.do_nothing()
                lines.append([l1, l2])
                starts.append(lst[0])
        gr.Draw('[]')
        gr.Draw('p')
        self._add_buckets()
        save_name = 'PedestalRegions' if ped else 'SignalRegions'
        self.save_plots(save_name, sub_dir=self.ana_save_dir, ch=None)
        self.histos.append([h, c, gr, lines, titles])

    def _add_buckets(self, ymin, ymax, xmin, xmax):
        axis = []
        labels = []
        arrows = []
        start = self.run.signal_regions['b'][0] % 40
        stop = int(.8 * xmax) if xmax > 500 else int(xmax)
        bucket0 = self.run.signal_regions['b'][0] / 40
        x_range = xmax - xmin
        y_range = ymax - ymin
        l = self.make_tlatex(xmin - .015 * x_range, ymin - 0.1 * y_range, 'Bucket:', align=30, color=kGreen + 2, size=0.03)
        l.Draw()
        labels.append(l)

        # peak_fit = self.fit_peak_values(draw=False) if hasattr(self, 'fit_peak_values') else 0
        peak_fit = self.run.signal_regions['a'][0] / 2.
        for i, x in enumerate(xrange(start, stop, 20), -bucket0):
            a = self.make_tgaxis(x, ymin - 0.12 * y_range, ymin - 0.05 * y_range, '', kGreen + 2)
            if x <= stop - 20:
                l = self.make_tlatex(x + 10, ymin - 0.1 * y_range, str(i), align=20, color=kGreen + 2, size=0.03)
                labels.append(l)
                l.Draw()
                if peak_fit:
                    pos = peak_fit % 20
                    if i == -2:
                        l1 = self.make_tlatex(x + pos, ymin + 0.05 * y_range, 'Average Peak Position', color=kOrange + 7, size=0.03)
                        l1.Draw()
                        labels.append(l1)
                    ar = TArrow(x + pos, ymin + 1, x + pos, ymin + 0.04 * y_range, .005, '<|')
                    ar.SetLineWidth(2)
                    ar.SetFillColor(kOrange + 7)
                    ar.SetLineColor(kOrange + 7)
                    ar.Draw()
                    arrows.append(ar)
            a.Draw()
            axis.append(a)
        self.histos.append([axis, labels, arrows])

    def draw_peak_integrals(self, event=None, add_buckets=True, show=True):
        h = self.__draw_single_wf(event=event, show=False)
        self.format_histo(h, title='Waveform', name='wf', x_tit='Time [ns]', y_tit='Signal [mV]', markersize=.8, y_off=.4, stats=0, tit_size=.05)
        xmin, xmax = self.run.signal_regions['e'][0] / 2 - 20, self.run.signal_regions['e'][1] / 2
        h.GetXaxis().SetRangeUser(xmin, xmax)
        stuff = self.draw_histo(h, show=show, lm=.06, rm=.045, bm=.2, x=3000, y=1000, grid=True)
        gROOT.SetBatch(1) if not show else self.do_nothing()
        sleep(.5)
        # draw line at found peak and pedestal region
        ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
        peak_pos, ped_pos = self.__draw_peak_pos(event, ymin, ymax)
        # draw error bars
        gr1 = self.make_tgrapherrors('gr1', '', color=kGreen + 2, marker_size=0, asym_err=True, width=2)
        gr2 = self.make_tgrapherrors('gr2', '', color=kCyan - 3, marker_size=0, asym_err=True, width=2)
        gStyle.SetEndErrorSize(5)
        i = 0
        y = ymax - ymin
        for int_, lst in self.run.peak_integrals.iteritems():
            if len(int_) < 3:
                gr1.SetPoint(i, peak_pos, ymax - y * ((i + 1) / 6. + 1 / 3.))
                gr2.SetPoint(i, ped_pos, ymax - y * ((i + 1) / 6. + 1 / 3.))
                gr1.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr1.SetPointError(i, .5, .5, 0, 0)
                gr2.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr2.SetPointError(i, .5, .5, 0, 0)
                l1 = self.make_tlatex(gr1.GetX()[i], gr1.GetY()[i] + 5, ' ' + int_, color=kGreen + 2, align=10)
                gr1.GetListOfFunctions().Add(l1)
                i += 1
        for gr in [gr1, gr2]:
            gr.Draw('[]')
            gr.Draw('p')
        self._add_buckets(ymin, ymax, xmin, xmax) if add_buckets else self.do_nothing()
        self.save_plots('IntegralPeaks', sub_dir=self.ana_save_dir, ch=None)
        gROOT.SetBatch(0)
        self.histos.append([stuff, gr1, gr2])

    def __draw_peak_pos(self, event, ymin, ymax):
        peak_pos = self.get_peak_position(event) / 2. if hasattr(self, 'get_peak_position') else self.run.signal_regions['a'][0] / 2.
        ped_region = self.PedestalRegion if hasattr(self, 'PedestalRegion') else 'ab'
        ped_pos = self.run.pedestal_regions[ped_region][1] / 2.
        y = ymax - ymin
        l = self.make_tgaxis(peak_pos, ymin, ymax - y / 3., '', 4, 2)
        l2 = self.make_tgaxis(ped_pos, ymin, ymax - y / 3, '', kViolet + 3, 2)
        l2.Draw()
        l.Draw()
        t1 = self.make_tlatex(peak_pos, ymax - y / 3.1, 'found peak', color=4)
        t2 = self.make_tlatex(ped_pos, ymax - y / 3.1, 'ab', color=kViolet + 3)
        t1.Draw()
        t2.Draw()
        self.RootObjects.append([t1, l, l2, t2])
        return peak_pos, ped_pos

    # endregion

    # ============================================================================================
    # region TRACKS
    def show_chi2(self, mode=None, show=True):
        gROOT.SetBatch(1)
        assert mode in ['x', 'y', None], 'mode has to be in {lst}!'.format(lst=['x', 'y', None])
        n_bins = 500 if mode is None else 1000
        mode = 'tracks' if mode is None else mode
        h = TH1F('h', '#chi^{2} in ' + mode, n_bins, 0, 100)
        self.tree.Draw('chi2_{mod}>>h'.format(mod=mode), '', 'goff')
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Chi2 in ' + mode, 1000, 1000)
        c.SetLeftMargin(.13)
        if show or mode == 'tracks':
            yq = zeros(1)
            h.GetQuantiles(1, yq, array([.9]))
            # h.GetXaxis().SetRangeUser(0, yq[0])
        self.format_histo(h, x_tit='#chi^{2}', y_tit='Entries', y_off=1.8)
        h.Draw()
        self.histos.append([h, c])
        gROOT.SetBatch(0)
        return h

    def show_all_chi2(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        histos = [self.show_chi2(mode, show=False) for mode in [None, 'x', 'y']]
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        max_chi2 = int(max([h.GetMaximum() for h in histos])) / 1000 * 1000 + 1000
        histos[0].GetYaxis().SetRangeUser(0, max_chi2)
        histos[0].SetTitle('All #chi^{2}')
        legend = TLegend(.7, .7, .9, .9)
        leg_names = ['#chi^{2} ' + mode for mode in ['of both', 'in x', 'in y']]
        for i, h in enumerate(histos):
            h.SetStats(0)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, leg_names[i], 'l')
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.RootObjects.append([legend, histos, c])
        self.save_plots('Chi2', canvas=c, sub_dir=self.ana_save_dir, ch=None)

    def show_angle(self, mode='x', show=True):
        """
        Displays the angle distribution of the tracks.
        :param mode: has to be eiher 'x' or 'y'
        :param show:
        :return: histogram
        """
        assert mode in ['x', 'y']
        gROOT.SetBatch(1)
        h = TH1F('h', 'Track Angle Distribution in ' + mode, 320, -4, 4)
        self.tree.Draw('slope_{mod}>>h'.format(mod=mode), '', 'goff')
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Angle in ' + mode, 1000, 1000)
        c.SetLeftMargin(.13)
        self.format_histo(h, x_tit='Track Angle [deg]', y_tit='Entries', y_off=1.8, lw=2)
        h.Draw()
        self.RootObjects.append([h, c])
        gROOT.SetBatch(0)
        # a = gROOT.GetListOfCanvases()
        # print a[0]
        self.save_plots('TrackAngle{mod}'.format(mod=mode.upper()), sub_dir=self.ana_save_dir, ch=None)
        return h

    def calc_angle_fit(self, mode='x', show=True):
        pickle_path = self.PickleDir + 'Tracks/AngleFit_{tc}_{run}_{mod}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run_number, mod=mode)

        def func():
            h = self.show_angle(mode, show=show)
            return self.fit_fwhm(h, draw=show)

        fit = func() if show else 0
        return self.do_pickle(pickle_path, func, fit)

    def show_both_angles(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        self.get_color()
        histos = [self.show_angle(mode, show=False) for mode in ['x', 'y']]
        c = TCanvas('c', 'Chi2', 1000, 1000)
        c.SetLeftMargin(.13)
        max_angle = int(max([h.GetMaximum() for h in histos])) / 1000 * 1000 + 1000
        histos[0].GetYaxis().SetRangeUser(0, max_angle)
        legend = TLegend(.7, .7, .9, .9)
        leg_names = ['Angle in ' + mode for mode in ['x', 'y']]
        for i, h in enumerate(histos):
            h.SetStats(0)
            h.SetTitle('Track Angle Distributions')
            h.SetLineColor(self.get_color())
            h.Draw() if not i else h.Draw('same')
            legend.AddEntry(h, leg_names[i], 'l')
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.RootObjects.append([legend, c, histos])
        self.save_plots('TrackAngles', sub_dir=self.ana_save_dir, ch=None)
    # endregion

    # ============================================================================================
    # region PIXEL
    def draw_hitmap(self, show=True, cut=None, plane=None):
        planes = [0, 1, 2, 3] if plane is None else [int(plane)]
        cut = self.Cut.all_cut if cut is None else cut
        histos = [TH2F('h_hm{0}'.format(i_pl), 'Hitmap Plane {0}'.format(i_pl), 52, 0, 52, 80, 0, 80) for i_pl in planes]
        for plane in planes:
            cut_string = cut + TCut('plane == {0}'.format(plane))
            self.tree.Draw('row:col>>h_hm{0}'.format(plane), cut_string, 'goff')
            self.format_histo(histos[plane], x_tit='col', y_tit='row')
            self.RootObjects.append(self.save_histo(histos[plane], 'HitMap{0}'.format(plane), False, self.ana_save_dir, draw_opt='colz'))
        gROOT.SetBatch(1) if not show else self.do_nothing()
        c = TCanvas('c_hm', 'Hitmaps', 2000, 2000)
        c.Divide(2, 2)
        for i, h in enumerate(histos, 1):
            h.SetStats(0)
            pad = c.cd(i)
            pad.SetBottomMargin(.15)
            h.Draw('colz')
        self.save_plots('HitMap', sub_dir=self.ana_save_dir, ch=None)
        gROOT.SetBatch(1)

    # endregion
    # ==============================================
    # region ALIGNMENT
    def check_alignment(self, binning=5000, draw=True, save_plot=True):
        pickle_path = 'Configuration/Individual_Configs/Alignment/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run.run_number)

        def func():
            if self.DUTType == 'pad':
                nbins = self.run.n_entries / binning
                h = TProfile('h', 'Pulser Rate', nbins, 0, self.run.n_entries)
                self.tree.Draw('(@col.size()>1)*100:Entry$>>h', 'pulser', 'goff')
                self.format_histo(h, name='align', title='Event Alignment', x_tit='Event Number', y_tit='Hits per Event @ Pulser Events [%]', y_off=1.3, stats=0, fill_color=821)
                h.GetYaxis().SetRangeUser(0, 105)
                if save_plot:
                    self.RootObjects.append(self.save_histo(h, 'EventAlignment', draw, self.ana_save_dir, draw_opt='hist'))
                align = self.__check_alignment_histo(h)
                return align
            else:
                # todo put some function for the pixel here!
                pass

        aligned = func() if draw else None
        aligned = self.do_pickle(pickle_path, func, aligned)
        if not aligned:
            msg = 'The events of RUN {run} are not aligned!'.format(run=self.run_number)
            print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(msg)) * '!', msg=msg)
        return aligned

    def __check_alignment_histo(self, histo):
        h = histo
        for bin_ in xrange(h.FindBin(self.StartEvent), h.GetNbinsX()):
            if h.GetBinContent(bin_) > 40:
                return False
        return True
    # endregion

    # ==============================================
    # region SHOW & PRINT

    def draw_pix_map(self, n=1, start=None, plane=1):
        start_event = self.StartEvent if start is None else start
        h = TH2F('h', 'Pixel Map', 52, 0, 51, 80, 0, 79)
        self.tree.GetEntry(start_event)
        for pln, col, row, adc in zip(self.tree.plane, self.tree.col, self.tree.row, self.tree.adc):
            if pln == plane:
                h.SetBinContent(col + 1, row + 1, -adc)
        c = TCanvas('c', 'Pixel Map', 1000, 1000)
        c.SetBottomMargin(.15)
        c.SetRightMargin(.14)
        h.SetStats(0)
        h.GetZaxis().SetTitle('adc [au]')
        h.GetZaxis().SetTitleOffset(1.3)
        self.format_histo(h, x_tit='col', y_tit='row')
        h.Draw('colz')
        self.histos[0] = [c, h]
        self.save_plots('PixMapPlane{pln}{evts}'.format(pln=plane, evts=n), sub_dir=self.ana_save_dir, ch=None)

    def draw_preliminary(self):
        c = gROOT.GetListOfCanvases()[-1]
        text = TText((c.GetUxmax() - c.GetUxmin()) / 2., (c.GetUymax() - c.GetUymin()) / 2., "Preliminary")
        text.SetTextColor(19)
        text.SetTextSize(.2)
        text.SetTextAngle(30)
        text.SetTextAlign(20)
        h = None
        for obj in c.GetListOfPrimitives():
            print obj.IsA().GetName()
            if obj.IsA().GetName() in ['TH1F', 'TH2F', 'TGraph', 'TGraphErrors']:
                h = obj
        text.Draw()
        h.Draw('same')
        c.RedrawAxis()
        self.histos[1] = text

    # ==============================================
    # region RUN FUNCTIONS

    def get_event_at_time(self, time_sec):
        return self.run.get_event_at_time(time_sec)

    def get_flux(self):
        return self.run.get_flux()

    # endregion


if __name__ == "__main__":
    ana_parser = ArgumentParser()
    ana_parser.add_argument('run', nargs='?', default=392, type=int)
    args = ana_parser.parse_args()
    this_run = args.run
    print '\nAnalysing run', this_run, '\n'
    z = Analysis(this_run)
