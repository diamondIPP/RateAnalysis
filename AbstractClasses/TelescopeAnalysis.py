from argparse import ArgumentParser
from copy import deepcopy

from ROOT import TCanvas, TH2F, gROOT, TH1F, TLegend, gStyle, kGreen, TCut, TF1, TGraph, TH1I
from numpy import log, zeros

from Elementary import Elementary
from Run import Run
from Cut import Cut
from Utils import *
from Plots import Plots
from Langaus import Langau
from InfoLegend import InfoLegend


class Analysis(Elementary):
    """ Class for the analysis of the non-channel specific stuff of a single run. """

    def __init__(self, run, high_low_rate=None, binning=5000):
        """
        Parent class for all analyses, which contains all the basic stuff about the Telescope.
        :param run:             run object of type "Run" or integer run number
        :param high_low_rate:   list of highest and lowest rate runs for an analysis collection
        """
        self.run = self.init_run(run)
        self.RunNumber = self.run.RunNumber
        Elementary.__init__(self, verbose=run.verbose)
        self.histos = []
        self.RootObjects = []

        # basics
        self.run.analysis = self
        self.RunInfo = deepcopy(self.run.RunInfo)
        self.lowest_rate_run = high_low_rate['min'] if high_low_rate is not None else self.run.RunNumber
        self.highest_rate_run = high_low_rate['max'] if high_low_rate is not None else self.run.RunNumber
        self.TelSaveDir = '{run}'.format(run=self.run.RunNumber)
        self.set_titles()

        # DUT
        self.DUTType = self.run.DUTType

        # tree
        self.tree = self.run.tree

        # miscellaneous
        self.channel = self.channel if hasattr(self, 'channel') else None

        # general for pads and pixels
        self.Cut = Cut(self, skip=run.tree is None)

        if run.tree:
            self.NRocs = self.run.NPlanes
            self.StartEvent = self.Cut.CutConfig['EventRange'][0]
            self.EndEvent = self.Cut.CutConfig['EventRange'][1]
            self.Plots = Plots(self.run)

            # binning TODO: move to plots class
            self.BinSize = binning
            self.binning = self.__get_binning()
            self.time_binning = self.get_time_binning()
            self.n_bins = len(self.binning)

        self.InfoLegend = InfoLegend(self)

        # save histograms // canvases
        self.signal_canvas = None

    # ============================================================================================
    # region INIT

    @staticmethod
    def init_run(run):
        """ Make some assertions if the correct run instance was provided """
        if not isinstance(run, Run):
            raise TypeError('You have to pass a Run instance!')
            # return Run(run, 3, load_tree=load_tree, verbose=verbose)
        if run.RunNumber is None:
            raise ValueError('No run selected, choose run.SetRun(run_nr) before you pass the run object')
        return run
    # endregion

    # ============================================================================================
    # region REGIONS AND PEAK INTEGRAL

    def __draw_single_wf(self, event=None, show=True, tcorr=False):
        start = self.StartEvent if event is None else event
        if hasattr(self, 'draw_waveforms') and self.run.wf_exists(self.channel):
            return self.draw_waveforms(n=1, show=show, start_event=start, t_corr=tcorr)
        else:
            h = TH2F('regions', '', 1024, 0, 511, 1000, -200, 50)
            if self.run.wf_exists(0):
                self.tree.Draw('wf0:Iteration$/2>>regions', self.Cut.all_cut, 'goff', 1, start)
        h.GetXaxis().SetNdivisions(520)
        self.format_histo(h, markersize=0.3, x_tit='Time [ns]', y_tit='Signal [au]', stats=0)
        self.save_histo(h, 'Regions', show, self.TelSaveDir, lm=.075, rm=.045, x_fac=1.5, y_fac=.5)
        return h

    def draw_regions(self, ped=True, event=None, show=True):
        h = self.__draw_single_wf(event=event, show=False)
        self.draw_histo(h, show=show, lm=.07, rm=.045, bm=.2, x=1.5, y=.75)
        ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
        ydiff = ymax - ymin
        h.SetTitle('Pedestal Regions' if ped else 'Signal Regions')
        starts = []
        regions = self.run.pedestal_regions if ped else self.run.signal_regions
        gr = self.make_tgrapherrors('gr', '', color=2, marker_size=0, width=3)
        i = 0
        gStyle.SetEndErrorSize(4)
        sleep(.5)
        for reg, lst in regions.iteritems():
            if len(reg) < 3:
                offset = ydiff * .4 if not all(i[0] < lst[0] < i[1] for i in starts) or not all(i[0] < lst[1] < i[1] for i in starts) else ydiff * .2
                if lst[1] - lst[0] > 1:
                    gr.SetPoint(i, (lst[1] + lst[0]) / 4., ymax - offset)
                    gr.SetPointError(i, (lst[1] - lst[0]) / 4., 0)
                    l = self.draw_tlatex(gr.GetX()[i], gr.GetY()[i] + 3, '{sig}{reg}'.format(reg=reg, sig='p' if ped else 's'), color=2, size=.04)
                    gr.GetListOfFunctions().Add(l)
                    i += 1
                if lst[1] - lst[0] > 1:
                    self.draw_vertical_line(lst[0] / 2, ymin, ymax - offset, color=2)
                    starts.append(lst)
                else:
                    self.draw_vertical_line(lst[0] / 2, ymin, ymax - ydiff * .5, color=4, w=2)
                    self.draw_tlatex(lst[0] / 2, ymax - ydiff * .49, '{sig}{reg}'.format(reg=reg, sig='p' if ped else 's'), size=.04, color=4)
                self.draw_vertical_line(lst[1] / 2, ymin, ymax - offset, color=2) if lst[1] - lst[0] > 1 else do_nothing()
        gr.Draw('[]')
        gr.Draw('p')
        self._add_buckets(ymin, ymax, avr_pos=6, full_line=True)
        save_name = 'PedestalRegions' if ped else 'SignalRegions'
        self.save_plots(save_name, both_dias=True)
        self.histos.append([h, gr])

    def _add_buckets(self, ymin, ymax, xmin=0, xmax=512, avr_pos=-2, full_line=False, size=.03):
        start = self.run.signal_regions['b'][0] % 40 / 2
        stop = int(.8 * xmax) if xmax > 500 else int(xmax)
        bucket0 = self.run.signal_regions['b'][0] / 40
        x_range = xmax - xmin
        y_range = ymax - ymin
        self.draw_tlatex(xmin - .015 * x_range, ymin - 0.18 * y_range, 'Bucket:', align=30, color=418, size=size)
        peak_fit = self.run.signal_regions['a'][0] / 2.
        for i, x in enumerate(xrange(start, stop, 20), -bucket0):
            y2 = ymax if full_line else ymin - 0.1 * y_range
            self.draw_vertical_line(x, ymin - 0.22 * y_range, y2, 418, style=3 if full_line else 1, tline=True)
            if x <= stop - 20:
                self.draw_tlatex(x + 10, ymin - 0.18 * y_range, str(i), align=20, color=418, size=size)
                if peak_fit:
                    pos = peak_fit % 20
                    p_pos = round_down_to(x, 20) + pos
                    p_pos += 20 if p_pos < x else 0
                    if i == avr_pos:
                        self.draw_tlatex(p_pos, ymin + 0.05 * y_range, 'Average Peak Position', color=807, size=size)
                    self.draw_arrow(p_pos, p_pos, ymin + 1, ymin + 0.04 * y_range, col=807, width=2)

    def draw_peak_integrals(self, event=None, buckets=True, show=True, main=False):
        old_count = self.count
        event = self.StartEvent if event is None else event
        h, n = self.__draw_single_wf(event=event, show=False, tcorr=True)
        h.GetYaxis().SetNdivisions(305)
        self.format_histo(h, title='Peak Integrals', name='wf', x_tit='Time [ns]', y_tit='Signal [mV]', markersize=.8, y_off=.5, stats=0, tit_size=.07, lab_size=.06)
        xmin, xmax = self.run.signal_regions['e'][0] / 2 - 20, self.run.signal_regions['e'][1] / 2
        h.GetXaxis().SetRangeUser(xmin, xmax)
        self.draw_histo(h, show=show, lm=.07, rm=.045, bm=.24, x=1.5, y=.5, gridy=True, gridx=True)
        gROOT.SetBatch(1) if not show else do_nothing()
        sleep(.5)
        # draw line at found peak and pedestal region
        ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
        peak_pos, ped_pos = self.__draw_peak_pos(event + n - 1, ymin, ymax)
        # draw error bars
        gr1 = self.make_tgrapherrors('gr1', '', color=418, marker_size=0, asym_err=True, width=2)
        gr2 = self.make_tgrapherrors('gr2', '', color=429 - 3, marker_size=0, asym_err=True, width=2, style=2)
        gStyle.SetEndErrorSize(5)
        i = 0
        y = ymax - ymin
        spacing = 7.
        for int_, lst in self.run.peak_integrals.iteritems():
            if main:
                if hasattr(self, 'PeakIntegral') and int_ != self.PeakIntegral:
                    continue
            if len(int_) < 3:
                gr1.SetPoint(i, peak_pos, ymax - y * ((i + 1) / spacing + 1 / 3.))
                gr2.SetPoint(i, ped_pos, ymax - y * ((i + 1) / spacing + 1 / 3.))
                gr1.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr1.SetPointError(i, .5, .5, 0, 0)
                gr2.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr2.SetPointError(i, .5, .5, 0, 0)
                if not main:
                    l1 = self.draw_tlatex(gr1.GetX()[i], gr1.GetY()[i] + 5, ' ' + int_, color=kGreen + 2, align=10)
                    gr1.GetListOfFunctions().Add(l1)
                i += 1
        for gr in [gr1, gr2]:
            gr.Draw('[]')
            gr.Draw('p')
        self._add_buckets(ymin, ymax, xmin, xmax, full_line=True, size=.05) if buckets else do_nothing()
        self.save_plots('IntegralPeaks', both_dias=True)
        gROOT.SetBatch(0)
        self.count = old_count
        self.histos.append([gr1, gr2])

    def __draw_peak_pos(self, event, ymin, ymax):
        peak_pos = self.get_peak_position(event, tcorr=True) if hasattr(self, 'get_peak_position') else self.run.signal_regions['a'][0] / 2.
        ped_region = self.PedestalRegion if hasattr(self, 'PedestalRegion') else 'ab'
        ped_pos = self.run.pedestal_regions[ped_region][1] / 2.
        y = ymax - ymin
        self.draw_vertical_line(peak_pos, ymin, ymax - y / 3., color=418, w=2, name='peak')
        self.draw_vertical_line(ped_pos, ymin, ymax - y / 3, color=429, w=2, name='ped')
        t1 = self.draw_tlatex(peak_pos, ymax - y / 3.1, 'peak', color=418, size=.07)
        t2 = self.draw_tlatex(ped_pos, ymax - y / 3.1, 'pedestal', color=429, size=.07)
        t1.Draw()
        t2.Draw()
        self.RootObjects.append([t1, t2])
        return peak_pos, ped_pos

    # endregion

    # ============================================================================================
    # region TRACKS
    def show_chi2(self, mode=None, show=True, save=True, fit=False, prnt=True):
        assert mode in ['x', 'y', None], 'mode has to be in {lst}!'.format(lst=['x', 'y', None])
        n_bins = 500 if mode is None else 250
        mode = 'tracks' if mode is None else mode
        set_root_output(False)
        h = TH1F('h', '#chi^{2} in ' + mode, n_bins, 0, 100 if mode == 'tracks' else 50)
        self.tree.Draw('chi2_{mod}>>h'.format(mod=mode), '', 'goff')
        if show:
            yq = zeros(1)
            h.GetQuantiles(1, yq, array([.99]))
            h.GetXaxis().SetRangeUser(0, yq[0])
        if fit:
            f = TF1('f', '[0]*TMath::GammaDist(x, {ndf}/2, 0, 2)'.format(ndf=4 if mode == 'tracks' else 2))
            h.Fit(f, 'q')
            h.SetStats(1)
            set_statbox(only_fit=True)
        self.format_histo(h, x_tit='#chi^{2}', y_tit='Number of Entries', y_off=1.8, stats=0)
        self.save_histo(h, 'Chi2{0}'.format(mode.title() if mode is not 'tracks' else 'All'), show, save=save, lm=.13, prnt=prnt)
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
        self.save_plots('Chi2', canvas=c, sub_dir=self.TelSaveDir, both_dias=True)

    def draw_angle_distribution(self, mode='x', show=True, print_msg=True, cut=None):
        """ Displays the angle distribution of the tracks. """
        assert mode in ['x', 'y']
        cut = cut if cut is not None else TCut('')
        set_root_output(False)
        h = TH1F('had', 'Track Angle Distribution in ' + mode, 320, -4, 4)
        self.tree.Draw('{v}_{mod}>>had'.format(v='slope' if self.run.has_branch('slope') else 'angle', mod=mode), cut, 'goff')
        self.format_histo(h, x_tit='Track Angle [deg]', y_tit='Entries', y_off=1.8, lw=2, stats=0)
        self.save_histo(h, 'TrackAngle{mod}'.format(mod=mode.upper()), show, lm=.13, prnt=print_msg)
        return h

    def draw_track_length(self, show=True, save=True, t_dia=500):
        h = TH1F('htd', 'Track Distance in Diamond', 200, t_dia, t_dia + 1)
        draw_var = 'slope' if self.run.has_branch('slope_x') else 'angle'
        length = '{t}*TMath::Sqrt(TMath::Power(TMath::Tan(TMath::DegToRad()*{v}_x), 2) + TMath::Power(TMath::Tan(TMath::DegToRad()*{v}_y), 2) + 1)'.format(t=t_dia, v=draw_var)
        self.tree.Draw('l>>hdd'.format(l=length), 'n_tracks', 'goff')
        self.format_histo(h, x_tit='Distance [#mum]', y_tit='Entries', y_off=2, lw=2, stats=0, fill_color=self.FillColor)
        h.GetXaxis().SetNdivisions(405)
        self.save_histo(h, 'DistanceInDia', show, lm=.16, save=save)
        return h

    def calc_angle_fit(self, mode='x', show=True):
        pickle_path = self.PickleDir + 'Tracks/AngleFit_{tc}_{run}_{mod}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.RunNumber, mod=mode)

        def func():
            h = self.draw_angle_distribution(mode, show=show)
            return self.fit_fwhm(h, draw=show)

        fit = func() if show else None
        return do_pickle(pickle_path, func, fit)

    def show_both_angles(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        self.get_color()
        histos = [self.draw_angle_distribution(mode, show=False) for mode in ['x', 'y']]
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
        self.save_plots('TrackAngles', sub_dir=self.TelSaveDir, both_dias=True)

    def _draw_residuals(self, roc, mode=None, cut=None, x_range=None, fit=False, show=True):
        mode = '' if mode is None else mode.lower()
        cut = TCut(cut) if cut is not None else TCut('')
        set_statbox(only_entries=True, only_fit=True, w=0.3, entries=14) if fit else set_statbox(only_entries=True)
        h = TH1F('htr', '{m} Residuals for Plane {n}'.format(n=roc, m=mode.title()), 800, -.4, .4)
        self.tree.Draw('residuals{m}[{r}]>>htr'.format(m='_{m}'.format(m=mode) if mode else '', r=roc), cut, 'goff')
        self.format_histo(h, name='Fit Result', y_off=2.0, y_tit='Number of Entries', x_tit='Distance [cm]', fill_color=self.FillColor, x_range=x_range)
        self.draw_histo(h, '', show, lm=.16)
        if fit:
            fit = TF1('f', 'gaus(0) + gaus(3)', -.4, .4)
            sigma = get_fwhm(h) / (2 * sqrt(2 * log(2)))
            fit.SetParameters(h.GetMaximum() / 10, 0, sigma * 5, h.GetMaximum(), 0, sigma)
            fit.SetParName(2, '#sigma1')
            fit.SetParName(5, '#sigma2')
            fit.SetNpx(500)
            h.Fit(fit, 'q')
            f2 = TF1('f2', 'gaus', -1, 1)
            f2.SetParameters(fit.GetParameters())
            f2.SetLineStyle(2)
            f2.Draw('same')
            self.RootObjects.append(f2)
        self.save_plots('{m}ResidualsRoc{n}'.format(m=mode.title(), n=roc))
        return h

    def _draw_cluster_size(self, roc, name=None, cut='', show=True):
        h = TH1I('h_cs', 'Cluster Size {d}'.format(d='ROC {n}'.format(n=roc) if name is None else name), 10, 0, 10)
        self.tree.Draw('cluster_size[{d}]>>h_cs'.format(d=roc), TCut(cut), 'goff')
        set_statbox(only_entries=True)
        self.format_histo(h, x_tit='Cluster Size', y_tit='Number of Entries', y_off=1.3, fill_color=self.FillColor)
        self.save_histo(h, 'ClusterSize', show, logy=True)
        return h

    def draw_event(self, event, plane, show=True):
        cut = 'plane == {r}'.format(r=plane)
        h = TH2F('h_ed{i}'.format(i=plane), 'Event Hits for Plane {r}'.format(r=plane), *self.Plots.Settings['2DBins'])
        self.tree.Draw('row:col>>h_ed{i}'.format(i=plane), cut, 'goff', 1, event)
        self.format_histo(h, x_tit='col', y_tit='row', y_off=1.3, stats=0)
        self.save_histo(h, 'EventDisplay{e}_{p}'.format(e=event, p=plane), draw_opt='col', show=show)

    def get_events(self, cut='', prnt=False):
        n = self.tree.Draw('event_number', TCut(cut), 'goff')
        events = [self.tree.GetV1()[i] for i in xrange(n)]
        if prnt:
            print events[:20]
        return events

    # endregion

    # ============================================================================================
    # region PIXEL
    def _draw_occupancy(self, plane, name=None, cluster=True, tel_coods=False, cut='', show=True):
        name = 'ROC {i}'.format(i=plane) if name is None else name
        bins = self.Plots.get_global_bins(sqrt(12)) if tel_coods else self.Plots.Settings['2DBins']
        h = TH2F('h_hm{i}'.format(i=plane), '{h} Occupancy {n}'.format(n=name, h='Hit' if not cluster else 'Cluster'), *bins)
        cut_string = self.Cut.all_cut if cut is None else TCut(cut)
        cut_string += 'plane == {0}'.format(plane) if not cluster else ''
        draw_string = 'cluster_row[{i}]:cluster_col[{i}]' if cluster else 'row:col'
        draw_string = 'cluster_ypos_local[{i}]:cluster_xpos_local[{i}]' if tel_coods else draw_string
        set_statbox(only_entries=True, x=.83)
        self.tree.Draw('{ds}>>h_hm{i}'.format(ds=draw_string.format(i=plane), i=plane), cut_string, 'goff')
        self.format_histo(h, x_tit='col', y_tit='row', y_off=1.2)
        self.save_histo(h, 'HitMap{0}'.format(plane), show, draw_opt='colz', rm=.15)
        return h

    def _draw_occupancies(self, planes=None, cut='', cluster=True, show=True):
        planes = range(4) if planes is None else list(planes)
        histos = [self._draw_occupancy(plane, cluster=cluster, cut=cut, show=False) for plane in planes]
        c = TCanvas('c_hm', 'Hitmaps', 2000, 2000)
        c.Divide(2, 2)
        for i, h in enumerate(histos, 1):
            h.SetStats(0)
            pad = c.cd(i)
            pad.SetBottomMargin(.15)
            h.Draw('colz')
        self.save_plots('HitMap', sub_dir=self.TelSaveDir, both_dias=True, show=show)

    # endregion

    def draw_trigger_phase(self, dut=True, show=True, cut=None):
        cut_string = self.Cut.generate_special_cut(excluded=['trigger_phase']) if cut is None else TCut(cut)
        h = TH1I('h_tp', 'Trigger Phase', 10, 0, 10)
        self.tree.Draw('trigger_phase[{r}]>>h_tp'.format(r=1 if dut else 0), cut_string, 'goff')
        set_statbox(only_entries=True)
        self.format_histo(h, x_tit='Trigger Phase', y_tit='Number of Entries', y_off=1.95, fill_color=self.FillColor)
        self.save_histo(h, '{m}TriggerPhase'.format(m='DUT' if dut else 'Tel'), show, lm=.145)

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
        self.save_plots('PixMapPlane{pln}{evts}'.format(pln=plane, evts=n), sub_dir=self.TelSaveDir, both_dias=True)

    # ==============================================
    # region TIME AND BINNING

    def set_bin_size(self, value):
        if value is None:
            return
        self.BinSize = value
        self.binning = self.__get_binning()
        self.time_binning = self.get_time_binning()
        self.n_bins = len(self.binning)
        return value

    def __get_binning(self):
        jumps = self.Cut.Interruptions
        n_jumps = len(jumps)
        bins = [self.Cut.get_min_event()]
        ind = 0
        for dic in jumps:
            start = dic['i']
            stop = dic['f']
            gap = stop - start
            # continue if first start and stop outside min event
            if stop < bins[-1]:
                ind += 1
                continue
            # if there is a jump from the start
            if start < bins[-1] < stop:
                bins[-1] = stop
                ind += 1
                continue
            # add bins until hit interrupt
            while bins[-1] + self.BinSize < start:
                bins.append(bins[-1] + self.BinSize)
            # two jumps shortly after one another
            if ind < n_jumps - 2:
                next_start = jumps[ind + 1]['i']
                next_stop = jumps[ind + 1]['f']
                if bins[-1] + self.BinSize + gap > next_start:
                    gap2 = next_stop - next_start
                    bins.append(bins[-1] + self.BinSize + gap + gap2)
                else:
                    bins.append(bins[-1] + self.BinSize + gap)
            else:
                bins.append(bins[-1] + self.BinSize + gap)
            ind += 1
        # fill up the end
        if ind == n_jumps - 1 and bins[-1] >= jumps[-1]['f'] or ind == n_jumps:
            while bins[-1] + self.BinSize < self.run.n_entries:
                bins.append(bins[-1] + self.BinSize)
        return bins

    def get_time_binning(self):
        time_bins = []
        for event in self.binning:
            time_bins.append(self.run.get_time_at_event(event))
        return time_bins

    def draw_time(self, show=True):
        entries = self.tree.Draw('time', '', 'goff')
        t = [self.tree.GetV1()[i] for i in xrange(entries)]
        t = [(i - t[0]) / 1000 for i in t]
        gr = TGraph(len(t), array(xrange(len(t)), 'd'), array(t, 'd'))
        gr.SetNameTitle('g_t', 'Time vs Events')
        fit = gr.Fit('pol1', 'qs0')
        self.log_info('Average data taking rate: {r:5.1f} Hz'.format(r=1 / fit.Parameter(1)))
        self.format_histo(gr, x_tit='Entry Number', y_tit='Time [s]', y_off=1.5)
        self.draw_histo(gr, show=show, draw_opt='al', lm=.13, rm=.08)

    def get_time_bins(self, evts_per_bin=None):
        self.set_bin_size(evts_per_bin)
        return [len(self.time_binning), array([self.run.StartTime] + [t / 1000. for t in self.time_binning], 'd')]

    def get_bins(self, binning):
        self.set_bin_size(binning)
        return [len(self.binning) - 1, array(self.binning, 'd')]

    def get_event_at_time(self, time_sec):
        return self.run.get_event_at_time(time_sec)

    # endregion

    # ==============================================
    # region RUN METHODS
    def get_flux(self):
        return self.run.get_flux()

    def has_branch(self, branch):
        return self.run.has_branch(branch)

    # endregion

    def draw_beam_current(self, rel_t=True, show=True):
        if not self.has_branch('beam_current'):
            log_warning('Branch "beam_current" does not exist!')
            return
        n = self.tree.Draw('beam_current:time/1000.>>test', 'beam_current<10000', 'goff')
        current = [self.tree.GetV1()[i] for i in xrange(n)]
        t = [self.tree.GetV2()[i] for i in xrange(n)]
        g = self.make_tgrapherrors('gbc', 'Beam Current', x=t + [t[-1]], y=current + [0])
        self.format_histo(g, x_tit='Time [hh:mm]', y_tit='Beam Current [mA]', fill_color=self.FillColor, markersize=.4, t_ax_off=self.run.StartTime if rel_t else 0,
                          x_range=[g.GetX()[0], g.GetX()[n]])
        self.save_histo(g, 'BeamCurrent', draw_opt='afp', lm=.08, x_fac=1.5, y_fac=.75, ind=None, show=show)
        return g

    def fit_langau(self, h=None, nconv=30, show=True, chi_thresh=8):
        h = self.draw_signal_distribution(show=show) if h is None and hasattr(self, 'draw_signal_distribution') else h
        h = self.draw_pulse_height_disto(show=show) if h is None and hasattr(self, 'draw_pulse_height_disto') else h
        fit = Langau(h, nconv)
        fit.langaufit()
        fit.Fit.Draw('lsame')
        c = get_last_canvas()
        c.Modified()
        c.Update()
        if fit.Chi2 / fit.NDF > chi_thresh:
            self.count += 5
            self.log_info('Chi2 too large ({c:2.2f}) -> increasing number of convolutions by 5'.format(c=fit.Chi2 / fit.NDF))
            fit = self.fit_langau(h, nconv + self.count, chi_thresh=chi_thresh)
        print 'MPV:', fit.Parameters[1]
        self.count = 0
        self.RootObjects.append(fit)
        return fit

    # endregion


if __name__ == "__main__":
    ana_parser = ArgumentParser()
    ana_parser.add_argument('run', nargs='?', default=392, type=int)
    args = ana_parser.parse_args()
    this_run = args.run
    print '\nAnalysing run', this_run, '\n'
    z = Analysis(this_run)
