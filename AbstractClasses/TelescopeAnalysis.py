from argparse import ArgumentParser
from collections import OrderedDict
from copy import deepcopy
from time import sleep

from ROOT import TCanvas, TH2F, gROOT, TProfile, TH1F, TLegend, gStyle, kGreen, TArrow, kOrange, kViolet, kCyan, TText
from numpy import array, zeros

from Elementary import Elementary
from RunClass import Run
from Cut import Cut


class Analysis(Elementary):
    """ Class for the analysis of the non-channel specific stuff of a single run. """

    def __init__(self, run, diamonds=3, verbose=False, high_low_rate=None):
        """
        An Analysis Object collects all the information and data for the analysis of one single run.
        The most important fields are:
            .run        - instance of Run class, containing run info and data
            .cut        - dict containing two instances of Cut class
            .RunInfo    - dict containing run infos
            .Pads       - dict containing two intances of
        The most important methods are:
            .CreateMeanSignalHistogram(channel)
            .Draw(varexp)
            .DrawRunInfo()
            .GetCut(channel)
            .GetEventAtTime(dt)
            .MakePreAnalysis()
            .ShowPulserRate()
        :param run:         run object of type "Run" or integer run number
        :param diamonds:    An integer number defining the diamonds activated for analysis: 0x1=ch0 (diamond 1) 0x2=ch3 (diamond 2)
        :param verbose:     if True, verbose printing is activated
        """
        Elementary.__init__(self, verbose=verbose)

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
        self.ana_save_dir = '{tc}_{run}'.format(tc=self.TESTCAMPAIGN[2:], run=self.run.run_number)

        # DUT
        self.DUTType = self.run.DUTType
        
        # tree
        self.tree = self.run.tree

        # miscellaneous
        self.channel = self.channel if hasattr(self, 'channel') else None

        # regions // ranges // for PAD
        if self.DUTType == 'pad':
            self.IntegralNames = self.get_integral_names()
            self.SignalRegion = self.ana_config_parser.get('BASIC', 'signal_region')
            self.PedestalRegion = self.ana_config_parser.get('BASIC', 'pedestal_region')
            self.PeakIntegral = self.ana_config_parser.get('BASIC', 'peak_integral')

        # general for pads and pixels
        if self.DUTType == 'pad':
            self.Cut = Cut(self)
            self.StartEvent = self.Cut.CutConfig['EventRange'][0]
            self.EndEvent = self.Cut.CutConfig['EventRange'][1]
        else:
            self.StartEvent = 1  # DA: for now... TODO pixel cuts!


        # save histograms // canvases
        self.signal_canvas = None
        self.histos = []
        self.canvases = {}
        self.lines = {}

        # alignment
        self.IsAligned = self.check_alignment(draw=False)

    # ============================================================================================
    # region INIT

    def get_integral_names(self):
        names = OrderedDict()
        self.tree.GetEntry(0)
        for i, name in enumerate(self.tree.IntegralNames):
            names[name] = i
        return names

    def init_run(self, run):
        if not isinstance(run, Run):
            assert type(run) is int, 'run has to be either a Run instance or an integer run number'
            return Run(run, self.diamonds)
        else:
            assert run.run_number is not None, 'No run selected, choose run.SetRun(run_nr) before you pass the run object'
            return run
    # endregion

    # ============================================================================================
    # region REGIONS AND PEAK INTEGRAL

    def __draw_single_wf(self, event=None, show=True):
        start = self.StartEvent if event is None else event
        if hasattr(self, 'draw_waveforms') and self.run.wf_exists(self.channel):
            h = self.draw_waveforms(n=1, show=show, start_event=start)
        else:
            h = TH2F('regions', '', 1024, 0, 511, 1000, -200, 50)
            if self.run.wf_exists(0):
                self.tree.Draw('wf0:Iteration$/2>>regions', self.Cut.all_cut, 'goff', 1, start)
        if not show:
            gROOT.SetBatch(1)
        c = TCanvas('c2', 'Regions', 1000, 500)
        c.SetMargin(.075, .045, .1, .1)
        c.SetGrid()
        h.SetStats(0)
        h.GetXaxis().SetNdivisions(26)
        self.format_histo(h, markersize=0.3, x_tit='Time [ns]', y_tit='Signal [au]')
        h.Draw()
        gROOT.SetBatch(0)
        self.histos[0] = [c, h]
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
        self.histos[0] = [h, c, gr, lines, titles]

    def _add_buckets(self, canvas=None):
        c = gROOT.GetSelectedPad() if canvas is None else canvas
        axis = []
        labels = []
        arrows = []
        start = self.run.signal_regions['b'][0] % 40
        stop = int(.8 * c.GetUxmax()) if c.GetUxmax() > 500 else int(c.GetUxmax())
        print 'sta-sto:', start, stop
        bucket0 = self.run.signal_regions['b'][0] / 40
        x_range = c.GetUxmax() - c.GetUxmin()
        y_range = c.GetUymax() - c.GetUymin()
        l = self.make_tlatex(c.GetUxmin() - .015 * x_range, c.GetUymin() - 0.1 * y_range, 'Bucket:', align=30, color=kGreen + 2, size=0.03)
        l.Draw()
        labels.append(l)
        l1 = self.make_tlatex(c.GetUxmin() - .015 * x_range, c.GetUymin() + 0.04 * y_range, 'Peak:', align=30, color=kOrange + 7, size=0.03)
        l1.Draw()
        labels.append(l1)
        peak_fit = self.fit_peak_values(draw=False) if hasattr(self, 'fit_peak_values') else 0
        for i, x in enumerate(xrange(start, stop, 20), -bucket0):
            a = self.make_tgaxis(x, c.GetUymin() - 0.12 * y_range, c.GetUymin() - 0.05 * y_range, '', kGreen + 2)
            if x <= stop - 20:
                l = self.make_tlatex(x + 10, c.GetUymin() - 0.1 * y_range, str(i), align=20, color=kGreen + 2, size=0.03)
                labels.append(l)
                l.Draw()
                if peak_fit:
                    pos = peak_fit.Parameter(1) % 20
                    ar = TArrow(x + pos, c.GetUymin() + 1, x + pos, c.GetUymin() + 0.04 * y_range, .005, '<|')
                    ar.SetLineWidth(2)
                    ar.SetFillColor(kOrange + 7)
                    ar.SetLineColor(kOrange + 7)
                    ar.Draw()
                    arrows.append(ar)
            a.Draw()
            axis.append(a)
        self.histos[1] = [axis, labels, arrows]

    def draw_peak_integrals(self, event=None):
        h = self.__draw_single_wf(event=event, show=False)
        c = TCanvas('c1', 'Regions', 1000, 500)
        c.SetMargin(.075, .045, .2, .1)
        c.SetGrid()
        self.format_histo(h, title='Peak Integrals', markersize=.5)
        h.GetXaxis().SetRangeUser(self.run.signal_regions['e'][0] / 2 - 20, self.run.signal_regions['e'][1] / 2)
        h.Draw()
        sleep(.5)
        # draw line at found peak and pedestal 'ab'
        peak_pos = self.get_peak_position(event) / 2. if hasattr(self, 'get_peak_position') else self.run.signal_regions['a'][0] / 2.
        ped_pos = self.run.pedestal_regions['ab'][1] / 2.
        l = self.make_tgaxis(peak_pos, c.GetUymin(), c.GetUymax() - 100, '', 4, 2)
        l2 = self.make_tgaxis(ped_pos, c.GetUymin(), c.GetUymax() - 100, '', kViolet + 3, 2)
        l2.Draw()
        l.Draw()
        t1 = self.make_tlatex(peak_pos, c.GetUymax() - 97, 'found peak', color=4)
        t2 = self.make_tlatex(ped_pos, c.GetUymax() - 97, 'ab', color=kViolet + 3)
        t1.Draw()
        t2.Draw()
        # draw error bars
        gr1 = self.make_tgrapherrors('gr1', '', color=kGreen + 2, marker_size=0, asym_err=True, width=3)
        gr2 = self.make_tgrapherrors('gr2', '', color=kCyan - 3, marker_size=0, asym_err=True, width=3)
        gStyle.SetEndErrorSize(4)
        i = 0
        for int_, lst in self.run.peak_integrals.iteritems():
            if len(int_) < 3:
                gr1.SetPoint(i, peak_pos, c.GetUymax() - 30 * (i + 1) - 100)
                gr2.SetPoint(i, ped_pos, c.GetUymax() - 33 * (i + 1) - 100)
                gr1.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr1.SetPointError(i, .5, .5, 0, 0)
                gr2.SetPointError(i, lst[0] / 2., lst[1] / 2., 0, 0) if lst[1] - lst[0] > 1 else gr2.SetPointError(i, .5, .5, 0, 0)
                l1 = self.make_tlatex(gr1.GetX()[i], gr1.GetY()[i] + 5, ' ' + int_, color=kGreen + 2, align=10)
                gr1.GetListOfFunctions().Add(l1)
                i += 1
        for gr in [gr1, gr2]:
            gr.Draw('[]')
            gr.Draw('p')
        self._add_buckets()
        self.save_plots('IntegralPeaks', sub_dir=self.ana_save_dir, ch=None)
        self.histos[0] = [gr1, gr2, c, l, t1, h, l2, t2]
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
            h.GetXaxis().SetRangeUser(0, yq[0])
        self.format_histo(h, x_tit='#chi^{2}', y_tit='Entries', y_off=1.8)
        h.Draw()
        self.histos[0] = h
        self.canvases[0] = c
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
            self.histos[i] = h
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.histos.append(legend)
        self.save_plots('Chi2', canvas=c, sub_dir=self.ana_save_dir, ch=None)
        self.canvases[0] = c

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
        self.histos[0] = h
        self.canvases[0] = c
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

        fit = self.do_pickle(pickle_path, func)
        return fit

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
            self.histos[i] = h
        legend.Draw()
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
        self.canvases[0] = c
        self.histos.append(legend)
        self.save_plots('TrackAngles', sub_dir=self.ana_save_dir, ch=None)
    # endregion

    # ==============================================
    # region ALIGNMENT
    def check_alignment(self, binning=5000, draw=True):
        pickle_path = 'Configuration/Individual_Configs/Alignment/{tc}_{run}.pickle'.format(tc=self.TESTCAMPAIGN, run=self.run.run_number)

        def func():
            if self.DUTType == 'pad':
                gROOT.SetBatch(1)
                nbins = self.run.n_entries / binning
                h = TProfile('h', 'Pulser Rate', nbins, 0, self.run.n_entries)
                self.tree.Draw('(@col.size()>1)*100:Entry$>>h', 'pulser', 'goff')
                self.format_histo(h, name='align', title='Event Alignment', x_tit='Event Number', y_tit='Hits per Event @ Pulser Events [%]', y_off=1.3)
                h.GetYaxis().SetRangeUser(0, 100)
                if draw:
                    gROOT.SetBatch(0)
                c = TCanvas('c', 'Pulser Rate Canvas', 1000, 1000)
                h.SetStats(0)
                h.Draw('hist')
                self.save_plots('EventAlignment', sub_dir=self.ana_save_dir, ch=None)
                gROOT.SetBatch(0)
                self.histos[0] = [h, c]
                align = self.__check_alignment_histo(h)
                return align
            else:
                return True  # DA: pixel doesn't have pulser, todo make correlations between planes to test misalignment

        aligned = func() if draw else self.do_pickle(pickle_path, func)
        if not aligned:
            msg = 'The events of RUN {run} are not aligned!'.format(run=self.run_number)
            print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(msg)) * '!', msg=msg)
        return aligned

    def find_alignment_offset(self):
        offsets = [i for i in xrange(-5, 5) if i]
        h = TH1F('h', 'Pixel Hits @ Pulser Events', 20, 0, 20)
        right_offset = None
        for offset in offsets:
            pulser_events = 0
            for event in xrange(self.StartEvent, self.run.n_entries):
                print '\rpulser events: {0:04d}'.format(pulser_events),
                if pulser_events >= 1000:
                    break
                self.tree.GetEntry(event)
                if self.tree.pulser:
                    pulser_events += 1
                    self.tree.GetEntry(event + offset)
                    hits = len(self.tree.col)
                    h.Fill(hits)
            h.SetName(str(offset))
            h.Draw()
            sleep(.2)
            c = gROOT.GetSelectedPad()
            c.Update()
            sleep(2)
            # find offset
            glob_max = [h.GetMaximumBin(), h.GetMaximum()]
            h.GetXaxis().SetRangeUser(1, 20)
            hit_max = [h.GetMaximumBin(), h.GetMaximum()]
            if glob_max[0] == 1 and glob_max[1] > 2 * hit_max[1]:
                right_offset = offset
                break
            h.Reset()
        h.Draw()
        self.histos[0] = h
        print '\nThe event offset is {off}'.format(off=right_offset)

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

    def ShowSignalSpread(self, channel=None, cut=""):
        if channel == None:
            channels = self.run.get_active_channels()
            namesuffix = ""
        else:
            channels = [channel]
            namesuffix = "_ch{ch}".format(ch=channel)

        canvas = TCanvas("{run}signalspreadcanvas".format(run=self.run.run_number), "{run}signalspreadcanvas".format(run=self.run.run_number), 800, len(channels) * 300)
        canvas.Divide(1, len(channels))

        for ch in channels:
            canvas.cd(channels.index(ch) + 1)
            if cut == "":
                thiscut = self.GetCut(ch)
            else:
                thiscut = cut.format(channel=ch)

            self.Draw(("sig_spread[{channel}]>>signalspread{run}{channel}(400, 0, 400)").format(channel=ch, run=self.run.run_number), thiscut)
            hist = gROOT.FindObject("signalspread{run}{channel}".format(channel=ch, run=self.run.run_number))
            if hist: hist.SetStats(0)

        canvas.Update()
        self.save_plots("Run{run}_SignalSpread{ns}.png".format(run=self.run.run_number, ns=namesuffix), canvas=canvas, sub_dir="Cuts")
        self.if_wait("Peak Position shown")

    def AnalyzeMultiHitContribution(self):

        single_particle = self.run.tree.Draw("1",
                                             "!pulser&&clusters_per_plane[0]>=1&&clusters_per_plane[1]>=1&&clusters_per_plane[2]>=1&&clusters_per_plane[3]>=1&&(clusters_per_plane[0]>1||clusters_per_plane[1]>1||clusters_per_plane[2]>1||clusters_per_plane[3]>1)")
        multiparticle = self.run.tree.Draw("1", "!pulser&&clusters_per_plane[0]<=1&&clusters_per_plane[1]<=1&&clusters_per_plane[2]<=1&&clusters_per_plane[3]<=1")

        total = single_particle + multiparticle

        return 1. * multiparticle / total

if __name__ == "__main__":
    ana_parser = ArgumentParser()
    ana_parser.add_argument('run', nargs='?', default=392, type=int)
    args = ana_parser.parse_args()
    this_run = args.run
    print '\nAnalysing run', this_run, '\n'
    z = Analysis(this_run)
