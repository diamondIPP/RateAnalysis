#!/usr/bin/env python
# --------------------------------------------------------
#       Main class for Rate Pixel Analysis
# created some time in 2016 by D. Sanz (sandiego@phys.ethz.ch), maintained by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from __future__ import print_function

from collections import Counter

from ROOT import TFormula, THStack, TProfile2D, TPie, gRandom, TH3F, TMultiGraph
from numpy import corrcoef, ceil

from dut_analysis import *
from pix_cut import CutPix


class PixAnalysis(DUTAnalysis):
    def __init__(self, run_number, dut, test_campaign=None, tree=True, t_vec=None, verbose=False, prnt=True):

        DUTAnalysis.__init__(self, run_number, dut, test_campaign, tree, t_vec, verbose, prnt)

        # Main
        self.Dut = dut + 3
        self.PX, self.PY = self.Run.PixelSize

        if self.Tree:
            self.Cut = CutPix.from_parent(self.Cut)
            self.IsAligned = self.check_alignment()

            # Pulse Height Calibrations
            self.Fit = TF1('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -500, 255 * 7)
            if self.check_calibration_files():
                self.Parameters = self.load_calibration_fitpars()
                self.Vcals = self.load_vcals()
                self.Points = self.load_calibration_points()

        self.print_finished(prnt=prnt)

    def __del__(self):
        for c in gROOT.GetListOfCanvases():
            c.Close()

    # ----------------------------------------
    # region INIT
    def update_config(self):
        self.Config.read(join(self.Dir, 'config', self.TCString, 'PixelConfig.ini'))

    def load_calibration_files(self, fits=False):
        calibration_list_dir = join(self.Run.Converter.TrackingDir, 'calibration_lists', 'GKCalibrationList_Telescope{n}.txt'.format(n=self.Run.Converter.TelescopeID))
        with open(calibration_list_dir) as f:
            calibration_dir = join(self.Run.Converter.TrackingDir, f.readline().strip('./\n'))
            return [join(calibration_dir, line.strip('\n ') if fits else 'phCalibration_C{}.dat'.format(i)) for i, line in enumerate(f.readlines())]

    def check_calibration_files(self):
        files = concatenate([self.load_calibration_files(), self.load_calibration_files(fits=True)])
        for f in files:
            if not file_exists(f):
                log_warning('Calibration file {} does not exist...'.format(basename(f)))
                return False
        return True

    def load_calibration_fitpars(self, redo=False):
        def f():
            self.check_calibration_files()
            split_at = arange(self.Bins.NRows, self.Bins.NCols * self.Bins.NRows, self.Bins.NRows)  # split at every new column (after n_rows)
            return array([split(genfromtxt(filename, skip_header=3, usecols=arange(4)), split_at) for filename in self.load_calibration_files(fits=True)])
        return do_pickle(self.make_pickle_path('Calibration', 'Pars', run=self.Run.Converter.TelescopeID), f, redo=redo)

    def load_vcals(self, redo=False):
        def func():
            return array([concatenate([genfromtxt(f, 'i2', skip_header=1, max_rows=1)[2:], 7 * genfromtxt(f, 'i2', skip_header=2, max_rows=1)[2:]]) for f in self.load_calibration_files()])
        return do_pickle(self.make_pickle_path('Calibration', 'Vcal', run=self.Run.Converter.TelescopeID), func, redo=redo)

    def load_calibration_points(self, redo=False):
        def func():
            split_at = arange(self.Bins.NRows, self.Bins.NCols * self.Bins.NRows, self.Bins.NRows)
            return [array(split(genfromtxt(f, 'i2', skip_header=4, usecols=arange(self.Vcals[i].size)), split_at)) for i, f in enumerate(self.load_calibration_files())]
        return do_pickle(self.make_pickle_path('Calibration', 'Points', run=self.Run.Converter.TelescopeID), func, redo=redo)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_ph_str(self):
        return 'cluster_charge[{}]'.format(self.Dut)

    def get_pulse_height(self, bin_size=None, cut=None, redo=False, corr=None):
        suffix = '{bins}_{c}'.format(bins=self.Bins.BinSize if bin_size is None else bin_size, c=self.Cut(cut).GetName())
        picklepath = self.make_pickle_path('Ph_fit', 'Fit', self.RunNumber, self.DUT.Number, suf=suffix)

        def f():
            p, fit_pars = self.draw_pulse_height(bin_size=bin_size, cut=self.Cut(cut), show=False)
            return fit_pars

        return make_ufloat(do_pickle(picklepath, f, redo=redo), par=0)

    def get_thresholds(self, cols=None, pix=None, vcal=True):
        columns, rows = split(array(self.Cut.CutConfig['local_fiducial']), 2) if self.Cut.CutConfig['local_fiducial'] is not None else [0, self.Bins.NCols - 1], [0, self.Bins.NRows - 1]
        columns = array([cols]).flatten() if cols is not None else columns
        columns, rows = (full(2, pix[0]), full(2, pix[1])) if pix is not None else (columns, rows)
        dic = {}
        for col in xrange(columns[0], columns[1] + 1):
            for row in xrange(rows[0], rows[1] + 1):
                self.Fit.SetParameters(*self.Parameters[self.Dut][col][row])
                dic[(col, row)] = self.Fit.GetX(0) * (self.Bins.VcalToEl if not vcal else 1)
        return dic

    def get_vcal(self, redo=False):
        h = self.draw_vcal_distribution(show=False, redo=redo)
        return ufloat(h.GetMean(), h.GetMeanError())
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region OCCUPANCY
    def draw_occupancy(self, roc=None, name=None, cluster=True, tel_coods=False, cut='', show=True):
        """ draw hitmap or cluster map """
        name = self.DUT.Name if roc is None else name
        roc = self.Dut if roc is None else roc
        return self._draw_occupancy(roc, name, cluster, tel_coods, cut, show)

    def draw_time_occupancy(self, cut=None, fid=False, binning=None, show=True, show_vid=False):
        self.Bins.set_bin_size(binning)
        cut_string = self.Cut.generate_custom(exclude='fiducial' if not fid else [], cluster=False) if cut is None else TCut(cut)
        h = TH3F('h_to', 'to', *(self.Bins.get_time() + self.Bins.get_pixel()))
        self.Tree.Draw('row:col:{}>>h_to'.format(self.get_t_var()), cut_string, 'goff')
        format_histo(h, y_tit='col', z_tit='row')
        gx, gy = [self.make_tgrapherrors('g_to{}'.format(i), 'Mean {}'.format(i.title())) for i in ['x', 'y']]
        for ibin in xrange(h.GetNbinsX() - 1):
            h.GetXaxis().SetRange(ibin, ibin + 1)
            t = h.GetXaxis().GetBinCenter(ibin)
            p = h.Project3D('zy')
            px, py = p.ProfileX(), p.ProfileY()
            for g, p in zip([gx, gy], [px, py]):
                g.SetPoint(ibin, t, p.GetMean())
                g.SetPointError(ibin, 0, p.GetMeanError())
            if show_vid:
                self.draw_histo(p, draw_opt='colz') if not ibin else p.Draw('samecolz')
                c = gROOT.GetListOfCanvases()[-1]
                c.Update()
                c.Modified()
                sleep(.5)
        set_root_output(show)
        c = self.make_canvas('cto', 'Hit Map Evolution', divide=(2, 1), x=2)
        for i, g in enumerate([gx, gy], 1):
            p = c.cd(i)
            p.SetLeftMargin(.13)
            format_histo(g, x_tit='Time [hh:mm]', y_tit='Pixel Position', t_ax_off=0, y_off=1.8)
            g.Draw('ap')
        return h

    # endregion OCCUPANCY
    # ----------------------------------------

    # ----------------------------------------
    # region DISTRIBUTIONS
    def draw_adc_distribution(self, cut=None, show=True, col=None, pix=None):
        h = TH1I('h_adc', 'ADC Distribution {d}'.format(d=self.DUT.Name), *self.Bins.get_adc())
        cut_string = self.Cut(cut) + self.Cut.generate_masks(col=col, pixels=pix, exclude=False)()
        self.format_statbox(entries=True)
        self.Tree.Draw('adc>>h_adc', cut_string, 'goff')
        format_histo(h, x_tit='adc', y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor)
        self.save_histo(h, 'ADCDisto', show, lm=.13, logy=True)
        return h

    def draw_vcal_distribution(self, cut=None, mcol=None, mpix=None, vcal=True, show=True, redo=False):
        def f():
            h1 = TH1F('h_vcal', 'vcal Distribution {d}'.format(d=self.DUT.Name), *self.Bins.get_ph(vcal))
            cut_string = self.Cut(cut) + self.Cut.generate_masks(col=mcol, pixels=mpix, exclude=False)() + TCut('n_hits[{}] == 1'.format(self.Dut))
            n = self.Tree.Draw('col:row:adc', cut_string, 'goff')
            cols, rows, adcs = self.Run.get_root_vecs(n, 3, dtype=int)
            for i, (col, row, adc) in enumerate(zip(cols, rows, adcs)):
                self.Fit.SetParameters(*self.Parameters[self.Dut][col][row])
                h1.Fill(self.Fit.GetX(adc) * (self.Bins.VcalToEl if not vcal else 1))
            return h1
        h = do_pickle(self.make_simple_pickle_path(sub_dir='VCAL'), f, redo=redo)
        if show:
            self.format_statbox(all_stat=True)
            format_histo(h, x_tit='Pulse Height [{u}]'.format(u='vcal' if vcal else 'e'), y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor)
            self.save_histo(h, '{p}Disto'.format(p='Ph' if not vcal else 'Vcal'), show, lm=0.13)
        return h

    def draw_signal_distribution(self, cut=None, show=True, prnt=True, roc=None, vcal=False, redo=False, draw_thresh=False):
        roc = self.Dut if roc is None else roc
        cut_string = self.Cut.generate_custom(exclude='masks') if cut is None else TCut(cut)
        pickle_path = self.make_pickle_path('PulseHeight', run=self.RunNumber, suf='{}_{}_{}'.format(roc, cut_string.GetName(), int(vcal)))

        def f():
            set_root_output(False)
            h1 = TH1F('h_phd', 'Pulse Height Distribution - {d}'.format(d=self.DUT.Name), *self.Bins.get_ph(vcal))
            self.Tree.Draw('cluster_charge[{d}]{v}>>h_phd'.format(d=self.Dut, v='/{}'.format(self.Bins.VcalToEl) if vcal else ''), cut_string, 'goff')
            return h1

        h = do_pickle(pickle_path, f, redo=redo)
        x_range = [h.GetXaxis().GetXmin(), h.GetBinCenter(h.FindLastBinAbove(2)) * 1.2]
        self.format_statbox(all_stat=True, x=.92, w=.25)
        format_histo(h, x_tit='Pulse Height [{u}]'.format(u='vcal' if vcal else 'e'), y_tit='Number of Entries', y_off=1.8, fill_color=self.FillColor, x_range=x_range)
        self.draw_histo(h, show=show, lm=.13, rm=.06)
        if draw_thresh:
            self.draw_y_axis(1500, h.GetYaxis().GetXmin(), h.GetMaximum(), 'threshold #approx {}e  '.format(1500), off=.3, line=True, opt='-L')
        self.save_plots('PulseHeightDisto', prnt=prnt)
        return h

    def draw_cluster_disto(self, n=1, cut=None, redo=False, show=True, prnt=True):
        cut_string = TCut(str(n), 'n_clusters[{r}] == {n}'.format(r=self.Dut, n=n)) + self.Cut(cut)
        return self.draw_signal_distribution(cut=cut_string, redo=redo, show=show, prnt=prnt)

    def draw_cluster_ph_distos(self, cut=None, redo=False, show=True):
        hs = OrderedDict([(str(n), self.draw_cluster_disto(n, redo=redo, show=False, prnt=False)) for n in xrange(1, 4)])
        hs['>3'] = self.draw_signal_distribution(cut=TCut('g3', 'n_clusters[{r}] > 3'.format(r=self.Dut)) + self.Cut(cut), show=False, prnt=False, redo=redo)
        stack = THStack('h_cph', 'Pulser Height per Cluster Size')
        l1 = self.make_legend(y2=.88, nentries=5, x1=.7)
        for i, (name, h) in enumerate(hs.iteritems()):
            format_histo(h, color=self.Colors[i], normalise=True, fill_color=self.Colors[i], fill_style=3002, rebin=5)
            l1.AddEntry(h, name, 'fl')
            stack.Add(h)
        x_max = max(h.GetBinCenter(h.FindLastBinAbove(0)) for h in hs.itervalues())
        format_histo(stack, x_tit='Pulse Height [e]', y_tit='Number of Entries', y_off=1.4, draw_first=True, x_range=[self.Bins.MinPH, x_max])
        self.save_histo(stack, 'ClusterPulseHeight', show, lm=.13, leg=l1, draw_opt='nostackhist', gridy=True)
        return hs
    # endregion DISTRIBUTIONS
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def draw_pulse_height(self, cut=None, bin_size=30000, show=True, adc=False):
        """ Pulse height analysis vs event for a given cut. If no cut is provided it will take all. """
        h = TProfile('hpht', '{} -  {}'.format('ADC' if adc else 'Pulse Height', self.DUT.Name), *self.Bins.get_time(bin_size))
        self.format_statbox(only_fit=True, w=.3)
        self.Tree.Draw('{}:{} >> hpht'.format(self.get_ph_str() if not adc else 'adc', self.get_t_var()), self.Cut(cut), 'goff')
        self.draw_histo(h, show=show)
        fit_par = h.Fit('pol0', 'qs')
        format_histo(h, name='Fit Result', y_off=1.9, y_tit='Pulse Height [e]', x_tit='Time [hh:mm]', markersize=.6, t_ax_off=self.Run.StartTime)
        self.save_histo(h, '{}Time'.format('ADC' if adc else 'PulseHeight'), show, lm=.16, draw_opt='e1', save=show, canvas=get_last_canvas())
        return h, FitRes(fit_par)

    def draw_ph_pull(self, event_bin_width=None, fit=True, bin_width=100, save=True, show=True, adc=False):
        return self._draw_ph_pull(event_bin_width, fit, bin_width, bins=self.Bins.get_ph(adc=adc, bin_width=bin_width), show=show, save=save)

    def draw_adc_vs_event(self, cut=None, show=True):
        return self.draw_pulse_height(cut, show=show, adc=True)

    def draw_adc_pull(self, event_bin_width=None, fit=True, save=True, show=True):
        return self.draw_ph_pull(event_bin_width, fit, show=show, save=save, adc=True)

    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region 2D DISTRIBUTIONS
    def get_2d_coods(self, roc=None, global_=False, native=False, mm=True):
        fac = '*10' if mm else ''
        coods = 'cluster_ypos_tel[{r}-4]{f}:cluster_xpos_tel[{n}-4]{f}' if global_ else 'row:col' if native else 'dia_track_y_local[{r}]{f}:dia_track_x_local[{r}]{f}'
        return coods.format(r=self.Dut if roc is None else roc, f=fac)

    def draw_adc_map(self, show=True, cut=None, zero_ratio=False):
        cut_string = self.Cut(cut) + self.Cut.generate_masks(cluster=False)()
        set_root_output(False)
        h = TProfile2D('p_am', 'ADC Map', *self.Bins.get_pixel())
        self.Tree.Draw('adc{}:row:col >> p_am'.format('<0' if zero_ratio else ''), cut_string, 'goff')
        self.format_statbox(x=.81, entries=True)
        format_histo(h, x_tit='col', y_tit='row', z_tit='ADC', y_off=1.3, z_off=1.5, z_range=[0, 1] if zero_ratio else None)
        self.save_histo(h, 'ADCMap', show, rm=.17, lm=.13, draw_opt='colz')

    def draw_threshold_map(self, vcal=True, cols=None, show=True):
        h = TProfile2D('p_tm', 'Artificial Threshold Map', *self.Bins.get_pixel())
        for (col, row), thresh in self.get_thresholds(cols, vcal=vcal).iteritems():
            h.Fill(col, row, thresh if vcal else thresh / 1000.)
        format_histo(h, x_tit='col', y_tit='row', z_tit='Artificial Treshold [{u}]'.format(u='vcal' if vcal else 'ke'), y_off=1.3, z_off=1.7, stats=0)
        self.save_histo(h, 'ThresholdMap', show, rm=.17, lm=.13, draw_opt='colz')
        return h

    def draw_adc_fixed_vcal_map(self, roc=None, vcal=200, show=True):
        roc = self.Dut if roc is None else roc
        h = TProfile2D('p_pm', 'ADC Map for Vcal {v}'.format(v=vcal), *self.Bins.get_pixel())
        cols, rows = split(array(self.Cut.CutConfig['local_fiducial']), 2) if self.Cut.CutConfig['local_fiducial'] is not None else [0, self.Bins.NCols - 1], [0, self.Bins.NRows - 1]
        for col in xrange(cols[0], cols[1] + 1):
            for row in xrange(rows[0], rows[1] + 1):
                self.Fit.SetParameters(*self.Parameters[roc][col][row])
                h.Fill(col, row, self.Fit(vcal))
        format_histo(h, x_tit='col', y_tit='row', z_tit='Pulse Height [adc]', y_off=1.3, z_off=1.5, stats=0)
        self.save_histo(h, 'ADCMap{v}'.format(v=vcal), show, rm=.17, lm=.13, draw_opt='colz')

    def draw_sig_map_disto(self, res=None, cut=None, fid=True, x_range=None, redo=False, normalise=False, ret_value=False, show=True, save=True):
        return self._draw_sig_map_disto(res, cut, fid, x_range, redo, normalise, ret_value, ph_bins=self.Bins.get_ph(), show=show, save=save)
    # endregion 2D DISTRIBUTIONS
    # ----------------------------------------

    # ----------------------------------------
    # region CALIBRATION
    def draw_calibration_fit(self, col, row, show=True, roc=None):
        roc = self.Dut if roc is None else roc
        self.Fit.SetParameters(*self.Parameters[roc][col][row])
        format_histo(self.Fit, title='Calibration Fit for Pix {c} {r}'.format(c=col, r=row), x_tit='vcal', y_tit='adc', y_off=1.4, color=632, lw=2)
        gr = TGraph(len(self.Vcals[roc]), array(self.Vcals[roc], 'd'), array(self.Points[roc][col][row], 'd'))
        format_histo(gr, marker=20, name='gr_fp', title='Calibration Fit for Pix {c} {r}'.format(c=col, r=row), x_tit='vcal', y_tit='adc', y_off=1.4)
        self.draw_histo(gr, lm=.12, draw_opt='ap')
        self.Fit.Draw('same')
        self.save_plots('CalFit{c}{r}'.format(c=col, r=row), show=show)

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
            format_histo(gr, marker=20, x_tit='vcal', y_tit='adc', y_off=1.3, title='Calibration Fit for Pix {c} {r}'.format(c=col, r=row))
            self.draw_histo(gr, draw_opt='ap', show=show, lm=.12)
            return fit
    # endregion CALIBRATION
    # ----------------------------------------

    # ----------------------------------------
    # region EFFICIENCY
    def get_efficiency_cut(self, trig_phase=True):
        return self.Cut.generate_custom(include=['fiducial', 'rhit', 'tracks', 'chi2_x', 'chi2_y', 'aligned', 'event_range', 'beam_interruptions'] + (['trigger_phase'] if trig_phase else []))

    def get_hit_efficiency(self, roc=None, cut=None):
        cut_string = self.get_efficiency_cut() if cut is None else TCut(cut)
        roc = self.Dut if roc is None else roc
        n = self.Tree.Draw('n_hits[{r}]>0'.format(r=roc), cut_string, 'goff')
        values = self.Run.get_root_vec(n, dtype=bool)
        return calc_eff(values=values)

    def draw_eff_vs_chi2(self, step=10, show=True):
        g = self.make_tgrapherrors('gec2', 'Efficiency vs Chi2')
        for i, chi2 in enumerate(arange(step, 100 + step / 100., step)):
            self.Cut.set_chi2(chi2)
            eff = self.get_hit_efficiency()
            g.SetPoint(i, chi2, eff.n)
            g.SetPointError(i, 0, eff.s)
        format_histo(g, x_tit='#chi^{2} [%quantile]', y_tit='Efficiency [%]', y_off=1.3)
        self.draw_histo(g, draw_opt='ap', lm=.12, show=show)

    def draw_hit_efficiency(self, save=True, cut=None, vs_time=True, bin_width=5000, n=1e9, start=0, show=True):
        set_root_output(False)
        h = TProfile('h_he', 'Hit Efficiency {s}'.format(s=self.Dut), *self.Bins.get(bin_width, vs_time))
        cut_string = self.Cut.generate_custom(exclude=['masks']) if cut is None else TCut(cut)
        x_var = self.get_t_var() if vs_time else 'event_number'
        self.Tree.Draw('(n_hits[{r}]>0)*100:{x} >> h_he'.format(r=self.Dut, x=x_var), cut_string, 'goff', int(n), start)
        g = self.make_graph_from_profile(h)
        fit = fix_chi2(g, .01, show)
        format_histo(g, x_tit='Time [hh:mm]' if vs_time else 'Event Number', y_tit='Efficiency [%]', y_off=1.4, y_range=[-5, 115], stats=0,
                     t_ax_off=self.Run.StartTime if vs_time else None, markersize=1.2, draw_first=True)
        self.draw_histo(g, show=show, lm=.13, gridy=True, draw_opt='apz', bm=.2)
        self.draw_stats(fit, width=.35, y2=.35, names=['Efficiency'])
        self.draw_preliminary()
        self.save_plots('HitEfficiencyROC{n}'.format(n=self.Dut), save=save, show=show)
        return fit if fit.Parameter(0) is not None else 0

    def draw_efficiency_map(self, res=None, cut='all', show=True):
        cut_string = TCut(cut) + self.Cut.CutStrings['tracks']
        cut_string = self.Cut.generate_custom(exclude=['masks', 'fiducial']) if cut == 'all' else cut_string
        p = TProfile2D('p_em', 'Efficiency Map {d}'.format(d=self.DUT.Name), *self.Bins.get_global(res_fac=res, mm=True))
        self.Tree.Draw('(n_hits[{r}]>0)*100:{}:{}>>p_em'.format(r=self.Dut, *self.Cut.get_track_vars(self.Dut - 4, mm=True)), cut_string, 'goff')
        self.format_statbox(entries=True, x=.81)
        format_histo(p, x_tit='Track Position X [mm]', y_tit='Track Position Y [mm]', z_tit='Efficiency [%]', y_off=1.4, z_off=1.5)
        self.draw_histo(p, 'Efficiency Map', show, lm=.13, rm=.17, draw_opt='colz')
        self.draw_fid_cut(scale=10)
        self.draw_preliminary()
        self.save_plots('Efficiency Map')

    def get_fiducial_cell(self, n):
        x1, x2, y1, y2 = self.Cut.CutConfig['fiducial']
        nx = int(round((x2 - x1) / self.PX))
        # ny = int(round((y2 - y1) / .010))
        return round(x1 + self.PX * (n % nx), 4), round(y1 + self.PY * (n / nx), 4)

    def draw_cell_efficiency(self, cell=0, res=2, show=True):
        # x, y = self.get_fiducial_cell(cell)
        cut_string = self.Cut.generate_custom(exclude=['masks'])
        p = TProfile2D('pce', 'Efficiency for Fiducial Cell {}'.format(cell), res, 0, self.PX, res, 0, self.PY)
        n = self.Tree.Draw('(n_hits[{r}]>0)*100:dia_track_y_local[{r1}]:dia_track_x_local[{r1}]>>pce'.format(r=self.Dut, r1=self.Dut - 4), cut_string, 'goff')
        effs = [self.Tree.GetV1()[i] for i in xrange(n)]
        x_vals = [self.Tree.GetV3()[i] for i in xrange(n)]
        y_vals = [self.Tree.GetV2()[i] for i in xrange(n)]
        for e, x, y in zip(effs, x_vals, y_vals):
            p.Fill(x % self.PX, y % self.PY, e)
        self.format_statbox(entries=True, x=.81)
        format_histo(p, x_tit='Track x [cm]', y_tit='Track y [cm]', z_tit='Efficiency [%]', y_off=1.4, z_off=1.5)
        self.draw_histo(p, show=show, lm=.13, rm=.17, draw_opt='colz')

    def draw_efficiency_vs_trigphase(self, show=True):
        n = self.Tree.Draw('n_hits[{}]>0:trigger_phase[1]'.format(self.Dut), self.get_efficiency_cut(trig_phase=False), 'goff')
        hits, tp = self.Run.get_root_vecs(n, 2)
        y = [calc_eff(values=hits[where(tp == i)]) for i in arange(10)]
        g = self.make_tgrapherrors('getp', 'Efficiency per Trigger Phase', x=arange(10), y=y)
        format_histo(g, fill_color=self.FillColor, x_tit='Trigger Phase', y_tit='Efficiency [%]', y_off=1.4, x_range=[-1, 10])
        self.save_histo(g, 'EffVsTrigPhase', show, draw_opt='ba', lm=.13)
    # endregion EFFICIENCY
    # ----------------------------------------

    # ----------------------------------------
    # region TRIGGER PHASE
    def draw_trigger_phase(self, cut=None, show=True):
        return self._draw_trigger_phase(dut=True, cut=self.Cut.generate_custom(exclude='trigger_phase') if cut is None else cut, show=show)

    def draw_trigger_phase_time(self, bin_width=30000, cut=None, show=True):
        return self._draw_trigger_phase_time(dut=True, bin_width=bin_width, cut=cut, show=show)

    def draw_trigphase_offset(self, cut=None, show=True):
        h = TH1I('h_tp', 'Trigger Phase Offset', 18, -9.5, 8.5)
        self.Tree.Draw('trigger_phase[1] - trigger_phase[0]>>h_tp', self.Cut.generate_custom(exclude='trigger_phase') if cut is None else cut, 'goff')
        self.format_statbox(entries=True, y=0.88)
        format_histo(h, x_tit='Trigger Phase', y_tit='Number of Entries', y_off=1.8, fill_color=self.FillColor, ndivx=20)
        self.save_histo(h, 'TPOff', show, lm=.16)

    def draw_trigphase_off_time(self, bin_width=30000, cut=None, show=True):
        h = TProfile('htpot', 'Trigger Phase vs Time', *self.Bins.get(bin_width, vs_time=True))
        self.Tree.Draw('(trigger_phase[1] - trigger_phase[0]):{}>>htpot'.format(self.get_t_var()), self.Cut.generate_custom(exclude='trigger_phase') if cut is None else cut, 'goff')
        self.format_statbox(entries=True, y=0.88)
        format_histo(h, x_tit='Time [hh:mm]', y_tit='Trigger Phase', y_off=1.8, fill_color=self.FillColor, t_ax_off=self.Run.StartTime)
        self.save_histo(h, 'TPOffTime', show, lm=.16)
    # endregion TRIGGER PHASE
    # ----------------------------------------

    # ----------------------------------------
    # region ARTIFICIAL THRESHOLD
    def find_landau(self, aver=10, m1=2500, m2=5000, s1=500, s2=1600):
        seed = self.draw_signal_distribution(show=False)
        h = deepcopy(seed)
        m_range = range(m1, m2 + 1, 100)
        s_range = range(s1, s2 + 1, 50)
        p = TProfile2D('g_fl', 'Find Landau', len(m_range) - 1, m_range[0], m_range[-1], len(s_range) - 1, s_range[0], s_range[-1])
        self.PBar.start(len(m_range) * len(s_range) * aver)
        i = 0
        r_min, r_max = .5, 1.5
        for _ in xrange(aver):

            for m in m_range:
                for s in s_range:
                    i += 1
                    self.PBar.update(i)
                    if r_min < m / 4. / s < r_max:
                        diff = self.model_landau(seed, h, m, s, show=False, thresh=True)
                        p.Fill(m, s, diff)
        self.PBar.finish()
        self.format_statbox(entries=True, x=.82)
        format_histo(p, x_tit='MPV [e]', y_tit='Sigma [e]', z_tit='#chi^{2} to Seed Function', y_off=1.7, z_off=1.3)
        self.draw_histo(p, draw_opt='colz', lm=.13, rm=0.16)
        self.draw_ms_ratios(r_min, r_max, m1, m2, s1, s2)
        self.save_plots('FindLandau')
        self.find_working_point(p)

    def draw_ms_ratios(self, r_min, r_max, m1, m2, s1, s2, step=.1):
        ratios = arange(r_min, r_max + step, step)
        off = .01
        for i, ratio in enumerate(ratios):
            cut = TCutG('ms{n}'.format(n=i), 2, array([0, 20000], 'd'), array([0, 20000 / ratio / 4]))
            x_pos = m1 if m1 / ratio / 4 > s1 else 4 * ratio * s1
            y_pos = m1 / ratio / 4 if m1 / ratio / 4 > s1 else s1
            self.draw_tlatex(x_pos + off * (m2 - m1), y_pos + off * (s2 - s1), text='{0:3.1f}'.format(ratio), size=.02, align=11)
            cut.Draw('same')
            self.Objects.append(cut)

    @staticmethod
    def find_working_point(h):
        ps = [h.ProfileY(), h.ProfileX()]
        fits = [TF1('f{n}'.format(n=i), 'pol2', 0, 10000) for i in xrange(2)]
        for fit, p in zip(fits, ps):
            p.Fit(fit, 'qs0')
        mins = [fit.GetMinimumX() for fit in fits]
        print(mins, mins[1] / mins[0] / 4)

    def model_landau(self, seed=None, h=None, m=10000, s=1000, show=True, thresh=False):
        seed = self.draw_signal_distribution(show=False) if seed is None else seed
        h = deepcopy(seed) if h is None else h
        h.SetName('h_ml')
        n = seed.GetEntries()
        h.Reset()
        thresholds = [47 * i for i in [150, 160, 170]]
        for _ in xrange(int(n)):
            v = gRandom.Landau(m, s)
            threshold = thresholds[int(gRandom.Rndm() * len(thresholds))]
            h.Fill(v if v > threshold else 0 if thresh else v)
            # h.Fill(v if v > 11421 else 0 if thresh else v)
        diff = mean([(h.GetBinContent(i) - seed.GetBinContent(i)) ** 2 for i in xrange(h.GetNbinsX()) if i is not h.FindBin(0)])
        seed.SetFillColor(2)
        h.SetFillColor(self.FillColor)
        seed.SetFillStyle(3017)
        h.SetFillStyle(3018)
        l1 = self.make_legend(y2=.76)
        l1.AddEntry(h, 'Simulation', 'f')
        l1.AddEntry(seed, 'Original', 'f')
        self.draw_histo(h, show=show, leg=l1)
        self.draw_histo(seed, show=show, draw_opt='same', canvas=gROOT.GetListOfCanvases()[-1])
        seed.Draw('same')
        return diff

    def landau_vid(self, save=False, mpv=5000, sigma=820):
        h = self.draw_signal_distribution()
        h.GetYaxis().SetRangeUser(0, 2500)
        zero_bin = h.FindBin(0)
        zero_entries = int(h.GetBinContent(zero_bin))
        entries = int(h.GetEntries())
        print(entries)
        c = gROOT.GetListOfCanvases()[-1]
        thresholds = self.get_thresholds(vcal=False)
        for i in xrange(entries):
            h.SetBinContent(zero_bin, zero_entries)
            v = gRandom.Landau(mpv, sigma)
            threshold = thresholds[int(rand() * self.Bins.NCols), int(rand() * self.Bins.NRows)]
            if v < threshold:
                h.Fill(v)
                zero_entries -= 1
            if i % 100 == 0:
                c.Update()
                c.Modified()
            if i % 100 == 0 and save:
                self.save_canvas(c, name='l{i:04d}'.format(i=i), show=False, print_names=False)
    # endregion ARTIFICIAL THRESHOLD
    # ----------------------------------------

    # ----------------------------------------
    # region ALIGNMENT
    def draw_alignment(self, plane1=2, plane2=None, mode='y', binning=5000, chi2=1, show=True, vs_time=True, redo=False):
        plane2 = self.Dut if plane2 is None else plane2
        picklepath = self.make_pickle_path('Alignment', run=self.RunNumber, suf='{m}_{p1}{p2}_{b}_{t}'.format(m=mode, p1=plane1, p2=plane2, b=binning, t='Time' if vs_time else 'EvtNr'))

        def func():
            start = self.info('Checking for alignment between plane {p1} and {p2} ... '.format(p1=plane1, p2=plane2), next_line=False)
            h = TH3F('hpa', 'pa', *(self.Bins.get(binning, vs_time=vs_time) + 2 * self.Bins.get_global_cood(mode, res_fac=sqrt(12))))
            y, z_ = ('cluster_{m}pos_tel[{r}]'.format(m=mode, r=plane) for plane in [plane1, plane2])
            cut_string = TCut('n_clusters[{p1}]==1 && n_clusters[{p2}]==1'.format(p1=plane1, p2=plane2)) + self.Cut.generate_chi2(mode, chi2)()
            self.Tree.Draw('{z}:{y}:{x}>>hpa'.format(z=z_, y=y, x=self.get_t_var() if vs_time else 'Entry$'), cut_string, 'goff')
            g = self.make_tgrapherrors('g_pa', 'Plane Correlation {p1} {p2}'.format(p1=plane1, p2=plane2), marker_size=.5)
            for ibin in xrange(h.GetNbinsX() - 1):
                h.GetXaxis().SetRange(ibin, ibin + 1)
                p = h.Project3D('yz')
                g.SetPoint(ibin, h.GetXaxis().GetBinCenter(ibin), p.GetCorrelationFactor())
            format_histo(g, x_tit='Time [hh::mm]' if vs_time else 'Event Number', y_tit='Correlation Factor', y_off=1.5, y_range=[0, 1], t_ax_off=self.Run.StartTime if vs_time else None)
            self.add_to_info(start)
            return g

        gr = do_pickle(picklepath, func, redo=redo)
        self.save_histo(gr, 'PixelAligment', show, draw_opt='alp', lm=.13, prnt=show)
        return gr

    def check_alignment(self):
        pickle_path = self.make_pickle_path('Alignment', run=self.RunNumber)

        def f():
            g = self.draw_alignment(show=False)
            x, y = get_graph_vecs(g)
            return mean_sigma(y)

        m, s = do_pickle(pickle_path, f)
        if m < .4:
            log_warning('Planes are not correlated!')
        elif s > .05:
            log_warning('Large fluctuations in correlation!')
        return m > .4

    def draw_correlation(self, plane1=2, plane2=None, mode='y', chi2=None, res=.7, start=0, evts=int(1e10), cut=None, show=True):
        old_chi2 = self.Cut.CutConfig['chi2_x']
        self.Cut.set_chi2(old_chi2 if chi2 is None else chi2)
        plane2 = self.Dut if plane2 is None else plane2
        h = TH2F('h_pc', 'Plane Correlation', *self.Bins.get_global(res_fac=res))
        cut = self.Cut(cut)
        self.Tree.Draw('cluster_{m}pos_tel[{p1}]:cluster_{m}pos_tel[{p2}]>>h_pc'.format(m=mode, p1=plane1, p2=plane2), cut, 'goff', evts, start)
        self.info('Correlation Factor: {f:4.3f}'.format(f=h.GetCorrelationFactor()))
        format_histo(h, x_tit='{m} Plane {p}'.format(p=plane1, m=mode), y_tit='{m} Plane {p}'.format(p=plane2, m=mode), y_off=1.5, stats=0, z_tit='Number of Entries', z_off=1.5)
        self.save_histo(h, 'PlaneCorrelation{m}{p1}{p2}'.format(m=mode.title(), p1=plane1, p2=plane2), show, lm=.13, draw_opt='colz', rm=.17)
        self.Cut.set_chi2(old_chi2)

    def draw_event_offsets(self, evnts=1000, start=0, pnts=None, rnge=1, show=True):
        pnts = int(ceil(self.Run.NEntries / evnts)) if pnts is None else pnts
        n_graphs = 2 * rnge + 1
        graphs = [self.make_tgrapherrors('g_to{i}'.format(i=i), '', color=self.get_color(), marker_size=.5) for i in xrange(n_graphs)]
        l1 = self.make_legend(x1=.8, y2=.5, nentries=n_graphs - 1)
        self.PBar.start(pnts * evnts - start)
        p1 = self.Dut
        p2 = 2
        for k in xrange(pnts):
            x, y = OrderedDict(), OrderedDict()
            n = self.Tree.Draw('plane:row:event_number', 'n_hits[2]==1 || n_hits[{p1}]==1'.format(p1=p1), 'goff', evnts, start + k * evnts)
            planes = [int(self.Tree.GetV1()[i]) for i in xrange(n)]
            rows = [int(self.Tree.GetV2()[i]) for i in xrange(n)]
            nrs = Counter([int(self.Tree.GetV3()[i]) for i in xrange(n)])
            for ev, size in sorted(nrs.iteritems()):
                self.PBar.update(ev + 1 - start)
                plane = [planes.pop(0) for _ in xrange(size)]
                row = [rows.pop(0) for _ in xrange(size)]
                if plane.count(p1) == 1:
                    x[ev] = row[plane.index(p1)]
                if plane.count(p2) == 1:
                    y[ev] = row[plane.index(p2)]
            xts = [[] for _ in xrange(n_graphs)]
            yts = [[] for _ in xrange(n_graphs)]
            for i in xrange(-rnge, rnge + 1):
                j = i + rnge
                for ev, row in x.iteritems():
                    if ev + i in y:
                        xts[j].append(row)
                        yts[j].append(y[ev + i])
                corr = corrcoef(xts[j], yts[j])[0][1]
                graphs[j].SetPoint(k, k * evnts + start, corr)
        self.PBar.finish()
        mg = TMultiGraph('m_eo', 'Correlations')
        for i, gr in enumerate(graphs):
            l1.AddEntry(gr, 'offset: {i}'.format(i=i - rnge), 'pl')
            mg.Add(gr, 'pl')
        format_histo(mg, x_tit='Event Number', y_tit='Correlation Factor', y_off=1.5, y_range=[-.2, 1], draw_first=True)
        self.save_histo(mg, 'EventOffsets', show, draw_opt='ap', leg=l1, lm=.13)
        self.reset_colors()
    # endregion ALIGNMENT
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_dut_hits(self, dut1=4, dut2=5, cut=None, show=True):
        h = TH2F('h_hd', 'Number of Hits', 40, 0, 40, 40, 0, 40)
        self.Tree.Draw('n_hits[{}]:n_hits[{}]>>h_hd'.format(dut2, dut1), self.Cut(cut), 'goff')
        self.format_statbox(entries=True, x=.81)
        x_tit, y_tit = ('Number of Hits in {}'.format(self.Run.DUT.Names[dut - 4]) for dut in [dut1, dut2])
        format_histo(h, x_tit=x_tit, y_tit=y_tit, y_off=1.3, z_tit='Number of Entries', z_off=1.4, x_range=[0, h.FindLastBinAbove(0, 1) + 1], y_range=[0, h.FindLastBinAbove(0, 2) + 1])
        self.save_histo(h, 'HitsDiaSil', show, draw_opt='colz', rm=0.17, lm=.13, logz=True)

    def draw_hit_pie(self, d2=5):
        zero, dia, dut2, both = (self.Tree.GetEntries(s.format(d2)) for s in ['n_hits[4]==0&&n_hits[{}]==0', 'n_hits[4]>0&&n_hits[{}]==0', 'n_hits[4]==0&&n_hits[{}]>0', 'n_hits[4]>0&&n_hits[{}]>0'])
        names = ['No Hits', '{} Hit'.format(self.DUT.Name), '{} Hit'.format(self.Run.DUT.Names[d2 - 4]), 'Both Hits']
        pie = TPie('pie', 'Hit Contributions', 4, array([zero, dia, dut2, both], 'f'), array(self.get_colors(4), 'i'))
        for i, label in enumerate(names):
            pie.SetEntryRadiusOffset(i, .05)
            pie.SetEntryLabel(i, label.title())
        format_pie(pie, h=.04, r=.2, text_size=.025, angle3d=70, angle_off=250)
        self.draw_histo(pie, draw_opt='3drsc')
        self.reset_colors()

    def draw_cluster_size(self, cut='', show=True):
        return self._draw_cluster_size(self.Dut, self.DUT.Name, cut, show)

    def draw_residuals(self, mode=None, cut=None, show=True, x_range=None):
        return self._draw_residuals(self.Dut, mode=mode, cut=cut, show=show, x_range=x_range)

    def draw_slope(self, show=True):
        with open(join(self.Dir, 'data', 'vcalCalibration.txt')) as f:
            d = load(f)
            h = TH1F('ho', 'Vcal Calibration Slopes', 30, 35, 65)
            for value in d['slope']:
                h.Fill(value)
            self.format_statbox(fit=True, all_stat=True, w=.3)
            format_histo(h, x_tit='Slope [e]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor)
            self.draw_histo(h, show=show)
            h.Fit('gaus', 'q')

    def draw_offset(self, show=True):
        with open(join(self.Dir, 'data', 'vcalCalibration.txt')) as f:
            d = load(f)
            h = TH1F('ho', 'Vcal Calibration Offsets', 15, -200, 400)
            for value in d['offset']:
                h.Fill(value)
            self.format_statbox(fit=True, all_stat=True, w=.3)
            format_histo(h, x_tit='Offset [e]', y_tit='Number of Entries', y_off=1.2, fill_color=self.FillColor)
            self.format_statbox(fit=True, all_stat=True)
            self.draw_histo(h, show=show)
            h.Fit('gaus', 'q')
    # endregion DRAW
    # ----------------------------------------

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


if __name__ == '__main__':

    pargs = init_argparser(run=139, tc='201810', dut=1, has_verbose=True, tree=True)
    z = PixAnalysis(pargs.run, pargs.dut, pargs.testcampaign, pargs.tree, pargs.verbose)
