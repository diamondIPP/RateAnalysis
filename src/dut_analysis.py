#!/usr/bin/env python
# --------------------------------------------------------
#       parent class for the analysis of a single device under test
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TProfile2D, TH2I
from numpy import vectorize
from numpy.random import rand
from uncertainties import umath

from src.currents import Currents
from src.cut import Cut
from src.run import Run
from src.tracks import Tracks
from src.telescope import Telescope
from helpers.save_plots import *
from src.analysis import Analysis
from helpers.fit import Langau


class DUTAnalysis(Analysis):
    def __init__(self, run_number, diamond_nr=1, testcampaign=None, tree=None, t_vec=None, verbose=False, prnt=True):

        self.Run = self.init_run(run_number, testcampaign, tree, t_vec, verbose)
        self.DUT = self.Run.DUTs[diamond_nr - 1]
        super(DUTAnalysis, self).__init__(testcampaign, results_dir=join(self.DUT.Name, str(self.Run.Number)), verbose=verbose)

        self.Tree = self.Run.Tree
        self.StartTime = self.Run.StartTime if self.Tree.Hash() else time_stamp(self.Run.LogStart)

        self.print_start(run_number, prnt, dut=self.DUT.Name)
        self.update_config()

        # Sub-Analyses
        self.Currents = Currents(self)

        if self.Tree.Hash():
            self.Cut = Cut(self)
            self.NRocs = self.Run.NPlanes
            self.StartEvent = self.Cut.get_min_event()
            self.EndEvent = self.Cut.get_max_event()
            self.Bins = Bins(self.Run, cut=self.Cut)
            self.Tracks = Tracks(self)
            self.Tel = Telescope(self)

    @staticmethod
    def init_run(run_number, testcampaign, tree, t_vec, verbose):
        return Run(run_number, testcampaign, tree, t_vec, verbose)

    def update_config(self):
        pass

    @staticmethod
    def get_info_header():
        return ['Run', 'Type', 'Diamond', 'Flux [kHz/cm2]', 'HV [V]']

    def show_information(self, header=True, prnt=True):
        rows = [[self.Run.Number, self.Run.Info['runtype'], self.DUT.Name, '{:14.1f}'.format(self.Run.Flux.n), '{:+6.0f}'.format(self.DUT.Bias)]]
        return print_table(rows, self.get_info_header() if header else None, prnt=prnt)

    # ----------------------------------------
    # region GET
    def has_branch(self, branch):
        return self.Run.has_branch(branch)

    def get_time(self):
        return self.Run.get_time()

    def get_t_var(self):
        return 'time / 1000.' if self.Run.TimeOffset is None else '(time - {}) / 1000.'.format(self.Run.TimeOffset)

    def get_t_off(self, rel_time):
        return self.Run.StartTime if rel_time else 0

    def get_root_vec(self, n=0, ind=0, dtype=None, var=None, cut='', nentries=None, firstentry=0):
        return self.Run.get_root_vec(n, ind, dtype, var, cut, nentries, firstentry)

    def get_events(self, cut=None, redo=False):
        cut = self.Cut(cut)
        return do_hdf5(self.make_simple_hdf5_path('Events', suf=cut.GetName()), self.get_root_vec, redo, dtype='i4', var='Entry$', cut=cut)

    def get_event_at_time(self, seconds, rel=False):
        return self.Run.get_event_at_time(seconds, rel)

    def get_n_entries(self, cut=None):
        return self.Tree.GetEntries(self.Cut(cut).GetTitle())

    def get_current(self):
        return self.Currents.get_current()

    def get_irradiation(self):
        return self.DUT.get_irradiation(self.TCString)

    def get_attenuator(self):
        return False

    def get_ph_str(self):
        return ''

    def get_pulse_height(self, *args, **kwargs):
        pass

    def get_sm_data(self, cut=None, fid=False):
        """ :return: signal map data as numpy array [[x], [y], [ph]] with units [[mm], [mm], [mV]]
            :param cut: applies all cuts if None is provided.
            :param fid: return only values within the fiducial region set in the AnalysisConfig.ini"""
        y, x = self.Cut.get_track_vars(self.DUT.Number - 1, mm=True)
        cut = self.Cut.generate_custom(exclude=['fiducial'], prnt=False) if not fid and cut is None else self.Cut(cut)
        n = self.Tree.Draw('{x}:{y}:{z}'.format(z=self.get_ph_str(), x=x, y=y), cut, 'goff')  # *10 to get values in mm
        return self.Run.get_root_vecs(n, 3)

    def get_uniformity(self, use_fwc=True, redo=False):
        return do_pickle(self.make_simple_pickle_path('Uniformity', int(use_fwc), 'Signal'), self.draw_uniformity, redo=redo, show=False, use_fwc=use_fwc)

    def get_track_length_var(self):
        dx2, dy2 = ['TMath::Power(TMath::Tan(TMath::DegToRad() * {}_{}), 2)'.format('slope' if self.Run.has_branch('slope_x') else 'angle', direction) for direction in ['x', 'y']]
        return '{} * TMath::Sqrt({} + {} + 1)'.format(self.DUT.Thickness, dx2, dy2)

    def get_flux(self, corr=True, rel_error=0, show=False):
        return (self.get_flux(show=show) if self.Tree and self.has_branch('rate') else self.Run.get_flux(rel_error)) * self.get_flux_correction(corr)

    def get_flux_correction(self, correct=True):
        def f():
            if not self.Cut.has('fiducial'):
                return 1
            n1, n2 = self.get_n_entries(self.Cut.generate_custom(include=['tracks', 'fiducial'], prnt=False)), self.get_n_entries(self.Cut.get('tracks'))
            a1, a2 = self.Cut.get_fiducial_area(), min(self.Run.get_unmasked_area().values())
            return n1 / a1 * a2 / n2
        return do_pickle(self.make_pickle_path('Flux', 'Corr', self.Run.Number, self.DUT.Number), f) if correct else 1

    def get_additional_peak_height(self):
        pass

    def get_peak_flux(self):
        pass

    def get_ph_values(self, *args, **kwargs):
        """ :returns: all pulse height values for a given cut. [numpy.ndarray] """

    def get_signal_name(self, *args, **kwargs):
        """ :returns: the pulse height variable in the tree. [str] """

    def generate_signal_name(self, *args, **kwargs):
        """ :returns: the pulse height variable in the tree + corrections. [str] """
    # endregion GET
    # ----------------------------------------
    def draw_chi2(self, *args, **kwargs):
        return self.Tracks.draw_chi2(*args, **kwargs)

    def draw_chi2s(self, *args, **kwargs):
        return self.Tracks.draw_chi2s(*args, **kwargs)

    def draw_angle(self, *args, **kwargs):
        return self.Tracks.draw_angle(*args, **kwargs)

    def draw_angles(self, *args, **kwargs):
        return self.Tracks.draw_angles(*args, **kwargs)

    def draw_occupancies(self, *args, **kwargs):
        return self.Tel.draw_occupancies(*args, **kwargs)

    def get_mean_angle(self, mode):
        return self.Tracks.get_mean_angle(mode)

    def draw_current(self, *args, **kwargs):
        return self.Currents.draw(*args, **kwargs)

    def draw_flux(self, *args, **kwargs):
        return self.Tel.draw_flux(*args, **kwargs)

    def draw_pulse_height(self, *args, **kwargs):
        return TProfile()

    def draw_ph_pull(self, *args, **kwargs):
        return self._draw_ph_pull(*args, **kwargs)

    def _draw_ph_pull(self, event_bin_width=None, fit=True, bin_width=.5, binning=None, show=True, save=True):
        p = self.draw_pulse_height(event_bin_width, show=False, save=False)[0]
        format_statbox(all_stat=True, fit=fit)
        h = get_pull(p, 'Signal Bin{0} Distribution'.format(Bins.Size), binning=self.Bins.get_pad_ph(bin_width=bin_width) if binning is None else binning, fit=fit)
        format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.5, fill_color=Draw.FillColor, draw_first=True)
        self.Draw(h, 'SignalBin{0}Disto'.format(Bins.Size), save=save, lm=.12, show=show)
        return h

    def draw_signal_distribution(self, *args, **kwargs):
        pass

    def draw_fid_cut(self, scale=10):
        self.Cut.draw_fid_cut(scale)

    def draw_track_length(self, show=True):
        h = TH1F('htd', 'Track Distance in Diamond', 200, self.DUT.Thickness, self.DUT.Thickness + 1)
        self.Tree.Draw('{}>>htd'.format(self.get_track_length_var()), 'n_tracks', 'goff')
        format_histo(h, x_tit='Distance [#mum]', y_tit='Entries', y_off=2, lw=2, stats=0, fill_color=Draw.FillColor, ndivx=405)
        self.Draw(h, show, lm=.16)
        return h

    def draw_signal_vs_angle(self, mode='x', bin_size=.1, show=True):
        p = TProfile('psa{}'.format(mode), 'Pulse Height vs. Angle in {}'.format(mode.title()), *self.Bins.get_angle(bin_size))
        self.Tree.Draw('{}:angle_{m}>>psa{m}'.format(self.generate_signal_name(), m=mode), self.Cut(), 'goff')
        format_statbox(entries=True)
        format_histo(p, x_tit='Track Angle [deg]', y_tit='Pulse Height [mV]', y_off=1.5)
        self.Draw(p, show=show, lm=.12)

    # ----------------------------------------
    # region SIGNAL MAP
    def draw_signal_map(self, res=None, cut=None, fid=False, hitmap=False, redo=False, binning=None, z_range=None, size=None, show=True, prnt=True):

        cut = self.Cut.generate_custom(exclude=['fiducial'], prnt=prnt) if not fid and cut is None else self.Cut(cut)
        suf = '{c}_{ch}_{res}'.format(c=cut.GetName(), ch=self.Cut.CutConfig['chi2_x'], res=res if binning is None else '{}x{}'.format(binning[0], binning[2]))
        pickle_path = self.make_pickle_path('SignalMaps', 'Hit' if hitmap else 'Signal', run=self.Run.Number, ch=self.DUT.Number, suf=suf)

        def func():
            set_root_output(0)
            name = 'h_hm' if hitmap else 'h_sm'
            atts = [name, 'Track Hit Map' if hitmap else 'Signal Map'] + (self.Bins.get_global(res, mm=True) if binning is None else binning)
            h1 = TH2I(*atts) if hitmap else TProfile2D(*atts)
            self.info('drawing {mode}map of {dia} for Run {run}...'.format(dia=self.DUT.Name, run=self.Run.Number, mode='hit' if hitmap else 'signal '), prnt=prnt)
            y, x = self.Cut.get_track_vars(self.DUT.Number - 1, mm=True)
            self.Tree.Draw('{z}{y}:{x}>>{h}'.format(z=self.get_ph_str() + ':' if not hitmap else '', x=x, y=y, h=name), cut, 'goff')
            set_2d_ranges(h1, *([3, 3] if size is None else size))
            adapt_z_range(h1) if not hitmap else do_nothing()
            return h1

        set_palette(pal=1 if hitmap else 53)
        format_statbox(entries=True, x=0.82)
        h = do_pickle(pickle_path, func, redo=redo)
        z_tit = 'Number of Entries' if hitmap else 'Pulse Height [mV]'
        format_histo(h, x_tit='Track Position X [mm]', y_tit='Track Position Y [mm]', y_off=1.4, z_off=1.5, z_tit=z_tit, ncont=50, ndivy=510, ndivx=510, z_range=z_range)
        self.Draw(h, show=show, lm=.12, rm=.16, draw_opt='colzsame')
        self.draw_fid_cut(scale=10)
        # self.draw_detector_size(scale=10)
        self.Draw.save_plots('HitMap' if hitmap else 'SignalMap2D', prnt=prnt)
        return h

    def draw_hitmap(self, res=None, cut=None, fid=False, redo=False, z_range=None, size=None, show=True, prnt=True):
        cut = self.Cut.get('tracks') if cut is None else self.Cut(cut)
        return self.draw_signal_map(res, cut, fid, hitmap=True, redo=redo, binning=None, z_range=z_range, size=size, show=show, prnt=prnt)

    def split_signal_map(self, m=2, n=2, grid=True, redo=False, show=True):
        fid_cut = array(self.Cut.CutConfig['fiducial']) * 10
        if not fid_cut.size:
            critical('fiducial cut not defined for {}'.format(self.DUT.Name))
        x_bins = linspace(fid_cut[0], fid_cut[1], m + 1)
        y_bins = linspace(fid_cut[2], fid_cut[3], n + 1)
        binning = [m, x_bins, n, y_bins]
        h = self.draw_signal_map(binning=binning, show=False, fid=True, redo=redo)
        format_histo(h, x_range=[fid_cut[0], fid_cut[1] - .01], y_range=[fid_cut[2], fid_cut[3] - .01], name='hssm', stats=0)
        self.Draw(h, show=show, lm=.12, rm=.16, draw_opt='colzsame')
        Draw.grid(x_bins, y_bins, width=2) if grid else do_nothing()
        self.Draw.save_plots('SplitSigMap')
        return h, x_bins, y_bins

    def draw_beam_profile(self, mode='x', fit=True, fit_range=.8, res=.7, show=True, prnt=True):
        h = self.draw_hitmap(res, show=False, prnt=prnt)
        p = h.ProjectionX() if mode.lower() == 'x' else h.ProjectionY()
        format_statbox(all_stat=True)
        format_histo(p, title='Profile {}'.format(mode.title()), name='pbp{}'.format(self.Run.Number), y_off=1.3, y_tit='Number of Hits', fill_color=Draw.FillColor)
        self.Draw(p, lm=.13, show=show)
        if fit:
            fit = self.fit_beam_profile(p, fit_range)
        self.Draw.save_plots('BeamProfile{}'.format(mode.title()), prnt=prnt)
        return FitRes(fit) if fit else p

    @staticmethod
    def fit_beam_profile(p, fit_range):
        x_peak = p.GetBinCenter(p.GetMaximumBin())
        x_min, x_max = [p.GetBinCenter(i) for i in [p.FindFirstBinAbove(0), p.FindLastBinAbove(0)]]
        return p.Fit('gaus', 'qs', '', x_peak - (x_peak - x_min) * fit_range, x_peak + (x_max - x_peak) * fit_range)

    def draw_sig_map_disto(self, *args, **kwargs):
        return self._draw_sig_map_disto(*args, **kwargs)

    def _draw_sig_map_disto(self, res=None, cut=None, fid=True, x_range=None, redo=False, normalise=False, ret_value=False, ph_bins=None, show=True, save=True):
        source = self.draw_signal_map(res, cut, fid, redo=redo, show=False)
        h = TH1F('h_smd', 'Signal Map Distribution', *([400, 0, 4] if normalise else self.Bins.get_pad_ph(bin_width=.2) if ph_bins is None else ph_bins))
        normalisation = 1 if not normalise else self.get_pulse_height()
        values = get_2d_hist_vec(source) / normalisation
        [h.Fill(v.n) for v in values]
        format_statbox(all_stat=True)
        format_histo(h, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=2, fill_color=Draw.FillColor, x_range=choose(x_range, ax_range(5, 5, .3, .3, h)))
        self.Draw(h, 'SignalMapDistribution', lm=.15, show=show, save=save)
        return mean_sigma(values) if ret_value else h

    def get_sm_std(self, res=None, redo=False):
        return self.draw_sig_map_disto(show=False, res=res, redo=redo, normalise=True, ret_value=True, save=False)[1]

    def draw_sm_profile(self, mode='x', factor=1.5, cut=None, fid=False, hitmap=False, redo=False, show=True):
        s = self.draw_signal_map(factor, cut, fid, hitmap=hitmap, redo=redo, show=False)
        g = Draw.make_tgrapherrors('g_smp', 'Signal Map Profile')
        values = [[] for _ in range(s.GetNbinsX() if mode == 'x' else s.GetNbinsY())]
        for xbin in range(s.GetNbinsX()):
            for ybin in range(s.GetNbinsY()):
                value = s.GetBinContent(xbin, ybin)
                if value:
                    values[(xbin if mode == 'x' else ybin)].append(value)
        for i, lst in enumerate(values):
            m, sigma = mean_sigma(lst) if lst else (0., 0.)
            xval = s.GetXaxis().GetBinCenter(i) if mode == 'x' else s.GetYaxis().GetBinCenter(i)
            g.SetPoint(i, xval, m)
            g.SetPointError(i, 0, sigma)

        format_histo(g, x_tit='Track in {m} [cm]'.format(m=mode), y_tit='Pulse Height [au]', y_off=1.5, ndivx=515)
        self.Draw(g, 'SignalMapProfile', draw_opt='ap', lm=.14, show=show, gridx=True)

    def draw_error_signal_map(self, show=True):
        h = self.draw_signal_map(show=False, fid=True).ProjectionXY('', 'c=e')
        format_histo(h, name='hsme', title='Signal Map Errors', y_off=1.6)
        self.Draw(h, 'SignalMapErrors', lm=.12, rm=.11, show=show, draw_opt='colz')
        return h

    def get_signal_spread(self, min_percent=5, max_percent=99, prnt=True):
        """ Calculates the relative spread of mean signal response from the 2D signal response map. """
        pickle_path = self.make_pickle_path('SignalMaps', 'Spread', self.Run.Number, self.DUT.Number)

        def f():
            h = self.draw_sig_map_disto(show=False)
            q = array([min_percent, max_percent]) / 100.
            y = zeros(2, 'd')
            mean_error = mean([v.n for v in get_2d_hist_vec(self.draw_error_signal_map(show=False))])
            h.GetQuantiles(2, y, q)
            return make_ufloat([y[1], mean_error]) / make_ufloat([y[0], mean_error]) - 1
        ratio = do_pickle(pickle_path, f)
        self.info('Relative Signal Spread is: {:2.2f} %'.format(ratio * 100), prnt=prnt)
        return ratio
    # endregion SIGNAL MAP
    # ----------------------------------------

    def draw_uniformity(self, histo=None, use_fwc=False, show=True):
        noise = self.Pedestal.get_fwhm() if histo is None and hasattr(self, 'Pedestal') else 0
        h = self.draw_signal_distribution(show=show, normalise=True, sumw2=True) if histo is None else histo
        (low, high), m = get_fwhm(h, ret_edges=True), get_fw_center(h) if use_fwc else ufloat(h.GetMean(), h.GetMeanError())
        Draw.vertical_line(m.n, 0, 1e7, style=7, w=2)
        fwhm, half_max = high - low, h.GetMaximum() / 2
        Draw.tlatex(m.n + 5, .1 * half_max, 'FWC' if use_fwc else 'Mean', align=10)
        Draw.arrow(low.n, m.n, half_max, half_max, col=2, width=3, opt='<', size=.02)
        Draw.arrow(high.n, m.n, half_max, half_max, col=2, width=3, opt='<', size=.02)
        Draw.tlatex(high.n + 5, half_max, 'FWHM', align=12, color=2)
        fwhm = umath.sqrt(fwhm ** 2 - noise ** 2)
        value = fwhm / m
        legend = Draw.make_legend(w=.3, y2=.768, nentries=1, margin=.1, cols=2, scale=1.1)
        legend.AddEntry('', 'FWHM/{}'.format('FWC' if use_fwc else 'Mean'), '')
        legend.AddEntry('', '{:.2f} ({:.2f})'.format(value.n, value.s), '')
        legend.Draw()
        self.info('FWHM / MPV: {}'.format(value))
        h.Sumw2(False)
        update_canvas()
        return m, fwhm, value

    def model_trap_number(self, f=1000, t=1, max_traps=10000, steps=20, show=True):
        filled_traps = zeros(steps, dtype=int)
        decay = vectorize(self.decay)
        n_traps = []
        for i in range(steps):
            filled_traps[i] = f
            filled_traps = decay(filled_traps, t)
            n_traps.append(min(sum(filled_traps), max_traps))
        g = Draw.make_tgrapherrors(x=arange(steps), y=n_traps)
        format_histo(g, title='Number of Filled Traps', x_tit='Time [s]', y_tit='Number of Traps', y_off=1.7)
        self.Draw(g, draw_opt='ap', lm=.13, show=show)
        return n_traps[-1]

    def draw_n_traps(self, t=1, max_traps=1e5, nbins=20):
        x, y = log_bins(nbins, 100, 1e6)[1], []
        self.PBar.start(x.size)
        for f in x:
            y.append(self.model_trap_number(f, t, max_traps, show=False))
            self.PBar.update()
        g = Draw.make_tgrapherrors(x=x / 1000, y=y)
        format_histo(g, title='Number of Filled Traps vs Flux', x_tit='Flux [kHz/cm^{2}]', y_tit='Number of Filled Traps', y_off=1.7)
        self.Draw(g, draw_opt='ap', lm=.13, logx=True)

    @staticmethod
    def decay(n, t):
        return count_nonzero(rand(n) <= exp(-1. / t))

    def save_tree(self, cut=None):
        f = TFile('test.root', 'RECREATE')
        t = self.Tree.CloneTree(0)
        n = self.Tree.Draw('Entry$', self.Cut(cut), 'goff')
        good_events = self.Run.get_root_vec(n, dtype='i4')
        self.PBar.start(n)
        for i, ev in enumerate(good_events):
            self.Tree.GetEntry(ev)
            t.Fill()
            self.PBar.update(i)
        f.cd()
        t.Write()
        macro = self.Run.RootFile.Get('region_information')
        if macro:
            macro.Write()
        f.Write()
        f.Close()
        self.info('successfully saved tree with only cut events.')

    def fit_langau(self, h=None, nconv=30, show=True, chi_thresh=8, fit_range=None):
        h = self.draw_signal_distribution(show=show) if h is None and hasattr(self, 'draw_signal_distribution') else h
        fit = Langau(h, nconv, fit_range)
        fit.get_parameters()
        fit(show=show)
        get_last_canvas().Modified()
        get_last_canvas().Update()
        if fit.get_chi2() > chi_thresh and nconv < 80:
            Draw.Count += 5
            self.info('Chi2 too large ({c:2.2f}) -> increasing number of convolutions by 5'.format(c=fit.get_chi2()))
            fit = self.fit_langau(h, nconv + Draw.Count, chi_thresh=chi_thresh, show=show)
        print('MPV: {:1.1f}'.format(fit.get_mpv()))
        self.Draw.Count = 0
        self.Draw.add(fit)
        return fit


if __name__ == '__main__':
    pargs = init_argparser(run=88, has_verbose=True, tree=True, dut=True)
    z = DUTAnalysis(pargs.run, pargs.dut, pargs.testcampaign, pargs.tree, verbose=True)
