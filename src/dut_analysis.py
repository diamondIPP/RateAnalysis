#!/usr/bin/env python
# --------------------------------------------------------
#       parent class for the analysis of a single device under test
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from Extrema import Extrema2D
from telescope_analysis import *
from currents import Currents
from ROOT import TProfile2D
from numpy import linspace, exp, vectorize
from numpy.random import rand
from uncertainties import umath
from dut import DUT


class DUTAnalysis(TelecopeAnalysis):
    def __init__(self, run_number, diamond_nr, test_campaign=None, tree=None, t_vec=None, verbose=False, prnt=True):

        TelecopeAnalysis.__init__(self, run_number, test_campaign, tree, t_vec, verbose)

        self.DUT = DUT(diamond_nr, self.Run.RunInfo)

        self.print_start(run_number, prnt, dut=self.DUT.Name)

        self.update_config()
        self.set_save_directory(join(self.DUT.Name, str(self.Run.Number).zfill(3)))

        self.Currents = Currents(self)

    def update_config(self):
        pass

    def load_dut_nr(self, dut_nr):
        if dut_nr not in self.Run.DUT.Numbers:
            critical('wrong diamond number "{}". The following diamond numbers are valid: {}'.format(dut_nr, self.Run.DUT.Numbers))
        return dut_nr

    def draw_current(self, relative_time=False, averaging=1, show=True, v_range=None, draw_opt='al'):
        self.Currents.draw_indep_graphs(rel_time=relative_time, averaging=averaging, show=show, v_range=v_range, draw_opt=draw_opt)

    @staticmethod
    def get_info_header():
        return ['Run', 'Type', 'Diamond', 'Flux [kHz/cm2]', 'HV [V]']

    def show_information(self, header=True, prnt=True):
        rows = [[self.Run.Number, self.Run.RunInfo['runtype'], self.DUT.Name, '{:14.1f}'.format(self.Run.Flux.n), '{:+6.0f}'.format(self.DUT.Bias)]]
        return print_table(rows, self.get_info_header() if header else None, prnt=prnt)

    # ----------------------------------------
    # region GET
    def get_events(self, cut=None, redo=False):
        cut = self.Cut(cut)
        return do_hdf5(self.make_hdf5_path('Events', run=self.Run.Number, ch=self.DUT.Number, suf=cut.GetName()), self.Run.get_root_vec, redo, dtype='i4', var='Entry$', cut=cut)

    def get_n_entries(self, cut=None):
        return self.Tree.GetEntries(self.Cut(cut).GetTitle())

    def get_current(self):
        return self.Currents.get_current()

    def get_irradiation(self):
        return self.DUT.get_irradiation(self.TCString)

    def get_attenuator(self):
        return False

    def get_ph_str(self):
        pass

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

    def get_uniformity(self, bins=10, redo=False):
        pickle_path = self.make_pickle_path('Signal', 'Uniformity', self.Run.Number, ch=self.DUT.Number, suf=bins)

        def f():
            return self.draw_uniformity(bins=bins, show=False)

        return do_pickle(pickle_path, f, redo=redo)

    def get_track_length_var(self):
        dx2, dy2 = ['TMath::Power(TMath::Tan(TMath::DegToRad() * {}_{}), 2)'.format('slope' if self.Run.has_branch('slope_x') else 'angle', direction) for direction in ['x', 'y']]
        return '{} * TMath::Sqrt({} + {} + 1)'.format(self.DUT.Thickness, dx2, dy2)

    def get_flux(self, corr=True, rel_error=0, show=False):
        return (self._get_flux(prnt=False, show=show) if self.Tree and self.has_branch('rate') else self.Run.get_flux(rel_error)) * (self.get_flux_correction() if corr else 1)

    def get_flux_correction(self):
        def f():
            if not self.Cut.has('fiducial'):
                return 1
            n1, n2 = self.get_n_entries(self.Cut.generate_custom(include=['tracks', 'fiducial'], prnt=False)), self.get_n_entries(self.Cut.get('tracks'))
            a1, a2 = self.Cut.get_fiducial_area(), min(self.Run.get_unmasked_area().values())
            return n1 / a1 * a2 / n2
        return do_pickle(self.make_pickle_path('Flux', 'Corr', self.Run.Number, self.DUT.Number), f)

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

    def draw_pulse_height(self, *args, **kwargs):
        pass

    def draw_ph_pull(self, *args, **kwargs):
        return self._draw_ph_pull(*args, **kwargs)

    def _draw_ph_pull(self, event_bin_width=None, fit=True, bin_width=.5, bins=None, show=True, save=True):
        p = self.draw_pulse_height(event_bin_width, show=False, save=False)[0]
        self.format_statbox(all_stat=True, fit=fit)
        h = get_pull(p, 'Signal Bin{0} Distribution'.format(self.Bins.BinSize), bins=self.Bins.get_pad_ph(bin_width=bin_width) if bins is None else bins, fit=fit)
        format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.5, fill_color=self.FillColor, draw_first=True)
        self.save_histo(h, 'SignalBin{0}Disto'.format(self.Bins.BinSize), save=save, lm=.12, show=show)
        return h

    def draw_signal_distribution(self, *args, **kwargs):
        pass

    def draw_fid_cut(self, scale=10):
        self.Cut.draw_fid_cut(scale)

    def draw_track_length(self, show=True, save=True):
        h = TH1F('htd', 'Track Distance in Diamond', 200, self.DUT.Thickness, self.DUT.Thickness + 1)
        self.Tree.Draw('{}>>htd'.format(self.get_track_length_var()), 'n_tracks', 'goff')
        format_histo(h, x_tit='Distance [#mum]', y_tit='Entries', y_off=2, lw=2, stats=0, fill_color=self.FillColor, ndivx=405)
        self.save_histo(h, 'DistanceInDia', show, lm=.16, save=save)
        return h

    def draw_signal_vs_angle(self, mode='x', bin_size=.1, show=True):
        p = TProfile('psa{}'.format(mode), 'Pulse Height vs. Angle in {}'.format(mode.title()), *self.Bins.get_angle(bin_size))
        self.Tree.Draw('{}:angle_{m}>>psa{m}'.format(self.generate_signal_name(), m=mode), self.Cut(), 'goff')
        self.format_statbox(entries=True)
        format_histo(p, x_tit='Track Angle [deg]', y_tit='Pulse Height [mV]', y_off=1.5)
        self.draw_histo(p, show, lm=.12)

    # ----------------------------------------
    # region SIGNAL MAP
    def draw_signal_map(self, res=None, cut=None, fid=False, hitmap=False, redo=False, bins=None, z_range=None, size=None, show=True, save=True, prnt=True):

        cut = self.Cut.generate_custom(exclude=['fiducial'], prnt=prnt) if not fid and cut is None else self.Cut(cut)
        suf = '{c}_{ch}_{res}'.format(c=cut.GetName(), ch=self.Cut.CutConfig['chi2_x'], res=res if bins is None else '{}x{}'.format(bins[0], bins[2]))
        pickle_path = self.make_pickle_path('SignalMaps', 'Hit' if hitmap else 'Signal', run=self.Run.Number, ch=self.DUT.Number, suf=suf)

        def func():
            set_root_output(0)
            name = 'h_hm' if hitmap else 'h_sm'
            atts = [name, 'Track Hit Map' if hitmap else 'Signal Map'] + (self.Bins.get_global(res, mm=True) if bins is None else bins)
            h1 = TH2I(*atts) if hitmap else TProfile2D(*atts)
            self.info('drawing {mode}map of {dia} for Run {run}...'.format(dia=self.DUT.Name, run=self.Run.Number, mode='hit' if hitmap else 'signal '), prnt=prnt)
            y, x = self.Cut.get_track_vars(self.DUT.Number - 1, mm=True)
            self.Tree.Draw('{z}{y}:{x}>>{h}'.format(z=self.get_ph_str() + ':' if not hitmap else '', x=x, y=y, h=name), cut, 'goff')
            set_2d_ranges(h1, *([3, 3] if size is None else size))
            adapt_z_range(h1) if not hitmap else do_nothing()
            return h1

        set_palette(pal=1 if hitmap else 53)
        self.format_statbox(entries=True, x=0.82)
        h = do_pickle(pickle_path, func, redo=redo)
        z_tit = 'Number of Entries' if hitmap else 'Pulse Height [mV]'
        format_histo(h, x_tit='Track Position X [mm]', y_tit='Track Position Y [mm]', y_off=1.4, z_off=1.5, z_tit=z_tit, ncont=50, ndivy=510, ndivx=510, z_range=z_range)
        self.draw_histo(h, '', show, lm=.12, rm=.16, draw_opt='colzsame')
        self.draw_fid_cut(scale=10)
        # self.draw_detector_size(scale=10)
        self.save_plots('HitMap' if hitmap else 'SignalMap2D', prnt=prnt, save=save)
        return h

    def draw_hitmap(self, res=None, cut=None, fid=False, redo=False, z_range=None, size=None, show=True, save=True, prnt=True):
        cut = self.Cut.get('tracks') if cut is None else self.Cut(cut)
        return self.draw_signal_map(res, cut, fid, hitmap=True, redo=redo, bins=None, z_range=z_range, size=size, show=show, save=save, prnt=prnt)

    def split_signal_map(self, m=2, n=2, grid=True, redo=False, show=True):
        fid_cut = array(self.Cut.CutConfig['fiducial']) * 10
        if not fid_cut.size:
            log_critical('fiducial cut not defined for {}'.format(self.DUT.Name))
        x_bins = linspace(fid_cut[0], fid_cut[1], m + 1)
        y_bins = linspace(fid_cut[2], fid_cut[3], n + 1)
        bins = [m, x_bins, n, y_bins]
        h = self.draw_signal_map(bins=bins, show=False, fid=True, redo=redo)
        format_histo(h, x_range=[fid_cut[0], fid_cut[1] - .01], y_range=[fid_cut[2], fid_cut[3] - .01], name='hssm', stats=0)
        self.draw_histo(h, show=show, lm=.12, rm=.16, draw_opt='colzsame')
        self.draw_grid(x_bins, y_bins, width=2) if grid else do_nothing()
        self.save_plots('SplitSigMap')
        return h, x_bins, y_bins

    def draw_sig_map_disto(self, *args, **kwargs):
        return self._draw_sig_map_disto(*args, **kwargs)

    def _draw_sig_map_disto(self, res=None, cut=None, fid=True, x_range=None, redo=False, normalise=False, ret_value=False, ph_bins=None, show=True, save=True):
        source = self.draw_signal_map(res, cut, fid, redo=redo, show=False, save=False)
        h = TH1F('h_smd', 'Signal Map Distribution', *([400, 0, 4] if normalise else self.Bins.get_pad_ph(bin_width=.2) if ph_bins is None else ph_bins))
        normalisation = 1 if not normalise else self.get_pulse_height()
        values = get_2d_hist_vec(source) / normalisation
        [h.Fill(v.n) for v in values]
        x_range = increased_range([h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(5), h.FindLastBinAbove(5)]], .3, .3) if x_range is None else x_range
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=2, fill_color=self.FillColor, x_range=x_range)
        self.save_histo(h, 'SignalMapDistribution', lm=.15, show=show, save=save)
        return mean_sigma(values) if ret_value else h

    def get_sm_std(self, res=None, redo=False):
        return self.draw_sig_map_disto(show=False, res=res, redo=redo, normalise=True, ret_value=True, save=False)[1]

    def draw_sm_profile(self, mode='x', factor=1.5, cut=None, fid=False, hitmap=False, redo=False, show=True):
        s = self.draw_signal_map(factor, cut, fid, hitmap=hitmap, redo=redo, show=False)
        g = self.make_tgrapherrors('g_smp', 'Signal Map Profile')
        values = [[] for _ in xrange(s.GetNbinsX() if mode == 'x' else s.GetNbinsY())]
        for xbin in xrange(s.GetNbinsX()):
            for ybin in xrange(s.GetNbinsY()):
                value = s.GetBinContent(xbin, ybin)
                if value:
                    values[(xbin if mode == 'x' else ybin)].append(value)
        for i, lst in enumerate(values):
            m, sigma = mean_sigma(lst) if lst else (0., 0.)
            xval = s.GetXaxis().GetBinCenter(i) if mode == 'x' else s.GetYaxis().GetBinCenter(i)
            g.SetPoint(i, xval, m)
            g.SetPointError(i, 0, sigma)

        format_histo(g, x_tit='Track in {m} [cm]'.format(m=mode), y_tit='Pulse Height [au]', y_off=1.5, ndivx=515)
        self.save_histo(g, 'SignalMapProfile', draw_opt='ap', lm=.14, show=show, gridx=True)

    def draw_error_signal_map(self, show=True):
        h = self.draw_signal_map(show=False, fid=True).ProjectionXY('', 'c=e')
        format_histo(h, name='hsme', title='Signal Map Errors', y_off=1.6)
        self.save_histo(h, 'SignalMapErrors', lm=.12, rm=.11, show=show, draw_opt='colz')
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

    # TODO check these methods
    def make_region_cut(self):
        return self.Cut.generate_region(self.draw_signal_map(show=False), self.draw_sig_map_disto(show=False))

    def find_2d_regions(self):
        extrema = Extrema2D(self.draw_signal_map(show=False), self.draw_sig_map_disto(show=False))
        extrema.clear_voting_histos()
        extrema.region_scan()
        extrema.show_voting_histos()
        self.save_plots('Regions2D')
        return extrema

    def find_2d_extrema(self, size=1, histo=None, show=True):
        extrema = Extrema2D(self.draw_signal_map(show=False), self.draw_sig_map_disto(show=False))
        extrema.clear_voting_histos()
        extrema.square_scan(size, histo)
        if show:
            extrema.show_voting_histos()
        self.save_plots('Extrema2D')
        return extrema
    # endregion SIGNAL MAP
    # ----------------------------------------

    def draw_uniformity(self, histo=None, bins=10, show=True, redo=False):

        fwhm_gauss = self.Pedestal.get_fwhm() if histo is None and hasattr(self, 'Pedestal') else 0
        h = self.draw_signal_distribution(show=show, normalise=True, sumw2=False, redo=redo) if histo is None else histo
        if histo is not None:
            self.draw_histo(h)
        sleep(.1)
        fit = TF1('fit', 'gaus', 0, 10000)
        h.Fit('fit', 'qs', '', h.GetBinCenter(h.GetMaximumBin() - bins), h.GetBinCenter(h.GetMaximumBin() + bins))
        mpv = make_ufloat((fit.GetParameter(1), fit.GetParError(1) + fit.GetParameter(1) * .02))
        half_max = fit(mpv.n) / 2
        x_fwhm_min, x_fwhm_max = (make_ufloat((h.GetBinCenter(i), h.GetBinWidth(1))) for i in [h.FindFirstBinAbove(half_max), h.FindLastBinAbove(half_max)])
        fwhm_total = x_fwhm_max - x_fwhm_min
        self.draw_vertical_line(mpv.n, 0, 1e7, style=7, w=2)
        self.draw_tlatex(mpv.n + 5, .1 * half_max, 'MPV', align=10)
        self.draw_arrow(x_fwhm_min.n, mpv.n, half_max, half_max, col=2, width=3, opt='<', size=.02)
        self.draw_arrow(x_fwhm_max.n, mpv.n, half_max, half_max, col=2, width=3, opt='<', size=.02)
        self.draw_tlatex(x_fwhm_max.n + 5, half_max, 'FWHM', align=12, color=2)
        if mpv < 2 or fwhm_total < 1:
            log_warning('Could not determine fwhm or mpv')
            return None, None, None
        fwhm = umath.sqrt(fwhm_total ** 2 - fwhm_gauss ** 2)
        value = fwhm / mpv
        legend = self.make_legend(w=.3, y2=.78, nentries=1, margin=.1, cols=2, scale=1.1)
        legend.AddEntry('', 'FWHM/MPV', '')
        legend.AddEntry('', '{:.2f} ({:.2f})'.format(value.n, value.s), '')
        legend.Draw()
        self.info('FWHM / MPV: {}'.format(fwhm / mpv))
        h.Sumw2(False)
        get_last_canvas().Update()
        return mpv, fwhm, value

    def model_trap_number(self, f=1000, t=1, max_traps=10000, steps=20, show=True):
        filled_traps = zeros(steps, dtype=int)
        decay = vectorize(self.decay)
        n_traps = []
        for i in xrange(steps):
            filled_traps[i] = f
            filled_traps = decay(filled_traps, t)
            n_traps.append(min(sum(filled_traps), max_traps))
            print filled_traps
        g = self.make_tgrapherrors('gt', 'Number of Filled Traps', x=arange(steps), y=n_traps)
        format_histo(g, x_tit='Time [s]', y_tit='Number of Traps', y_off=1.7)
        self.draw_histo(g, draw_opt='ap', lm=.13, show=show)
        return n_traps[-1]

    def draw_n_traps(self, t, max_traps=1e5, bins=20):
        x, y = log_bins(bins, 100, 1e6)[1], []
        self.PBar.start(x.size)
        for f in x:
            y.append(self.model_trap_number(f, t, max_traps, show=False))
            self.PBar.update()
        g = self.make_tgrapherrors('gnt', 'Number of Filled Traps vs Flux', x=x / 1000, y=y)
        format_histo(g, x_tit='Flux [kHz/cm^{2}]', y_tit='Number of Filled Traps', y_off=1.7)
        self.draw_histo(g, draw_opt='ap', lm=.13, logx=True)

    @staticmethod
    def decay(n, t):
        return count_nonzero(rand(n) <= exp(-1. / t))
