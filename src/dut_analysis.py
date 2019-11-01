#!/usr/bin/env python
# --------------------------------------------------------
#       parent class for the analysis of a single device under test
# created on Oct 30th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from Extrema import Extrema2D
from telescope_analysis import *
from CurrentInfo import Currents
from ROOT import TProfile2D
from numpy import linspace


class DUTAnalyis(TelecopeAnalysis):
    def __init__(self, run_number, diamond_nr, test_campaign=None, tree=None, t_vec=None, verbose=False, prnt=True):

        TelecopeAnalysis.__init__(self, run_number, test_campaign, tree, t_vec, verbose, prnt)

        self.DUTNumber = self.load_dut_nr(diamond_nr)
        self.DUTName = self.Run.DUTNames[diamond_nr - 1]
        self.Bias = self.Run.Bias[diamond_nr - 1]

        self.update_config()
        self.set_save_directory(join(self.DUTName, str(self.RunNumber).zfill(3)))

        self.Currents = Currents(self)

    def update_config(self):
        pass

    def load_dut_nr(self, dut_nr):
        if dut_nr not in self.Run.DUTNumbers:
            critical('wrong diamond number "{}". The following diamond numbers are valid: {}'.format(dut_nr, self.Run.DUTNumbers))
        return dut_nr

    def draw_current(self, relative_time=False, averaging=1, show=True, v_range=None):
        self.Currents.draw_indep_graphs(rel_time=relative_time, averaging=averaging, show=show, v_range=v_range)

    # ----------------------------------------
    # region GET
    def get_current(self):
        return self.Currents.get_current()

    def get_irradiation(self):
        return self.Run.get_irradiations()[self.DUTNumber - 1]

    def get_attenuator(self):
        return False

    def get_ph_str(self):
        pass

    def get_pulse_height(self, bin_size=None, cut=None, redo=False):
        pass

    def get_sm_data(self, cut=None, fid=False):
        """ :return: signal map data as numpy array [[x], [y], [ph]] with units [[mm], [mm], [mV]]
            :param cut: applies all cuts if None is provided.
            :param fid: return only values within the fiducial region set in the AnalysisConfig.ini"""
        y, x = self.Cut.get_track_vars(self.DUTNumber - 1, mm=True)
        cut = self.Cut.generate_special_cut(excluded=['fiducial'], prnt=False) if not fid and cut is None else self.Cut(cut)
        n = self.Tree.Draw('{x}:{y}:{z}'.format(z=self.get_ph_str(), x=x, y=y), cut, 'goff')  # *10 to get values in mm
        return self.Run.get_root_vecs(n, 3)
    # endregion GET
    # ----------------------------------------

    def draw_pulse_height(self, *args, **kwargs):
        pass

    def draw_ph_pull(self, *args, **kwargs):
        return self._draw_ph_pull(*args, **kwargs)

    def _draw_ph_pull(self, event_bin_width=None, fit=True, bin_width=.5, bins=None, show=True, save=True):
        p = self.draw_pulse_height(event_bin_width, show=False)[0]
        self.format_statbox(all_stat=True, fit=fit)
        h = get_pull(p, 'Signal Bin{0} Distribution'.format(self.Bins.BinSize), bins=self.Bins.get_pad_ph(bin_width=bin_width) if bins is None else bins, fit=fit)
        format_histo(h, x_tit='Pulse Height [au]', y_tit='Entries', y_off=1.5, fill_color=self.FillColor, draw_first=True)
        self.save_histo(h, 'SignalBin{0}Disto'.format(self.Bins.BinSize), save=save, lm=.12, show=show)
        return h

    def draw_signal_distribution(self, *args, **kwargs):
        pass

    def draw_fid_cut(self, scale=1):
        self.Cut.draw_fid_cut(scale)

    def draw_detector_size(self, scale=1):
        split_runs = self.Cut.get_fiducial_splits()
        first_cut_name = 'detector size' if self.Config.has_option('CUT', 'detector size') else 'detector size 1'
        values = next(self.Cut.load_dia_config('detector size {n}'.format(n=i + 1) if i else first_cut_name) for i in xrange(len(split_runs)) if self.RunNumber <= split_runs[i])
        if values is None:
            return
        x, y, lx, ly = values
        cut = TCutG('det{}'.format(scale), 5, array([x, x, x + lx, x + lx, x], 'd') * scale, array([y, y + ly, y + ly, y, y], 'd') * scale)
        cut.SetVarX(self.Cut.get_track_var(self.DUTNumber - 1, 'x'))
        cut.SetVarY(self.Cut.get_track_var(self.DUTNumber - 1, 'y'))
        self.Objects.append(cut)
        cut.SetLineWidth(3)
        cut.Draw()

    # ----------------------------------------
    # region SIGNAL MAP
    def draw_signal_map(self, res=None, cut=None, fid=False, hitmap=False, redo=False, bins=None, z_range=None, size=None, show=True, save=True, prnt=True):

        cut = self.Cut.generate_special_cut(excluded=['fiducial'], prnt=prnt) if not fid and cut is None else self.Cut(cut)
        suf = '{c}_{ch}_{res}'.format(c=cut.GetName(), ch=self.Cut.CutConfig['chi2X'], res=res if bins is None else '{}x{}'.format(bins[0], bins[2]))
        pickle_path = self.make_pickle_path('SignalMaps', 'Hit' if hitmap else 'Signal', run=self.RunNumber, ch=self.DUTNumber, suf=suf)

        def func():
            set_root_output(0)
            name = 'h_hm' if hitmap else 'h_sm'
            atts = [name, 'Track Hit Map' if hitmap else 'Signal Map'] + (self.Bins.get_global(res, mm=True) if bins is None else bins)
            h1 = TH2I(*atts) if hitmap else TProfile2D(*atts)
            self.info('drawing {mode}map of {dia} for Run {run}...'.format(dia=self.DUTName, run=self.RunNumber, mode='hit' if hitmap else 'signal '), prnt=prnt)
            y, x = self.Cut.get_track_vars(self.DUTNumber - 1, mm=True)
            self.Tree.Draw('{z}{y}:{x}>>{h}'.format(z=self.get_ph_str() + ':' if not hitmap else '', x=x, y=y, h=name), cut, 'goff')
            set_2d_ranges(h1, *([3, 3] if size is None else size))
            adapt_z_range(h1) if not hitmap else do_nothing()
            z_tit = 'Number of Entries' if hitmap else 'Pulse Height [mV]'
            format_histo(h1, x_tit='Track Position X [mm]', y_tit='Track Position Y [mm]', y_off=1.4, z_off=1.5, z_tit=z_tit, ncont=50, ndivy=510, ndivx=510, z_range=z_range)
            return h1

        set_palette(pal=1 if hitmap else 53)
        self.format_statbox(entries=True, x=0.82)
        h = do_pickle(pickle_path, func, redo=redo)
        self.draw_histo(h, '', show, lm=.12, rm=.16, draw_opt='colzsame')
        self.draw_fid_cut(scale=10)
        # self.draw_detector_size(scale=10)
        self.save_plots('HitMap' if hitmap else 'SignalMap2D', prnt=prnt, save=save)
        return h

    def draw_hitmap(self, res=None, cut=None, fid=False, redo=False, z_range=None, size=None, show=True, save=True, prnt=True):
        cut = self.Cut.CutStrings['tracks'] if cut is None else self.Cut(cut)
        self.draw_signal_map(res, cut, fid, hitmap=True, redo=redo, bins=None, z_range=z_range, size=size, show=show, save=save, prnt=prnt)

    def split_signal_map(self, m=2, n=2, grid=True, redo=False, show=True):
        fid_cut = array(self.Cut.CutConfig['fiducial']) * 10
        if not fid_cut.size:
            log_critical('fiducial cut not defined for {}'.format(self.DUTName))
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

    def calc_signal_spread(self, min_percent=5, max_percent=99):
        """ Calculates the relative spread of mean signal response from the 2D signal response map. """
        h = self.draw_sig_map_disto(show=False)
        q = array([min_percent, max_percent]) / 100.
        y = zeros(2, 'd')
        mean_error = mean([v.n for v in get_2d_hist_vec(self.draw_error_signal_map(show=False))])
        h.GetQuantiles(2, y, q)
        ratio = make_ufloat([y[1], mean_error]) / make_ufloat([y[0], mean_error]) - 1
        self.info('Relative Signal Spread is: {:2.2f} %'.format(ratio * 100))
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
