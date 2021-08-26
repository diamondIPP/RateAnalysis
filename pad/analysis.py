from ROOT import gRandom, TCut
from numpy import insert, sum, round, in1d, max, unravel_index, argmax

from pad.pedestal import PedestalAnalysis
from pad.timing import TimingAnalysis
from src.dut_analysis import *
from pad.cut import PadCut
from pad.run import PadRun
from pad.peaks import PeakAnalysis
from pad.pulser import PulserAnalysis
from pad.waveform import Waveform
from src.mc_signal import MCSignal


class PadAnalysis(DUTAnalysis):
    def __init__(self, run_number, diamond_nr, test_campaign=None, load_tree=True, verbose=False, prnt=True):

        DUTAnalysis.__init__(self, run_number, diamond_nr, test_campaign, load_tree, verbose, prnt)

        # MAIN
        self.Channel = self.get_channel()
        self.DigitiserBinWidth = self.Run.DigitiserBinWidth

        if self.Tree.Hash():
            self.Polarity = self.get_polarity()
            # Regions & Ranges
            self.SignalRegionName = self.load_region_name()
            self.SignalRegion = self.get_region()
            self.PeakIntegralName = self.load_peak_integral()
            self.PeakIntegral = self.Run.PeakIntegrals[self.DUT.Number - 1][self.PeakIntegralName]

            # Signal Names
            self.SignalName = self.get_signal_name()
            self.RawName = self.get_signal_name(peak_int=1)
            self.PedestalName = self.get_pedestal_name()
            self.PeakName = self.get_peak_name()

            # Modules
            self.Timing = TimingAnalysis(self)
            self.Pedestal = PedestalAnalysis(self)
            self.Waveform = Waveform(self)
            self.Peaks = PeakAnalysis(self)
            self.Pulser = PulserAnalysis(self)
            self.MC = MCSignal(self)

        self.print_finished(prnt=prnt)

    @staticmethod
    def get_info_header():
        return ['Run', 'Type', 'Diamond', 'Flux [kHz/cm2]', 'HV [V]', 'Region', 'Integral']

    @quiet
    def save_plots(self):
        super(PadAnalysis, self).save_plots()
        self.Pedestal.draw_disto_fit(show=False)
        self.Pulser.draw_pulse_height(show=False)
        self.Pulser.draw_distribution_fit(show=False)
        self.Pulser.draw_pedestal_fit(show=False)

    @reload_tree
    def get_data(self):
        return [self.get_flux(), self.get_current(), self.get_pulse_height(), self.get_pedestal(), self.get_pedestal(par=2), self.Pulser.get_pulse_height(),
                self.Pulser.get_sigma(), self.get_pedestal(pulser=True), self.get_pedestal(pulser=True, par=2), ufloat(self.get_n_entries(), 0)]

    # ----------------------------------------
    # region INFO
    def show_integral_names(self):
        rows = []
        for i, name in enumerate(self.get_integral_names()):
            if name.startswith('ch{}'.format(self.Channel)):
                v = array(name.split('_'))
                v = v if v.size == 4 else insert(v, 2, '')
                r, i = array([self.Run.IntegralRegions[self.DUT.Number - 1]['_'.join(v[1:3]).strip('_')], self.Run.PeakIntegrals[self.DUT.Number - 1][v[3]]]) * self.DigitiserBinWidth
                rows.append([str(i).rjust(3), v[1].title(), v[2], str(r), remove_letters(v[3]), str(i), name])
        print_table(rows, ['Nr.', 'Type', 'Region', 'Value [ns]', 'Int', 'Value [ns]', 'Name'])

    def show_information(self, header=True, prnt=True, ret_row=False):
        peak_int = f'{self.PeakIntegral} ({remove_letters(self.PeakIntegralName)})'
        region = f'{self.SignalRegion} ({self.SignalRegionName.split("_")[-1]})'
        rows = [[self.Run.Number, self.Run.Info['runtype'], self.DUT.Name, f'{self.Run.Flux.n:14.1f}', f'{self.DUT.Bias:+5.0f}', region, peak_int]]
        if ret_row:
            return rows[0]
        print_table(rows, self.get_info_header() if header else None, prnt=prnt)
    # endregion INFO
    # ----------------------------------------

    # ----------------------------------------
    # region INIT
    @staticmethod
    def init_run(run_number, testcampaign, load_tree, verbose):
        return PadRun(run_number, testcampaign, load_tree, verbose)

    def init_cut(self):
        return PadCut(self)

    def update_config(self):
        self.Config.read(join(self.Dir, 'config', self.TCString, 'PadConfig.ini'))

    def load_region_name(self, sig_type='signal', region=None):
        short = choose(region, self.Config.get_value('SIGNAL', '{} region'.format(sig_type), default=''))
        return '{}_{}'.format(sig_type, short) if short else [r for r in self.Run.IntegralRegions[0] if r.startswith(sig_type)][0]

    def load_peak_integral(self):
        peak_int = 'PeakIntegral{}'.format(self.Config.get('SIGNAL', 'peak integral'))
        return peak_int if peak_int in self.Run.PeakIntegrals[self.DUT.Number - 1] else self.Run.PeakIntegrals[self.DUT.Number - 1].keys()[0]

    def get_signal_number(self, region=None, peak_integral=None, sig_type='signal'):
        region = self.load_region_name(sig_type) if region is None else region if sig_type in region else '_'.join([sig_type] + ([region] if region else []))
        peak_integral = self.load_peak_integral() if peak_integral is None else peak_integral if 'Peak' in str(peak_integral) else 'PeakIntegral{}'.format(peak_integral)
        int_name = 'ch{ch}_{reg}_{int}'.format(ch=self.get_channel(), reg=region, int=peak_integral)
        return self.get_integral_names().index(int_name)

    def get_signal_name(self, region=None, peak_int=None, sig_type='signal', t_corr=True):
        return '{}IntegralValues[{}]'.format('Time' if t_corr else '', self.get_signal_number(region, peak_int, sig_type))

    def get_pedestal_name(self, region=None, peak_int=None):
        return self.get_signal_name(region=region, peak_int=peak_int, sig_type='pedestal')

    def get_peak_name(self, region=None, type_='signal', t_corr=True):
        peak_name = 'IntegralPeakTime' if t_corr else 'IntegralPeaks'
        return '{name}[{num}]'.format(name=peak_name, num=self.get_signal_number(region, None, type_))

    def get_all_signal_names(self, sig_type='signal'):
        reg_ints = [(r, i) for r in self.Run.IntegralRegions[self.DUT.Number - 1] for i in self.Run.PeakIntegrals[self.DUT.Number - 1] if sig_type in r]
        return {self.get_signal_name(r, i, sig_type): r.replace(sig_type, '').strip('_') + i.replace('PeakIntegral', '') for r, i in reg_ints}
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_polarity(self):
        return int(self.Run.TreeConfig.get('General', 'polarities').split()[self.get_channel()])

    def get_channel(self):
        return self.Run.Channels[self.DUT.Number - 1]

    def get_signal_var(self, signal=None, evnt_corr=True, off_corr=False, cut=None, region=None, peak_int=None, sig_type='signal', t_corr=True):
        sig_name = choose(signal, self.get_signal_name(region, peak_int, sig_type, t_corr))
        ped = self.Pedestal.get_mean(cut=cut, raw=True).n if off_corr else self.get_pedestal_name() if evnt_corr else None
        return f'{self.get_polarity()} * {sig_name if ped is None else f"({sig_name} - {ped})"}'

    def get_raw_signal_var(self, region=None, peak_int=None, sig_type='signal'):
        return self.get_signal_var(None, False, False, None, region, peak_int, sig_type)

    def get_eff_var(self, n_sigma=3):
        return '{} > {} * {}'.format(self.get_signal_var(evnt_corr=False, off_corr=False), n_sigma, self.Pedestal.get_signal_var())

    def get_ph_var(self, ped=False):
        return self.Pedestal.get_signal_var() if ped else self.get_signal_var()

    def get_integral_names(self):
        return self.Run.TreeConfig.get_value('Integral Names', 'Names', dtype=list)

    def get_region(self, sig_type='signal', region=None):
        return array(self.Run.IntegralRegions[self.DUT.Number - 1][region if region is not None and '_' in region else self.load_region_name(sig_type, region)])

    def get_signal_range(self, lfac=0, rfac=0, t_corr=True, region=None):
        return ax_range(self.get_region(region=region) * (self.DigitiserBinWidth if t_corr else 1), None, lfac, rfac)

    def get_peak_integral(self, name=None):
        name = choose(name, self.PeakIntegralName)
        return array(self.Run.PeakIntegrals[self.DUT.Number - 1]['PeakIntegral{}'.format(name) if 'Peak' not in str(name) else name])

    def get_short_name(self, signal=None, signal_type='signal'):
        return self.get_all_signal_names(signal_type)[signal if signal is not None else self.SignalName]

    def get_attenuator(self):
        return self.DUT.Attenuator

    def get_ph_data(self, cut=None):
        """ :return: pulse height data as numpy array [[time] [ph]] with units [[s], [mV]]
            :param cut: applies all cuts if None is provided."""
        return self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var()], cut=self.Cut(cut))

    def get_ph_values(self, region=None, name=None, evnt_corr=True, cut=None):
        """ :returns: all pulse height values for a given cut. """
        return self.Run.get_tree_vec(var=self.get_signal_var(name, evnt_corr, cut=cut, region=region), cut=self.Cut(cut))

    @save_pickle('Fit', sub_dir='PH', suf_args='[0, 1, 2, 3]')
    def _get_pulse_height(self, bin_size=None, cut=None, corr=True, sig=None, _redo=False):
        """ :returns: fitted (pol0) pulse height over time """
        return self.draw_pulse_height(sig=sig, bin_size=bin_size, cut=cut, corr=corr, show=False, save=False, redo=_redo)[1][0]

    @update_pbar
    def get_pulse_height(self, bin_size=None, cut=None, corr=True, sig=None, sys_err=0, peaks=False, corr_ph=True, redo=False):
        ph = self.Peaks.get_bunch_height() if peaks else self.correct_ph(self._get_pulse_height(bin_size, cut, corr, sig, _redo=redo), cut, corr_ph)
        return ufloat(ph.n, sqrt(ph.s ** 2 + (sys_err * ph.n) ** 2))

    def correct_ph(self, ph=None, cut=None, corr=True):
        ph = choose(ph, self._get_pulse_height, cut=self.Cut.get_tp())
        if (not self.Cut.has('bucket') or self.Cut(cut).GetName() == 'tp') and corr and not self.Cut(cut).GetName() == 'bful':
            r = self.estimate_bucket_ratio() / 100  # % -> value
            return ph if r == 0 else (ph - r * self.Pedestal.get_under_signal()[0]) / (1 - r)
        return ph

    def get_min_signal(self, name=None):
        h = self.draw_signal_distribution(show=False, save=False, sig=name)
        return h.GetBinCenter(h.FindFirstBinAbove(h.GetMaximum() * .01))

    def get_t_bins(self, bin_size=None):
        xmin, xmax = self.SignalRegion * self.DigitiserBinWidth
        return Bins.make(xmin, xmax, choose(bin_size, default=self.DigitiserBinWidth))

    def get_alignment(self):
        from pad.alignment import PadAlignment
        return PadAlignment

    def print_results(self, prnt=True):
        rows = [[u_to_str(v, prec=2) for v in [self.get_pulse_height(), self.Pedestal.get_mean(), self.Pulser.get_pulse_height()]]]
        print_table(header=['Signal [mV]', 'Pedestal [mV]', 'Pulser [mV]'], rows=rows, prnt=prnt)

    def get_efficiency(self, n_sigma=3, redo=False):
        def f():
            return calc_eff(values=self.get_tree_vec(self.get_eff_var(n_sigma), self.Cut()))
        return do_pickle(self.make_simple_pickle_path(suf=str(n_sigma), sub_dir='Efficiency'), f, redo=redo)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region ALIASES
    def get_additional_peak_height(self):
        return self.Peaks.draw_additional(show=False)

    def get_peak_flux(self, corr=True, prnt=True):
        return self.Peaks.get_flux(corr=corr, prnt=prnt)

    def get_n_peaks(self, start_bunch=None, end_bunch=None):
        return self.Peaks.get_n_additional(start_bunch, end_bunch)

    def get_pedestal(self, pulser=False, par=1, redo=False):
        return self.Pulser.get_pedestal(par, redo=redo) if pulser else self.Pedestal.get(par, redo=redo)

    def get_peak_timing(self, par=1, redo=False):
        return self.Timing.get(par, redo)

    def get_peak_times(self, cut=None):
        return self.Peaks.get_from_tree(cut=cut)

    def draw_timing(self):
        self.Timing.draw_all()

    def draw_pulser_rate(self, show=True, prnt=True):
        self.Pulser.draw_rate(show=show, prnt=prnt)

    def draw_region(self, region=None):
        self.Waveform.draw_region(region)
    # endregion ALIASES
    # ----------------------------------------

    # ----------------------------------------
    # region 2D MAPS
    def draw_efficiency_vs_threshold(self, thresh=None, bin_width=.5, show=True):
        thresh = self.Pedestal.get_noise().n * 4 if thresh is None else thresh
        h = self.draw_signal_distribution(show=False, bin_width=bin_width)
        bins = range(h.FindBin(0), h.FindBin(thresh) + 2)
        err = zeros(1)
        efficiencies = [ufloat(h.IntegralAndError(ibin, h.GetNbinsX(), err), err[0]) / h.Integral(0, h.GetNbinsX()) * 100 for ibin in bins]
        g = Draw.make_tgrapherrors([h.GetBinCenter(ibin) for ibin in bins], efficiencies, title='Detector Efficiency')
        format_histo(g, x_tit='Threshold [mV]', y_tit='Efficiency [%]', y_off=1.3)
        self.Draw(g, draw_opt='ap', lm=.12, show=show)
        Draw.vertical_line(x=self.Pedestal.get_noise().n * 3, ymin=0, ymax=110, w=2)
        Draw.tlatex(x=self.Pedestal.get_noise().n * 3, y=95, text=' 3 #times noise', align=10)
        self.Draw.save_plots('EffThresh')

    def draw_low_sig_map(self, high=10, low=None, fid=False, cut=None, hit_map=True):
        low = '&&{}>{}'.format(self.get_signal_var(), low) if low is not None else ''
        fid_cut = self.Cut.generate_custom(exclude='fiducial') if cut is None else self.Cut(cut)
        kwargs = {'redo': True, 'cut': TCut('{}<{}{}'.format(self.get_signal_var(), high, low)) + (self.Cut(cut) if fid else fid_cut)}
        self.draw_hitmap(**kwargs) if hit_map else self.draw_signal_map(**kwargs)

    @staticmethod
    def _correlate_sm(sm1, sm2, col=True):
        """ optimise correlation of two signal maps by translating one map. """
        a1, a2 = [get_2d_hist_vec(sm, err=False, flat=False) for sm in [sm1, sm2]]
        n1, n2 = [get_2d_bin_entries(sm) for sm in [sm1, sm2]]
        a1[n1 < .1 * max(n1)] = 0  # set bins with low stats to 0
        a2[n2 < .1 * max(n2)] = 0
        c = array([[correlate(a1, roll(a2, [x, y], axis=[0, 1])) for x in range(a1.shape[0])] for y in range(a1.shape[1])])
        return c if col else c.T

    def correlate_sm(self, run, col=True, **kwargs):
        sm1, sm2 = [ana.get_signal_map(fid=True, **kwargs) for ana in [self, PadAnalysis(run, self.DUT.Number, self.TCString, verbose=False, prnt=False) if isint(run) else run]]
        return self._correlate_sm(sm1, sm2, col)

    def find_best_sm_correlation(self, run):
        c = self._correlate_sm(*run) if is_iter(run) else self.correlate_sm(run, res=1)
        self.info(f'best correlation: {max(c):.2f} at shift {unravel_index(argmax(c), c.shape)}')
        return max(c)

    def draw_sm_correlation_vs_shift(self, run, col=True, res=1, n=4, **kwargs):
        c = self.correlate_sm(run, col, res=res)
        iy, ix = unravel_index(argmax(c), c.shape)
        x = (arange(2 * n + 1) - n + ix)
        x = x % c.shape[1] if max(x) > c.shape[1] else x
        self.Draw.graph(x, c[iy][x], **prep_kw(kwargs, x_tit=f'Shift in {"X" if col else "Y"}', y_tit='Correlation Factor'))

    def draw_sm_correlation(self, run, m=10, show=True):
        x0, x1 = [get_2d_hist_vec(f.split_signal_map(m, show=False)[0], err=False) for f in [self, PadAnalysis(run, self.DUT.Number, self.TCString)]]
        g = self.Draw.histo_2d(x0, x1, self.Bins.get_pad_ph(2) * 2, 'Signal Map Correlation', show=show, x_tit='Pulse Height {} [mV]'.format(self.Run.Number), y_tit='Pulse Height {} [mV]'.format(run),
                               x_range=ax_range(x0, 0, .1, .1), y_range=ax_range(x1, 0, .1, .1))
        Draw.info('Correlation Factor: {:.2f}'.format(g.GetCorrelationFactor()))

    def draw_corr_coeff(self, run, show=True):
        x = [5, 10, 25, 50, 100, 200]
        ana = PadAnalysis(run, self.DUT.Number, self.TCString)
        y = [correlate(*[get_2d_hist_vec(f.split_signal_map(m, show=False)[0], err=False, zero_supp=0) for f in [self, ana]]) for m in x]
        self.Draw.graph(x, y, 'Signal Map Correlation', x_tit='Number of Bins', y_tit='Correlation Coefficient', show=show)

    # endregion 2D MAPS
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def draw_disto_vs_time(self, bin_width=.2, signal=None, evnt_corr=False, off_corr=False, rel_t=False, show=True):
        t, ph = self.get_tree_vec(var=[self.get_t_var(), self.get_signal_var(signal, evnt_corr, off_corr)], cut=self.Cut())
        bins = self.Bins.get_time() + self.Bins.get_pad_ph(bin_width)
        h = self.Draw.histo_2d(t, ph, bins, 'Signal vs. Time', pal=53, x_tit='Time [min]', y_tit='Pulse Height [au]', t_ax_off=self.get_t_off(rel_t), show=False)
        self.Draw(h, 'SignalTime', show, draw_opt='colz', rm=.15)
        return h

    def draw_pulse_height_vs_binsize(self, show=True):
        bin_sizes = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]
        pulse_heights = [fit2u(self.draw_pulse_height(bin_size=bin_size, show=False)[1], par=0) for bin_size in bin_sizes]
        g = Draw.make_tgrapherrors(bin_sizes, pulse_heights, title='Pulse Height vs Number of Events per Bin', x_tit='Number of Events per Bin', y_tit='Pulse Height [mV]', y_off=1.2, x_off=1.2)
        self.Draw(g, lm=.12, show=show, gridy=True, logx=True)

    def draw_mean_vs_binsize(self, steps=10, show=True):
        values = self.get_ph_values()
        splits = 2 ** arange(min(steps, int(log2(values.size))) + 1)
        x = values.size / splits
        means = [mean_sigma([mean(ar) for ar in split(values[:-(values.size % s)] if values.size % s else values, s)])[0] for s in splits]
        g = Draw.make_tgrapherrors(x, means, title='Mean vs Bin Size', x_tit='Bin Size', y_tit='Pulse Height [mV]', x_range=[min(x) / 2, max(x) * 2])
        self.Draw(g, lm=.14, show=show, gridy=True, logx=True)

    @save_pickle('Trend', sub_dir='PH', suf_args='[0, 1, 2, 3]')
    def get_pulse_height_trend(self, bin_size=None, sig=None, cut=None, corr=True, _redo=False):
        ph, t = self.get_tree_vec(var=[self.get_signal_var(sig, corr), self.get_t_var()], cut=self.Cut(cut))
        return self.Draw.profile(t, ph, self.Bins.get_time(bin_size, cut), 'Pulse Height Trend', y_tit='Mean Pulse Height [mV]', **self.get_t_args(), graph=True, show=False)

    def draw_pulse_height(self, bin_size=None, sig=None, cut=None, corr=True, redo=False, rel_t=True, **kwargs):
        g = self.get_pulse_height_trend(bin_size, sig, cut, corr, _redo=redo)
        fit = FitRes(g.Fit('pol0', 'qs', '', 0, self.__get_max_fit_pos(g)))
        kwargs = prep_kw(kwargs, y_range=ax_range(get_graph_y(g, err=False), 0, .6, .8), ndivx=505, stats=set_statbox(fit=True, form='.2f'))
        g = self.Draw(g, file_name=f'PulseHeight{Bins.w(bin_size)}', **self.get_t_args(rel_t), **kwargs)
        return g, fit

    def __get_max_fit_pos(self, g, thresh=.8):
        """ look for huge fluctiations in ph graph and return last stable point"""
        x, y = get_graph_vecs(g, err=False)
        if mean(y) > 10:  # if the pulse height is very low there will be always big fluctuations!
            averages = cumsum(y) / (arange(y.size) + 1)
            j = next((i for i, is_bad in enumerate(y[2:] < thresh * averages[1:-1], 2) if is_bad), None)  # find next entry that is below the average of the previous
            if j is not None:
                warning(f'Found PH fluctiation in {self.Run}! Stopping fit after {100 * j / y.size:2.0f}%')
                return x[j - 2]
        return x[-1] + 1000

    @save_pickle('Disto', sub_dir='PH', print_dur=True, suf_args='[0, 1, 2, 3, 4]')
    def get_signal_distribution(self, sig=None, cut=None, evnt_corr=True, off_corr=False, bin_width=None, _redo=False):
        x = self.get_tree_vec(var=self.get_signal_var(sig, evnt_corr, off_corr, cut), cut=self.Cut(cut))
        bin_width = choose(bin_width, Bins.find_width, x=x)
        return self.Draw.distribution(x, self.Bins.get_pad_ph(bin_width), 'Pulse Height Distribution', show=False, x_tit='Pulse Height [mV]', y_off=1.65)

    def draw_signal_distribution(self, sig=None, cut=None, evnt_corr=True, off_corr=False, bin_width=None, redo=False, **kwargs):
        h = self.get_signal_distribution(sig, cut, evnt_corr, off_corr, bin_width, _redo=redo)
        return self.Draw(h, **prep_kw(kwargs, lm=.15, y_off=1.65, file_name='SignalDistribution'))

    def draw_raw_signal_distribution(self, cut=None, bin_width=None, redo=False, **kwargs):
        return self.draw_signal_distribution(self.RawName, cut, False, True, bin_width, redo, file_name='RawSignalDistribution', **kwargs)

    def find_bunch_region(self, n=1):
        w = self.BunchSpacing / self.DigitiserBinWidth
        m0 = mean(self.SignalRegion) + (n - 1) * w
        for key, r in self.Run.IntegralRegions[self.DUT.Number - 1].items():
            if 'signal' in key and r[0] < m0 < r[1] and r[1] - r[0] < 1.1 * w:
                return key

    def draw_bunch_distribution(self, n=1, region=None, corr=False, bin_width=5, cut=True, show=True):
        cut = self.get_pulser_cut() & invert(self.Peaks.get_ebunch_cut(n - 1)) & invert(self.Peaks.get_ebunch_cut(n + 1)) & cut
        x = self.get_tree_vec(self.get_signal_var(region=choose(region, self.find_bunch_region(n)), evnt_corr=corr))
        h = self.Draw.distribution(x[cut], self.Bins.get_pad_ph(bin_width), 'Bunch {} Distribution'.format(n), x_tit='Pulse Height [mV]', show=show)
        if h.GetBinCenter(h.GetMaximumBin()) < 10:
            m, s = self.Pedestal()
            x, y = get_hist_vecs(h, err=False)
            format_histo(h, y_range=ax_range(0, max(y[x > m + 6 * s]), 0, .05))

    def draw_bunch_comparison(self, n, regions, bin_width=5, show=True):
        histos = [self.draw_bunch_distribution(i, r, bin_width, show=False) for i, r in zip(n, regions)]
        self.Draw.stack(histos, 'Bunch Distributions', ['Bucket {}'.format(i + 2) for i in n], scale=True, show=show)

    def draw_ph_peaktime(self, bin_size=None, fine_corr=False, cut=None, region=None, x=None, y=None, y_range=None, xbins=None, prof=True, normalise=False, logz=False, show=True):
        xvar, yvar = self.Timing.get_peak_var(corr=True, fine_corr=fine_corr, region=region), self.get_signal_var()
        x, y = self.get_tree_vec(var=[xvar, yvar], cut=self.Cut(cut)) if x is None and y is None else [choose(i, self.get_tree_vec(var=var, cut=self.Cut(cut))) for i, var in [(x, xvar), (y, yvar)]]
        bins = choose(xbins, self.get_t_bins(bin_size)) + ([] if prof else self.Bins.get_pad_ph(1))
        h = (self.Draw.profile if prof else self.Draw.histo_2d)(x, y, bins, 'Signal vs Peak Position', show=show, logz=logz)
        format_histo(h, x_tit='Signal Peak Position [ns]', y_tit='Pulse Height [mV]', stats=0, y_range=y_range)
        if normalise:
            normalise_bins(h)
        return h

    def draw_signal_vs_tot(self, bin_size=.5, show=True):
        x, y = self.Peaks.get_tot(), self.Run.get_tree_vec(var=self.get_signal_var(), cut=self.Cut())
        cut = (x > 0) & (x < 3 * mean(x))
        m, s = mean_sigma(x[cut], err=False)
        return self.Draw.histo_2d(x[cut], y[cut], Bins.make(m - 4 * s, m + 4 * s, bin_size) + self.Bins.get_pad_ph(1), 'Pulse Height vs. ToT', y_tit='Pulse Height [mV]', x_tit='ToT [ns]', show=show)

    def draw_ph_tot(self, bin_size=.5, show=True):
        h = self.draw_signal_vs_tot(bin_size, show=False).ProfileX()
        format_histo(h, stats=0, y_tit='Pulse Height [mV]')
        self.Draw(h, show=show)

    def draw_tot_vs_peaktime(self, corr=True, fine_corr=False, bin_size=None, show=True):
        x, y = self.Run.get_tree_vec(var=self.Timing.get_peak_var(corr, fine_corr), cut=self.Cut()), self.Peaks.get_tot()
        cut = (y > 0) & (y < 3 * mean(y))
        self.Draw.profile(x[cut], y[cut], self.get_t_bins(bin_size), 'ToT vs. Peaktime', x_tit='Signal Peak Position [ns]', y_tit='ToT [ns]', stats=0, show=show)

    def draw_signal_vs_cft(self, bin_size=None, x=None, y=None, show=True):
        x, y = choose(x, self.Peaks.get_cft()), choose(y, self.get_ph_values())
        title = 'Signal vs {:.0f}% Constant Fraction Time'
        return self.Draw.profile(x, y, self.get_t_bins(bin_size), title, x_tit='Constant Fraction Time [ns]', y_tit='Pulse Height [mV]', show=show)

    def draw_signal_vs_times(self, bin_size=.2, show=True):
        x, y, zz = array(self.Peaks.get_from_tree()), self.Peaks.get_cft(), self.get_ph_values()
        self.Draw.prof2d(x, y, zz, self.get_t_bins(bin_size) * 2, 'Signal vs. CFD and Peak Time', y_tit='Constant Fraction Time [ns]', x_tit='Peak Time [ns]', z_tit='Pulse Height [mV]', show=show)

    def draw_ph_triggercell(self, bin_width=10, t_corr=True, cut=None, show=True):
        x, y = self.get_tree_vec(var=['trigger_cell', self.get_signal_var(t_corr=t_corr)], cut=self.Cut(cut))
        return self.Draw.profile(x, y, Bins.make(0, self.Run.NSamples, bin_width), 'Signal vs. Trigger Cell', x_tit='Trigger Cell', y_tit='Pulse Height [mV]', show=show)

    def draw_ph_pathlength(self, bin_size=.1, show=True):
        x, y = self.get_tree_vec(var=[self.get_track_length_var(), self.get_signal_var()], cut=self.Cut.generate_custom(exclude=['slope_x', 'slope_y']))
        self.Draw.profile(x, y, Bins.make(*array([0, 1]) + self.DUT.Thickness, bin_size), 'Pulse Height vs. Path Length', x_tit='Distance [#mum]', y_tit='Pulse Height [mV]', ndivx=405, show=show)

    def draw_ph_peakint(self, show=True):
        x = sum(array([ip for ip in self.Run.PeakIntegrals[self.DUT.Number - 1].values()]) / 2, axis=1)
        y = [self.get_pulse_height(sig=self.get_signal_name(peak_int=name)) for name in self.Run.PeakIntegrals[self.DUT.Number - 1]]
        self.Draw.graph(x, y, title='Signal vs. Peak Integral', x_tit='Integralwidth [ns]', y_tit='Pulse Height [mV]', show=show, x_range=ax_range(x, 0, .1, .1))
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region BUCKET
    def get_b2_var(self):
        return f'b2_integral[{self.DUT.Number - 1}]'

    def get_bucket_ph(self, cut=None):
        return self.get_tree_vec(self.get_b2_var(), cut=self.Cut(cut))

    def get_bucket_region(self, w=.5):
        """ width [w] in bunch spacings"""
        return round(self.Peaks.get_bunch_centre(2) + self.BunchSpacing * array([-w / 2, w / 2]), 2)

    def get_bucket_pulse_heights(self, redo=False):
        cuts = [self.Cut.generate_custom(exclude='bucket', name='nobucket', prnt=0), self.Cut.generate_custom(exclude='bucket', name='prebucket', prnt=0) + self.Cut.generate_pre_bucket()(), None]
        return [self.get_pulse_height(cut=cut, redo=redo) for cut in cuts]

    @update_pbar
    @save_pickle('Ratio', suf_args='[0, 1]', sub_dir='Bucket')
    def get_bucket_ratio(self, cut=None, all_cuts=False, _redo=False):
        """ return fraction of bucket events with no signal in the signal region. """
        default = self.Cut.generate_custom(include=['pulser', 'event range', 'beam stops', 'saturated'], prnt=False)
        c = choose(cut, self.Cut.generate_custom(exclude=['bucket', 'bucket2'], prnt=False) if all_cuts else default)
        n = self.get_n_entries(c, _redo=_redo)
        n_bucket = n - self.get_n_entries(self.Cut['bucket'] + c, _redo=_redo)
        return ufloat(n_bucket, sqrt(n_bucket)) / n

    @save_pickle('TPRatio', suf_args=0, sub_dir='Bucket')
    def get_bucket_tp_ratio(self, all_cuts=False, _redo=False):
        """ :returns the ratio of correctly identified bucket events with the trigger phase cut """
        x = self.get_tree_vec('trigger_phase', cut=self.Cut.get_bucket(all_cuts))
        return calc_eff(values=in1d(x, self.Cut.get_trigger_phases())) if x.size else [50] * 3

    def get_tp_ratio(self, redo=False):
        """ :returns the ratio of bucket events after the trigger phase cut. """
        return self.get_bucket_ratio(cut=self.Cut.get_tp(), _redo=redo)

    def calc_bucket_ratio(self):
        br = self.get_bucket_ratio(all_cuts=True)  # ratio of bucket events
        tpr = self.get_bucket_tp_ratio(all_cuts=True) / 100  # ratio of the bucket trigger phase
        return br * (1 - tpr)

    def estimate_bucket_ratio(self):
        if not self.Run.Config.has_option('BASIC', 'bucket scale') and not self.Run.Config.has_option('BASIC', 'bucket tpr'):
            warning('bucket parameters not defined -> applying no correction!')
            return 0
        return self.Run.Config.get_ufloat('BASIC', 'bucket scale') * self.get_flux() * (1 - self.Run.Config.get_ufloat('BASIC', 'bucket tpr'))

    def draw_bucket_ph(self, cut=None, fid=False, bin_width=2, logz=True, draw_cut=True, use_wf_int=False, redo=False, **kwargs):
        if use_wf_int:
            cut = choose(cut, self.get_event_cut(self.Cut.generate_custom(include=['pulser', 'ped sigma', 'event range'] + (['fiducial'] if fid else []), prnt=False, name=f'bph{fid}')))
            x, y = [self.Waveform.get_integrals(r, redo=redo)[cut] for r in [None, self.get_bucket_region()]]
        else:
            cut = choose(cut, self.Cut.generate_custom(include=['pulser', 'ped sigma', 'event range'] + (['fiducial'] if fid else []), prnt=False))
            x, y = self.get_ph_values(cut=cut), self.get_bucket_ph(cut=cut)
        h = self.Draw.histo_2d(x, y, Bins.get_pad_ph(bin_width) * 2, x_tit='Signal Pulse Height [mV]', y_tit='Bucket 2 Pulse Height [mV]', logz=logz, **kwargs)
        if draw_cut:
            m, s = self.Pedestal.get_under_signal()
            v = m.n + 3 * s.n
            lv, lh, b = Draw.vertical_line(v, color=2, w=2), Draw.horizontal_line(v, color=2, w=2), Draw.box(-100, v, v, 1000, line_color=2, opacity=.2, fillcolor=2)
            f = self.draw_b2_fit(color=1)
            Draw.fypolygon(f, -100, 600, 1000, fill_color=2, opacity=.2, line_color=1)
            Draw.legend([b, f], ['excluded', 'pedestal fit'], ['lf', 'l'], y2=.822)
            format_statbox(h, entries=True, w=.2)  # draw above cut

    def draw_b2_fit(self, draw_opt='same', **kwargs):
        fit, (m, s) = self.Cut.get_b2_fit(), self.Pedestal()
        return self.Draw.function(self.Draw.make_f(None, f'{fit[0].n} + {fit[1].n} * x + {fit[2].n} * pow(x, 2) + {m + 4 * s}', -50, 500, **kwargs), draw_opt=draw_opt)

    def draw_bucket_pedestal(self, bin_width=1, cut=None, **kwargs):
        self.Draw.distribution(self.get_bucket_ph(cut), Bins.get_pad_ph(bin_width), x_tit='Pulse Height [mV]', **kwargs)

    def draw_bucket_profile(self, cut=None, bin_width=2, logz=True, **kwargs):
        cut = self.get_event_cut(Cut.to_string(self.Cut() if cut is None else cut))
        x, y = self.get_tree_vec(self.get_raw_signal_var())[cut], self.get_tree_vec(self.get_b2_var())[cut]
        return self.Draw.profile(x, y, Bins.get_pad_ph(bin_width), x_tit='Signal Pulse Height [mV]', y_tit='Bucket 2 Pulse Height [mV]', logz=logz, **kwargs)

    def draw_bucket_fraction(self, redo=False, **kwargs):
        cuts = {key: value for key, value in list(self.Cut.ConsecutiveCuts.items()) if 'bucket' not in key}
        cuts['trigger phase'] = Cut.make(f'{len(cuts) + 1}', list(cuts.values())[-1] + self.Cut.generate_trigger_phase()())
        self.PBar.start(len(cuts), counter=True) if redo or not file_exists(self.make_simple_pickle_path('BucketRatio', len(cuts) - 1)) else do_nothing()
        y = array([self.get_bucket_ratio(cut, _redo=redo) for cut in cuts.values()])
        self.Draw.graph(arange(y.size), y * 100, y_tit='Fraction of Bucket Events [%]', bin_labels=cuts.keys(), **kwargs, y_range=[0, max(y).n * 110])

    def draw_bucket_map(self, **kwargs):
        cut = self.Cut.get('bucket', invert=True) + self.Cut.generate_custom(exclude=['bucket', 'bucket2', 'fiducial'])
        self.draw_hitmap(cut=cut, **kwargs)
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region SNR
    @update_pbar
    def calc_snr(self, region=None, peak_int=None, redo=False):
        return self.get_pulse_height(sig=self.get_signal_name(region=region, peak_int=peak_int), redo=redo) / self.Pedestal.get_noise(self.Pedestal.get_signal_name(peak_int=peak_int), redo=redo)

    def get_snrs(self, region=None, redo=False):
        peak_ints = self.Run.PeakIntegrals[self.DUT.Number - 1]
        self.PBar.start(len(peak_ints), counter=True, t='min')
        values = array([self.calc_snr(region, peak_int=name, redo=redo) for name in peak_ints])
        return values, array(list(peak_ints.values())) / 2

    def draw_snrs(self, show=True, redo=False):
        values, peak_ints = self.get_snrs(redo=redo)
        self.Draw.graph(arange(values.size), values, title='Signal To Noise Ratios', y_tit='SNR', show=show, x_range=ax_range(0, len(values) - 1, .1, .1))
        for i, ip in enumerate(peak_ints):
            Draw.tlatex(i, values[i].n + (max(values) - min(values)).n * .05, str(ip))

    def draw_2d_snrs(self, show=True, redo=False):
        values, peak_ints = self.get_snrs(redo=redo)
        x, y = peak_ints.T
        gStyle.SetPaintTextFormat('5.4g')
        return self.Draw.prof2d(x, y, values, Bins.make2d(x, y, off=-.25), 'Signal to Noise Ratios', x_tit='Left Width [ns]', y_tit='Right Width [ns]', draw_opt='coltext',
                                show=show, z_range=ax_range(values, fl=-.5, fh=.1), stats=False, rm=.03)

    def draw_snr_left(self, right=False):
        p = getattr(self.draw_2d_snrs(show=False), 'Profile{}'.format('Y' if right else 'X'))('{} Length'.format('Right' if right else 'Left'))
        format_histo(p, x_tit='Integral Length [ns]', y_tit='SNR', y_off=1.35, stats=0)
        self.Draw.histo(p, draw_opt='Pe')

    def draw_snr_right(self):
        self.draw_snr_left(right=True)

    @staticmethod
    def create_peak_ints(lmin=0, lmax=13, rmin=0, rmax=20, step=1):
        i = 1
        for lw in range(lmin, lmax, step):
            for rw in range(rmin, rmax, step):
                print(f'PeakIntegral{i}_range = [{lw},{rw}]')
                i += 1
    # endregion SNR
    # ----------------------------------------

    @staticmethod
    def sim_landau(n=1e6, m=80, s=5, cut_off=None):
        gRandom.Landau(m, s)
        v = array([gRandom.Landau(m, s) for _ in range(int(n))])
        v = v[v < cut_off] if cut_off else v
        while v.size < int(n):
            new = array([gRandom.Landau(m, s) for _ in range(int(n) - v.size)])
            v = append(v, new[new < cut_off])
        return v

    def draw_landau_stats(self, n=40, m=80, s=5, cut_off=None, use_quantile=False):
        x, y, = 2 ** arange(0, n // 2, .5), []
        self.PBar.start(n)
        for i in x:
            values = self.sim_landau(i, m, s, cut_off)
            y.append(quantile(values, .9) if use_quantile else mean_sigma(values)[0])
            self.PBar.update()
        tit = 'Landau Statistics{}'.format('' if cut_off is None else ' with cutoff {}'.format(cut_off))
        self.Draw.graph(x, y, title=tit, x_tit='Number of Events', y_tit='90 % Quantile' if use_quantile else 'Mean Pulse Height', logx=True)


if __name__ == '__main__':

    args = init_argparser(run=88, tc='201908', dut=1, tree=True, has_verbose=True)
    z = PadAnalysis(args.run, args.dut, args.testcampaign, args.tree, verbose=args.verbose)
