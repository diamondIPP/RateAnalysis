from ROOT import gRandom, TCut
from numpy import quantile, insert, sum

from src.pedestal import PedestalAnalysis
from src.timing import TimingAnalysis
from src.dut_analysis import *
from src.pad_cut import PadCut
from src.pad_run import PadRun
from src.peaks import PeakAnalysis
from src.pulser import PulserAnalysis
from src.waveform import Waveform
from src.mc_signal import MCSignal
from src.analysis import update_pbar


class PadAnalysis(DUTAnalysis):
    def __init__(self, run_number, diamond_nr, test_campaign=None, load_tree=True, verbose=False, prnt=True):

        DUTAnalysis.__init__(self, run_number, diamond_nr, test_campaign, load_tree, verbose, prnt)

        # MAIN
        self.Channel = self.get_channel()
        self.DigitiserBinWidth = .5 if self.Run.Digitiser == 'drs4' else .4

        if self.Tree.Hash():
            # Polarities
            self.Polarity = self.get_polarity()
            self.PulserPolarity = self.get_polarity(pulser=True)

            # Regions & Ranges
            self.SignalRegionName = self.load_region_name()
            self.SignalRegion = self.get_region()
            self.PeakIntegralName = self.load_peak_integral()
            self.PeakIntegral = self.Run.PeakIntegrals[self.DUT.Number - 1][self.PeakIntegralName]

            # Signal Names
            self.SignalName = self.get_signal_name()
            self.PedestalName = self.get_pedestal_name()
            self.PeakName = self.get_peak_name()

            # cuts
            self.Timing = TimingAnalysis(self)
            self.Pedestal = PedestalAnalysis(self)
            self.Waveform = Waveform(self)
            self.Peaks = PeakAnalysis(self)
            self.Pulser = PulserAnalysis(self)
            self.MC = MCSignal(self)

            # alignment
            self.IsAligned = self.check_alignment()

        self.print_finished(prnt=prnt)

    @staticmethod
    def get_info_header():
        return ['Run', 'Type', 'Diamond', 'Flux [kHz/cm2]', 'HV [V]', 'Region', 'Integral']

    def make_all(self, redo=False):
        self.draw_signal_distribution(redo=redo, show=False)
        self.draw_pulse_height(redo=redo, show=False)
        self.draw_signal_map(redo=redo, show=False)
        self.draw_hitmap(redo=redo, show=False)

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
        return rows[0] if ret_row else print_table(rows, self.get_info_header() if header else None, prnt=prnt)
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

    def check_alignment(self, redo=False):
        """ check if the events from telescope and digitiser are aligned"""
        def f():
            from src.pad_alignment import PadAlignment
            return PadAlignment(self.Run.Converter).IsAligned
        is_aligned = do_pickle(self.make_simple_pickle_path(sub_dir='Alignment'), f, redo=redo)
        warning('Run {r} is misaligned :-('.format(r=self.Run.Number)) if not is_aligned else do_nothing()
        return is_aligned
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_polarity(self, pulser=False):
        return int(self.Run.TreeConfig.get('General', '{}polarities'.format('pulser ' if pulser else '')).split()[self.get_channel()])

    def get_channel(self):
        return self.Run.Channels[self.DUT.Number - 1]

    def get_signal_var(self, signal=None, evnt_corr=True, off_corr=False, cut=None, region=None, peak_int=None, sig_type='signal', t_corr=True):
        sig_name = choose(signal, self.get_signal_name(region, peak_int, sig_type, t_corr))
        pol = self.get_polarity('pulser' in sig_type)
        if not any([off_corr, evnt_corr]):
            return '{} * {}'.format(pol, sig_name)
        return '{} * ({} - {})'.format(pol, sig_name, self.Pedestal.get_mean(cut=cut, raw=True).n if off_corr else self.get_pedestal_name())

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

    def get_short_regint(self, signal=None, signal_type='signal'):
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

    # @reload_tree
    def get_pulse_height(self, bin_size=None, cut=None, redo=False, corr=True, sig=None, sys_err=0, peaks=False):
        """ :returns: fitted (pol0) pulse height over time """
        def f():
            return self.draw_pulse_height(sig=sig, bin_size=bin_size, cut=self.Cut(cut), corr=corr, show=False, save=False, redo=redo)[1]
        suffix = '{}_{}_{}_{}'.format(choose(bin_size, Bins.get_size(bin_size)), int(corr), self.get_short_regint(sig), self.Cut(cut).GetName())
        ph = do_pickle(self.make_simple_pickle_path('Fit', suffix, sub_dir='Ph_fit'), f, redo=redo)[0]
        return self.Peaks.get_bunch_height() if peaks else ufloat(ph.n, ph.s + sys_err)

    def get_bucket_pulse_heights(self, redo=False):
        cuts = [self.Cut.generate_custom(exclude='bucket', name='nobucket', prnt=0), self.Cut.generate_custom(exclude='bucket', name='prebucket', prnt=0) + self.Cut.generate_pre_bucket()(), None]
        return [self.get_pulse_height(cut=cut, redo=redo) for cut in cuts]

    def get_bucket_ratio(self, fid=False, redo=False):
        def f():
            cut = self.Cut.generate_custom(include=['pulser', 'ped sigma', 'event range'] + (['fiducial' if fid else []]), prnt=False)
            return 1 - self.get_n_entries(cut + self.Cut['bucket']) / self.get_n_entries(cut)
        return do_pickle(self.make_simple_pickle_path('BucketRatio', int(fid)), f, redo=redo)

    def get_bucket_tp_ratio(self, fid=False, show=False, redo=False):
        def f():
            x = get_hist_vec(self.Tel.draw_trigger_phase(cut=self.Cut.get_bucket(fid), show=show))
            return max(x) / sum(x) if sum(x) else ufloat(.5, .5)
        return do_pickle(self.make_simple_pickle_path('BucketTPRatio', int(fid)), f, redo=redo or show)

    def get_min_signal(self, name=None):
        h = self.draw_signal_distribution(show=False, save=False, sig=name)
        return h.GetBinCenter(h.FindFirstBinAbove(h.GetMaximum() * .01))

    def get_t_bins(self, bin_size=None):
        xmin, xmax = self.SignalRegion * self.DigitiserBinWidth
        return Bins.make(xmin, xmax, choose(bin_size, default=self.DigitiserBinWidth))

    def get_alignment(self, bin_size=200, thresh=2):
        x, y = get_graph_vecs(self.Pulser.draw_hit_efficiency(bin_size, show=False), err=False)
        return x, y > thresh

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
        return self.Pulser.get_pedestal(par, redo) if pulser else self.Pedestal.get_par(par, redo=redo)

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

    def draw_pulse_height(self, bin_size=None, cut=None, y_range=None, redo=False, corr=True, sig=None, rel_t=True, show=True, save=True, prnt=True):
        def f():
            ph, t = self.get_tree_vec(var=[self.get_signal_var(sig, corr), self.get_t_var()], cut=self.Cut(cut))
            return self.Draw.profile(t, ph, self.Bins.get_time(bin_size, cut), 'Pulse Height Evolution', x_tit='Time [hh:mm]', y_tit='Mean Pulse Height [mV]', y_off=1.6, show=False)
        pickle_path = self.make_simple_pickle_path('', '{}{}_{}{}'.format(Bins.w(bin_size), int(corr), self.get_short_regint(sig), self.Cut(cut).GetName()), 'Ph_fit')
        p = do_pickle(pickle_path, f, redo=redo)
        y = get_hist_vec(p, err=False)
        y_range = ax_range(min(y), max(y), .5, .5) if y_range is None else y_range
        format_histo(p, name='Fit Result', x_range=[self.Run.StartTime, self.Bins.get_time()[1][-1]], t_ax_off=self.get_t_off(rel_t), y_range=y_range, ndivx=505)
        self.Draw(p, show=show, lm=.14, prnt=save)
        fit = self.fit_pulse_height(p, pickle_path)
        format_statbox(p, fit=True)
        self.Draw.save_plots('PulseHeight{}'.format(Bins.w(bin_size)), show=show, save=save, prnt=prnt)
        return p, fit

    def fit_pulse_height(self, p, picklepath):
        fit = p.Fit('pol0', 'qs', '', 0, self.__get_max_fit_pos(p))
        SaveDraw.server_pickle(picklepath, FitRes(fit))
        return FitRes(fit)

    def __get_max_fit_pos(self, h):
        """ look for huge fluctiations in ph graph and return last stable point"""
        if mean([h.GetBinContent(i) for i in range(h.GetNbinsX())]) < 10:  # if the pulse height is very low there will be always big fluctuations!
            return h.GetBinCenter(h.GetNbinsX()) + 1000
        sum_ph = h.GetBinContent(1)
        for i in range(2, h.GetNbinsX() + 1):
            sum_ph += h.GetBinContent(i)
            if h.GetBinContent(i) < .8 * sum_ph / (i + 1):
                if not h.GetBinEntries(i):
                    continue  # if the bin is empty
                warning('Found PH fluctiation in run {}! Stopping fit after {:2.0f}%'.format(self.Run.Number, 100. * (i - 1) / h.GetNbinsX()))
                return h.GetBinCenter(i - 1)
        return h.GetBinCenter(h.GetNbinsX()) + 1000

    def draw_signal_distribution(self, cut=None, evnt_corr=True, off_corr=False, show=True, sig=None, bin_width=None, events=None,
                                 start=0, x_range=None, redo=False, prnt=True, save=True, normalise=None):
        def func():
            self.info('Drawing signal distribution for run {} and {}...'.format(self.Run.Number, self.DUT.Name), prnt=prnt)
            nentries = self.Run.NEvents if events is None else self.Run.find_n_events(n=events, cut=str(cut), start=start)
            values = self.get_tree_vec(var=self.get_signal_var(sig, evnt_corr, off_corr, cut), cut=self.Cut(cut), nentries=nentries, firstentry=start)
            return self.Draw.distribution(values, self.Bins.get_pad_ph(bin_width, mean(values)), show=False, x_tit='Pulse Height [mV]', y_off=2)

        suffix = f'{bin_width}_{evnt_corr:d}_{self.Cut(cut).GetName()}_{self.get_short_regint(sig)}'
        h = do_pickle(self.make_simple_pickle_path('Histo', suffix, 'PulseHeight'), func, redo=redo)
        format_histo(h, x_range=choose(x_range, ax_range(0, 3, .1, h=h)), normalise=normalise)
        self.Draw(h, 'SignalDistribution', lm=.15, show=show, prnt=prnt, save=save, stats=None)
        return h

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

    def draw_bucket_ph(self, cut=None, fid=False, bin_width=2, logz=True, draw_cut=True, redo=False, **kwargs):
        cut = choose(cut, self.get_event_cut(self.Cut.generate_custom(include=['pulser', 'ped sigma', 'event range'] + (['fiducial'] if fid else []), prnt=False, name=f'bph{fid}')))
        x, y = [self.Waveform.get_integrals(r, redo=redo)[cut] for r in [None, self.SignalRegion * self.DigitiserBinWidth + self.BunchSpacing]]
        self.Draw.histo_2d(x, y, Bins.get_pad_ph(bin_width) * 2, x_tit='Signal Pulse Height [mV]', y_tit='Bucket 2 Pulse Height [mV]', logz=logz, **kwargs)
        if draw_cut:
            m, s = self.Pedestal.get_under_signal()
            v = m.n + 3 * s.n
            lv, lh, b = Draw.vertical_line(v, color=2, w=2), Draw.horizontal_line(v, color=2, w=2), Draw.box(-100, v, v, 1000, line_color=2, opacity=.2, fillcolor=2)
            Draw.legend([b], ['excluded'], y2=.822)
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
