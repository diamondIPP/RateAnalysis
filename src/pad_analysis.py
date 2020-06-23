from __future__ import print_function

from dut_analysis import *
from json import loads
from ROOT import gRandom, TProfile2D, THStack, Double, Long

from pad_cut import PadCut
from peaks import PeakAnalysis
from Pedestal import PedestalAnalysis
from pulser import PulserAnalysis
from Timing import TimingAnalysis
from converter import Converter
from waveform import Waveform
from pad_alignment import PadAlignment
from numpy import ones


class PadAnalysis(DUTAnalysis):
    def __init__(self, run_number, diamond_nr, test_campaign=None, tree=True, t_vec=None, verbose=False, prnt=True):

        DUTAnalysis.__init__(self, run_number, diamond_nr, test_campaign, tree, t_vec, verbose, prnt)

        # Main
        self.Channel = self.Run.Channels[diamond_nr - 1]
        self.DigitiserBinWidth = .5 if self.Run.Digitiser == 'drs4' else .4

        if self.Tree:
            self.Converter = Converter(self.Run)
            # Polarities
            self.Polarity = self.get_polarity()
            self.PulserPolarity = self.get_pulser_polarity()

            # Regions & Ranges
            self.IntegralNames = self.get_integral_names()
            self.IntegralRegions = self.load_regions()
            self.SignalRegionName = self.IntegralRegions['signal']
            self.SignalRegion = array(self.Run.IntegralRegions[self.DUT.Number - 1][self.SignalRegionName])
            self.PedestalRegion = array(self.IntegralRegions['pedestal'])
            self.PeakIntegralName = self.load_peak_integral()
            self.PeakIntegral = self.Run.PeakIntegrals[self.DUT.Number - 1][self.PeakIntegralName]

            # Signal Names
            self.SignalDefinition = '({pol}*TimeIntegralValues[{num}])'
            self.SignalNumber = self.get_signal_number()
            self.SignalName = self.get_signal_name()
            self.PedestalName = self.get_pedestal_name()
            self.PeakName = self.get_peak_name()

            # cuts
            self.Timing = TimingAnalysis(self)
            self.Cut = PadCut.from_parent(self.Cut)

            # subclasses
            self.Pedestal = PedestalAnalysis(self)
            self.Pulser = PulserAnalysis(self)
            self.Waveform = Waveform(self)
            self.Peaks = PeakAnalysis(self)

            # alignment
            self.IsAligned = self.check_alignment()

            self.Timing.reload_cut()

        self.print_finished(prnt=prnt)

    def draw_timing(self):
        self.Timing.draw_all()

    def draw_pulser_rate(self, show=True, prnt=True):
        self.Pulser.draw_rate(show=show, prnt=prnt)

    @staticmethod
    def get_info_header():
        return ['Run', 'Type', 'Diamond', 'Flux [kHz/cm2]', 'HV [V]', 'Region', 'Integral']

    # ----------------------------------------
    # region INIT
    def update_config(self):
        self.Config.read(join(self.Dir, 'config', self.TCString, 'PadConfig.ini'))

    def get_short_regint(self, signal=None, signal_type='signal'):
        return self.get_all_signal_names(signal_type)[signal if signal is not None else self.SignalName]

    def get_integral_names(self):
        if self.Run.TreeConfig.has_section('Integral Names'):
            return [str(name) for name in loads(self.Run.TreeConfig.get('Integral Names', 'Names'))]
        self.Tree.GetEntry(0)
        return [str(name) for name in self.Tree.IntegralNames]

    def get_polarity(self):
        if not self.Run.TreeConfig.has_option('General', 'polarities'):
            warning('OLD DATA! Take polarities from config...')
            return loads(self.Run.Converter.load_polarities())[self.Channel]
        return int(self.Run.TreeConfig.get('General', 'polarities').split()[self.Channel])

    def get_pulser_polarity(self):
        if not self.Run.TreeConfig.has_option('General', 'pulser polarities'):
            warning('OLD DATA! Take pulser polarities from config...')
            return loads(self.Run.Converter.load_polarities(pulser=True))[self.Channel]
        return int(self.Run.TreeConfig.get('General', 'pulser polarities').split()[self.Channel])

    def load_regions(self):
        all_regions = {}
        for name in ['signal', 'pedestal', 'pulser']:
            option = '{} region'.format(name)
            region = '{name}_{region}'.format(name=name, region=self.Config.get('SIGNAL', option)) if option in self.Config.options('SIGNAL') else ''
            regions = [reg for reg in self.Run.IntegralRegions[self.DUT.Number - 1] if reg.startswith(name)]
            all_regions[name] = region if region in regions else regions[0]
        return all_regions

    def load_peak_integral(self):
        peak_int = 'PeakIntegral{}'.format(self.Config.get('SIGNAL', 'peak integral'))
        return peak_int if peak_int in self.Run.PeakIntegrals[self.DUT.Number - 1] else self.Run.PeakIntegrals[self.DUT.Number - 1].keys()[0]

    def get_signal_number(self, region=None, peak_integral=None, sig_type='signal'):
        region = self.IntegralRegions[sig_type] if region is None else self.make_region(sig_type, region)
        peak_integral = self.PeakIntegralName if peak_integral is None else 'PeakIntegral{}'.format(peak_integral)
        int_name = 'ch{ch}_{reg}_{int}'.format(ch=self.Channel, reg=region, int=peak_integral)
        return self.IntegralNames.index(int_name)

    def get_signal_name(self, region=None, peak_integral=None, sig_type='signal'):
        num = self.get_signal_number(region, peak_integral, sig_type)
        return self.SignalDefinition.format(pol=self.Polarity, num=num)

    def get_signal_region(self, name=None):
        return self.Run.IntegralRegions[self.DUT.Number - 1][self.SignalRegionName if name is None else 'signal_{}'.format(name)]

    def get_signal_range(self, lfac=0, rfac=.3):
        return increased_range(self.SignalRegion * self.DigitiserBinWidth, lfac, rfac)

    def set_signal_definitions(self, use_time=True, sig_region=None, peak_int=None):
        signal = 'TimeIntegralValues' if use_time else 'IntegralValues'
        signal = '({{pol}}*{sig}[{{num}}])'.format(sig=signal)
        print('changed SignalDefinition to:', signal)
        self.SignalDefinition = signal
        self.update_signal_definitions(sig_region, peak_int)

    def update_signal_definitions(self, sig_region=None, peak_int=None):
        self.SignalNumber = self.get_signal_number(sig_region, peak_int)
        self.SignalName = self.get_signal_name(sig_region, peak_int)

    def get_pedestal_name(self, region=None, peak_int=None):
        return self.get_signal_name(region=region, peak_integral=peak_int, sig_type='pedestal')

    def get_peak_name(self, region=None, type_='signal', t_corr=True):
        peak_name = 'IntegralPeakTime' if t_corr else 'IntegralPeaks'
        return '{name}[{num}]'.format(name=peak_name, num=self.get_signal_number(region, None, type_))
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_ph_str(self):
        return self.generate_signal_name()

    def get_attenuator(self):
        return self.DUT.Attenuator

    def get_ph_data(self, cut=None):
        """ :return: pulse height data as numpy array [[time] [ph]] with units [[s], [mV]]
            :param cut: applies all cuts if None is provided."""
        n = self.Tree.Draw('{}:{}'.format(self.get_t_var(), self.generate_signal_name()), self.Cut(cut), 'goff')
        return self.Run.get_root_vecs(n, 2)

    def get_ph_values(self, region=None, name=None, cut=None):
        """ :returns: all pulse height values for a given cut. """
        return self.Run.get_root_vec(var=self.generate_signal_name(name, cut=cut, region=region), cut=self.Cut(cut))

    def get_pulse_height(self, bin_size=None, cut=None, redo=False, corr=True, sig=None, sys_err=0, peaks=False):
        if peaks:
            return self.Peaks.get_signal_ph()
        correction = '' if not corr else '_eventwise'
        suffix = '{bins}{cor}_{c}'.format(bins=self.Bins.BinSize if bin_size is None else bin_size, cor=correction, reg=self.get_short_regint(sig), c=self.Cut(cut).GetName())
        picklepath = self.make_simple_pickle_path('Fit', suffix, sub_dir='Ph_fit')

        def f():
            p, fit_pars = self.draw_pulse_height(bin_size=bin_size, cut=self.Cut(cut), corr=corr, show=False, save=False, redo=redo)
            return fit_pars

        ph = make_ufloat(do_pickle(picklepath, f, redo=redo), par=0)
        return ufloat(ph.n, ph.s + sys_err)

    def get_pedestal(self, pulser=False, par=1, redo=False):
        return self.Pulser.get_pedestal(par, redo) if pulser else self.Pedestal.get_par(par, redo=redo)

    def get_results(self):
        return [self.get_pulse_height(), self.Pedestal.get_mean(), self.Pulser.get_pulse_height()]

    def get_peak_timing(self, par=1, redo=False):
        return self.Timing.get(par, redo)

    def get_min_signal(self, name=None):
        h = self.draw_signal_distribution(show=False, save=False, sig=name)
        return h.GetBinCenter(h.FindFirstBinAbove(h.GetMaximum() * .01))

    def get_additional_peak_height(self):
        return self.Peaks.find_additional(scale=True, show=False)

    def get_peak_flux(self, prnt=True, redo=False):
        return self.Peaks.get_flux(prnt=prnt, redo=redo)

    def get_n_peaks(self, start_bunch=None, end_bunch=None):
        return self.Peaks.get_n_additional(start_bunch, end_bunch)

    def get_t_bins(self, bin_size=None):
        xmin, xmax = self.SignalRegion * self.DigitiserBinWidth
        return make_bins(xmin, xmax, choose(bin_size, default=self.DigitiserBinWidth))

    def print_results(self, prnt=True):
        rows = [[u_to_str(v, prec=2) for v in self.get_results()]]
        print_table(header=['Signal [mV]', 'Pedestal [mV]', 'Pulser [mV]'], rows=rows, prnt=prnt)
    # endregion GET
    # ----------------------------------------

    def make_all(self, redo=False):
        self.draw_signal_distribution(redo=redo, show=False)
        self.draw_pulse_height(redo=redo, show=False)
        self.draw_signal_map(redo=redo, show=False)
        self.draw_hitmap(redo=redo, show=False)

    # ----------------------------------------
    # region 2D SIGNAL DISTRIBUTION
    def draw_efficiency_map(self, res=None, cut='all', show=True):
        cut_string = TCut(cut) + self.Cut.CutStrings['tracks']
        cut_string = self.Cut.generate_custom(exclude=['fiducial']) if cut == 'all' else cut_string
        p = TProfile2D('p_em', 'Efficiency Map {d}'.format(d=self.DUT.Name), *self.Bins.get_global(res, mm=True))
        y, x = self.Cut.get_track_vars(self.DUT.Number - 1, mm=True)
        thresh = self.Pedestal.get_mean() * 4
        self.Tree.Draw('({s}>{t})*100:{y}:{x}>>p_em'.format(s=self.generate_signal_name(), x=x, y=y, t=thresh), cut_string, 'goff')
        self.format_statbox(entries=True, x=.81)
        set_2d_ranges(p, dx=3, dy=3)
        format_histo(p, x_tit='Track x [cm]', y_tit='Track y [cm]', z_tit='Efficiency [%]', y_off=1.4, z_off=1.5, ncont=100, z_range=[0, 100])
        self.draw_histo(p, show=show, lm=.13, rm=.17, draw_opt='colz', x=1.15 if self.Title else 1)
        self.draw_fid_cut(scale=10)
        self.draw_detector_size()
        self.save_plots('EffMap')

    def draw_efficiency_vs_threshold(self, thresh=None, bin_width=.5, show=True):
        thresh = self.Pedestal.get_noise().n * 4 if thresh is None else thresh
        h = self.draw_signal_distribution(show=False, bin_width=bin_width)
        bins = range(h.FindBin(0), h.FindBin(thresh) + 2)
        err = Double()
        efficiencies = [make_ufloat(h.IntegralAndError(ibin, h.GetNbinsX(), err), err) / h.Integral() * 100 for ibin in bins]
        g = self.make_tgrapherrors('get', 'Detector Efficiency', x=[h.GetBinCenter(ibin) for ibin in bins], y=efficiencies)
        format_histo(g, x_tit='Threshold [mV]', y_tit='Efficiency [%]', y_off=1.3)
        self.draw_histo(g, draw_opt='ap', lm=.12, show=show)
        self.draw_vertical_line(x=self.Pedestal.get_noise().n * 3, ymin=0, ymax=110, w=2)
        self.draw_tlatex(x=self.Pedestal.get_noise().n * 3, y=95, text=' 3 #times noise', align=10)
        self.save_plots('EffThresh')

    def draw_pedestal_map(self, high=10, low=None, fid=False, cut=None, hit_map=True):
        low = '&&{}>{}'.format(self.generate_signal_name(), low) if low is not None else ''
        fid_cut = self.Cut.generate_custom(exclude='fiducial') if cut is None else self.Cut(cut)
        kwargs = {'redo': True, 'cut': TCut('{}<{}{}'.format(self.generate_signal_name(), high, low)) + (self.Cut(cut) if fid else fid_cut)}
        self.draw_hitmap(**kwargs) if hit_map else self.draw_signal_map(**kwargs)

    def find_center(self, redo=False):
        def f():
            h = self.draw_signal_map(cut='', show=False)
            px, py = h.ProjectionX(), h.ProjectionY()
            mx = mean([px.GetBinCenter(b) for b in [px.FindFirstBinAbove(px.GetMaximum() / 2), px.FindLastBinAbove(px.GetMaximum() / 2)]])
            my = mean([py.GetBinCenter(b) for b in [py.FindFirstBinAbove(py.GetMaximum() / 2), py.FindLastBinAbove(py.GetMaximum() / 2)]])
            return array([mx, my])
        return do_pickle(self.make_simple_pickle_path('Center', sub_dir='Center'), f, redo=redo)

    # endregion 2D SIGNAL DISTRIBUTION
    # ----------------------------------------

    # ----------------------------------------
    # region PULSE HEIGHT
    def generate_signal_name(self, signal=None, evnt_corr=True, off_corr=False, cut=None, region=None):
        sig_name = signal if signal is not None else self.get_signal_name(region)
        # pedestal polarity is always the same as signal polarity
        ped_pol = '1'
        # change polarity if pulser has opposite polarity to signal
        if hasattr(self, 'Pulser') and signal == self.Pulser.SignalName:
            ped_pol = '-1' if self.PulserPolarity != self.Polarity else ped_pol
        if off_corr:
            sig_name += '-{pol}*{ped}'.format(ped=self.Pedestal.get_mean(cut=cut).n, pol=ped_pol)
        elif evnt_corr:
            sig_name += '-{pol}*{ped}'.format(ped=self.PedestalName, pol=ped_pol)
        return sig_name

    def make_signal_time_histos(self, bin_width=.2, signal_name=None, evnt_corr=False, off_corr=False, bin_corr=False, rel_t=False, show=True):
        signal_name = self.generate_signal_name(self.SignalName if signal_name is None else signal_name, evnt_corr, off_corr, bin_corr)
        h = TH2F('h_st', 'Signal vs. Time', *(self.Bins.get_time() + self.Bins.get_pad_ph(bin_width)))
        self.format_statbox(entries=True, x=.83)
        self.Tree.Draw('{}:{}>>h_st'.format(signal_name, self.get_t_var()), self.Cut(), 'goff')
        format_histo(h, x_tit='Time [min]', y_tit='Pulse Height [au]', y_off=1.4, t_ax_off=self.Run.StartTime if rel_t else 0, pal=53)
        self.save_histo(h, 'SignalTime', show, lm=.12, draw_opt='colz', rm=.15)
        return h

    def draw_pulse_height_vs_binsize(self, show=True):
        bin_sizes = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000]
        # pulse_heights = [make_ufloat(self.draw_pulse_height(bin_width=bin_size, show=False)[1], par=0) for bin_size in bin_sizes]
        pulse_heights = [make_ufloat(self.draw_ph(bin_size=bin_size, show=False)[1], par=0) for bin_size in bin_sizes]
        g = self.make_tgrapherrors('gdbs', 'Pulse Height vs Number of Events per Bin', x=bin_sizes, y=pulse_heights)
        format_histo(g, x_tit='Number of Events per Bin', y_tit='Pulse Height [mV]', y_off=1.2)
        self.draw_histo(g, lm=.12, show=show, gridy=True, logx=True)

    def draw_pulse_height(self, bin_size=None, cut=None, y_range=None, redo=False, corr=True, sig=None, rel_t=True, show=True, save=True, prnt=True):

        # TODO fix errors or extract from mean

        sig = self.SignalName if sig is None else sig
        correction = '' if not corr else '_eventwise'
        cut_str = self.Cut(cut)
        bin_size = self.Bins.BinSize if bin_size is None else bin_size
        suffix = '{bins}{cor}_{reg}{c}'.format(bins=bin_size, cor=correction, reg=self.get_short_regint(sig), c=cut_str.GetName())
        picklepath = self.make_pickle_path('Ph_fit', None, self.Run.Number, self.DUT.Number, suf=suffix)

        def func():
            signal = self.generate_signal_name(sig, corr)
            prof = TProfile('pph', 'Pulse Height Evolution', *self.Bins.get_time(bin_size))
            self.Tree.Draw('{}:{}>>pph'.format(signal, self.get_t_var()), cut_str, 'goff')
            self.PulseHeight = prof
            return prof

        p = do_pickle(picklepath, func, redo=redo)
        self.format_statbox(only_fit=True, w=.3)
        y = get_hist_vec(p)
        format_histo(p, name='Fit Result', x_tit='Time [hh:mm]', y_tit='Mean Pulse Height [mV]', y_off=1.6, x_range=[self.Run.StartTime, self.Bins.get_time()[1][-1]],
                     t_ax_off=self.Run.StartTime if rel_t else 0, y_range=increased_range([y.min().n, y.max().n], .5, .5) if y_range is None else y_range, ndivx=505)
        self.draw_histo(p, show=show, lm=.14, prnt=save)
        fit = self.fit_pulse_height(p, picklepath)
        self.save_plots('PulseHeight{}'.format(bin_size), show=show, save=save, prnt=prnt)
        return p, fit

    def fit_pulse_height(self, p, picklepath):
        fit = p.Fit('pol0', 'qs', '', 0, self.__get_max_fit_pos(p))
        self.server_pickle(picklepath, FitRes(fit))
        return FitRes(fit)

    def __get_max_fit_pos(self, h):
        """ look for huge fluctiations in ph graph and return last stable point"""
        if mean([h.GetBinContent(i) for i in xrange(h.GetNbinsX())]) < 10:  # if the pulse height is very low there will be always big fluctuations!
            return h.GetBinCenter(h.GetNbinsX()) + 1000
        sum_ph = h.GetBinContent(1)
        for i in xrange(2, h.GetNbinsX() + 1):
            sum_ph += h.GetBinContent(i)
            if h.GetBinContent(i) < .8 * sum_ph / (i + 1):
                if not h.GetBinEntries(i):
                    continue  # if the bin is empty
                log_warning('Found PH fluctiation in run {}! Stopping fit after {:2.0f}%'.format(self.Run.Number, 100. * (i - 1) / h.GetNbinsX()))
                return h.GetBinCenter(i - 1)
        return h.GetBinCenter(h.GetNbinsX()) + 1000

    def draw_ph(self, bin_size=10000, y_range=None, rel_t=False, show=True):
        """ get pulse height by fitting every time bin disto with a Landau and then extrapolate with a pol0 """
        gr = self.make_tgrapherrors('hphl', 'Pulser Height Evolution')
        h = TH2F('tempph', '', *(self.Bins.get_time(bin_size) + self.Bins.get_pad_ph(bin_width=20)))
        self.Tree.Draw('{}:{}>>tempph'.format(self.SignalName, self.get_t_var()), self.Cut(), 'goff')
        i = 0
        for xbin in xrange(2, h.GetNbinsX() + 1):  # first bin is always empty
            py = h.ProjectionY('_py{}'.format(xbin), xbin, xbin)
            self.draw_histo(py)
            try:
                fit = self.fit_langau(py, nconv=50, show=True)
                raw_input()
                if fit.ParErrors[1] < .5:
                    continue
                gr.SetPoint(i, h.GetXaxis().GetBinCenter(xbin), fit.Parameters[1])
                gr.SetPointError(i, h.GetXaxis().GetBinWidth(xbin) / 2., fit.ParErrors[1])
                i += 1
            except ZeroDivisionError:
                pass
        self.format_statbox(only_fit=True)
        y_vals = [gr.GetY()[i] for i in xrange(gr.GetN())]
        format_histo(gr, x_tit='Time [min]', y_tit='Mean Pulse Height [au]', y_off=1.6, x_range=[self.Run.StartTime, self.Bins.get_time()[1][-1]],
                     t_ax_off=self.Run.StartTime if rel_t else 0, y_range=increased_range([min(y_vals), max(y_vals)], .5, .5) if y_range is None else y_range, ndivx=505)
        fit = gr.Fit('pol0', 'qs')
        self.draw_histo(gr, draw_opt='ap', show=show)
        return gr, FitRes(fit)

    def show_ph_overview(self, binning=None):
        self.draw_pulse_height(bin_size=binning, show=False)
        h1 = self.draw_pulse_height(show=False)[0]
        format_histo(h1, y_off=1.4)
        h2 = self.draw_ph_pull(event_bin_width=binning, show=False)
        print(h1, h2)
        c = TCanvas('c', 'Pulse Height Distribution', 1500, 750)
        c.Divide(2, 1)
        for i, h in enumerate([h1, h2], 1):
            pad = c.cd(i)
            pad.SetBottomMargin(.15)
            h.Draw()
        self.save_plots('PHEvolutionOverview{0}'.format(self.Bins.BinSize))
        self.add(h2, c)

    def draw_signal_distribution(self, cut=None, evnt_corr=True, off_corr=False, show=True, sig=None, bin_width=.5, events=None,
                                 start=None, x_range=None, redo=False, prnt=True, save=True, normalise=None, sumw2=False):
        cut = self.Cut(cut)
        suffix = '{b}_{c}_{cut}_{n}'.format(b=bin_width, c=int(evnt_corr), cut=cut.GetName(), n=self.get_short_regint(sig))
        pickle_path = self.make_pickle_path('PulseHeight', 'Histo', run=self.Run.Number, ch=self.DUT.Number, suf=suffix)

        def func():
            self.info('Drawing signal distribution for run {run} and {dia}...'.format(run=self.Run.Number, dia=self.DUT.Name), prnt=prnt)
            set_root_output(False)
            h1 = TH1F('h_sd', 'Pulse Height {s}'.format(s='with Pedestal Correction' if evnt_corr else ''), *self.Bins.get_pad_ph(bin_width))
            sig_name = self.generate_signal_name(sig, evnt_corr, off_corr, cut)
            start_event = int(float(start)) if start is not None else 0
            n_events = self.Run.find_n_events(n=events, cut=str(cut), start=start_event) if events is not None else self.Run.NEntries
            self.Tree.Draw('{name}>>h_sd'.format(name=sig_name), str(cut), 'goff', n_events, start_event)
            h1.Rebin(max(1, int(h1.GetMean() / 30)))
            return h1

        self.format_statbox(all_stat=1, w=.3)
        h = do_pickle(pickle_path, func, redo=redo)
        x_range = increased_range([h.GetBinCenter(i) for i in [h.FindFirstBinAbove(0), h.FindLastBinAbove(3)]], .1) if x_range is None else x_range
        format_histo(h, x_tit='Pulse Height [mV]', y_tit='Number of Entries', y_off=2, fill_color=self.FillColor, x_range=x_range, normalise=normalise)
        self.save_histo(h, 'SignalDistribution', lm=.15, show=show, prnt=prnt, save=save, sumw2=sumw2)
        return h

    def draw_signal_vs_peaktime(self, region=None, cut=None, bin_size=None, fine_corr=False, prof=True, x=None, y=None, y_range=None, xbins=None, show=True):
        h = (TProfile if prof else TH2F)('hspt', 'Signal vs Peak Position', *(choose(xbins, self.get_t_bins(bin_size)) + ([] if prof else self.Bins.get_pad_ph())))
        n = 0
        if x is None or y is None:
            n = self.Tree.Draw('{}:{}'.format(self.generate_signal_name(), self.Timing.get_peak_name(corr=True, fine_corr=fine_corr, region=region)), self.Cut(cut), 'goff')
        fill_hist(h, choose(x, self.Run.get_root_vec, n=n, ind=1), choose(y, self.Run.get_root_vec, n=n, ind=0))
        if not prof:
            px = h.ProjectionX()
            for xbin in range(h.GetNbinsX()):
                for ybin in range(h.GetNbinsY()):
                    h.SetBinContent(xbin, ybin, h.GetBinContent(xbin, ybin) / (px.GetBinContent(xbin) if px.GetBinContent(xbin) else 1))
        format_histo(h, x_tit='Signal Peak Position [ns]', y_tit='Pulse Height [mV]', y_off=1.4, stats=0, y_range=y_range)
        self.save_histo(h, 'SignalVsPeakPos', show, lm=.11, draw_opt='' if prof else 'colz', rm=.03 if prof else .18)
        return h

    def draw_signal_vs_tot(self, show=True):
        values = self.Peaks.get_all_tot()
        x = self.Run.get_root_vec(var=self.generate_signal_name(), cut=self.Cut())
        v = values[(values > 0) & (values < 3 * mean(values))]
        x = x[(values > 0) & (values < 3 * mean(values))]
        m, s = mean_sigma(v)
        h = TH2F('hstot', 'Pulse Height vs. Time over Threshold', 200, *(increased_range([m - 3 * s, m + 3 * s], .3, .3) + self.Bins.get_pad_ph()))
        h.FillN(x.size, v.astype('d'), x, ones(x.size))
        self.format_statbox(entries=True, x=.83)
        format_histo(h, y_tit='Pulse Height [mV]', x_tit='Time over Threshold [ns]', y_off=1.4, z_tit='Number of Entries', z_off=1.1)
        self.draw_histo(h, show=show, lm=.11, draw_opt='colz', rm=.15)

    def draw_tot_vs_peaktime(self, corr=True, fine_corr=False, bin_size=None, show=True):
        values = self.Peaks.get_all_tot()
        x = self.Run.get_root_vec(var=self.Timing.get_peak_name(corr, fine_corr), cut=self.Cut())
        v = values[(values > 0) & (values < 3 * mean(values))]
        x = x[(values > 0) & (values < 3 * mean(values))]
        p = TProfile('ptotpt', 'ToT vs. Peaktime', *self.get_t_bins(bin_size))
        p.FillN(x.size, x, v.astype('d'), ones(x.size))
        format_histo(p, x_tit='Signal Peak Position [ns]', y_tit='Time over Threshold [ns]', y_off=1.4, stats=0)
        self.save_histo(p, 'ToTVsPeakPos{}{}'.format(int(corr), int(fine_corr)), show, lm=.11)

    def draw_signal_vs_cft(self, bin_size=None, thresh=.5, x=None, y=None, show=True):
        p = TProfile('pscfd', 'Signal vs {:.0f}% Constant Fraction Time'.format(thresh * 100), *self.get_t_bins(bin_size))
        self.format_statbox(all_stat=True)
        fill_hist(p, choose(x, self.Peaks.get_all_cft(thresh)), choose(y, self.get_ph_values()))
        format_histo(p, x_tit='Constant Fraction Time [ns]', y_tit='Pulse Height [mV]', y_off=1.4)
        self.draw_histo(p, show=show, lm=.12)
        return p

    def draw_signal_vs_times(self, bin_size=.2, show=True):
        p = TProfile2D('hscfp', 'Signal vs. CFD and Peak Time', *(self.get_t_bins(bin_size) + self.get_t_bins(bin_size)))
        x = array(self.Peaks.get_from_tree())
        y = self.Peaks.get_all_cft()
        zz = self.Run.get_root_vec(var=self.generate_signal_name(), cut=self.Cut())
        for i in range(x.size):
            p.Fill(x[i], y[i], zz[i])
        format_histo(p, y_tit='Constant Fraction Time [ns]', x_tit='Peak Time [ns]', z_tit='Pulse Height [mV]', y_off=1.4, stats=0, z_off=1.2)
        self.format_statbox(entries=True)
        self.draw_histo(p, show=show, lm=.12, rm=.13, draw_opt='colz')

    def draw_signal_vs_triggercell(self, bin_width=10, cut=None, show=True):
        p = TProfile('pstc', 'Signal vs. Trigger Cell', self.Run.NSamples / bin_width, 0, self.Run.NSamples)
        self.Tree.Draw('{}:trigger_cell>>pstc'.format(self.generate_signal_name()), self.Cut(cut), 'goff')
        format_histo(p, x_tit='Trigger Cell', y_tit='Pulse Height [au]', y_off=1.2, stats=0)
        self.save_histo(p, 'SignalVsTriggerCell', show, lm=.11)
    # endregion PULSE HEIGHT
    # ----------------------------------------

    # ----------------------------------------
    # region CUTS
    def show_bucket_histos(self):
        h = TH1F('h', 'Bucket Cut Histograms', 350, -50, 300)
        self.Tree.Draw('{name}>>h'.format(name=self.SignalName), self.Cut.get_bucket(), 'goff')
        self.draw_histo(h)
        h1 = deepcopy(h)
        fit = fit_bucket(h1, show=False)
        sig_fit = TF1('f1', 'gaus', -50, 300)
        sig_fit.SetParameters(fit.GetParameters())
        ped1_fit = TF1('f2', 'gaus', -50, 300)
        ped2_fit = TF1('f2', 'gaus', -50, 300)
        ped1_fit.SetParameters(*[fit.GetParameter(i) for i in xrange(3, 6)])
        ped2_fit.SetParameters(*[fit.GetParameter(i) for i in xrange(6, 9)])
        h_sig = deepcopy(h)
        h_ped1 = deepcopy(h)
        h_ped2 = deepcopy(h)
        h_sig.Add(ped1_fit, -1)
        h_sig.Add(ped2_fit, -1)
        h_ped1.Add(ped2_fit, -1)
        h_ped2.Add(ped1_fit, -1)
        h_ped1.Add(h_sig, -1)
        h_ped2.Add(h_sig, -1)
        c = TCanvas('c', 'Bucket Histos', 1000, 1000)
        for i, h in enumerate([h_ped1, h_ped2, h_sig]):
            h.SetStats(0)
            h.SetLineColor(self.get_color())
            h.SetLineWidth(2)
            h.Draw('same') if i else h.Draw()
        self.save_plots('BucketHistos')
        self.add(h, h_sig, h_ped1, h_ped2, c)

    def show_bucket_numbers(self, show=True):
        pickle_path = self.PickleDir + 'Cuts/BucketEvents_{tc}_{run}_{dia}.pickle'.format(tc=self.TestCampaign, run=self.Run.Number, dia=self.DUT.Name)

        def func():
            print('getting number of bucket events for run {run} and {dia}...'.format(run=self.Run.Number, dia=self.DUT.Name))
            n_new = self.Tree.Draw('1', '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['bucket'], pul=self.Cut.CutStrings['pulser']), 'goff')
            n_old = self.Tree.Draw('1', '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['old_bucket'], pul=self.Cut.CutStrings['pulser']), 'goff')
            if show:
                print('New Bucket: {0} / {1} = {2:4.2f}%'.format(n_new, self.Run.NEntries, n_new / float(self.Run.NEntries) * 100))
                print('Old Bucket: {0} / {1} = {2:4.2f}%'.format(n_old, self.Run.NEntries, n_old / float(self.Run.NEntries) * 100))
            return {'old': n_old, 'new': n_new, 'all': float(self.Run.NEntries)}

        return do_pickle(pickle_path, func)

    def show_bucket_hits(self, show=True):
        # hit position
        h = TH2F('h', 'Diamond Margins', 80, -.3, .3, 52, -.3, .3)
        nr = 1 if not self.Channel else 2
        cut = '!({buc})&&{pul}'.format(buc=self.Cut.CutStrings['old_bucket'], pul=self.Cut.CutStrings['pulser'])
        self.Tree.Draw('dia_track_x[{nr}]:dia_track_y[{nr}]>>h'.format(nr=nr), cut, 'goff')
        projections = [h.ProjectionX(), h.ProjectionY()]
        zero_bins = [[], []]
        for i, proj in enumerate(projections):
            last_bin = None
            for bin_ in xrange(proj.GetNbinsX()):
                efficiency = proj.GetBinContent(bin_) / float(proj.GetMaximum())
                if bin_ > 1:
                    if efficiency > .05 and last_bin < 5:
                        zero_bins[i].append(proj.GetBinCenter(bin_ - 1))
                    elif efficiency < .05 and last_bin > 5:
                        zero_bins[i].append((proj.GetBinCenter(bin_)))
                last_bin = proj.GetBinContent(bin_)
        if show:
            print(zero_bins)
            c = TCanvas('c', 'Diamond Hit Map', 1000, 1000)
            h.GetXaxis().SetRangeUser(zero_bins[0][0], zero_bins[0][-1])
            h.GetYaxis().SetRangeUser(zero_bins[1][0], zero_bins[1][-1])
            h.Draw('colz')
            self.add(h, c)
        return h

    def draw_bucket_pedestal(self, show=True, corr=True, additional_cut=''):
        set_palette(55)
        cut_string = self.Cut.generate_custom(exclude=['bucket', 'timing'])
        cut_string += additional_cut
        self.draw_signal_vs_peaktime('e', cut_string, fine_corr=corr, prof=False, show=show)
        self.save_plots('BucketPedestal')

    def draw_bucket_waveforms(self, show=True, t_corr=True, start=100000):
        good = self.Waveform.draw(1, show=False, start_event=None, t_corr=t_corr)[0]
        cut = self.Cut.generate_custom(exclude=['bucket', 'timing']) + TCut('!({0})'.format(self.Cut.CutStrings['bucket']))
        bucket = self.Waveform.draw(1, cut=cut, show=False, start_event=start, t_corr=t_corr)[0]
        cut = self.Cut.generate_custom(exclude=['bucket', 'timing']) + TCut('{buc}&&!({old})'.format(buc=self.Cut.CutStrings['bucket'], old=self.Cut.CutStrings['old_bucket']))
        bad_bucket = self.Waveform.draw(1, cut=cut, show=False, t_corr=t_corr, start_event=None)[0]
        self.reset_colors()
        mg = TMultiGraph('mg_bw', 'Bucket Waveforms')
        leg = self.make_legend(.85, .4, nentries=3)
        names = ['good wf', 'bucket wf', 'both wf']
        for i, gr in enumerate([good, bucket, bad_bucket]):
            format_histo(gr, color=self.get_color(), markersize=.5)
            mg.Add(gr, 'lp')
            leg.AddEntry(gr, names[i], 'lp')
        format_histo(mg, draw_first=True, x_tit='Time [ns]', y_tit='Signal [mV]')
        x = [self.Run.signal_regions['e'][0] / 2, self.Run.signal_regions['e'][1] / 2 + 20]
        format_histo(mg, x_range=x, y_off=.7)
        self.draw_histo(mg, show=show, draw_opt='A', x=1.5, y=0.75, lm=.07, rm=.045, bm=.2, leg=leg)
        # y = mg.GetYaxis().GetXmin(), mg.GetYaxis().GetXmax()
        # self._add_buckets(y[0], y[1], x[0], x[1], avr_pos=-1, full_line=True)
        self.save_plots('BucketWaveforms')
        self.reset_colors()

    def show_bucket_means(self, show=True, plot_histos=True):
        pickle_path = self.PickleDir + 'Cuts/BucketMeans_{tc}_{run}_{dia}.pickle'.format(tc=self.TestCampaign, run=self.Run.Number, dia=self.DUT.Name)

        def func():
            gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
            cuts_nobucket = TCut('no_bucket', '')
            cuts_oldbucket = TCut('old_bucket', '')
            for key, value in self.Cut.CutStrings.iteritems():
                if not key.startswith('old') and key not in ['AllCuts', 'bucket']:
                    cuts_nobucket += value
                if key not in ['AllCuts', 'bucket']:
                    cuts_oldbucket += value
            h1 = self.draw_signal_distribution(show=False, evnt_corr=True)
            h2 = self.draw_signal_distribution(show=False, evnt_corr=True, cut=cuts_nobucket)
            h3 = self.draw_signal_distribution(show=False, evnt_corr=True, cut=cuts_oldbucket)
            if plot_histos:
                c = TCanvas('c', 'Bucket Histos', 1000, 1000)
                format_histo(h1, color=self.get_color(), lw=1, x_tit='Pulse Height [au]', y_tit='Entries')
                h1.Draw()
                format_histo(h2, color=self.get_color(), lw=1)
                h2.Draw('same')
                format_histo(h3, color=self.get_color(), lw=1)
                h3.Draw('same')
                self.Objects.append([h1, h2, h3, c])
            result = {name: [h.GetMean(), h.GetMeanError()] for name, h in zip(['new', 'no', 'old'], [h1, h2, h3])}
            gROOT.ProcessLine('gErrorIgnoreLevel = 0;')
            if show:
                print(result)
            return result

        res = func() if plot_histos else None
        return do_pickle(pickle_path, func, res)

    def compare_single_cuts(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        c1 = TCanvas('single', '', 1000, 1000)
        c2 = TCanvas('all', '', 1000, 1000)
        c2.SetLeftMargin(0.15)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        histos = []
        drawn_first = False
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                print('saving plot', key)
                save_name = 'signal_distribution_{cut}'.format(cut=key)
                histo_name = 'signal {range}{peakint}'.format(range=self.SignalRegionName, peakint=self.PeakIntegralName)
                histo_title = 'signal with cut ' + key
                histo = TH1F(histo_name, histo_title, 350, -50, 300)
                # safe single plots
                c1.cd()
                self.Tree.Draw("{name}>>{histo}".format(name=self.SignalName, histo=histo_name), value)
                self.save_plots(save_name)
                # draw all single plots into c2
                c2.cd()
                histo.SetLineColor(self.get_color())
                if not drawn_first:
                    format_histo(histo, title='Signal Distribution of Different Single Cuts', x_tit='Pulse Height [au]', y_tit='Entries', y_off=2)
                    histo.SetStats(0)
                    histo.Draw()
                    drawn_first = True
                else:
                    if key == 'AllCuts':
                        histo.SetLineWidth(2)
                    histo.Draw('same')
                histos.append(histo)
                legend.AddEntry(histo, key, 'l')
        # save c2
        legend.Draw()
        self.save_plots('all')
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        gROOT.SetBatch(0)

    def compare_normalised_cuts(self, scale=False, show=True):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        gROOT.SetBatch(1)
        self.reset_colors()
        c1 = TCanvas('single', '', 1000, 1000)
        name = 'sCutComparison'
        if scale:
            name += "_scaled"
        else:
            name += "_noarmalized"
        if scale:
            title = 'Scaled Signal Distribution with Single Cuts'
        else:
            title = 'Normalised Signal Distribution with Single Cuts'
        title += ';Pulse Height [au];Normalised Entries'

        stack = THStack(name, title)

        entries = 0
        for value in self.Cut.CutStrings.itervalues():
            if str(value):
                entries += 1
        legend = self.make_legend(x1=.57, nentries=entries - 2)
        histos = []
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                save_name = 'signal_distribution_normalised_{cut}'.format(cut=key)
                histo_name = 'signal {range}{peakint}'.format(range=self.SignalRegionName, peakint=self.PeakIntegralName)
                histo_title = 'normalized' if not scale else 'scaled'
                histo_title += ' signal with cut ' + key
                histo = TH1F(histo_name, histo_title, 350, -50, 300)
                # safe single plots
                c1.cd()
                self.Tree.Draw("{name}>>{histo}".format(name=self.SignalName, histo=histo_name), value)
                if scale:
                    histo = scale_histo(histo, to_max=True)
                else:
                    histo = normalise_histo(histo, from_min=True, x_range=[0, 30])
                histo.Draw()
                c1.Update()
                self.save_plots(save_name)
                # draw all single plots into c2
                histo.SetLineColor(self.get_color())

                if key == 'AllCuts':
                    histo.SetLineWidth(2)
                stack.Add(histo)
                histos.append(histo)
                legend.AddEntry(histo, key, 'l')
        stack.Draw()
        gROOT.SetBatch(0)

        for h in histos:
            h.SetStats(False)
        name = '{0}Cuts'.format('Normalised' if not scale else 'Scaled')
        format_histo(stack, y_off=1.4, x_off=1.1)
        self.save_histo(stack, name, show, lm=.15, leg=legend, draw_opt='nostack')
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        gROOT.SetBatch(0)

    def compare_consecutive_cuts(self, scale=False, show=True, save_single=True, short=False, x_range=None, redo=False):
        short_cuts = ['raw', 'saturated', 'timing', 'pulser', 'tracks', 'bucket', 'fiducial']
        legend = self.make_legend(.75 if short else .71, .88, nentries=len(self.Cut.ConsecutiveCuts) + 1 if not short else len(short_cuts) + 1, scale=.7)
        stack = THStack('scc', 'Signal Distribution with Consecutive Cuts')
        leg_style = 'l' if scale else 'f'
        for i, (key, cut) in enumerate(self.Cut.ConsecutiveCuts.iteritems()):
            if short and key not in short_cuts:
                continue
            self.info('adding cut {0}'.format(key))
            h = self.draw_signal_distribution(cut=cut, show=False, redo=redo)
            if scale:
                scale_histo(h, to_max=True, x_range=[30, 500])
            self.save_histo(h, 'signal_distribution_{n}cuts'.format(n=i), show=False, save=save_single)
            color = self.get_color()
            format_histo(h, color=color, stats=0, fill_color=color if not scale else None)
            stack.Add(h)
            leg_entry = '+ {0}'.format(key) if i else key
            legend.AddEntry(h, leg_entry, leg_style)
        if short:
            h = self.draw_signal_distribution(show=False, x_range=x_range)
            color = self.get_color()
            format_histo(h, color=color, stats=0, fill_color=color if not scale else None)
            stack.Add(h)
            legend.AddEntry(h, '+ other', leg_style)
        format_histo(stack, x_tit='Pulse Height [au]', y_tit='Number of Entries', y_off=1.9, draw_first=True)
        save_name = 'Consecutive{1}{0}'.format('Scaled' if scale else '', 'Short' if short else '')
        self.save_histo(stack, save_name, show, leg=legend, draw_opt='nostack', lm=0.14)
        stack.SetName(stack.GetName() + 'logy')
        # stack.SetMaximum(stack.GetMaximum() * 1.2)
        self.save_histo(stack, '{name}LogY'.format(name=save_name), show, logy=True, draw_opt='nostack', lm=0.14)
        self.reset_colors()

    def draw_cut_means(self, show=True, short=False):
        gr = self.make_tgrapherrors('gr_cm', 'Mean of Pulse Height for Consecutive Cuts')
        gr.SetPoint(0, 0, 0)
        short_keys = ['raw', 'saturated', 'timing', 'bucket', 'pulser', 'tracks', 'fiducial']
        cuts = OrderedDict((key, item) for key, item in self.Cut.ConsecutiveCuts.iteritems() if not short or key in short_keys)
        for i, (key, cut) in enumerate(cuts.iteritems(), 1):
            self.info('adding cut {0}'.format(key))
            h = self.draw_signal_distribution(cut=cut, show=False)
            self.info('{0}, {1}, {2}'.format(key, h.GetMean(), h.GetMeanError()))
            gr.SetPoint(i, i, h.GetMean())
            gr.SetPointError(i, 0, h.GetMeanError())
        if short:
            h = self.draw_signal_distribution(show=False)
            gr.SetPoint(gr.GetN(), gr.GetN(), h.GetMean())
            gr.SetPointError(gr.GetN(), 0, h.GetMeanError())
        format_histo(gr, markersize=.2, fill_color=17, y_tit='Mean Pulse Height [au]', y_off=1.4)
        y = [gr.GetY()[i] for i in xrange(1, gr.GetN())]
        gr.GetYaxis().SetRangeUser(min(y) - 1, max(y) + 1)
        gr.GetXaxis().SetLabelSize(.05)
        for i in xrange(1, gr.GetN()):
            bin_x = gr.GetXaxis().FindBin(i)
            gr.GetXaxis().SetBinLabel(bin_x, cuts.keys()[i - 1])
        self.Objects.append(self.save_histo(gr, 'CutMeans{s}'.format(s='Short' if short else ''), show, bm=.25, draw_opt='bap', lm=.12, x=1.5))
        gROOT.ProcessLine('gErrorIgnoreLevel = 0;')

    def draw_distance_vs_ph(self, show=True, bins=10):
        p = TProfile('pdph', 'Pulse Height vs. Track Distance', bins, self.DUT.Thickness, self.DUT.Thickness + 1)
        self.Tree.Draw('{}:{}>>pdph'.format(self.get_signal_name(), self.get_track_length_var()), self.Cut.generate_custom(exclude=['slope_x', 'slope_y']), 'goff')
        self.format_statbox(all_stat=True)
        format_histo(p, x_tit='Distance [#mum]', y_tit='Pulse Height [mV]', y_off=1.3, ndivx=405)
        self.save_histo(p, 'PHTrackLength', show, lm=.11)

    def test_landau_stats(self):
        gr = self.make_tgrapherrors('gr_ls', 'Landau Statistics')
        set_root_output(False)
        self.PBar.start(sum(int(pow(2, i / 2.)) for i in xrange(1, 40)))
        k = 0
        for j, i in enumerate(xrange(1, 40)):
            h = TH1F('h', 'h', 500, 0, 1000)
            for _ in xrange(int(pow(2, i / 2.))):
                k += 1
                h.Fill(gRandom.Landau(80, 5))
            self.PBar.update(k)
            gr.SetPoint(j, pow(2, i), h.GetMean())
            gr.SetPointError(j, 0, h.GetMeanError())
        self.PBar.finish()
        self.draw_histo(gr, draw_opt='alp', logx=True)

    def find_conv(self):
        gr = self.make_tgrapherrors('gr_c', 'chi2 vs nconv')
        for i, j in enumerate(xrange(10, 70, 5)):
            print(j)
            f = self.fit_langau(j, False)
            gr.SetPoint(i, j, f.Chi2 / f.NDF)
        self.draw_histo(gr)
    # endregion CUTS
    # ----------------------------------------

    # ----------------------------------------
    # region SHOW
    def draw_size(self, size=None, color=1, name=''):
        if size is None:
            warning('The {} size of "{}" was not specified in the dia info'.format(' {}'.format(name) if name else '', self.DUT.Name))
            return
        cx, cy = self.find_center()
        x, y = size
        self.draw_box(cx - x / 2, cy - y / 2, cx + x / 2, cy + y / 2, line_color=color, width=3, name=name, fillstyle=4000)

    def draw_detector_size(self):
        self.draw_size(self.DUT.Size, color=432, name='detector')

    def draw_metal_size(self):
        size = [self.DUT.PadSize.n] * 2 if self.DUT.PadSize is not None else None
        self.draw_size(size, color=923, name='metal')

    def draw_guard_ring(self):
        size = [self.DUT.GuardRing] * 2 if self.DUT.GuardRing is not None else None
        self.draw_size(size, color=417, name='guard ring')

    def draw_all_sizes(self):
        self.draw_fid_cut()
        self.draw_detector_size()
        self.draw_metal_size()
        self.draw_guard_ring()

    def draw_signal_vs_signale(self, show=True):
        gStyle.SetPalette(53)
        cut = self.Cut.generate_custom(exclude=['bucket'])
        num = self.get_signal_number(region='e')
        cut += TCut('IntegralPeakTime[{0}]<94&&IntegralPeakTime[{0}]>84'.format(num))
        h = TH2F('hsse', 'Signal b vs Signal e', 62, -50, 200, 50, 0, 200)
        self.Tree.Draw('{sige}:{sigb}>>hsse'.format(sigb=self.SignalName, sige=self.get_signal_name(region='e')), cut, 'goff')
        format_histo(h, x_tit='Signal s_b [au]', y_tit='Signal s_e [au]', z_tit='Number of Entries', z_off=1.1, y_off=1.5, stats=0)
        self.Objects.append(self.save_histo(h, 'SignalEvsSignalB', show, rm=.15, lm=.13, draw_opt='colz'))
        gStyle.SetPalette(1)

    def show_integral_names(self):
        for i, name in enumerate(self.IntegralNames):
            if name.startswith('ch{}'.format(self.Channel)):
                print(str(i).zfill(3), name)

    def show_information(self, header=True, prnt=True):
        peak_int = '{} ({})'.format(self.PeakIntegral, remove_letters(self.PeakIntegralName))
        region = '{} ({})'.format(self.SignalRegion, self.SignalRegionName.split('_')[-1])
        rows = [[self.Run.Number, self.Run.RunInfo['runtype'], self.DUT.Name, '{:14.1f}'.format(self.Run.Flux.n), '{:+5.0f}'.format(self.DUT.Bias), region, peak_int]]
        print_table(rows, self.get_info_header() if header else None, prnt=prnt)
    # endregion SHOW
    # ----------------------------------------

    def check_alignment(self):
        """ check if the events from telescope and digitiser are aligned"""
        pickle_path = self.make_pickle_path('Alignment', run=self.Run.Number)

        def f():
            alignment = PadAlignment(self.Run.Converter, verbose=False)
            return alignment.IsAligned

        is_aligned = do_pickle(pickle_path, f)
        log_warning('\nRun {r} is misaligned :-('.format(r=self.Run.Number)) if not is_aligned else do_nothing()
        return is_aligned

    def draw_alignment(self, n_pulser=200, thresh=40, show=True):
        """ draw the aligment of telescope and digitiser events """
        xbins = self.Bins.get_pulser(n_pulser)
        p = self.Pulser.draw_hit_efficiency(xbins, show=False)
        h = TH2F('ha{}'.format(self.Run.Number), 'Event Alignment', *(xbins + [3, 0, 3]))
        for ibin in xrange(1, xbins[0]):
            h.SetBinContent(ibin, 2, int(p.GetBinContent(ibin) <= thresh) + 1)
        format_histo(h, x_tit='Event Number', y_tit='Alignment', stats=False, l_off_y=99, center_y=True)
        gStyle.SetPalette(3, array([1, 2, 3], 'i'))
        leg = self.make_legend(nentries=2, x2=.9, margin=.2)
        leg.AddEntry(self.draw_box(0, 0, 0, 0, line_color=3, name='b1'), 'aligned', 'f')
        leg.AddEntry(self.draw_box(0, 0, 0, 0, line_color=2), 'misaligned', 'f')
        self.save_histo(h, 'EventAlignment', draw_opt='col', rm=.08, leg=leg, show=show, prnt=show)
        return h

    def analyse_signal_histograms(self):
        gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
        # gROOT.SetBatch(1)
        legend = TLegend(0.7, 0.3, 0.98, .7)
        gr1 = TGraphErrors()
        gr1.SetTitle('mean values')
        gr1.SetMarkerStyle(20)
        gr2 = TGraph()
        gr2.SetTitle('median values')
        gr2.SetMarkerStyle(21)
        gr2.SetMarkerColor(2)
        gr3 = TGraph()
        gr3.SetMarkerStyle(22)
        gr3.SetMarkerColor(3)
        histos = []
        i = 0
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                print('process cut ' + key)
                # h = TH1F('h', '', 600, -100, 500)
                # self.Tree.Draw("{name}>>h".format(name=self.signal_name), value)
                h = self.draw_signal_distribution(evnt_corr=True, cut=value, show=False)
                i_mean = self.__get_mean(h)
                median = self.__get_median(h)
                mpv = self.__get_mpv(h)
                # print(mean, median, mpv)
                gr1.SetPoint(i, i, i_mean[0])
                gr1.SetPointError(i, 0, i_mean[1])
                gr2.SetPoint(i, i, median)
                gr3.SetPoint(i, i, mpv)
                histos.append(h)
                i += 1
        # rename bins
        legend.AddEntry(gr1, 'mean', 'lp')
        legend.AddEntry(gr2, 'median', 'lp')
        legend.AddEntry(gr3, 'mpv', 'lp')
        xaxis = gr1.GetXaxis()
        i = 0
        for key, value in self.Cut.CutStrings.iteritems():
            if str(value) or key == 'raw':
                bin_x = xaxis.FindBin(i)
                gr1.GetXaxis().SetBinLabel(bin_x, key[:7])
                i += 1
        gROOT.ProcessLine("gErrorIgnoreLevel = 0;")
        # gROOT.SetBatch(0)
        c1 = TCanvas('c1', '', 1000, 1000)
        c1.cd()
        gr1.GetXaxis().SetRangeUser(-1, len(histos) + 1)
        gr1.Draw('alp')
        gr2.Draw('lp')
        gr3.Draw('lp')
        legend.Draw()
        self.Objects.append(legend)
        return [gr1, gr2, gr3]

    @staticmethod
    def __get_histo_without_pedestal(histo):
        h = histo
        h.GetXaxis().SetRangeUser(0, 30)
        min_bin = h.GetMinimumBin()
        min_x = h.GetBinCenter(min_bin)
        h.GetXaxis().SetRangeUser(min_x, 500)
        return h

    def __get_mean(self, histo):
        h = self.__get_histo_without_pedestal(histo)
        h.GetXaxis().SetRangeUser(0, 30)
        min_bin = h.GetMinimumBin()
        min_x = h.GetBinCenter(min_bin)
        h.GetXaxis().SetRangeUser(min_x, 500)
        return [h.GetMean(), h.GetMeanError()]

    def __get_median(self, histo):
        h = self.__get_histo_without_pedestal(histo)
        integral = h.GetIntegral()
        median_i = 0
        for j in range(h.GetNbinsX() - 1):
            if integral[j] < 0.5:
                median_i = j
            else:
                break
        weight = (0.5 - integral[median_i]) / (integral[median_i + 1] - integral[median_i])
        median_x = h.GetBinCenter(median_i) + (h.GetBinCenter(median_i + 1) - h.GetBinCenter(median_i)) * weight
        return median_x

    def __get_mpv(self, histo):
        h = self.__get_histo_without_pedestal(histo)
        max_bin = h.GetMaximumBin()
        return h.GetBinCenter(max_bin)

    def draw_snrs(self, show=True, lego=True, proj=False, draw_opt='lego2'):
        self.Verbose = False
        gStyle.SetPaintTextFormat('5.4g')
        lego = False if proj else lego
        gr = self.make_tgrapherrors('gr', 'Signal to Noise Ratios')
        h = TProfile2D('h_snr', 'Signal to Noise Ratios', 10, 0, 5, 10, 0, 5)
        i = 0
        for name, region in self.get_all_signal_names().iteritems():
            if self.SignalRegionName.split('_')[-1] in region:
                peak_integral = self.get_peak_integral(remove_letters(region))
                snr = self.calc_snr(name=name, reg=self.get_all_signal_names()[name])
                h.Fill(peak_integral[0] / 2., peak_integral[1] / 2., snr.n)
                gr.SetPoint(i, i + 1, snr.n)
                gr.SetPointError(i, 0, snr.s)
                gr.GetListOfFunctions().Add(self.draw_tlatex(i + 1, snr.n + snr.s * 1.5, str(peak_integral), align=22, size=.02))
                i += 1
        format_histo(gr, y_tit='SNR', y_off=1.2, color=self.get_color(), fill_color=1)
        vals = sorted([h.GetBinContent(i) for i in xrange(h.GetNbinsX() * h.GetNbinsY()) if h.GetBinContent(i)])
        x, y, z1 = Long(0), Long(0), Long(0)
        xmin, ymin = h.GetXaxis().GetXmin(), h.GetYaxis().GetXmin()
        h.GetBinXYZ(h.GetMaximumBin(), x, y, z1)
        x1, y1 = (x - 1) / 2. + xmin, (y - 1) / 2. + ymin
        self.__draw_profiles(h, x, y, proj)
        format_histo(h, x_tit='Left Length [ns]', x_off=1.45, y_tit='Right Length [ns]', y_off=1.6, z_tit='snr', z_off=1.6, stats=0, z_range=[vals[2], max(vals)])
        h.SetContour(50)
        gStyle.SetPalette(53)
        if draw_opt == 'coltext':
            self.show_best_snr(h, x1, y1, show)
        else:
            self.save_histo(h, 'SNRLego', show and lego, draw_opt=draw_opt, bm=.2, rm=.1, lm=.13, phi=-30, theta=40)
        gStyle.SetPalette(1)
        self.save_histo(gr, 'SNR', not (lego or proj) and show, draw_opt='bap')

    def show_best_snr(self, histo, x, y, show):
        h = histo
        format_histo(h, x_off=1, y_off=1.15, stats=0, z_tit='snr [au]', z_off=1.35)
        self.draw_histo(h, '', show, draw_opt='colztext', rm=.16)
        self.draw_vertical_line(x, -1e5, 1e5, color=2, style=2, name='a', w=2)
        self.draw_vertical_line(x + .5, -1e5, 1e5, color=2, style=2, name='b', w=2)
        self.draw_horizontal_line(y, 0, 10, color=418, style=2, w=2, name='c')
        self.draw_horizontal_line(y + .5, 0, 100, color=418, style=2, w=2, name='d')
        self.save_plots('SNRColText')

    def __draw_profiles(self, histo, x, y, show=True):
        h = histo
        py = h.ProfileY('Right Length', x, x)
        px = h.ProfileX('Left Length', y, y)
        vals = [py.GetBinContent(i) for i in xrange(py.GetNbinsX()) if py.GetBinContent(i)] + [px.GetBinContent(i) for i in xrange(px.GetNbinsX()) if px.GetBinContent(i)]
        format_histo(py, stats=0, lw=2)
        format_histo(px, stats=0, lw=2)
        py.SetLineColor(2)
        px.SetLineColor(418)
        leg = self.make_legend(.68, .95)
        [leg.AddEntry(p, p.GetName(), 'fp') for p in [py, px]]
        stack = THStack('s_sp', 'SNR Profiles')
        stack.Add(py, 'histe')
        stack.Add(px, 'histe')
        format_histo(stack, draw_first=True, x_tit='Integral Length [ns]', y_tit='snr [au]', y_off=1.35)
        stack.SetMinimum(increased_range([min(vals), max(vals)], .5, .5)[0])
        stack.SetMaximum(increased_range([min(vals), max(vals)], .5, .5)[1])
        self.save_histo(stack, 'SNRProfiles', show, draw_opt='nostack', leg=leg, lm=.13)

    def calc_snr(self, name=None, reg=''):
        signal_name = self.SignalName if name is None else name
        peak_int = remove_letters(self.get_all_signal_names()[signal_name])
        ped_sigma = make_ufloat(self.Pedestal.draw_disto_fit(save=False, name=self.Pedestal.get_signal_name(peak_int=peak_int), show=False), par=2)
        signal = make_ufloat(self.draw_pulse_height(corr=True, show=False, sig=signal_name)[1])
        snr = signal / ped_sigma
        print('{name} {0}\t| SNR is: {snr}\t {1} {2}'.format(self.get_peak_integral(peak_int), signal.n, ped_sigma.n, name=reg, snr=snr))
        return snr

    # ----------------------------------------
    # region PEAK INTEGRAL
    def find_best_snr(self, show=True, same_width=False):
        gROOT.SetBatch(1)
        gr = self.make_tgrapherrors('gr', 'Signal to Noise Ratios')
        peak_integrals = OrderedDict(sorted({key: value for key, value in self.Run.peak_integrals.iteritems() if len(key) < 3}.items()))
        i = 0
        for name, value in peak_integrals.iteritems():
            signal = self.get_signal_name('b', name)
            snr = self.calc_snr(signal)
            print(value)
            x = (value[1] + value[0]) / 2. if not same_width else value[0] / 2.
            gr.SetPoint(i, x, snr[0])
            gr.SetPointError(i, 0, snr[1])
            i += 1
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'SNR', 1000, 1000)
        format_histo(gr, x_tit='Integralwidth [ns]', y_tit='SNR')
        gr.Draw('ap')
        gROOT.SetBatch(0)
        self.save_plots('BestSNR')
        self.Objects.append([gr, c])

    def signal_vs_peakintegral(self, show=True, ped=False):
        gROOT.SetBatch(1)
        gr = self.make_tgrapherrors('gr', '{sig} vs Peak Integral'.format(sig='Signal' if not ped else 'Pedestal'))
        peak_integrals = OrderedDict(sorted({key: value for key, value in self.Run.peak_integrals.iteritems() if len(key) < 3}.items()))
        i = 0
        ratio = '{0}{1}'.format(self.Run.peak_integrals.values()[0][0], self.Run.peak_integrals.values()[0][1])
        for name, value in peak_integrals.iteritems():
            sig_name = self.get_signal_name(region='b', peak_integral=name)
            signal = self.draw_pulse_height(corr=True, show=False, sig=sig_name) if not ped else self.Pedestal.draw_disto(save=False, name=self.Pedestal.get_signal_name(peak_int=name))
            par = 2 if ped else 0
            gr.SetPoint(i, (value[1] + value[0]) / 2., signal.Parameter(par))
            gr.SetPointError(i, 0, signal.ParError(par))
            i += 1
        if show:
            gROOT.SetBatch(0)
        c = TCanvas('c', 'Signal vs Peak Integral', 1000, 1000)
        format_histo(gr, x_tit='Integralwidth [ns]', y_tit='Signal [au]', y_off=1.3)
        gr.Draw('ap')
        gROOT.SetBatch(0)
        self.save_plots('{sig}PeakInt_{rat}'.format(rat=ratio, sig='Ped' if ped else 'Sig'))
        self.Objects.append([gr, c])
    # endregion PEAK INTEGRAL
    # ----------------------------------------

    # ----------------------------------------
    # region MISCELLANEOUS
    def get_all_signal_names(self, sig_type='signal'):
        names = OrderedDict()
        for region in self.Run.IntegralRegions[self.DUT.Number - 1]:
            if sig_type in region:
                for integral in self.Run.PeakIntegrals[self.DUT.Number - 1]:
                    name = 'ch{ch}_{reg}_{int}'.format(ch=self.Channel, reg=region, int=integral)
                    num = self.IntegralNames.index(name)
                    reg = region.replace(sig_type, '').strip('_') + integral.replace('PeakIntegral', '')
                    names[self.SignalDefinition.format(pol=self.Polarity, num=num)] = reg
        return names
    # endregion MISCELLANEOUS
    # ----------------------------------------

    def get_peak_integral(self, name):
        return self.Run.PeakIntegrals[self.DUT.Number - 1]['PeakIntegral{}'.format(name) if 'Peak' not in str(name) else name]

    @staticmethod
    def make_region(signal, region=''):
        return '{}{}'.format(signal, '_' + region if region else '')


if __name__ == '__main__':

    args = init_argparser(run=23, tc='201908', dut=1, tree=True, has_verbose=True)
    z = PadAnalysis(args.run, args.dut, args.testcampaign, args.tree, args.verbose)
