#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TMath, TGraph
from ROOT.gRandom import Landau as rLandau
from numpy import polyfit, argmax, insert, array_split, invert, sum, column_stack, any, sort
from numpy.random import normal, rand
from scipy.signal import find_peaks
from scipy.stats import poisson
from uncertainties import ufloat_fromstr
from uncertainties.umath import log as ulog, sqrt as usqrt

from helpers.fit import PoissonI, Gauss, Landau
from src.sub_analysis import PadSubAnalysis
from helpers.draw import *
from src.analysis import update_pbar


class PeakAnalysis(PadSubAnalysis):

    Prominence = 5
    Distance = 35

    def __init__(self, pad_analysis):
        super().__init__(pad_analysis, pickle_dir='Peaks')
        self.WF = self.Ana.Waveform
        self.NoiseThreshold = self.calc_noise_threshold()
        self.Threshold = max(self.NoiseThreshold, self.Ana.get_min_signal(self.Ana.get_signal_name(peak_int=1)))
        self.BinWidth = self.DigitiserBinWidth
        self.StartAdditional = self.get_start_additional()
        self.NBunches = self.calc_n_bunches()
        self.MaxBunch = self.NBunches + 3
        self.BunchSpacing = self.Ana.BunchSpacing
        self.Fit = TF1('lan', 'landau', 0, 512)

    # ----------------------------------------
    # region INIT
    def calc_noise_threshold(self):
        """ return peak threshold, 5 times the raw noise + mean of the noise. """
        return abs(self.Ana.Pedestal.get_raw_mean() + 5 * self.Ana.Pedestal.get_raw_noise()).n

    def get_start_additional(self, b0=None):
        """Set the start of the additional peaks 2.5 bunches after the signal peak to avoid the biased bunches after the signal. Signalbunch = bunch 1 """
        return int(self.Run.IntegralRegions[self.DUT.Number - 1]['signal_a'][0] + self.Ana.BunchSpacing * (choose(b0, 4) - 1.5) / self.BinWidth)

    def get_bunch_range(self, b0=None, b_end=None, bins=False):
        rng = array([self.get_start_additional(b0), self.get_start_additional(choose(b_end, self.NBunches + 4))])
        return rng.tolist() if bins else rng * self.BinWidth

    def calc_n_bunches(self):
        return int((self.Run.NSamples - self.StartAdditional) * self.BinWidth / self.BunchSpacing)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region DATA
    def get_all(self, cut=None, thresh=None, fit=False, redo=False):
        t, h, n = self.find_all(thresh, fit, redo)
        cut = self.get_event_cut(cut)
        lcut = cut.repeat(n)  # make cut for the flat arrays
        return array(t)[lcut] / 100, array(h)[lcut] / 100, array(n)[cut]

    def get(self, flat=False, thresh=None, fit=False, cut=None, i=0):
        t, h, n = self.get_all(cut, thresh, fit)
        data = [t, h]
        return array(data[i]) if flat else array(split(data[i], cumsum(n)[:-1]), dtype=object)

    def get_heights(self, flat=False, thresh=None, fit=False, cut=None, corr=False):
        return self.get(flat, thresh, fit, self.p2ecut(self.get_bunch_cut(1)) if corr else cut, i=1)

    def get_times(self, flat=False, thresh=None, fit=False, cut=None, corr=False):
        return self.get_corrected_times(True, thresh, fit) if corr else self.get(flat, thresh, fit, cut, i=0)

    def get_npeaks(self, thresh=None, cut=...):
        hdf5_path = self.make_simple_hdf5_path(suf=f'{choose(thresh, self.Threshold):.1f}_1', dut=self.Channel)
        return self.get_all(..., thresh, fit=file_exists(hdf5_path))[-1][cut]

    def get_from_tree(self, redo=False, cut=...):
        return array(do_hdf5(self.make_simple_hdf5_path('Peaks'), self.Run.get_tree_vec, var=self.Ana.PeakName, dtype='f4', redo=redo))[cut]

    def get_pars(self, i=0, cut=...):
        return self.find_pars()[:, i][cut]

    def get_cft(self, cut=...):
        return self.get_pars(2, cut)

    def get_tot(self, cut=...):
        return self.get_pars(3, cut)
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_binning(self, bin_size=None, n=None):
        return self.WF.get_binning(bin_size) if n is None else make_bins(*self.get_bunch_range(n, n + 1), choose(bin_size, self.BinWidth), last=True)

    def get_n_events(self):
        return count_nonzero(self.Ana.get_pulser_cut())

    def get_bunch_nrs(self, start=None):
        return arange(choose(start, 3), self.MaxBunch + 1)

    def get_corrected_times(self, pcut=True, thresh=None, fit=False):
        pcut = self.get_bunch_cut(n=1) & pcut
        times = self.get_times(True, thresh, fit, cut=...)
        t1 = times[pcut]
        times = times[self.p2e2pcut(pcut)]  # select all peaks if the event has a peak in bunch 1
        return times + mean(t1) - t1.repeat(self.get_npeaks(cut=self.p2ecut(pcut)))

    def get_n_total(self, cut=None, b0=None, b_end=None, thresh=None, fit=True, corr=True):
        n = count_nonzero(self.get_bunch_cut(b0, b_end, thresh, fit) & choose(cut, self.get_pulser_cut())) * (self.get_pad_gr_ratio(thresh) if corr else 1)
        return n + ufloat(0, sqrt(n.n if corr else n))

    def get_lambda(self, ecut=None, thresh=None, fit=True, corr=True):
        """calculate average number of peaks for single bunch."""
        cut = choose(ecut, self.Ana.get_pulser_cut())
        n_events, n_peaks = count_nonzero(cut), self.get_n_total(self.e2pcut(cut, thresh), thresh=thresh, fit=fit, corr=corr)
        return -ulog(1 - n_peaks / n_events / self.NBunches)  # calc lambda based on the number of events with no peak

    def get_flux(self, lam=None, n_bunches=None, thresh=None, fit=True, prnt=True, corr=True):
        n_bunches = 1 if lam is None else choose(n_bunches, self.NBunches)
        flux = choose(lam, self.get_lambda(thresh=thresh, fit=fit, corr=corr)) / (self.Ana.BunchSpacing * n_bunches * self.DUT.get_area()) * 1e6  # ns -> ms (kHz)
        self.info('Estimated flux by number of peaks: {}'.format(make_flux_string(flux, term=True)), prnt=prnt)
        return flux

    def get_pad_gr_ratio(self, thresh=None):
        return ufloat_fromstr(self.Run.Config.get_list('BASIC', 'pad guard ring ratio')['0' if thresh is None else str(thresh)])

    def get_pad_count(self, thresh=None):
        hsig, hall = self.draw_bunch_height(thresh=thresh, fit=True, show=False), self.draw_bunch_height(thresh=thresh, all_=True, fit=True, show=False)
        lo, hi = [hsig.GetBinCenter(i) for i in [hsig.FindFirstBinAbove(.5 * hsig.GetMaximum()), hsig.FindLastBinAbove(.2 * hsig.GetMaximum())]]
        sfit = FitRes(hsig.Fit('landau', 'qs', '', lo, hi))
        fit = self.Draw.make_f('f', 'landau', 0, 400, pars=sfit.Pars, line_style=2)
        fit.SetParError(1, sfit[1].s)
        fit.FixParameter(1, sfit[1].n)
        fit.FixParameter(2, sfit[2].n)
        hall.Fit(fit, 'qsb', '', hall.GetBinCenter(hall.GetMaximumBin() + 3), hi)
        return (hall.Integral(0, hall.FindBin(hi)) - ufloat(fit.Integral(0, hi), fit.IntegralError(0, hi) * 10) / hall.GetBinWidth(1)) / hall.Integral()

    def get_p(self, corr=False):
        return self.get_n_total(corr=corr) / self.get_n_events() / self.NBunches

    def get_eff_b0(self):
        return (self.get_n_total(b0=0, corr=False) - self.get_n_consecutive(corr=False, n_peaks=True)) / (self.get_p() * self.get_n_events())

    def get_p_extra(self):
        scale = self.get_n_total(b0=1) / self.get_n_events()  # assume the additional miss as many events as the signal
        return (1 - poisson.cdf(0, self.get_lambda().n)) / scale

    def get_n_consecutive(self, n=2, corr=True, n_peaks=False):
        if n_peaks:  # just calculate the number of observed peak
            return self.get_p(corr) ** n * self.get_n_events()
        lam = self.get_lambda(corr=corr).n  # lambda for a single bunch
        p = 1 - poisson.cdf(0, lam)  # probability to have at least one peak in a single bunch
        return sum(p ** (n - 1) * poisson.pmf(i, lam) * i for i in range(50)) * self.get_n_events()

    def get_peak_fac(self, redo=False):
        def f():
            h = self.draw(show=False)
            self.PBar.start(self.MaxBunch - 2)
            facs = [self.get_n_bunch(i) / h.Integral(*self.get_bunch_range(i, i + 1, bins=True)) for i in self.get_bunch_nrs()]
            return mean_sigma(facs)[0]
        return do_pickle(self.make_simple_pickle_path('PeakFac'), f, redo=redo)

    def get_bunch_eff(self, n=0):
        """return ratio of bunch events to total events """
        return self.get_n_bunch(n) / self.get_peak_fac() / count_nonzero(self.get_pulser_cut())

    @update_pbar
    def get_n_bunch(self, n=0, bin_width=.2, h=None, show=False):
        h = choose(h, self.draw(fit=True, show=show, bin_size=bin_width, logy=False, x_range=self.get_bunch_range(n, n + 1), cut=self.get_pulser_cut() & self.get_bunch_cut(n)))
        thresh = FitRes(h.Fit('gaus', 'qs0'))[0].n * .4
        f = Gauss(h, [h.GetBinCenter(i) for i in [h.FindFirstBinAbove(thresh), h.FindLastBinAbove(thresh)]]).fit(show=show)
        return f.Fit.Integral(0, self.Run.NSamples) / bin_width

    def get_n_additional(self, n=None, end=None, thresh=None, fit=False, pcut=None):
        cut, n = self.get_bunch_cut(n, end, thresh, fit), self.get_npeaks(thresh, cut=...)
        filled = insert(cut, insert(cumsum(n)[:-1], 0, 0).repeat(max(n) - n).astype('i'), False)  # make all events the same length and fill with False
        npeaks = sum(filled.reshape(filled.size // max(n), max(n)), axis=1)
        if pcut is not None:
            npeaks[invert(self.p2ecut(pcut))] = 0
        return npeaks

    def get_bunch_values(self, n, thresh=None, excl=0, cut=None, i=1, fit=False, isolated=True, all_=False, pars=False):
        pcut = self.get_isolated_cut(n, thresh, fit, all_) if isolated else self.get_pulser_cut() & self.get_bunch_cut(n, thresh=thresh)
        values = (self.get_pars(i) if pars else self.get(True, thresh, fit, cut=..., i=i))[pcut if cut is None else pcut & cut]
        return values[int(excl * values.size):]

    @update_pbar
    def get_bunch_heights(self, n, thresh=None, excl=0, cut=None, fit=False, isolated=True, all_=False):
        return self.get_bunch_values(n, thresh, excl, cut, 1, fit, isolated, all_)

    def get_bunch_height(self, n=1, thresh=80, fit=False, cut=None, all_=False, isolated=True):
        x = self.get_bunch_heights(n, cut=cut, fit=fit, isolated=isolated, all_=all_).astype('d')
        return mean_sigma(x[(x > thresh) & (x < 500)])[0]

    def get_bunch_mpv(self, n=1, fit=True, cut=None, all_=False, show=False, par=None):
        h = self.draw_bunch_height(n, bin_width=1, fit=fit, cut=cut, all_=all_, show=show, draw_opt='')
        landau = Landau(h, ax_range(*[i.n for i in get_fwhm(h, ret_edges=True)], -.2, .2))
        return landau.get_mpv(show) if par is None else landau.fit(show=show, minuit=False)[par]

    def get_bunch_times(self, n=1, thresh=None, excl=0, cut=None, fit=False, isolated=True, all_=False, cft=False):
        return self.get_bunch_cfts(n, excl, cut) if cft else self.get_bunch_values(n, thresh, excl, cut, 0, fit, isolated, all_)

    def get_bunch_time(self, n=1, fit=True, cut=None, show=False, par=1, redo=False):
        def f():
            return fit_fwhm(self.draw_bunch_times(n, fit=fit, cut=cut, show=show, draw_opt='', bin_width=.1), show=show)
        return do_pickle(self.make_simple_pickle_path(f'BT{n}', int(fit)), f, redo=redo)[par]

    def get_bunch_cfts(self, n, excl=0, cut=None):
        return self.get_bunch_values(n, None, excl, cut, 2, fit=True, pars=True)

    def get_bunch_tot(self, n, excl=0, cut=None):
        return self.get_bunch_values(n, None, excl, cut, 3, fit=True, pars=True)

    def get_spacings(self, n, cut=None, ecut=True):
        cut = choose(cut, self.get_spacing_cut(n, ecut))
        t, p = self.get_times(True, fit=True, cut=cut), self.get_heights(True, cut=cut, fit=True)
        return (t[1::2] - t[::2]) / (n - 1)

    def get_spacing(self, redo=False):
        return do_pickle(self.make_simple_pickle_path('BunchSpacing'), self.draw_all_spacings, redo=redo, show=False)[0]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region CUT
    def p2ecut(self, pcut, thresh=None, ecut=True):
        e = arange(self.Run.NEvents).repeat(self.get_npeaks(thresh=thresh))[pcut]
        return self.Ana.make_event_cut(e) & ecut

    def e2pcut(self, ecut, thresh=None):
        return ecut.repeat(self.get_npeaks(thresh=thresh))

    def p2e2pcut(self, pcut, thresh=None):
        return self.e2pcut(self.p2ecut(pcut, thresh), thresh)

    def get_event_cut(self, cut=None):
        return ones(self.Run.NEvents, '?') if cut is None or cut is Ellipsis else cut

    def get_ebunch_cut(self, n):
        return self.p2ecut(self.get_bunch_cut(n))

    def get_pre_bunch_cut(self, n, delta=-1):
        return None if delta is None else self.e2pcut(self.get_ebunch_cut(n + delta))

    def make_bunch_cut(self, n=1):
        return self.Ana.get_pulser_cut() & invert(self.get_ebunch_cut(n - 1)) & invert(self.get_ebunch_cut(n + 1))

    def get_pulser_cut(self, thresh=None):
        return self.e2pcut(self.Ana.get_pulser_cut(), thresh)

    def get_bunch_cut(self, n=None, end=None, thresh=None, fit=False):
        t = self.get_times(True, thresh, fit, cut=...)
        t0, t1 = self.get_bunch_range(n, n + 1 if end is None and n is not None else end)
        return (t > t0) & (t < t1)

    def get_time_cut(self, t0, t1, fit=True):
        t = self.get_times(True, None, fit, cut=...)
        return (t > t0) & (t < t1)

    def get_isolated_cut(self, n, thresh=None, fit=False, all_=False):
        # find events with peaks in buckets before and after, invert and repeat cut again
        if all_:
            return self.get_all_isolated(thresh, fit)
        t1, t2, np = self.get_bunch_cut(n - 1, fit=fit, thresh=thresh), self.get_bunch_cut(n + 1, fit=fit, thresh=thresh), self.get_npeaks(thresh=thresh)
        return invert(self.p2ecut(t1, thresh)).repeat(np) & invert(self.p2ecut(t2, thresh)).repeat(np) & self.get_pulser_cut(thresh) & self.get_bunch_cut(n, fit=fit, thresh=thresh)

    def get_all_isolated(self, thresh=None, fit=False):
        def f():
            return any([self.get_isolated_cut(i, thresh, fit) for i in arange(3, self.NBunches + 4)], axis=0)
        return do_pickle(self.make_simple_pickle_path('AllIso', f'{thresh}{int(fit)}'), f)

    def get_spacing_cut(self, n, ecut=True):
        return ecut & (self.get_npeaks() == 2) & self.p2ecut(self.get_isolated_cut(1, fit=1)) & self.p2ecut(self.get_isolated_cut(n, fit=1))
    # endregion CUT
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_thresh(self, thresh=None):
        self.Draw.info('Threshold: {:.1f}'.format(choose(thresh, self.Threshold)))

    def draw(self, corr=False, thresh=None, bin_size=None, y_range=None, cut=..., fit=False, logy=True, **kwargs):
        x = self.get_times(flat=True, cut=..., thresh=thresh, fit=fit, corr=corr)[cut]
        h = self.Draw.distribution(x, self.get_binning(bin_size), 'Peak Times', w=1.5, h=.75, logy=logy, stats=set_statbox(entries=True, w=.25), bm=.2, lm=.082, **kwargs)
        format_histo(h, x_tit='Time [ns]', y_tit='Number of Peaks', y_range=choose(y_range, [.9, h.GetMaximum() * 2]), tit_size=.06, lab_size=.05, y_off=.7)
        return h

    def draw_n(self, n=None, fit=False, show=True):
        """draw disto of additional number of peaks."""
        h = self.Draw.distribution(self.get_n_additional(n), make_bins(self.NBunches), 'Number of Peaks', x_tit='Number of Peaks', show=False, lw=2)
        PoissonI(h, p1=1).fit(show=fit)
        format_statbox(h, fit=fit, m=True, entries=True, w=.4 if fit else .2, form='.2f')
        self.Draw(h, f'PeakNumbers{"Fit" if fit else ""}', show, logy=True, lm=.11, draw_opt='')

    @update_pbar
    def draw_spacing(self, n=3, bin_size=None, show=True, redo=False, ecut=True):
        def f():
            x, bs = self.get_spacings(n, ecut=ecut), choose(bin_size, get_y(3, 22, .02, .002, n))
            h = self.Draw.distribution(x, make_bins(*(array([-3, 3]) + self.BunchSpacing), bs), f'Bunch {n} Spacing', x_tit='Time Delta [ns]', show=show, draw_opt='')
            format_statbox(h, fit=1, all_stat=True)
            return fit_fwhm(h, show=show, fit_range=.5)
        return do_pickle(self.make_simple_pickle_path(f'Spacing{n}'), f, redo=redo or show)

    def draw_spacing_vs_peaktime(self, n=3, bin_size=.2, ecut=True, **kwargs):
        cut = self.get_spacing_cut(n, ecut)
        y, t = self.get_spacings(n, cut), self.get_times(flat=True, cut=cut)[::2]
        self.Draw.profile(t, y, self.get_binning(bin_size, n=1), x_tit='Signal Peak Time [ns]', y_tit='Bunch Spacing [ns]', **kwargs)

    def draw_all_spacings(self, show=True, ecut=True, redo=False):
        self.PBar.start(self.NBunches) if redo or not file_exists(self.make_simple_pickle_path(f'Spacing{self.NBunches + 3}')) else do_nothing()
        x, y = self.get_bunch_nrs(), [self.draw_spacing(i, show=False, ecut=ecut, redo=redo)[1] for i in self.get_bunch_nrs()]
        self.Draw.graph(x, y, 'Bunch Spacings', x_tit='Bunch Number', y_tit='Spacing [ns]', **Draw.mode(2, lm=.11, y_off=.95), show=show, y_range=ax_range(y, 0, .5, .2))
        a, b, xa, c = [mean_sigma(y[-10:])[0]] * 2, [1e9 / ufloat(self.BeamFrequency, 1e4)] * 2, self.MaxBunch + array([-10, 0]), get_last_canvas()
        g = [self.Draw.graph(x0, v, '', c, fill_color=col, color=col, draw_opt='le3', lw=2, opacity=.2) for v, col, x0 in [(a, 2, xa), (b, 4, [0, self.MaxBunch + 2])]]
        self.Draw.legend(g, ['fit', 'PSI bunch spacing'], w=.3, scale=1.2)
        return mean_sigma(y[-10:])

    def draw_spacings(self, bin_size=.5, show=True):
        t = self.get_times(flat=True)
        c1, c2 = self.get_bunch_cut(1, fit=1), self.get_bunch_cut(fit=1)
        ecut = (self.get_npeaks() > 1) & (self.get_n_additional(1, fit=1) == 1) & self.p2ecut(c1) & self.p2ecut(c2)
        t = t[c2 & self.e2pcut(ecut)] - t[c1 & self.e2pcut(ecut)].repeat(self.get_n_additional(fit=1)[ecut])
        self.Draw.distribution(t, self.get_binning(bin_size), 'Peak Spacing', x_tit='Peak Spacing [ns]', **Draw.mode(2), show=show, stats=set_statbox(entries=True))

    def draw_overlayed_times(self, bin_size=.1, off=0, cut=True, thresh=None, fit=True, **kwargs):
        t = self.get_corrected_times(cut, thresh, fit)
        t = (t[t > self.get_bunch_time(n=1).n + 5] + off) % self.get_spacing().n
        return self.Draw.distribution(t, make_bins(0, self.BunchSpacing, bin_size), 'Overlayed Times', x_tit='Peak Time [ns]', y_off=1.8, lm=.12, stats=set_statbox(entries=True), **kwargs)

    def draw_additional(self, h=None, show=True):
        """draw heights of the additinoal peaks in the peak time distribution"""
        h = choose(h, self.draw(show=False, bin_size=.2, fit=True))
        y = get_hist_vec(h, err=False)[int(self.StartAdditional * 2.5):]
        peaks = find_peaks(y, height=max(y) / 2, distance=self.BunchSpacing * 2.5)[0] + int(self.StartAdditional * 2.5)
        fits = array([FitRes(h.Fit('gaus', 'qs0', '', h.GetBinCenter(i - 10), h.GetBinCenter(i + 10))) for i in peaks.tolist()])
        g = self.Draw.graph(fits[:, 1], fits[:, 0], title='Additional Peak Heights', show=show, gridy=True, **Draw.mode(2))
        g.Fit('pol0', 'qs')
        format_statbox(g, fit=True)
        return fits[:, 0]

    def draw_n_per_bunch(self, thresh=None, fit=True, show=True):
        n = [count_nonzero(self.get_bunch_cut(i, thresh=thresh, fit=fit) & self.get_pulser_cut(thresh)) for i in self.get_bunch_nrs()]
        self.Draw.graph(self.get_bunch_nrs(), [ufloat(v, sqrt(v)) for v in n], 'Peaks per Bunch', x_tit='Bunch Number', y_tit='Number of Peaks', **Draw.mode(2), show=show)

    def draw_bunch_height(self, n=1, bin_width=5, cut=None, thresh=None, fit=False, all_=False, **kwargs):
        x = self.get_bunch_heights(n, cut=cut, fit=fit, thresh=thresh, all_=all_)
        return self.Draw.distribution(x, self.Bins.get_pad_ph(bin_width), f'Peak Height B{n}', x_tit='Peak Height [mV]', **kwargs)

    def draw_bunch_times(self, n=1, bin_width=.2, cut=None, fit=False, **kwargs):
        x, bins = self.get_bunch_times(n, cut=cut, fit=fit), self.get_binning(bin_width, n)
        return self.Draw.distribution(x, bins, f'Peak Times B{n}', x_tit='Peak Time [ns]', **kwargs, lm=.12, y_off=1.7, stats=set_statbox(.43, all_stat=True))

    def _add_ph(self, n=1, c=None, bin_size=.2, fit=False, y_range=None, cft=False, show=True):
        if show:
            c = choose(c, get_last_canvas())
            g = self.draw_bunch_height_vs_time(n, bin_size, y_range=y_range, fit=fit, show=False, cft=cft, x_range=self.get_bunch_range(n, n + 1), l_off_x=1)
            format_histo(g, x_tit='', title=' ')
            Draw.tpad(transparent=True, lm=.12, rm=.12, c=c)
            g.Draw('apy+')
            return update_canvas(c)

    def draw_bunch_height_vs_time(self, n=1, bin_width=.2, cut=None, fit=False, cft=False, **kwargs):
        x, y = self.get_bunch_times(n, cut=cut, fit=fit, cft=cft), self.get_bunch_heights(n, cut=cut, fit=True if cft else fit)
        p = self.Draw.profile(x, y, self.get_binning(bin_width, n), y_tit='Peak Height [mV]', x_tit='Peak Time [ns]', show=False)
        (x, y), n = get_hist_vecs(p), get_h_entries(p)
        cut = n > max(n) / 1000
        return self.Draw.graph(x[cut], y[cut], y_tit='Peak Height [mV]', x_tit='Peak Time [ns]', **kwargs)

    def compare_bunch_heights(self, n, bin_width=5, fit=True, thresh=None, fill=None, show=True):
        h = [self.draw_bunch_height(i, bin_width, thresh=thresh, fit=fit, show=False, all_=i < 0) for i in make_list(n)]
        self.Draw.stack(h, 'Peak Height Comparison', [f'B-{i if i > 0 else "All"}' for i in make_list(n)], scale=True, show=show, fill=fill, opacity=.3)

    def compare_times(self, t1, t2, n=1, fill=None, show=True):
        """draw two distributions with different time cuts."""
        h = [self.draw_bunch_height(n, cut=self.get_time_cut(*t, n), show=False) for t in [t1, t2]]
        self.Draw.stack(h, 'Peak Height Comparison', [f'{t[0]}ns - {t[1]}ns' for t in [t1, t2]], scale=True, show=show, fill=fill, opacity=.3)

    def draw_map(self, i=0, fit=True, res=None, cut=None, show=True, hitmap=False):
        ecut = self.Ana.get_event_cut(self.Cut['tracks'])
        v = self.get_bunch_values(1, excl=0, fit=fit, cut=self.e2pcut(ecut), i=i)
        tcut = (ecut & self.p2ecut(self.get_isolated_cut(1, fit=fit)) & (1 if cut is None else self.p2ecut(cut)))[ecut]
        x, y = self.get_tree_vec(self.Ana.get_track_vars())
        if hitmap:
            p = self.Draw.histo_2d(x[tcut], y[tcut], self.Ana.Bins.get_global(res), f'Peak {"Height" if i else "Time"} Map', show=show, draw_opt='colz')
        else:
            p = self.Draw.prof2d(x[tcut], y[tcut], v, self.Ana.Bins.get_global(res), f'Peak {"Height" if i else "Time"} Map', show=show, draw_opt='colz')
        format_histo(p, x_tit='track x [cm]', y_tit='track y [cm]', z_tit=f'Peak {"Height [mV]" if i else "Time [ns]"}', ncont=20, ndivy=510, ndivx=510)
        self.Ana.draw_fid_cut()
        return p

    def draw_time_map(self, fit=True, res=None, show=True):
        return self.draw_map(0, fit, res, show)

    def draw_height_map(self, fit=True, res=None, show=True):
        return self.draw_map(1, fit, res, show)

    def draw_bunch_heights(self, thresh=150, fit=True, cut=None, show=True):
        self.PBar.start(self.MaxBunch)
        x, heights = self.get_bunch_nrs(1), [self.get_bunch_height(i, thresh=thresh, cut=self.get_pre_bunch_cut(i, cut), fit=fit, isolated=cut is None) for i in self.get_bunch_nrs(1)]
        g = Draw.make_tgrapherrors(x, heights, x_tit='Bunch', y_tit='Peak Height [mV]')
        ga = Draw.make_tgrapherrors([ufloat(mean(x[3:]), (x[-1] - x[3]) / 2)], [mean_sigma(heights[3:])[0]], color=2)
        return self.Draw.multigraph([g, ga], 'Bunch Heights', ['single', 'average'], w=2, show=show, color=False, gridy=True)

    def draw_pre_bunch_heights(self, b=-1, thresh=150, y_range=None):
        g = [ig for cut in [None, b] for ig in self.draw_bunch_heights(thresh=thresh, show=False, cut=cut).GetListOfGraphs()]
        tits = ['single peak', 'single average', 'peak in bunch {}'.format(b), 'B{} average'.format(b)]
        self.Draw.multigraph(g, 'Bunch Peak Heights', tits, w=2, y_range=y_range, gridy=True)

    def draw_height_evolution(self, n=1, excl=0, bin_size=20000, thresh=0):
        x = self.get_tree_vec(self.get_t_var())[self.p2ecut(self.get_isolated_cut(n))]
        x, y = x[int(excl * x.size):], self.get_bunch_heights(n, excl=excl).astype('d')
        x = append(x, zeros(y.size - x.size))
        cut = (y > thresh) & (y < 500)
        bins = x[cut][::bin_size].astype('d')
        p = self.Draw.profile(x[cut], y[cut], [bins.size - 1, bins], x_tit='Time [hh:mm]', y_tit='Pulse Height [mV]', t_ax_off=0)
        format_histo(p, y_range=ax_range(get_hist_vec(p), fl=.2, fh=.5))

    def draw_tot(self, n=1, show=True):
        x = self.get_bunch_tot(n)
        self.Draw.distribution(x, title='Time Over 75% Peak Height', x_tit='ToT [ns]', show=show, thresh=.1)
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region FIND
    def find_bunches(self, center=False, bin_size=.5, show=False):
        h = self.draw(corr=False, show=show, bin_size=bin_size)
        values = get_hist_vec(self.draw(corr=False, show=False, bin_size=bin_size))[self.StartAdditional:]
        peaks, d = find_peaks([v.n for v in values], height=max(values).n / 2., distance=self.BunchSpacing)
        peaks += self.StartAdditional
        fit_peaks = []
        for p, ph in zip(peaks, d['peak_heights']):
            low, hi = int(p - 4 / bin_size), int(p + 4 / bin_size)
            f = h.Fit('gaus', 'qs', '', *[h.GetBinCenter(ibin) for ibin in [h.FindFirstBinAbove(ph / 3, 1, low, hi), h.FindLastBinAbove(ph / 2, 1, low, hi)]])
            fit_peaks.append(f.Parameter(1))
        fit_peaks = array(fit_peaks) - (self.BunchSpacing / 2. if not center else 0)
        return concatenate([fit_peaks, [fit_peaks[-1] + self.BunchSpacing]])

    def _find_all(self, ind, tc, thresh=None, fit=False):
        """find peaks for a subset. used for parallelising."""
        f = TF1('lan', 'landau', 0, 512)
        pbar = PBar(ind.size) if not ind[0] else None
        values = []
        for i in ind:
            values.append(self.find(i, tc[i], thresh, fit=f if fit else False))
            if i % 100 and not ind[0]:
                pbar.update(i)
        if pbar:
            pbar.finish()
        return values

    def find_all(self, thresh=None, fit=False, redo=False):
        hdf5_path = self.make_simple_hdf5_path(suf=f'{choose(thresh, self.Threshold):.1f}_{int(fit)}', dut=self.Channel)
        self.WF.get_all()
        if file_exists(hdf5_path) and not redo:
            f = h5py.File(hdf5_path, 'r')
            return f['times'], f['heights'], f['n_peaks']
        if file_exists(hdf5_path):
            remove_file(hdf5_path)

        with Pool() as pool:
            tc = self.WF.get_trigger_cells()
            result = pool.starmap(self._find_all, [(i, tc, thresh, fit) for i in array_split(arange(tc.size), cpu_count())])
            times = [tup[0] for lst in result for tup in lst]
            heights = [tup[1] for lst in result for tup in lst]
            f = h5py.File(hdf5_path, 'w')
            f.create_dataset('times', data=(concatenate(times) * 100).astype('u2'))  # save five digits of the value (reduce file size), max 51200
            f.create_dataset('heights', data=(concatenate(heights) * 100).astype('u2'))
            f.create_dataset('n_peaks', data=array([ar.size for ar in heights]).astype('u1'))
            return f['times'], f['heights'], f['n_peaks']

    def find(self, i, tc, thresh=None, fit=False):
        values = self.WF.get_all()[i]
        peaks = find_peaks(values, height=self.Threshold if thresh is None else thresh, distance=self.Distance, prominence=self.Prominence)
        return self.fit_landau(values, tc, peaks[0], fit) if fit else ([self.WF.get_calibrated_time(tc, value) for value in peaks[0]], peaks[1]['peak_heights'])

    def find_pars(self, redo=False):
        def f():
            t, h, n = self.get_all(cut=..., fit=True)
            peak_info, tc = split(column_stack([t, h]), cumsum(n)[:-1]), self.WF.get_trigger_cells()
            data = []
            self.PBar.start(n[n].size)
            rw, rwi = 10 / self.BinWidth, int(10 / self.BinWidth)
            delay = int(round(self.WF.get_average_rise_time()) / self.BinWidth)
            for i in where(n)[0]:
                v, t = self.WF.get_all()[i], self.WF.get_calibrated_times(tc[i])
                for pt, ph in peak_info[i]:
                    cut = (t > pt - rw) & (t < pt + 1)
                    iv, it = v[cut], t[cut]
                    j = where(iv > ph / 2)[0][0]  # find first value above half the peak height
                    k = where(t > pt)[0][0]  # get peak index
                    data.append([self.find_width(iv, it, pt, ph, j),
                                 self.find_slope(iv, it, j),
                                 self.find_cft(v, t, k, delay),
                                 self.find_tot(v[k - rwi:k + 2 * rwi], t[k - rwi:k + 2 * rwi], .75 * ph)])
                if i % 1000:
                    self.PBar.update(i)
            return array(data).astype('d')
        return do_hdf5(self.make_simple_hdf5_path('Pars'), f, redo)

    @staticmethod
    def find_width(values, times, pt, ph, i):
        return pt - get_x(*times[i - 1:i + 1], *values[i - 1:i + 1], ph / 2) if i else 0

    @staticmethod
    def find_slope(values, times, i):
        return polyfit(times[i - 2:i + 2].astype('d'), values[i - 2:i + 2].astype('d'), deg=1)[0] if i > 1 else 0

    @staticmethod
    def find_cft(values, times, i, delay, fac=.3):
        imin, imax = max(i - 2 * delay, 0), i + delay
        v, t = values[imin:imax], times[imin:imax]
        v -= roll(v, -delay) * fac  # subtract a shifted reduced waveform
        j = where(v < 0)[0]
        j = j[-1] if j.size else None
        return -999 if j is None or j >= 3 * delay - 1 else get_x(t[j], t[j + 1], v[j], v[j + 1], 0)  # interpolate the zero crossing

    @staticmethod
    def find_tot(values, times, thresh):
        i, j = where(values >= thresh)[0][[0, -1]] if any(values > thresh) else None, None
        return (times[-1] if j == times.size - 1 else get_x(*times[j:j + 2], *values[j:j + 2], thresh)) - get_x(*times[i:i + 2], *values[i:i + 2], thresh) if i else -999

    def find_particle_pos(self, h=None, off=0, bin_size=.1):
        h = self.draw_overlayed_times(bin_size, off, show=False) if h is None else h
        x0, w0 = fit_fwhm(h)[1:]
        spec = TSpectrum(6)
        x = [spec.GetPositionX()[i] for i in range(spec.Search(h, 2, '', .001))][1:]
        return append(x0, sorted([FitRes(h.Fit('gaus', 'qs0', '', ix - w0.n, ix + w0.n))[1] for ix in x]))

    def find_composition(self, off=0, bin_size=.1):
        h = self.draw_overlayed_times(bin_size, off, show=False)
        w = fit_fwhm(h)[2].n
        ints = array([Gauss(h, [ix.n - w, ix.n + w]).fit(show=1).get_integral(-10, 30) for ix in self.find_particle_pos(h)]) / bin_size
        ints = array([i + ufloat(0, sqrt(i.n)) for i in ints])
        n, nb, nc = h.Integral(), ints[1], mean_sigma(ints[3:])[0]
        na = n - 2 * nb - 2 * nc  # all events but the four additional peaks
        pa = usqrt(na + usqrt(na ** 2 - 4 * (nb ** 2 + nc ** 2))) / sqrt(2 * n)
        return pa, nb / n / pa, nc / n / pa

    def find_momentum(self, off=0, bin_size=.1, e_pl=1):
        h, s = self.draw_overlayed_times(bin_size, off, show=False), self.get_spacing()
        x, r = self.find_particle_pos(h), self.Momentum + array([-10, 10])
        x = sort((x[1:] - x[0]) % s)
        f = self.Draw.make_tf1(None, lambda p: sum(i.n / i.s for i in (x - self.get_time_differences(None, p, s)) ** 2), *r)
        y0, p0 = f.GetMinimum(), f.GetMinimumX()
        ey = sum(i / i.s for i in (x - self.get_time_differences(None, p0, s)) ** 2).s
        e_sys = self.Draw.make_tf1(None, lambda p: sum(i.n / i.s for i in (x - self.get_time_differences(self.PathLength + e_pl, p, s)) ** 2), *r).GetMinimumX()
        return p0, abs(f.GetX(y0 + ey) - p0), abs(e_sys - p0)
    # endregion FIND
    # ----------------------------------------

    # ----------------------------------------
    # region cft
    def draw_cft(self, n=1, bin_size=.2, smear=None, **kwargs):
        x, bins = self.smear_times(self.get_bunch_cfts(n), smear), self.get_binning(bin_size, n)
        return self.Draw.distribution(x, bins, f'{30}% Constant Fraction Times', lm=.12, x_tit='Constant Fraction Time [ns]', y_off=1.8, **kwargs)

    def draw_height_vs_cft(self, n=1, bin_width=.2, show=True):
        return self.draw_bunch_height_vs_time(n, bin_width, cft=True, show=show)

    def draw_cft_vs_time(self, n=1, bin_size=.2, show=True):
        x, y, bins = self.get_bunch_times(n, fit=True), self.get_bunch_cfts(n), self.get_binning(bin_size, n) * 2
        tit, xtit, ytit = 'Constant Fraction vs. Peak Time', 'Peak Time [ns]', 'Constrant Fraction Time [ns]'
        self.Draw.histo_2d(x, y, bins, tit, x_tit=xtit, y_tit=ytit, stats=set_statbox(entries=True, w=.2), show=show)
    # endregion cft
    # ----------------------------------------

    # ----------------------------------------
    # region MODEL
    def find_scale(self, n=100, redo=False, show=False):
        def f():
            cut = self.get_isolated_cut(n=1)
            ind, t, h = where(self.p2ecut(cut))[0][:n], self.get_bunch_times(n=1)[:n], self.get_bunch_heights(n=1)[:n]
            v = array([FitRes(self.WF.draw_single(ind=j, show=False).Fit('landau', 'qs0', '', t[i] - 3, t[i] + 4))[0] for i, j in enumerate(ind)])
            return FitRes(self.Draw.graph(h, v, 'Model Scale', show=show).Fit('pol1', 'qs'))[1]
        return do_pickle(self.make_simple_pickle_path('ModelScale'), f, redo=redo or show)

    @staticmethod
    def _signal0(x, height, peak_time, rise_time, rise_fac):
        x0, y0 = peak_time, height
        x1, x2 = x0 - rise_time, x0 + rise_fac * rise_time
        p1 = get_p1(x0, x1 if x < x0 else x2, y0, 0)
        return p1 * x + get_p0(x0, y0, p1) if x1 <= x <= x2 else 0

    def signal0(self, height, peak_time, rise_time, rise_fac=3, noise=5):
        x0, x1 = ax_range(peak_time - rise_time, peak_time + rise_time * rise_fac, .5, .5)
        x = arange(x0, x1, self.BinWidth, dtype='d')
        y = array([self._signal0(ix, height, peak_time, rise_time, rise_fac) for ix in x])
        return x, y + normal(scale=noise, size=x.size)

    def signal1(self, height, peak_time, scale, landau_width=3, noise=5):
        x0, x1 = ax_range(peak_time - 2 * landau_width, peak_time + 4 * landau_width, .5, .5)
        x = arange(x0, x1, self.BinWidth, dtype='d')
        y = array([height * scale * TMath.Landau(ix, peak_time, landau_width) for ix in x])
        return x, y + normal(scale=noise, size=x.size)

    def get_signal(self, n, *args, **kwargs):
        return self.signal0(*args, **kwargs) if not n else self.signal1(*args, **kwargs)

    def draw_model_signal(self, model=0, *args, **kwargs):
        x, y = self.get_signal(model, *args, **kwargs)
        return self.Draw.graph(x, y, 'Model Signal', **Draw.mode(3))

    def draw_model_signal1(self, height=None, peak_time=None, landau_width=3):
        self.draw_model_signal(1, choose(height, self.get_bunch_height()), choose(peak_time, self.get_bunch_time().n), None, landau_width)

    def _get_peak_time(self, i, model, *args, **kwargs):
        x, y = self.get_signal(model, *args, **kwargs)
        if i % 1000:
            self.PBar.update(i)
        return x[argmax(y)]

    def draw_raw_model(self, n=1e5, model=1, *args, **kwargs):
        pt, ph, noise = self.get_bunch_time().n, self.get_bunch_height().n, self.Ana.Pedestal.get_raw_noise().n
        self.PBar.start(int(n))
        times = array([self._get_peak_time(i, model, ph, pt, noise=noise, *args, **kwargs) for i in range(self.PBar.N)])
        return self.Draw.distribution(times, self.get_binning(n=1), 'Raw Model Peak Times', x_tit='Signal Peak Time [ns]', stats=set_statbox(entries=True, w=.2))

    def _get_model_data(self, i, model, cft, *args, **kwargs):
        x, y = self.get_signal(model, *args, **kwargs)
        if i % 100 == 0:
            self.PBar.update(i)
        return array([self.find_cft(y, x, argmax(y), delay=10) if cft else x[argmax(y)], max(y).n])

    def model(self, n=1e6, model=1, noise=None, cft=False, redo=False, *args, **kwargs):
        n = int(n)
        hdf5_path = self.make_simple_hdf5_path('M', f'{n}_{model}_{int(cft)}')
        if file_exists(hdf5_path) and not redo:
            f = h5py.File(hdf5_path, 'r')
            return f['times'], f['heights']
        if redo and file_exists(hdf5_path):
            remove_file(hdf5_path)
        rLandau(5, 2)  # first value is always crap
        noise = choose(noise, self.Ana.Pedestal.get_raw_noise().n)
        mpv, w = [v.n for v in self.get_bunch_mpv(par=[1, 2])]
        pt, pts = [v.n for v in self.get_bunch_time(par=[1, 2])]
        self.PBar.start(n)
        peak_times = normal(scale=pts, size=n) + pt
        values = array([self._get_model_data(i, model, cft, min(500, rLandau(mpv, w)), peak_times[i], noise=noise, *args, **kwargs) for i in range(self.PBar.N)])
        self.PBar.finish()
        f = h5py.File(hdf5_path, 'w')
        f.create_dataset('times', data=values[:, 0].astype('f2'))
        f.create_dataset('heights', data=values[:, 1].astype('f4'))
        return f['times'], f['heights']

    def model1(self, n=1e6, redo=False, scale=None, landau_width=3, cft=False):
        scale = self.find_scale() if scale is None else scale
        return self.model(n, model=1, cft=cft, redo=redo, scale=scale, landau_width=landau_width)

    def draw_model(self, n=1e5, bin_width=.2, cft=False, redo=False, **kwargs):
        x, y = self.model1(n, cft=cft, redo=redo)
        return self.Draw.profile(x, y, self.get_binning(bin_width, n=1), y_tit='Peak Height [mV]', x_tit='Peak Time [ns]', **kwargs)
    # endregion MODEL
    # ----------------------------------------

    def fit_landau(self, values, tc, peak_indices, f):
        t, h = [], []
        for ip in peak_indices:
            x, y = self.WF.get_calibrated_times(tc)[max(0, ip - 6):ip + 8], values[max(0, ip - 6):ip + 8]
            g = TGraph(len(x), x.astype('d'), y.astype('d'))
            g.Fit(f, 'q0')
            t.append(f.GetMaximumX())
            h.append(f.GetMaximum())
        return t, array(h)

    @staticmethod
    def smear_times(times, width=2.5, n=5, gaus=False):
        if width is not None:
            if gaus:  # only Gaussian smear
                times += normal(0, width, times.size) if width else 0  # gaussian width
            else:  # flat plus Gaussian smear
                times += rand(times.size) * width - width / 2 if width else 0
                times[::n] += normal(0, width / 2, times.size // 5 + 1)[:times[::n].size] if width else 0
        return times
