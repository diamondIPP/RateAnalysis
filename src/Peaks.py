#!/usr/bin/env python
# --------------------------------------------------------
#       Peak analysis of the high rate pad beam tests at PSI
# created on June 7th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import Analysis, format_histo
from ROOT import TH1F, TCut, gROOT
from utils import fit_poissoni, warning
from math import pi


class PeakAnalysis(Analysis):

    def __init__(self, pad_analysis):
        self.Ana = pad_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)
        self.Run = self.Ana.Run
        self.Channel = self.Ana.Channel
        self.DUT = self.Ana.DUT.Number
        self.Tree = self.Ana.Tree
        self.set_save_directory(self.Ana.SubDir)

    def draw_n_peaks(self, spec=False, show=True, fit=False):
        h = TH1F('h_pn', 'Number of Peaks', 10, 0, 10)
        draw_var = '@peaks{ch}_x.size()>>h_pn' if spec and self.Run.has_branch('peaks{ch}_x'.format(ch=self.Channel)) else 'n_peaks[{ch}] - 1>>h_pn'
        self.Tree.Draw(draw_var.format(ch=self.Channel), self.Ana.Cut(), 'goff')
        self.format_statbox(only_fit=True, w=.3) if fit else self.format_statbox(entries=True)
        format_histo(h, x_tit='Number of Peaks', y_tit='Number of Entries', y_off=1.4, fill_color=self.FillColor, lw=2)
        self.save_histo(h, 'PeakNumbers', show, logy=True, lm=.11)
        return h

    def draw_n_peaks_fit(self, show=True, spec=False):
        h = self.draw_n_peaks(spec, show=False, fit=True)
        format_histo(h, 'Fit Result')
        self.draw_histo(h, '', show, logy=True, lm=.11)
        fit = fit_poissoni(h, show=show)
        self.save_histo(h, 'PeakNumbersFit', show, logy=True, lm=.11, draw_opt='e1same', canvas=gROOT.GetListOfCanvases()[-1])
        self.get_flux(fit=fit)
        return fit

    def get_flux(self, fit=None):
        fit = self.draw_n_peaks_fit(show=False) if fit is None else fit
        lambda_ = fit.GetParameter(1)
        err = fit.GetParError(1)
        n_bunches = 25
        proc_frequency = 50.6e6
        flux = lambda_ * proc_frequency / (n_bunches * self.get_area()) / 1000.
        err *= proc_frequency / (n_bunches * self.get_area()) / 1000.
        print 'Estimated Flux by number of peaks: {f:0.4f} kHz/cm2'.format(f=flux)
        return flux, err

    def get_area(self):
        """ return the total area of the BCM' pad sizes """
        i = int(self.Ana.DUT.Name.split('-')[-1]) - 1
        base_length = 0.0928125  # [cm]
        spacing = 0.0025
        radius = 0.0049568
        rounded_edge = radius ** 2 * (4 - pi)
        area = base_length ** 2
        return 2 ** i * area + self.get_spacings(i, spacing, base_length) - rounded_edge

    def get_spacings(self, i, spacing, length):
        """ return the additional spacings for the BCM' pad sizes """
        if i > 0:
            j = 2 ** (i / 2)
            return spacing * (j * length + (j - 1) * spacing) + 2 * self.get_spacings(i - 1, spacing, length)
        else:
            return 0

    def draw_positions(self, cut='', corr=False, show=True):
        if not self.Ana.has_branch('peak_positions'):
            warning('The peak_positions branch does not exist!')
            return
        h = TH1F('h_pt', 'Peak {m}'.format(m='Timings' if corr else 'Positions'), 1024, 0, 512 if corr else 1024)
        self.Tree.Draw('peak_{p}[{c}]>>h_pt'.format(c=self.Channel, p='positions' if not corr else 'times'), TCut(cut) + TCut('!pulser'), 'goff')
        format_histo(h, x_tit='Time [ns]' if corr else 'Digitiser Bin', y_tit='Number of Entries', y_off=.4, fill_color=self.FillColor, lw=1, tit_size=.05, stats=0)
        self.save_histo(h, 'PeakTimings', show, logy=True, lm=.045, rm=.045, x=4, y=.5)
        return h

    def draw_timings(self, cut='', show=True):
        self.draw_positions(cut, corr=True, show=show)

    def draw_max_position(self, cut='', corr=False, show=True):
        h = TH1F('h_pt', 'Max Peak {m}'.format(m='Timings' if corr else 'Positions'), 1024, 0, 512 if corr else 1024)
        cut = TCut(cut) + TCut('!pulser') if 'pulser' not in cut else TCut(cut)
        self.Tree.Draw('max_peak_{p}[{c}]>>h_pt'.format(c=self.Channel, p='position' if not corr else 'time'), cut, 'goff')
        format_histo(h, x_tit='Time [ns]' if corr else 'Digitiser Bin', y_tit='Number of Entries', y_off=.4, fill_color=self.FillColor, lw=1, tit_size=.07, stats=0, lab_size=.06)
        self.save_histo(h, 'MaxPeak{m}'.format(m='Timings' if corr else 'Positions'), show, logy=True, lm=.073, rm=.045, x=2, y=.5)

    def draw_max_timing(self, cut='', show=True):
        self.draw_max_position(cut, corr=True, show=show)

    def check_peaks(self):
        h, n = self.Ana.draw_waveforms(t_corr=False)
        self.Tree.GetEntry(self.Ana.StartEvent + self.Count - 1)
        print n, self.Ana.StartEvent, self.Ana.count
        a = self.Tree.Draw('peak_positions[{ch}]:n_peaks[{ch}]'.format(ch=self.Channel), '', 'goff', 1, self.Ana.StartEvent + self.Ana.count - 1)
        print [self.Tree.GetV1()[j] / 2 for j in xrange(a)], [self.Tree.GetV2()[j] for j in xrange(a)][0]
        a = self.Tree.Draw('peaks3_x', '', 'goff', 1, self.Ana.StartEvent + self.Ana.count - 1)
        print [i / 2. for i in [self.Tree.GetV1()[j] for j in xrange(a)]]
