#!/usr/bin/env python
# --------------------------------------------------------
#       class for the pulse height calibration of the digital CMS pixel chips
# created on August 2nd 2021 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import gRandom
from helpers.utils import deepcopy, arange, array, mean
from helpers.draw import format_statbox, format_histo, Draw
from numpy.random import randint
from src.sub_analysis import SubAnalysis
from src.dut import Plane


class HighThresh(SubAnalysis):

    def __init__(self, pix_analysis):
        super().__init__(pix_analysis, pickle_dir='HighThresh')

    def find_landau(self, aver=10, m1=2500, m2=5000, s1=500, s2=1600):
        seed = self.Ana.draw_signal_distribution(show=False)
        h = deepcopy(seed)
        m_range = range(m1, m2 + 1, 100)
        s_range = range(s1, s2 + 1, 50)
        p = TProfile2D('g_fl', 'Find Landau', len(m_range) - 1, m_range[0], m_range[-1], len(s_range) - 1, s_range[0], s_range[-1])
        self.PBar.start(len(m_range) * len(s_range) * aver)
        i = 0
        r_min, r_max = .5, 1.5
        for _ in range(aver):
            for m in m_range:
                for s in s_range:
                    i += 1
                    self.PBar.update(i)
                    if r_min < m / 4. / s < r_max:
                        delta = self.model_landau(seed, h, m, s, show=False, thresh=True)
                        p.Fill(m, s, delta)
        self.PBar.finish()
        format_statbox(h, entries=True, x2=.82)
        format_histo(p, x_tit='MPV [e]', y_tit='Sigma [e]', z_tit='#chi^{2} to Seed Function', y_off=1.7, z_off=1.3)
        self.Draw(p, draw_opt='colz', lm=.13, rm=0.16)
        self.draw_ms_ratios(r_min, r_max, m1, m2, s1, s2)
        self.Draw.save_plots('FindLandau')
        self.find_working_point(p)

    @staticmethod
    def draw_ms_ratios(r_min, r_max, m1, m2, s1, s2, step=.1):
        ratios = arange(r_min, r_max + step, step)
        off = .01
        for i, ratio in enumerate(ratios):
            cut = TCutG('ms{n}'.format(n=i), 2, array([0, 20000], 'd'), array([0, 20000 / ratio / 4]))
            x_pos = m1 if m1 / ratio / 4 > s1 else 4 * ratio * s1
            y_pos = m1 / ratio / 4 if m1 / ratio / 4 > s1 else s1
            Draw.tlatex(x_pos + off * (m2 - m1), y_pos + off * (s2 - s1), text='{0:3.1f}'.format(ratio), size=.02, align=11)
            cut.Draw('same')
            Draw.add(cut)

    @staticmethod
    def find_working_point(h):
        ps = [h.ProfileY(), h.ProfileX()]
        fits = [TF1('f{n}'.format(n=i), 'pol2', 0, 10000) for i in range(2)]
        for fit, p in zip(fits, ps):
            p.Fit(fit, 'qs0')
        mins = [fit.GetMinimumX() for fit in fits]
        print(mins, mins[1] / mins[0] / 4)

    def model_landau(self, seed=None, h=None, m=10000, s=1000, show=True, thresh=False):
        seed = self.Ana.draw_signal_distribution(show=False) if seed is None else seed
        h = deepcopy(seed) if h is None else h
        h.SetName('h_ml')
        n = seed.GetEntries()
        h.Reset()
        thresholds = [47 * i for i in [150, 160, 170]]
        for _ in range(int(n)):
            v = gRandom.Landau(m, s)
            threshold = thresholds[int(gRandom.Rndm() * len(thresholds))]
            h.Fill(v if v > threshold else 0 if thresh else v)
            # h.Fill(v if v > 11421 else 0 if thresh else v)
        delta = mean([(h.GetBinContent(i) - seed.GetBinContent(i)) ** 2 for i in range(h.GetNbinsX()) if i is not h.FindBin(0)])
        seed.SetFillColor(2)
        h.SetFillColor(Draw.FillColor)
        seed.SetFillStyle(3017)
        h.SetFillStyle(3018)
        l1 = Draw.make_legend(y2=.76)
        l1.AddEntry(h, 'Simulation', 'f')
        l1.AddEntry(seed, 'Original', 'f')
        self.Draw(h, show=show, leg=l1)
        self.Draw(seed, show=show, draw_opt='same', canvas=gROOT.GetListOfCanvases()[-1])
        seed.Draw('same')
        return delta

    def landau_vid(self, save=False, mpv=5000, sigma=820):
        h = self.Ana.draw_signal_distribution()
        h.GetYaxis().SetRangeUser(0, 2500)
        zero_bin = h.FindBin(0)
        zero_entries = int(h.GetBinContent(zero_bin))
        entries = int(h.GetEntries())
        c = gROOT.GetListOfCanvases()[-1]
        thresholds = self.Ana.Calibration.get_thresholds(arange(0, Plane.NCols), arange(0, Plane.NRows), vcal=False).T[2].resphape((Plane.NCols, Plane.NRows))
        for i in range(entries):
            h.SetBinContent(zero_bin, zero_entries)
            v = gRandom.Landau(mpv, sigma)
            if v < thresholds[randint(0, Plane.NCols - 1), randint(0, Plane.NCols - 1)]:
                h.Fill(v)
                zero_entries -= 1
            if i % 100 == 0:
                c.Update()
                c.Modified()
            if i % 100 == 0 and save:
                self.Draw.save_plots(f'l{i:04d}', show=False, print_names=False)
