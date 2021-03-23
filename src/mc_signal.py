#!/usr/bin/env python
# --------------------------------------------------------
#       module to create MC signal distributions
# created on January 12th 2021 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy.random import normal, shuffle
from numpy import concatenate, arange, zeros, invert
from ROOT import gRandom
from helpers.draw import Draw, choose, get_base_dir, join
from scipy.stats import poisson


class MCSignal(object):

    def __init__(self):
        self.Noise = 7
        self.Pedestal = -2

        self.N = int(1e6)
        self.Rate = 1e7
        self.Area = .3 * .3  # cm²
        self.AreaScint = 1 ** 2  # cm²
        self.NOther = self.Rate * self.Area * 19.8e-9
        self.NScint = self.Rate * (1 - self.Area) * 19.8e-9
        # self.PBucket = (1 - self.Area / self.AreaScint) * (1 - poisson.cdf(1, self.NOther, 1))
        self.PBucket = 0.05246

        self.Draw = Draw(join(get_base_dir(), 'config', 'main.ini'))

        self.B0 = self.gen_noise()
        self.ESig, self.EBuc = self.gen_events()
        self.B1, self.B2 = self.gen_data()
        self.BCut = invert((self.B2 > 2 * self.Noise) & (self.B1 < 2 * self.Noise))

    def gen_noise(self, n=None):
        return normal(self.Pedestal, self.Noise, choose(n, self.N))

    def gen_signal(self, n=None, w0=20, w1=20, m0=70, m1=85, p=.6):
        g = concatenate([normal(choose(m, 80), choose(w, 20), round(choose(n, self.N) * ip)) for w, ip, m in [(w0, p, m0), (w1, 1 - p, m1)]])
        return [gRandom.Landau(ig, ig / 8) for ig in g]  # assuming R=FWHM/MPV=.5, FWHM=4*width

    def draw_signal(self, x=None, n=None, w0=20, w1=20, m0=80, m1=95, p=.6):
        x = choose(x, self.gen_signal(n, w0, w1, m0, m1, p))
        self.Draw.distribution(x, x_range=[-50, 400], x_tit='Pulse Height [mV]', lm=.15, y_off=2, w=1.2)

    def gen_events(self):
        events = arange(self.N)
        shuffle(events)
        return events[int(self.PBucket * self.N):], events[:int(self.PBucket * self.N)]

    def gen_data(self):
        b1, b2 = zeros((2, self.N))
        b1[self.ESig] = self.gen_signal(self.ESig.size)
        b1[self.EBuc] = self.gen_noise(self.EBuc.size)
        b2[self.EBuc] = self.gen_signal(self.EBuc.size)
        b2_events = poisson.rvs(self.NOther, size=self.ESig.size)
        b2_signals = self.ESig[b2_events > 0]
        b2_noise = self.ESig[b2_events == 0]
        b2[b2_noise] = self.gen_noise(b2_noise.size)
        b2[b2_signals] = self.gen_signal(b2_signals.size)
        return b1, b2

    def get_bcut_eff(self):
        return invert(self.BCut)[self.EBuc].nonzero()[0].size / self.EBuc.size


if __name__ == '__main__':
    z = MCSignal()
