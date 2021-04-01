#!/usr/bin/env python
# --------------------------------------------------------
#       module to create MC signal distributions
# created on January 12th 2021 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy.random import normal, randint
from numpy import concatenate, cumsum, linspace, where, polyfit, array, append, arange, invert, mean
from ROOT import gRandom
from helpers.draw import Draw, choose, get_base_dir, join, get_hist_vecs, PBar, calc_eff, mean_sigma, ax_range


class MCSignal(object):

    Pedestal = 2
    Noise = 5

    def __init__(self, ana=None, n=1e6):

        self.Ana = ana
        self.HasAna = self.Ana is not None
        self.Draw = Draw(join(get_base_dir(), 'config', 'main.ini'))
        self.PBar = PBar()

        self.Pedestal, self.Noise = ana.Pedestal.get_under_signal(err=False) if self.HasAna else (MCSignal.Pedestal, MCSignal.Noise)

        self.N = int(n)
        self.Rate = 1e7
        self.Area = .3 * .3  # cm²
        self.AreaScint = 1 ** 2  # cm²
        self.NOther = self.Rate * self.Area * 19.8e-9
        self.NScint = self.Rate * (1 - self.Area) * 19.8e-9
        # self.PBucket = (1 - self.Area / self.AreaScint) * (1 - poisson.cdf(1, self.NOther, 1))
        self.PBucket = 0.05246 if ana is None else ana.get_bucket_ratio()
        self.NBuc = round(self.N * self.PBucket)
        self.NCons = None

        # self.B1, self.B2 = self.gen_data()
        # self.BCut = invert((self.B2 > 2 * self.Noise) & (self.B1 < 2 * self.Noise))

    def gen_noise(self, n=None):
        return normal(self.Pedestal, self.Noise, choose(n, self.N))

    def gen_signal(self, n=None, w0=20, w1=20, m0=70, m1=85, p=.6):
        if self.Ana is not None:
            return MCSignal.gen_sig_from_dist(self.Ana.draw_signal_distribution(show=False, prnt=False, evnt_corr=False), n=n)
        g = concatenate([normal(choose(m, 80), choose(w, 20), round(choose(n, self.N) * ip)) for w, ip, m in [(w0, p, m0), (w1, 1 - p, m1)]])
        return [gRandom.Landau(ig, ig / 8) for ig in g]  # assuming R=FWHM/MPV=.5, FWHM=4*width

    @staticmethod
    def gen_sig_from_dist(h, n=None, p_steps=100000):
        """interpolate the cdf in equal prob steps and sample it"""
        n = int(choose(n, 1e6))
        x, y = get_hist_vecs(h, err=False)
        ints = cumsum(y) / h.GetEntries()
        p = linspace(0, 1, p_steps, endpoint=False)
        ip = array([where(ip >= ints)[0][-1] for ip in p])
        fits = polyfit(x[:2], [ints[ip], ints[ip + 1]], deg=1)
        xp = append((p - fits[1]) / fits[0] + (x[ip] - x[0] + (x[1] - x[0]) / 2), max(x))
        return xp[randint(0, p.size, size=n)]

    def draw_signal_distribution(self, cut=...):
        x = self.get_data()[0][cut]
        bins = self.Ana.Bins.get_pad_ph(2) if self.HasAna else None
        self.Draw.distribution(x, bins, x_range=[-50, 400], x_tit='Pulse Height [mV]', lm=.15, y_off=2, w=1.2)

    def get_data(self):
        b, s = self.get_bucket(), self.get_signal()
        return concatenate([b[0], s[0]]), concatenate([b[1], s[1]])

    def get_bucket(self):
        return self.gen_noise(self.NBuc), self.gen_signal(self.NBuc)

    def get_signal(self):
        n_sig = self.N - self.NBuc
        p_cons = 0.025 if self.HasAna is None else self.Ana.Peaks.get_p_extra()
        self.NCons = round(n_sig * p_cons)
        return self.gen_signal(n_sig), concatenate([self.gen_signal(self.NCons), self.gen_noise(n_sig - self.NCons)])

    def get_thresh(self, s):
        return self.Pedestal + s * self.Noise

    def get_bucket_cut(self, data=None, s1=2, s2=3):
        data = [self.get_bucket(), self.get_signal()] if data is None else data
        return [(i[0] < self.get_thresh(s1)) & (i[1] > self.get_thresh(s2)) for i in data]

    def get_bucket_stats(self, data=None, s1=2, s2=3, prnt=True):
        cb, cs = self.get_bucket_cut(data, s1, s2)
        eb, es = calc_eff(values=cb), calc_eff(values=invert(cs))
        if prnt:
            print(f'Sensitivity: {eb:.3f}')
            print(f'Specificity: {es:.3f}')
        return eb, es

    def draw_stats(self, step=.1, max_sig=10):
        s = arange(.5, max_sig + step, step)
        data = [self.get_bucket(), self.get_signal()]
        res = array([self.get_bucket_stats(data, i, i, False) for i in s])
        g = [Draw.make_tgrapherrors(s, res[:, i]) for i in [0, 1]]
        self.Draw.multigraph(g, 'Bucket Cut Stats', ['Sensitivity', 'Specificity'], x_tit='N Sigma', y_tit='Efficiency [%]', draw_opt='l', lw=2)

    def draw_signals(self, step=.1, max_sig=10, rel=True):
        s = arange(0, max_sig + step, step)
        data = [self.get_bucket(), self.get_signal()]
        res = array([mean_sigma(concatenate([data[0][0][invert(cb)], data[1][0][invert(cs)]]))[0] for cb, cs in [self.get_bucket_cut(data, i, i) for i in s]])
        real = mean(data[1][0])
        no_cut = mean(concatenate([data[0][0], data[1][0]])) / (real if rel else 1)
        res /= real if rel else 1
        g = self.Draw.graph(s, res, x_tit='N Sigma', y_tit=f'{"Relative " if rel else ""}Pulse Height {"" if rel else "[mV]"}', draw_opt='ae3', lw=2, fill_color=2, y_range=ax_range(no_cut, max(res).n, .1, .2))
        lc = Draw.horizontal_line(no_cut, w=2)
        lr = Draw.horizontal_line(1 if rel else real, w=2, color=3)
        self.Draw.legend([g, lr, lc], ['w/ cut', 'real signal', 'w/o cut'], ['f', 'l', 'l'])


if __name__ == '__main__':
    z = MCSignal()
