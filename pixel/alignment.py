#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import column_stack, polyfit, argmax
from src.event_alignment import *
from helpers.draw import Draw


class PixAlignment(EventAligment):
    def __init__(self, converter=None, tel_plane=None):

        # Info
        self.Threshold = .35
        self.NTelPlanes = converter.Run.NTelPlanes
        self.DUTPlane = self.find_dut_plane(converter.Run.DUTs)
        self.TelPlane = choose(tel_plane, 1 if self.DUTPlane > self.NTelPlanes else 2)
        self.MaxOffset = 1000

        self.X1, self.X2 = array([]), array([])  # col arrays
        self.Y1, self.Y2 = array([]), array([])  # row arrays
        self.C1, self.C2 = array([]), array([])  # cut arrays
        self.BinSize = None

        EventAligment.__init__(self, converter)

    # ----------------------------------------
    # region INIT
    def find_dut_plane(self, duts):
        return next((i for i, dut in enumerate(duts, self.NTelPlanes) if 'Si' in dut.Name or 'D' in dut.Name), 0)

    def load_variables(self):
        data = super().load_variables()
        t = self.Run.info('Loading pixel information from tree ... ', endl=False)
        self.init_data()
        self.FirstOffset = self.find_first_offset()
        self.FinalOffset = self.find_final_offset()
        self.Run.add_to_info(t)
        return data

    @staticmethod
    def get_e_cut(plane, plane_vec, n):
        """find all events with one hit in [plane]."""
        events = arange(n.size).repeat(n)[plane_vec == plane]  # events with hits in the plan
        dif = append(diff(events).astype('?'), True)
        dif[where(invert(dif))[0] + 1] = False  # select only events with a single hit
        cut = zeros(n.size, '?')
        cut[events[dif]] = True
        return cut

    def init_data(self, tree=None):
        if not self.Y1.size or tree is not None:
            self.X1, self.X2, self.Y1, self.Y2, self.C1, self.C2 = self.load_data(tree=tree)

    def load_data(self, firstentry=0, nentries=None, tree=None):
        tree = choose(tree, self.InTree)
        n = get_tree_vec(tree, self.HitVar, dtype='u1', nentries=nentries, firstentry=firstentry)
        tree.SetEstimate(sum(n))
        pl, col, row = get_tree_vec(tree, ['plane', 'col', 'row'], dtype='u1', firstentry=firstentry, nentries=nentries)
        c1, c2 = self.get_e_cut(self.TelPlane, pl, n), self.get_e_cut(self.DUTPlane, pl, n)
        cut1, cut2 = c1.repeat(n) & (pl == self.TelPlane), c2.repeat(n) & (pl == self.DUTPlane)
        return col[cut1], col[cut2], row[cut1], row[cut2], c1, c2

    def get_aligned(self, tree=None, bin_size=200):
        self.init_data(tree)
        x, y = self.correlate_all(n=bin_size)
        return (x > self.Threshold) & (y > self.Threshold)

    def set_aligned(self, bin_size=200):
        data = [self.get_data(off, start, end) for (start, off), end in zip(self.Offsets.items(), [*self.Offsets][1:] + [self.NEntries])]
        e, xydata = concatenate([d[2] for d in data]), concatenate([d[:2] for d in data], axis=1)
        aligned = invert([self.is_misaligned([ix, iy]) for ix, iy in zip(*PixAlignment.bin_data(*xydata, bin_size))])
        aligned = append(aligned, True if self.NEntries % bin_size else [])  # add last bin
        aligned[roll(invert(aligned), 1) & roll(invert(aligned), -1)] = False  # extend to neighbouring bins
        aligned = aligned.repeat(diff(concatenate([[0], e[::bin_size][1:], [self.NEntries]])))  # calculate how many events are in each bins'
        self.Aligned[:aligned.size] = aligned
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region DATA
    def get_data(self, offset=0, start=0, end=None):
        end = choose(end, self.NEntries)
        e1, e2 = [(e >= start) & (e < end) for e in [where(c)[0] for c in [self.C1, self.C2]]]
        x1, x2, y1, y2, c1, c2 = self.X1[e1], self.X2[e2], self.Y1[e1], self.Y2[e2], self.C1[start:end], self.C2[start:end]
        co1, co2 = self.find_events(offset, c1, c2)
        return column_stack([x1[co1], x2[co2]]), column_stack([y1[co1], y2[co2]]), where(c1)[0][co1] + start

    @staticmethod
    def bin_data(x, y, n):
        s = x.shape[0] // n
        return x[:s * n].reshape(s, n, 2), y[:s * n].reshape(s, n, 2)

    def fill_branches(self, ev, offset):
        dut_ev, tel_ev = ev, ev - offset
        dut_hits, tel_hits_p1 = self.NTot[dut_ev], self.NTot[tel_ev + 1]
        n_dut = where(self.Variables[0][dut_hits:self.NTot[dut_ev + 1]] < self.NTelPlanes)[0].size
        n_tel = where(self.Variables[0][self.NTot[tel_ev]:tel_hits_p1] >= self.NTelPlanes)[0].size
        self.Branches['n_hits_tot'][0][0] = n_dut + n_tel  # n hits
        ind = concatenate([arange(dut_hits, dut_hits + n_dut, dtype='i'), arange(tel_hits_p1 - n_tel, tel_hits_p1, dtype='i')])  # take ind from dut event and tel event
        for i, br in enumerate(self.get_tel_branches()):
            self.Branches[br][0][:n_dut + n_tel] = self.Variables[i][ind]
        self.Branches['trigger_phase'][0][:2] = self.Variables[-1][[2 * dut_ev, 2 * tel_ev + 1]]
        self.Branches['aligned'][0][0] = self.Aligned[ev]
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region OFFSETS
    def find_events(self, offset=0, c1=None, c2=None):
        """ find events (with offsets) which have both exactly one hit. """
        c1, c2 = choose(c1, self.C1), choose(c2, self.C2)
        o1, o2 = (0, offset) if offset < 0 else (-offset, 0)
        c = roll(c1, o1) & roll(c2, o2)
        if offset:
            c[-abs(offset):] = False
        return roll(c, -o1)[c1], roll(c, -o2)[c2]

    @update_pbar
    def find_next(self, x, y, n):
        """:returns: the next uncorrelated event number"""
        xs, ys = PixAlignment.bin_data(x, y, n)  # fill the data in buckets
        i_s = next((i for i in range(xs.shape[0]) if self.is_misaligned([xs[i], ys[i]])), None)  # first bucket below
        if i_s is None:  # all buckets are correlated
            return None
        i = next((i for i in range(max(0, (i_s - 1) * n), x.shape[0]) if self.is_misaligned([x, y], i, n + i)), None)  # first index
        if i == 0:  # very short jump
            c = [correlate(*x[i:n + i].T) < self.Threshold and correlate(*y[i:n + i].T) < self.Threshold for i in range(10)]
            if all(c):
                return 0
            i_min = argmax(invert(c))
            i = next(i for i in range(i_min, x.shape[0]) if correlate(*x[i:n + i].T) < self.Threshold or correlate(*y[i:n + i].T) < self.Threshold)
        return i

    def correlate(self, x, y, n):
        return [correlate(*x[i:n + i].T) < self.Threshold and correlate(*y[i:n + i].T) < self.Threshold for i in range(n)]

    def is_aligned(self, d, start=None, n=None):
        start, n = (0, choose(start, d[0].shape[0])) if n is None else (start, n)
        return correlate(*d[0][start:n].T) > self.Threshold and correlate(*d[1][start:n].T) > self.Threshold

    def is_misaligned(self, d, start=None, n=None):
        start, n = (0, choose(start, d[0].size)) if n is None else (start, n)
        return correlate(*d[0][start:n].T) < self.Threshold and correlate(*d[1][start:n].T) < self.Threshold

    def find_final_offset(self, n=100):
        start = where(self.C1)[0][-6 * n]
        for off in range(self.MaxOffset):
            if not all(self.correlate(*self.get_data(off, start)[:2], n)):
                return off

    def find_first_offset(self, n=50):
        start = where(self.C1[:20 * n] & self.C2[:20 * n])[0][n]
        for off in range(50):
            for sgn in [-1, 1]:
                if self.is_aligned(self.get_data(sgn * off, start, start + 10 * n), 0, n):
                    return off
        warning('could not determine starting offset! assuming 0 ...')
        return 0

    def find_all_offsets(self):
        return self.find_offsets()

    def find_off_event(self, off, start, n):
        x, y, e = self.get_data(off, start)
        i = self.find_next(x, y, n)
        if i is None:
            return
        if i <= 0:
            return start - int(n * i)
        r1, r2 = max(0, int(i - n // 3)), int(i + n // 3)
        xi, (c1, c2) = arange(r1, r2), array([[correlate(*x[i:i + n].T), correlate(*y[i:i + n].T)] for i in range(r1, r2)]).T
        f1, f2 = polyfit(xi, c1, deg=1), polyfit(xi, c2, deg=1)
        return e[int(round(mean([-f[1] / f[0] for f in [f1, f2]])))]

    def find_next_off(self, last_off, start, n=50):
        neg, pos = [self.get_data(last_off + off, start, start + 10 * n) for off in [-1, 1]]
        neg, pos = [[self.is_aligned(v, i, i + n) for i in range(n)] for v in [neg, pos]]
        return last_off + (1 if count_nonzero(pos) > 3 else -1 if count_nonzero(neg) > 3 else 0)

    def find_offsets(self, n=50):
        start, off = 0, self.FirstOffset
        info('STEP 1: Finding the offsets ...')
        self.PBar.start(self.NEntries, counter=False)
        while start is not None:
            self.Offsets[start] = off
            start = self.find_off_event(off, start, n)
            if start is not None:
                off = self.find_next_off(off, start)
                self.PBar.update(start)
        self.PBar.finish()
        info(f'found all offsets ({len(self.Offsets) - (1 if [*self.Offsets.values()][0] == 0 else 0)})! :-)')

    def correlate_all(self, offset=0, n=None, start=0):
        x, y, e = self.get_data(offset, start)
        s = x.shape[0] // n
        return array([correlate(*i.T) for i in x[:s * n].reshape(s, n, 2)]), array([correlate(*i.T) for i in y[:s * n].reshape(s, n, 2)])

    def find_all(self, offset=0, n=None, start=0):
        x, y = self.correlate_all(offset, n, start)
        return (x < self.Threshold) | (y < self.Threshold)
    # endregion OFFSETS
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def get_x(self, y, off=0, bin_size=50):
        return mean(where(self.C1)[0][self.find_events(off)[0]][:y.size * bin_size].reshape(y.size, bin_size), axis=1)

    def draw_hit_rate(self):
        y = self.get_tree_vec(self.HitVar)
        self.Draw.profile(arange(y.size), y, x_tit='Event Number', y_tit='Total Number of Hits', y_range=[0, 1.5 * mean(y)])

    def draw_correlation(self, off=0, bin_size=50, **kwargs):
        yx, yy = self.correlate_all(off, bin_size)
        g = [self.Draw.graph(self.get_x(y, off, bin_size), y, x_tit='Event Number', y_tit='Correlation Factor', show=False) for y in [yx, yy]]
        self.Draw.multigraph(g, 'Correlations', ['xx', 'yy'], **prep_kw(kwargs, draw_opt='l', **Draw.mode(2), y_range=[0, 1.18]))

    def draw(self, off=0, bin_size=50):
        y = invert(self.find_all(off, bin_size))
        self.Draw.graph(self.get_x(y, off, bin_size), y, f'{off} Alignment', x_tit='Event Number', y_tit='Aligned', draw_opt='al', **Draw.mode(2))
    # endregion DRAW
    # ----------------------------------------


if __name__ == '__main__':

    from pixel.run import PixelRun
    from src.converter import Converter

    pargs = init_argparser(run=489, tc='201610')
    zrun = PixelRun(pargs.run, testcampaign=pargs.testcampaign, load_tree=False, verbose=True)
    z = PixAlignment(Converter(zrun))
    z.reload()
