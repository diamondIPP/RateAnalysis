#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import append, column_stack, polyfit, argmax
from src.event_alignment import *
from helpers.draw import Draw


class PixAlignment(EventAligment):
    def __init__(self, converter=None):

        # Info
        self.Threshold = .35
        self.NDutPlanes = 4
        self.TelPlane = 2
        self.DUTPlane = 4
        self.MaxOffset = 1000

        self.Y1, self.Y2 = array([]), array([])  # row arrays
        self.X1, self.X2 = array([]), array([])  # row arrays
        self.C1, self.C2 = array([]), array([])  # cut arrays
        self.BinSize = None

        EventAligment.__init__(self, converter)

    # ----------------------------------------
    # region INIT
    def load_variables(self):
        data = super().load_variables()
        t = self.Run.info('Loading pixel information from tree ... ', endl=False)
        self.init_data()
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

    def init_data(self):
        if not self.Y1.size:
            self.X1, self.X2, self.Y1, self.Y2, self.C1, self.C2 = self.load_data()

    def load_data(self, firstentry=0, nentries=None):
        n = self.get_tree_vec('n_hits_tot', dtype='u1', nentries=nentries, firstentry=firstentry)
        self.InTree.SetEstimate(sum(n))
        pl, col, row = self.get_tree_vec(['plane', 'col', 'row'], dtype='u1', firstentry=firstentry, nentries=nentries)
        c1, c2 = self.get_e_cut(self.TelPlane, pl, n), self.get_e_cut(self.DUTPlane, pl, n)
        cut1, cut2 = c1.repeat(n) & (pl == self.TelPlane), c2.repeat(n) & (pl == self.DUTPlane)
        return col[cut1], col[cut2], row[cut1], row[cut2], c1, c2

    def get_aligned(self, n=100):
        self.init_data()
        x, y = self.correlate_all(n=n)
        return (x > self.Threshold) & (y > self.Threshold)
    # endregion INIT
    # ----------------------------------------

    def get_data(self, offset=0, start=0):
        e1, e2 = where(self.C1)[0] >= start, where(self.C2)[0] >= start
        x1, x2, y1, y2, c1, c2 = self.X1[e1], self.X2[e2], self.Y1[e1], self.Y2[e2], self.C1[start:], self.C2[start:]
        co1, co2 = self.find_events(offset, c1, c2)
        return column_stack([x1[co1], x2[co2]]), column_stack([y1[co1], y2[co2]]), where(c1)[0][co1] + start

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
        s = x.shape[0] // n
        xs, ys = x[:s * n].reshape(s, n, 2), y[:s * n].reshape(s, n, 2)
        i_s = next((i for i in range(xs.shape[0]) if correlate(*xs[i].T) < self.Threshold or correlate(*ys[i].T) < self.Threshold), None)  # first bucket below
        if i_s is None:  # all buckets uncorrelated
            return -.5
        i = next(i for i in range(max(0, (i_s - 1) * n), x.shape[0]) if correlate(*x[i:n + i].T) < self.Threshold and correlate(*y[i:n + i].T) < self.Threshold)  # first index
        if not i:  # very short jump
            c = [correlate(*x[i:n + i].T) < self.Threshold and correlate(*y[i:n + i].T) < self.Threshold for i in range(10)]
            if all(c):
                return 0
            i_min = argmax(invert(c))
            i = next(i for i in range(i_min, x.shape[0]) if correlate(*x[i:n + i].T) < self.Threshold or correlate(*y[i:n + i].T) < self.Threshold)
        return i

    def correlate(self, x, y, n):
        return [correlate(*x[i:n + i].T) < self.Threshold and correlate(*y[i:n + i].T) < self.Threshold for i in range(n)]

    def find_final_offset(self, n=100):
        start = where(self.C1)[0][-6 * n]
        for off in range(self.MaxOffset):
            if not all(self.correlate(*self.get_data(off, start)[:2], n)):
                return off

    def find_offsets(self, off=None, n=50, start=None):
        off = choose(off, self.FinalOffset)
        if off == self.FinalOffset:
            self.PBar.start(off)
        start = start if start is not None else self.find_offsets(off - 1, n) if off > 0 else 0
        if start > 0:
            self.Offsets[start] = off
        if off == self.FinalOffset:
            info(f'found all offsets ({len(self.Offsets) + 1}) :-) ')
            return
        x, y, e = self.get_data(off, start)
        i = self.find_next(x, y, n)
        if i <= 0:
            return start - int(n * i)
        r1, r2 = max(0, int(i - n // 3)), int(i + n // 3)
        xi, (c1, c2) = arange(r1, r2), array([[correlate(*x[i:i + n].T), correlate(*y[i:i + n].T)] for i in range(r1, r2)]).T
        f1, f2 = polyfit(xi, c1, deg=1), polyfit(xi, c2, deg=1)
        return e[int(round(mean([-f[1] / f[0] for f in [f1, f2]])))]

    def correlate_all(self, offset=0, n=None, start=0):
        x, y, e = self.get_data(offset, start)
        s = x.shape[0] // n
        return array([correlate(*i.T) for i in x[:s * n].reshape(s, n, 2)]), array([correlate(*i.T) for i in y[:s * n].reshape(s, n, 2)])

    def find_all(self, offset=0, n=None, start=0):
        x, y = self.correlate_all(offset, n, start)
        return (x < self.Threshold) | (y < self.Threshold)

    def fill_branches(self, ev, offset):
        if not offset:
            event = self.Variables[ev]
        else:
            tel_offset = abs(offset) if offset < 0 else 0
            dut_offset = offset if offset > 0 else 0
            # select the according events
            tel_event, dut_event = self.Variables[ev + tel_offset], self.Variables[ev + dut_offset]
            # merge the right parts of the events
            event = concatenate((tel_event[:, tel_event[0] < self.NDutPlanes], dut_event[:, dut_event[0] >= self.NDutPlanes]), axis=1)
        for i, name in enumerate(self.Branches.keys()):
            for value in event[i]:
                self.Branches[name].push_back(value)

    def draw(self, off=0, bin_size=50):
        y = invert(self.find_all(off, bin_size))
        x = mean(where(self.C1)[0][self.find_events(off)[0]][:y.size * bin_size].reshape(y.size, bin_size), axis=1)
        self.Draw.graph(x, y, f'{off} Alignment', x_tit='Event Number', y_tit='Aligned', draw_opt='al', **Draw.mode(2))


if __name__ == '__main__':

    from pixel.run import PixelRun
    from src.converter import Converter

    pargs = init_argparser(run=147, tc='201810')
    zrun = PixelRun(pargs.run, testcampaign=pargs.testcampaign, load_tree=False, verbose=True)
    z = PixAlignment(Converter(zrun))
