#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import histogram2d, sum
from src.event_alignment import *
from src.binning import make_bins
from helpers.draw import get_hist_vec


class PadAlignment(EventAligment):
    def __init__(self, converter):

        # Info
        self.Threshold = .4
        self.Pulser = array([])
        self.PulserEvents = array([])
        self.FirstOffset = None
        self.BinSize = 30
        self.NMaxHits = 3

        EventAligment.__init__(self, converter)

        if not self.IsAligned:
            self.PulserRate = self.get_pulser_rate()
            self.Threshold = .4 if self.PulserRate > .05 else .2

    # ----------------------------------------
    # region INIT
    def print_start(self):
        print_banner('STARTING PAD EVENT ALIGNMENT OF RUN {}'.format(self.Run.Number))

    @staticmethod
    def init_branches():
        return EventAligment.init_branches() + [('trigger_phase', zeros(1, 'u1'), 'trigger_phase/b')]

    def load_variables(self):
        data = super(PadAlignment, self).load_variables()
        t = self.Run.info('Loading pad information from tree ... ', endl=False)
        self.Pulser = self.get_pulser()
        self.PulserEvents = where(self.Pulser)[0]
        self.FirstOffset = self.find_start_offset()
        tp = get_tree_vec(self.InTree, 'trigger_phase', dtype='u1')
        self.Run.add_to_info(t)
        return data + [tp]

    def check_alignment_fast(self, bin_size=1000):
        """ just check the zero correlation """
        x, y = get_tree_vec(self.InTree, dtype='u4', var=['Entry$', self.get_hit_var()], cut='pulser')
        bins = histogram2d(x, y >= self.NMaxHits, bins=[self.NEntries // bin_size, [0, .5, 50]])[0]  # histogram the data to not over-count the empty events
        bin_average = bins[:, 1] / sum(bins, axis=1)
        misaligned = count_nonzero(bin_average > self.Threshold)
        self.Run.info(f'{100 * misaligned / bins.shape[0]:.1f}% of the events are misalignment :-(' if misaligned else f'Run {self.Run.Number} is perfectly aligned :-)')
        return misaligned == 0
    # endregion INIT
    # ----------------------------------------

    def fill_branches(self, ev):
        n = self.NHits[ev]
        hits = sum(self.NHits[:ev])
        self.Branches[0][1][0] = n  # n hits
        for i in range(1, 6):
            self.Branches[i][1][:n] = self.Variables[i - 1][hits:hits + n]
        self.Branches[-1][1][0] = self.Variables[-1][ev]  # trigger phase

    def get_pulser(self):
        return get_tree_vec(self.InTree, 'pulser', dtype='?')

    def get_pulser_rate(self):
        x = self.get_pulser()
        values = get_hist_vec(self.Run.Draw.profile(arange(x.size), x, make_bins(0, x.size, 100), show=False))
        return mean_sigma(values[values < mean(x) + .2])[0]

    def set_offset(self, pulser_event, offset):
        errors = self.Converter.DecodingErrors
        event = errors[(self.PulserEvents[pulser_event - 5] < errors) & (self.PulserEvents[pulser_event + 5] > errors)]
        self.Offsets[event[0] if event.size else self.PulserEvents[pulser_event]] = offset

    # ----------------------------------------
    # region OFFSETS
    def find_offset(self, off=-5, event=0, n=None):
        aligned = mean(self.NHits[roll(self.Pulser, off)][event:event + choose(n, self.BinSize)] > self.NMaxHits) < self.Threshold
        return off if aligned else self.find_offset(off + 1, event, n)

    def find_final_offset(self, off=-5):
        """return the offset at the end of the run"""
        return self.find_offset(off, event=self.PulserEvents.size - self.BinSize)

    def find_start_offset(self, off=-5):
        return self.find_offset(off)

    def find_offsets(self, off=0):
        start = self.find_offsets(off - 1) if off > self.FirstOffset else 0
        if start:
            self.set_offset(start, off)
        y = self.NHits[roll(self.Pulser, off)][start:] > self.NMaxHits
        v = mean(y[:y.size // self.BinSize * self.BinSize].reshape(y.size // self.BinSize, self.BinSize), axis=1)
        bad = where(v > self.Threshold)[0]
        if not bad.size:
            return info('found all offsets :-)')
        n = bad[0]
        y0 = y[n * self.BinSize:(n + 1) * self.BinSize]
        return where(y0 & (roll(y0, -1) == y0) & (roll(y0, -2) == y0))[0][0] + n * self.BinSize + start  # find first event with three misaligned in a row
    # endregion OFFSETS
    # ----------------------------------------


if __name__ == '__main__':

    from pad_run import PadRun
    from converter import Converter

    args = init_argparser(run=7, tc='201708')
    zrun = PadRun(args.run, testcampaign=args.testcampaign, load_tree=False, verbose=True)
    z = PadAlignment(Converter(zrun))
