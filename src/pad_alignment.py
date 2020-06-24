#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TProfile
from numpy import histogram2d, sum

from event_alignment import EventAligment
from utils import *


class PadAlignment(EventAligment):
    def __init__(self, converter, verbose=True):

        # Info
        self.Threshold = .4
        self.PulserEvents = array([])
        self.BinSize = 30

        EventAligment.__init__(self, converter, verbose)
        self.Offsets = sorted(range(-self.MaxOffset, self.MaxOffset + 1), key=abs)

        if not self.IsAligned:
            self.PulserRate = self.get_pulser_rate()
            self.Threshold = .4 if self.PulserRate > 5 else .2

    # ----------------------------------------
    # region INIT

    def print_start(self):
        print_banner('STARTING PAD EVENT ALIGNMENT OF RUN {}'.format(self.Run.Number))

    def update_variables(self):
        """ get additional vectors"""
        t = self.Run.info('Loading pad information from tree ... ', next_line=False)
        n = self.InTree.Draw('Entry$', 'pulser', 'goff')
        self.PulserEvents = get_root_vec(self.InTree, n, dtype=int)
        self.Run.add_to_info(t)

    def check_alignment_fast(self, bin_size=1000):
        """ just check the zero correlation """
        n = self.InTree.Draw('Entry$:@col.size()', 'pulser', 'goff')
        x, y = get_root_vecs(self.InTree, n, 2, dtype='u4')
        bins = histogram2d(x, y, bins=[self.NEntries / bin_size, [0, .5, 50]])[0]
        bin_average = bins[:, 1] / sum(bins, axis=1)
        misaligned = count_nonzero(bin_average > self.Threshold)
        self.Run.info('{:1.0f} % of the events are misalignment :-('.format(100. * misaligned / len(bins)) if misaligned else 'Run {} is perfectly aligned :-)'.format(self.Run.Number))
        return misaligned == 0
    # endregion INIT
    # ----------------------------------------

    def get_pulser_rate(self):
        h = TProfile('hprc', 'Pulser Rate', self.NEntries / 500, 0, self.NEntries)
        self.InTree.Draw('pulser:Entry$>>hprc', '', 'goff')
        values = [make_ufloat((h.GetBinContent(i), 1 / h.GetBinError(i))) for i in xrange(1, h.GetNbinsX() + 1) if h.GetBinContent(i) and h.GetBinError(i)]
        raw_mean = mean(values)
        values = [v for v in values if v < raw_mean + .1]
        m, s = mean_sigma(values)
        return m

    # ----------------------------------------
    # region OFFSETS
    def calc_hit_fraction(self, start, off=0, n_events=None):
        """" get the mean number of hits for a given offset in the pulser event interval [start: start + n] """
        n_events = self.BinSize if n_events is None else n_events
        return mean(self.NHits[self.PulserEvents[start:start + n_events] + off] > 3)

    def find_error_offset(self, pulser_event, offset):
        """check if the given event matches a decoding error"""
        for ev in self.Converter.DecodingErrors:
            if self.PulserEvents[pulser_event - self.BinSize / 2] < ev < self.PulserEvents[pulser_event + self.BinSize]:
                return ev, self.find_offset(self.find_next_pulser_event(ev), offset)
        return None, None

    def find_offset(self, start, initial_offset=0, n_events=None):
        n_events = self.BinSize if n_events is None else n_events
        for offset in self.Offsets:
            if self.calc_hit_fraction(start, offset + initial_offset, n_events) < self.Threshold:
                return offset
        return 0  # if nothing is found return no offset

    def find_start_offset(self):
        """return the offset at the beginning of the run"""
        return self.find_offset(0)

    def find_final_offset(self):
        """return the offset at the end of the run"""
        off = 0
        while self.calc_hit_fraction(len(self.PulserEvents) - self.BinSize - 1 - off, off) > self.Threshold:
            off += 1
        return off

    def find_next_pulser_event(self, event):
        return where(self.PulserEvents > event)[0][0]

    def find_shifting_offsets(self):
        t = self.Run.info('Scanning for precise offsets ... ', next_line=False)
        n = self.BinSize
        total_offset = self.find_start_offset()
        # add first offset
        offsets = OrderedDict([(0, total_offset)] if total_offset else [])
        rates = [self.calc_hit_fraction(0)]
        i = 1
        while i < len(self.PulserEvents) - abs(total_offset) - n:
            rate = self.calc_hit_fraction(i, total_offset)
            if rate > self.Threshold:
                # first check if the event is in the decoding errors
                off_event, this_offset = self.find_error_offset(i + n / 2, total_offset)
                if off_event is None:
                    # assume that the rate was good n/2 events before
                    good_rate = rates[-n / 2] if len(rates) > n / 2 else .1
                    for j, r in enumerate(rates[-n / 2:]):
                        if r > good_rate + .1:
                            # i + j - n/2 + n is the first misaligned event
                            off_event = self.PulserEvents[i + j - 1 + n / 2]
                            this_offset = self.find_offset(i + j - 1 + n / 2, total_offset)
                            if this_offset:
                                i += j - 1 + n / 2
                                break
                else:  # if error offset was found jump ahead to the next pulser event
                    i = self.find_next_pulser_event(off_event) - 1
                if this_offset:
                    offsets[off_event] = this_offset
                    total_offset += this_offset
            rates.append(rate)
            i += 1
            if len(rates) > n:
                del rates[0]
        self.Run.add_to_info(t)
        return offsets

    def find_manual_offsets(self):
        return self.find_shifting_offsets()
    # endregion OFFSETS
    # ----------------------------------------

    def fill_branches(self, offset):
        for i, name in enumerate(self.Branches.iterkeys()):
            for value in self.Variables[self.AtEntry + offset][i]:
                self.Branches[name].push_back(value)


if __name__ == '__main__':

    from pad_run import PadRun
    from converter import Converter

    args = init_argparser(run=23, tc='201908')
    zrun = PadRun(args.run, test_campaign=args.testcampaign, tree=False, verbose=True)
    z = PadAlignment(Converter(zrun))
