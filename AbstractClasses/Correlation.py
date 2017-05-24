#!/usr/bin/env python
# --------------------------------------------------------
#       Class to handle correlations of telescope planes
# created on March 28th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from collections import OrderedDict
from numpy import corrcoef
from ROOT import TProfile
from copy import deepcopy


class Correlation(object):

    def __init__(self, alignment, bucket_size=None, n_offsets=0):

        self.Alignment = alignment
        self.BucketSize = bucket_size if bucket_size is not None else self.Alignment.BucketSize
        self.NBuckets = None
        self.BucketWindow = 4  # number of bucket of the usual window
        self.MaxBuckets = 4  # we need to start at five buckets

        self.Offsets = [0] + [v for v in xrange(-n_offsets, n_offsets + 1) if v]
        self.Events = {off: [] for off in self.Offsets}
        self.TelRow = {off: [] for off in self.Offsets}
        self.DiaRow = {off: [] for off in self.Offsets}

    def set_bucket_size(self, value):
        self.BucketSize = int(value)

    def get_events(self, offset=0):
        return len(self.Events[offset])

    def fill(self, event, offset=0, tel_row=None, dia_row=None):
        # change the event number of the diamond row by the current offset
        dia_event = event + offset
        # fill the lists for n_off offsets around the current offset to cross check if these are aligned
        for off in self.Offsets:
            this_ev = dia_event + off
            if tel_row is not None and dia_row is not None:
                self.TelRow[off].append(tel_row)
                self.DiaRow[off].append(dia_row)
                self.Events[off].append(event)
            elif this_ev in self.Alignment.DiaRow:
                self.TelRow[off].append(self.Alignment.TelRow[event])
                self.DiaRow[off].append(self.Alignment.DiaRow[this_ev])
                self.Events[off].append(event)
        self.NBuckets = len(self.Events[0]) / self.BucketSize

    def del_first_bucket(self):
        for off in self.Offsets:
            while self.get_events(off) >= (self.NBuckets - 1) * self.BucketSize:
                del self.Events[off][0]
                del self.TelRow[off][0]
                del self.DiaRow[off][0]
        self.decrement_buckets()

    def reset_except_last(self, n):
        while len(self.Events[0]) >= self.BucketSize * n:
            self.del_first_bucket()

    def reshuffle(self, offset, n=2):
        tmp_ev = deepcopy(self.Events)
        tmp_tr = deepcopy(self.TelRow)
        tmp_dr = deepcopy(self.DiaRow)
        self.reset_except_last(2)
        # only keep the last n
        for off in self.Offsets:
            if off + offset in self.Offsets:
                self.Events[off] = tmp_ev[off + offset][(-n * self.BucketSize):]
                self.TelRow[off] = tmp_tr[off + offset][(-n * self.BucketSize):]
                self.DiaRow[off] = tmp_dr[off + offset][(-n * self.BucketSize):]

    def delete_events(self, start, stop):
        for off in self.Offsets:
            del self.Events[off][start:stop]
            del self.TelRow[off][start:stop]
            del self.DiaRow[off][start:stop]

    def reset(self):
        while self.get_events() >= self.BucketSize * self.BucketWindow:
            self.del_first_bucket()

    def start(self):
        # we can only start if we have at exactly N buckets filled
        return self.get_events() == self.MaxBuckets * self.BucketSize

    def increment_max_bucket(self):
        self.MaxBuckets += 1

    def decrement_max_bucket(self):
        self.MaxBuckets -= 1

    def decrement_buckets(self):
        self.NBuckets -= 1

    def get(self, offset=0, start_bucket=0, evt_offset=0, bucket_division=1., debug=False):
        s = int(start_bucket * self.BucketSize + evt_offset)
        e = int(s + self.BucketSize / bucket_division) if bucket_division else len(self.Events[offset])  # go one bucket further by default
        corr = correlate(self.TelRow[offset][s:e], self.DiaRow[offset][s:e])
        if debug:
            print start_bucket, len(self.Events[offset]), self.Events[offset][s], self.Events[offset][e if e < len(self.Events[offset]) else -1], corr
        return corr

    def get_inter_sliding(self):
        correlations = self.get_all_sliding()
        # need to sort the correlations into bins to be able to compare the different events
        profs = {}
        z_events = correlations[0].keys()
        nbins, s, e = len(z_events) / 2, z_events[0] / 50 * 50, z_events[-1] / 50 * 50 + 50
        for off, dic in correlations.iteritems():
            p = TProfile(str(off), 'p', nbins, s, e)
            for ev, corr in dic.iteritems():
                p.Fill(ev, corr)
            profs[off] = p
        # now get the correlation between the offset correlations
        correlations = {}
        z_list = [profs[0].GetBinContent(ibin) for ibin in xrange(profs[0].GetNbinsX())]
        for off, p in profs.iteritems():
            if off:
                off_list = [p.GetBinContent(ibin) for ibin in xrange(p.GetNbinsX())]
                # sort out statistical fluctuation below threshold
                if any(value > self.Alignment.Threshold for value in off_list):
                    # sort out the zeros
                    l1 = [v2 for v1, v2 in zip(off_list, z_list) if v1 and v2]
                    l2 = [v1 for v1, v2 in zip(off_list, z_list) if v1 and v2]
                    correlations[correlate(l1, l2)] = off
        return correlations if correlations else None

    def get_zero(self, start_bucket=-3, debug=False):
        div = 0 if start_bucket == -1 else 1
        return self.get(0, start_bucket, bucket_division=div, debug=debug)

    def get_shifted(self):
        """ get all the correlations for zero offset shifting through three buckets and the respective last events of the buckets"""
        n = self.BucketSize
        correlations = {self.Events[0][-n * 3 + k]: self.get(0, start_bucket=-4, evt_offset=k) for k in xrange(2 * n)}
        return OrderedDict(sorted(correlations.iteritems()))

    def get_sliding(self, offset=0):
        """ get all the correlations for zero offset sliding through all buckets and the respective last events of the buckets"""
        correlations = {}
        for k in xrange(self.NBuckets * self.BucketSize):  # try one bucket too much since not every list hast the same amount of entries
            try:
                ev = self.Events[offset][k + self.BucketSize - 1]
                correlations[ev] = self.get(offset=offset, evt_offset=k)
            except IndexError:
                pass
        return OrderedDict(sorted(correlations.iteritems()))

    def get_all_sliding(self):
        return {off: self.get_sliding(off) for off in self.Offsets}

    def get_shifted_long(self, last_offset):
        """ get all the correlations for zero offset shifting through three buckets """
        n = self.BucketSize
        correlations = {self.Events[0][-n * 2 + k]: self.get(0, start_bucket=-3, evt_offset=k) for k in xrange(n)}
        for k in xrange(n):
            correlations[self.Events[last_offset][-n * 3 + k]] = self.get(last_offset, start_bucket=-4, evt_offset=k)
        return OrderedDict(sorted(correlations.iteritems()))

    def get_off_all(self):
        # everything from -1.5 buckets till the end
        return {self.get(off, start_bucket=-3/2, bucket_division=0): off for off in self.Offsets if off}

    def get_all_sliced_offs(self, off_event):
        # sliced correlation from the off event until the end
        dic = {}
        old_size = self.BucketSize
        self.set_bucket_size(old_size / 2)
        n = self.BucketSize
        for offset in self.Offsets:
            start_event = next(ev for ev in self.Events[offset] if ev >= off_event)
            evt_offset = self.Events[offset].index(start_event)
            print offset, start_event, evt_offset
            correlations = {self.Events[offset][(b + 1) * n - 1]: self.get(offset=offset, start_bucket=b, evt_offset=evt_offset) for b in xrange((self.get_events() - evt_offset) / n)}
            dic[offset] = OrderedDict(sorted(correlations.iteritems()))
        self.set_bucket_size(old_size)
        return dic

    def get_all_zero(self):
        # all correlations for the zero offset
        return [self.get(0, start_bucket=j) for j in xrange(len(self.TelRow[0]) / self.BucketSize)]

    def get_detailed(self, division=5):
        n = self.BucketSize
        # start at -1.5 buckets and make sub buckets
        return {off: [self.get(off, start_bucket=-3/2, bucket_division=division, evt_offset=i) for i in xrange(0, n, n / division)] for off in self.Offsets}


def correlate(l1, l2):
    return corrcoef(l1, l2)[0][1]
