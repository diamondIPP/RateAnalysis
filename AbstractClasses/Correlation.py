#!/usr/bin/env python
# --------------------------------------------------------
#       Class to handle correlations of telescope planes
# created on March 28th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from collections import OrderedDict
from numpy import corrcoef


class Correlation(object):

    def __init__(self, alignment, bucket_size=None, n_offsets=0):

        self.Alignment = alignment
        self.BucketSize = bucket_size if bucket_size is not None else self.Alignment.BucketSize
        self.AtBucket = 4  # we need at least four buckets to start

        self.Offsets = [0] + [v for v in xrange(-n_offsets, n_offsets + 1) if v]
        self.Events = {off: [] for off in self.Offsets}
        self.TelRow = {off: [] for off in self.Offsets}
        self.DiaRow = {off: [] for off in self.Offsets}

    def set_bucket_size(self, value):
        self.BucketSize = int(value)

    def get_events(self):
        return len(self.Events[0])

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

    def del_first_bucket(self):
        for off in self.Offsets:
            del self.Events[off][0:self.BucketSize]
            del self.TelRow[off][0:self.BucketSize]
            del self.DiaRow[off][0:self.BucketSize]

    def delete_events(self, start, stop):
        for off in self.Offsets:
            del self.Events[off][start:stop]
            del self.TelRow[off][start:stop]
            del self.DiaRow[off][start:stop]

    def reset(self, do=True):
        if not do:
            pass
        while len(self.Events[0]) >= self.BucketSize * 4:
            self.del_first_bucket()
            self.decrement_bucket()

    def start(self):
        # we can only start if we have at least 4 buckets filled and the size is exactly
        if len(self.Events[0]) == self.AtBucket * self.BucketSize:
            self.increment_bucket()
            return True
        return False

    def increment_bucket(self):
        self.AtBucket += 1

    def decrement_bucket(self):
        self.AtBucket -= 1

    def get(self, offset=0, start_bucket=0, evt_offset=0, bucket_division=1., debug=False):
        s = int(start_bucket * self.BucketSize + evt_offset)
        e = int(s + self.BucketSize / bucket_division) if bucket_division else len(self.Events[offset])  # go one bucket further by default
        corr = corrcoef(self.TelRow[offset][s:e], self.DiaRow[offset][s:e])[0][1]
        if debug:
            print start_bucket, len(self.Events[offset]), self.Events[offset][s], self.Events[offset][e if e < len(self.Events[offset]) else -1], corr
        return corr

    def get_zero(self, start_bucket=-2, debug=False):
        div = 0 if start_bucket == -1 else 1
        return self.get(0, start_bucket, bucket_division=div, debug=debug)

    def get_shifted(self):
        """ get all the correlations for zero offset shifting through three buckets """
        n = self.BucketSize
        correlations = {self.Events[0][-n * 3 + k]: self.get(0, start_bucket=-4, evt_offset=k) for k in xrange(2 * n)}
        return OrderedDict(sorted(correlations.iteritems()))

    def get_off_all(self):
        # everything from -1.5 buckets till the end
        return {self.get(off, start_bucket=-3/2, bucket_division=0): off for off in self.Offsets if off}

    def get_all_zero(self):
        # all correlations for the zero offset
        return [self.get(0, start_bucket=j) for j in xrange(len(self.TelRow[0]) / self.BucketSize)]

    def get_detailed(self, division=5):
        n = self.BucketSize
        # start at -1.5 buckets and make sub buckets
        return {off: [self.get(off, start_bucket=-3/2, bucket_division=division, evt_offset=i) for i in xrange(0, n, n / division)] for off in self.Offsets}
