#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TH1F
from numpy import linspace
from src.correlation import Correlation
from src.event_alignment import *


class PixAlignment(EventAligment):
    def __init__(self, converter=None):

        # Info
        self.Threshold = .4
        self.NDutPlanes = 4
        self.TelPlane = 2
        self.DUTPlane = 4

        EventAligment.__init__(self, converter)

        self.TelRow = {}
        self.DUTRow = {}
        self.BinSize = None

    def print_start(self):
        print_banner('STARTING PIXEL EVENT ALIGNMENT OF RUN {}'.format(self.Converter.RunNumber))

    def update_variables(self):
        """Find all events with have both hits in the two checked planes"""
        t = self.Run.info('Loading pixel information from tree ... ', endl=False)
        for i, ev in enumerate(self.Variables):
            plane, row = ev[0], ev[2]
            if count_nonzero(plane == self.TelPlane) == 1:
                self.TelRow[i] = row[where(plane == self.TelPlane)][0]
            if count_nonzero(plane == self.DUTPlane) == 1:
                self.DUTRow[i] = row[where(plane == self.DUTPlane)][0]
        self.BinSize = self.find_bucket_size(show=False)
        self.Run.add_to_info(t)

    def check_alignment_fast(self):
        """ only check alignment for subsets of 10k events """
        t = self.Run.info('Fast check for event alignment ... ', endl=False)
        corrs = []
        for start_event in linspace(0, self.NEntries, 10, endpoint=False, dtype='i4'):
            correlation = Correlation(self, bucket_size=10000)
            n = self.InTree.Draw('plane:row', '', 'goff', 10000, start_event)
            planes, rows = get_tree_vec(self.InTree, n, 2, dtype='u1')
            n_hits = self.load_n_hits(10000, start_event)
            for i, (plane, row) in enumerate(split(array([planes, rows]), cumsum(n_hits), axis=1)):
                if count_nonzero(plane == self.TelPlane) == 1 and count_nonzero(plane == self.DUTPlane) == 1:  # if both planes have exactly one hit
                    correlation.fill(i, tel_row=row[where(plane == self.TelPlane)][0], dia_row=row[where(plane == self.DUTPlane)][0])
            corrs.append(correlation(debug=False))
        is_aligned = all(corr > self.Threshold for corr in corrs)
        self.Run.add_to_info(t)
        self.Run.info('Run {} is perfectly aligned :-)'.format(self.Run.RunNumber) if is_aligned else 'Fast check found misalignment :-(')
        return is_aligned

    def check_alignment(self):
        t = self.Run.info('Checking aligment ... ', endl=False)
        correlation = Correlation(self)
        for ev, row in self.TelRow.items():
            correlation.fill(ev)
        correlations = correlation.get_all_zero()
        h = TH1F('h_ee', 'Event Alignment', int(sqrt(len(correlations))), 0, 1)
        for cor in correlations:
            h.Fill(cor)
        set_root_output(0)
        fit = h.Fit('gaus', 'qs', '', .6, 1)
        # self.Run.format_histo(h, x_tit='Correlation Factor', y_tit='Number of Entries', y_off=1.4, stats=0)
        # self.Run.save_histo(h, 'EventAlignmentControl', show=False, lm=.13, prnt=False)
        mean_, sigma = fit.Parameter(1), fit.Parameter(2)
        if mean_ - 5 * sigma < .2:
            self.Run.add_to_info(t)
            info('run is very badly misaligned...')
            return False
        low_events = [cor for cor in correlations if cor < mean_ - 5 * sigma]
        misalignments = len(low_events) / float(len(correlations))
        self.Run.add_to_info(t)
        if misalignments > .02:
            info('found {v:5.2f} % misalignet events'.format(v=misalignments * 100))
            return False
        low_events = [cor for cor in correlations if cor < .3]
        misalignments = len(low_events) / float(len(correlations))
        if misalignments > .05:
            info('found {v:5.2f} % misalignet events'.format(v=misalignments * 100))
        else:
            self.Run.info('Everything is nicely aligned =)')
        return misalignments < .05

    def find_lose_corr_event(self, correlation, last_off_event, debug=False):
        # if all correlations are below threshold return the very first event
        if all(corr < self.Threshold for corr in correlation.get_all_zero()):
            return correlation.Events[0][0]
        correlations = correlation.get_sliding()
        if debug:
            self.draw_sliding(correlation)
        # find maximum and take mean around it
        start_index = next(correlations.keys().index(ev) + 20 for ev in correlations.keys() if ev > last_off_event) if last_off_event + 1000 > correlations.keys()[0] else 0
        max_index = correlations.values().index(max(correlations.values()[start_index:]))
        mean_ = self.get_mean(max_index, correlations.values())
        try:
            off_event = correlations.keys()[correlations.values().index(next(m for m in correlations.values()[max_index:] if m < mean_ - .1)) - 1]
        except StopIteration:
            off_event = correlations.keys()[max_index]
        if debug:
            print(mean_, off_event, next(m for m in correlations.values()[max_index:] if m < mean_ - .1))
        return off_event

    def find_offset(self, correlation, debug=False):
        inter_correlations = correlation.get_inter_sliding()
        if inter_correlations is None:
            return None
        min_anti_corr = min(inter_correlations.keys())
        if debug:
            print('anti correlation:', min_anti_corr)
        return inter_correlations[min_anti_corr] if min_anti_corr < -self.Threshold else None

    def find_start_offset(self):
        correlation = Correlation(self, bucket_size=200)
        correlation.fill_n(self.TelRow.keys()[:200])
        correlations = correlation.get_all()
        return max(correlations, key=correlations.get) if any(array(correlations.values()) > self.Threshold) else 0

    def find_final_offset(self):
        correlation = Correlation(self, bucket_size=400, n_offsets=200)
        correlation.fill_n(self.TelRow.keys()[-400:])
        correlations = correlation.get_all()
        return max(correlations, key=correlations.get)

    def draw_sliding(self, correlation):
        correlations = correlation.get_all_sliding()
        for off, corrs in correlations.items():
            self.Run.set_root_output(False)
            g = self.Run.make_tgrapherrors('g{o}'.format(o=off), 'Sliding correlation for offset {o}'.format(o=off), x=corrs.keys(), y=corrs.values())
            self.Run.histo(g, draw_opt='alp')

    def find_jump(self, correlation, debug=False):
        correlations = correlation.get_sliding()
        if debug:
            self.draw_sliding(correlation)
        # first find the minimum: it has to be in the central region since it has a maximum to each of his sides
        min_index = self.find_min_index(correlations.values(), start_index=len(correlations) / 3)
        # look left for a maximum
        l_max_index = self.get_max_index(correlations.values(), min_index)
        l_mean = self.get_mean(l_max_index, correlations.values())
        l_off_event = self.find_falling_event(correlations, l_max_index, l_mean)
        # now find the offset by checking the inter correlations and look when that correlation falls off again
        offset = self.find_offset(correlation, debug=debug)
        if offset is None:
            return None, None, None
        correlations = correlation.get_sliding(offset=offset)
        o_max_ind = self.get_max_index(correlations.values())
        o_mean = self.get_mean(o_max_ind, correlations.values())
        o_off_event = self.find_falling_event(correlations, o_max_ind, o_mean)
        if o_off_event is None:
            return None, None, None
        if debug:
            print(l_mean, l_off_event)
            print(o_off_event)
            print(offset)
        return l_off_event, o_off_event, offset

    @staticmethod
    def find_min_index(values, start_index=0):
        return values.index(min(values[start_index:]))

    @staticmethod
    def get_max_index(values, end_index=None):
        e = end_index if end_index is not None else len(values)
        return values.index(max(values[:e]))

    @staticmethod
    def get_mean(ind, values, window=5):
        s = ind - window if ind >= window else 0
        e = ind + window
        return mean(values[s:e])

    @staticmethod
    def find_falling_event(correlations, ind, mean_):
        try:
            return correlations.keys()[correlations.values().index(next(m for m in correlations.values()[ind:] if m < mean_ - .1)) - 1]
        # if the maximum is at the end we won't see the falloff
        except StopIteration:
            return

    def find_rising_event(self, correlations, ind):
        return correlations.keys()[correlations.values().index(next(c for c in correlations.values()[ind:] if c > self.Threshold)) - 3]

    def find_manual_offsets(self, debug=False):
        t = self.Run.info('Scanning for precise offsets ... ', endl=False)

        n = self.BinSize
        correlation = Correlation(self, n_offsets=2, bucket_size=n)

        offsets = OrderedDict()
        offset = 0
        for ev in self.TelRow.keys():
            # fill the correlation vectors
            correlation.fill(ev, offset)
            if correlation.start():
                # check zero offset correlation
                if correlation.get_zero(start_bucket=-2, debug=debug) < self.Threshold:
                    # there is a jump if the last bucket is already aligned again
                    if correlation.get_zero(start_bucket=-1, debug=debug) > self.Threshold:
                        # get the events when it starts losing correlation and when the offset correlation falls back into zero offset
                        l_off_event, o_off_event, off = self.find_jump(correlation, debug=debug)
                        if off is not None and o_off_event > l_off_event:
                            offsets[l_off_event] = off
                            offsets[o_off_event] = -off
                            if debug:
                                info('Found a jump of {v} between events {e1} and {e2}'.format(v=off, e1=l_off_event, e2=o_off_event))
                            # only keep the last two buckets since at least the last is correlated
                            correlation.reset_except_last(2)
                    # now we have lost correlation for at least three buckets
                    else:
                        last_off_event = list(offsets.keys())[-1] if offsets else 0
                        off_event = self.find_lose_corr_event(correlation, last_off_event, debug=debug)
                        off = self.find_offset(correlation)
                        if off is not None:
                            offset += off
                            offsets[off_event] = off
                            if debug:
                                info('Found an offset of {v} at event {e}'.format(v=off, e=off_event))
                            # delete the first two buckets and use the offset values for zero correlation
                            correlation.reshuffle(off)
                    if off is None and debug:
                        info('Could not find any offset!')
                # reset old lists for speed improvements
                correlation.reset()
        self.Run.add_to_info(t)
        return offsets

    def find_bucket_size(self, show=True):
        """ take first 10000 events and find a suitable bucket size to build the correlation """
        correlation = Correlation(self, bucket_size=10)
        max_ev = 10000
        for ev in self.TelRow.keys():
            if ev > max_ev:
                break
            correlation.fill(ev)
        sigmas = OrderedDict()
        size = 50
        while True:
            if correlation.get_events() < 700:
                break
            try:
                for i, n in enumerate(range(10, 100)):
                    correlation.set_bucket_size(n)
                    corrs = correlation.get_all_zero()
                    mean_, sigma = mean_sigma(corrs)
                    sigmas[sigma] = n
                # if show:
                #     g = self.Run.make_tgrapherrors('g_bs', 'Sigma of the Bucket Sizes', x=sigmas.values(), y=sigmas.keys())
                #     self.Run.draw_histo(g, draw_opt='alp')
                size = next(n for sig, n in sigmas.items() if sig < .09)
                break
            except StopIteration:
                correlation.delete_events(len(correlation.TelRow[0]) / 2, max_ev)
        return round_up_to(size, 5)

    def fill_branches(self, offset):
        if not offset:
            event = self.Variables[self.AtEntry]
        else:
            tel_offset = abs(offset) if offset < 0 else 0
            dut_offset = offset if offset > 0 else 0
            # select the according events
            tel_event, dut_event = self.Variables[self.AtEntry + tel_offset], self.Variables[self.AtEntry + dut_offset]
            # merge the right parts of the events
            event = concatenate((tel_event[:, tel_event[0] < self.NDutPlanes], dut_event[:, dut_event[0] >= self.NDutPlanes]), axis=1)
        for i, name in enumerate(self.Branches.keys()):
            for value in event[i]:
                self.Branches[name].push_back(value)


if __name__ == '__main__':

    from src.pix_run import PixelRun
    from src.converter import Converter

    args = init_argparser(run=147, tc='201810')
    zrun = PixelRun(args.run, testcampaign=args.testcampaign, load_tree=False, verbose=True)
    z = PixAlignment(Converter(zrun))
