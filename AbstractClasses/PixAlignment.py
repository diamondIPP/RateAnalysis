#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TFile, TH1F, vector
from collections import OrderedDict, Counter
from numpy import mean, sqrt
from Utils import set_root_output, log_message, time, print_elapsed_time, calc_mean, round_up_to
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from Correlation import Correlation


class PixAlignment:
    def __init__(self, converter):
        # main
        self.StartTime = time()
        self.Converter = converter
        self.Run = converter.Run
        self.NDutPlanes = 4
        self.Threshold = .4
        # progress bar
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
        self.ProgressBar = None
        # files/trees
        self.InFile = TFile(converter.get_root_file_path())
        self.InTree = self.InFile.Get(self.Run.treename)
        self.NewFile = None
        self.NewTree = None
        # alignment
        self.NEntries = int(self.InTree.GetEntries())
        self.AtEntry = 0
        self.IsAligned = self.check_alignment_fast()
        if not self.IsAligned:
            # branches
            self.Branches = self.init_branches()
            self.BranchLists = {name: [] for name in self.Branches}
            # info
            self.TelRow = {}
            self.DiaRow = {}
            self.load_variables()
            self.BucketSize = self.find_bucket_size(show=False)

    def __del__(self):
        self.InFile.Close()
        print_elapsed_time(self.StartTime, 'Pixel Alignment')

    @staticmethod
    def init_branches():
        dic = OrderedDict()
        dic['plane'] = vector('unsigned short')()
        dic['col'] = vector('unsigned short')()
        dic['row'] = vector('unsigned short')()
        dic['adc'] = vector('short')()
        dic['charge'] = vector('unsigned int')()
        return dic

    def load_variables(self):
        t = self.Run.log_info('Loading information from tree ... ', next_line=False)
        self.InTree.SetEstimate(self.InTree.Draw('plane', '', 'goff'))
        x, y = OrderedDict(), OrderedDict()
        p1, p2 = 2, 4
        dic = {name: None for name in self.BranchLists}
        n = self.InTree.Draw('plane:row:col:event_number', '', 'goff')
        dic['plane'] = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
        dic['row'] = [int(self.InTree.GetV2()[i]) for i in xrange(n)]
        dic['col'] = [int(self.InTree.GetV3()[i]) for i in xrange(n)]
        nrs = Counter([int(self.InTree.GetV4()[i]) for i in xrange(n)])
        n = self.InTree.Draw('adc:charge', '', 'goff')
        dic['adc'] = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
        dic['charge'] = [int(self.InTree.GetV2()[i]) for i in xrange(n)]
        n_ev = 0
        at_event = 0
        for ev, size in sorted(nrs.iteritems()):
            # account for events that have no pixel information
            while at_event != ev:
                for name in dic.iterkeys():
                    self.BranchLists[name].append([])
                at_event += 1
            for name, lst in dic.iteritems():
                self.BranchLists[name].append(lst[n_ev:size + n_ev])
            plane = dic['plane'][n_ev:size + n_ev]
            row = dic['row'][n_ev:size + n_ev]
            if plane.count(p1) == 1:
                x[ev] = row[plane.index(p1)]
            if plane.count(p2) == 1:
                y[ev] = row[plane.index(p2)]
            n_ev += size
            at_event += 1
        self.Run.add_info(t)
        self.TelRow = x
        self.DiaRow = y

    def check_alignment_fast(self):
        """ only check alignment of the first, last and middle 10k events """
        t = self.Run.log_info('Fast check for event alignment ... ', next_line=False)
        corrs = []
        for start_event in [self.NEntries / 10 * i for i in xrange(9)] + [self.NEntries - 10000]:
            correlation = Correlation(self, bucket_size=10000)
            n = self.InTree.Draw('plane:row:event_number', '', 'goff', 10000, start_event)
            planes = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
            rows = [int(self.InTree.GetV2()[i]) for i in xrange(n)]
            nrs = Counter([int(self.InTree.GetV3()[i]) for i in xrange(n)])
            n_ev = 0
            for ev, size in sorted(nrs.iteritems()):
                plane = planes[n_ev:size + n_ev]
                row = rows[n_ev:size + n_ev]
                if plane.count(2) == 1 and plane.count(4) == 1:
                    correlation.fill(ev, tel_row=row[plane.index(2)], dia_row=row[plane.index(4)])
                n_ev += size
            corrs.append(correlation.get(debug=False))
        is_aligned = all(corr > self.Threshold for corr in corrs)
        self.Run.add_info(t)
        if not is_aligned:
            self.Run.log_info('Fast check found misalignment :-(')
        return is_aligned

    def check_alignment(self):
        t = self.Run.log_info('Checking aligment ... ', next_line=False)
        correlation = Correlation(self)
        for ev, row in self.TelRow.iteritems():
            correlation.fill(ev)
        correlations = correlation.get_all_zero()
        h = TH1F('h_ee', 'Event Alignment', int(sqrt(len(correlations))), 0, 1)
        for cor in correlations:
            h.Fill(cor)
        set_root_output(0)
        fit = h.Fit('gaus', 'qs', '', .6, 1)
        self.Run.format_histo(h, x_tit='Correlation Factor', y_tit='Number of Entries', y_off=1.4, stats=0)
        self.Run.save_histo(h, 'EventAlignmentControl', show=False, lm=.13, prnt=False)
        mean_, sigma = fit.Parameter(1), fit.Parameter(2)
        if mean_ - 5 * sigma < .2:
            self.Run.add_info(t)
            log_message('run is very badly misaligned...')
            return False
        low_events = [cor for cor in correlations if cor < mean_ - 5 * sigma]
        misalignments = len(low_events) / float(len(correlations))
        self.Run.add_info(t)
        if misalignments > .02:
            log_message('found {v:5.2f} % misalignet events'.format(v=misalignments * 100))
            return False
        low_events = [cor for cor in correlations if cor < .3]
        misalignments = len(low_events) / float(len(correlations))
        if misalignments > .05:
            log_message('found {v:5.2f} % misalignet events'.format(v=misalignments * 100))
        else:
            self.Run.log_info('Everything is nicely aligned =)')
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
        off_event = correlations.keys()[correlations.values().index(next(m for m in correlations.values()[max_index:] if m < mean_ - .1)) - 1]
        if debug:
            print mean_, off_event, next(m for m in correlations.values()[max_index:] if m < mean_ - .1)
        return off_event

    def find_offset(self, correlation, debug=False):
        inter_correlations = correlation.get_inter_sliding()
        if inter_correlations is None:
            return None
        min_anti_corr = min(inter_correlations.keys())
        if debug:
            print 'anti correlation:', min_anti_corr
        return inter_correlations[min_anti_corr] if min_anti_corr < -self.Threshold else None

    def draw_sliding(self, correlation):
        correlations = correlation.get_all_sliding()
        for off, corrs in correlations.iteritems():
            self.Run.set_root_output(False)
            g = self.Run.make_tgrapherrors('g{o}'.format(o=off), 'Sliding correlation for offset {o}'.format(o=off), x=corrs.keys(), y=corrs.values())
            self.Run.draw_histo(g, draw_opt='alp')
        raw_input('Press any key to continue!')

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
            print l_mean, l_off_event
            print o_off_event
            print offset
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

    def find_offsets(self, debug=True):
        t = self.Run.log_info('Scanning for precise offsets ... ', next_line=False)

        n = self.BucketSize
        correlation = Correlation(self, n_offsets=2, bucket_size=n)

        offsets = OrderedDict()
        offset = 0
        for ev in self.TelRow.iterkeys():
            # fill the correlation vectors
            correlation.fill(ev, offset)
            if correlation.start():
                # check zero offset correlation
                if correlation.get_zero(start_bucket=-2, debug=debug) < self.Threshold:
                    # there is a jump if the last bucket is already aligned again
                    if correlation.get_zero(start_bucket=-1, debug=debug) > self.Threshold:
                        # get the events when it starts losing correlation and when the offset correlation falls back into zero offset
                        l_off_event, o_off_event, off = self.find_jump(correlation, debug=True)
                        if off is not None and o_off_event > l_off_event:
                            offsets[l_off_event] = off
                            offsets[o_off_event] = -off
                            if debug:
                                log_message('Found a jump of {v} between events {e1} and {e2}'.format(v=off, e1=l_off_event, e2=o_off_event))
                            # only keep the last two buckets since at least the last is correlated
                            correlation.reset_except_last(2)
                    # now we have lost correlation for at least three buckets
                    else:
                        last_off_event = offsets.keys()[-1] if offsets else 0
                        off_event = self.find_lose_corr_event(correlation, last_off_event, debug=True)
                        off = self.find_offset(correlation)
                        if off is not None:
                            offset += off
                            offsets[off_event] = off
                            if debug:
                                log_message('Found an offset of {v} at event {e}'.format(v=off, e=off_event))
                            # delete the first two buckets and use the offset values for zero correlation
                            correlation.reshuffle(off)
                    if off is None and debug:
                        log_message('Could not find any offset!')
                # reset old lists for speed improvements
                correlation.reset()
        self.Run.add_info(t)
        log_message('Found {n} offsets'.format(n=len(offsets)))
        return offsets

    def find_detailed_offset(self, corrs):
        # todo improve this... and maybe use it in general
        offs = {lst.index(corr): off for off, lst in corrs.iteritems() for corr in lst if corr > self.Threshold}
        print offs
        if not len(offs):
            for off, corr in corrs.iteritems():
                print off, corr
            return
        return offs[max(offs)]

    def find_bucket_size(self, show=True):
        """ take first 10000 events and find a suitable bucket size to build the correlation """
        correlation = Correlation(self, bucket_size=10)
        max_ev = 10000
        for ev in self.TelRow.iterkeys():
            if ev > max_ev:
                break
            correlation.fill(ev)
        sigmas = OrderedDict()
        size = 50
        while True:
            if correlation.get_events() < 700:
                break
            try:
                for i, n in enumerate(xrange(10, 100)):
                    correlation.set_bucket_size(n)
                    corrs = correlation.get_all_zero()
                    mean_, sigma = calc_mean(corrs)
                    sigmas[sigma] = n
                if show:
                    g = self.Run.make_tgrapherrors('g_bs', 'Sigma of the Bucket Sizes', x=sigmas.values(), y=sigmas.keys())
                    self.Run.draw_histo(g, draw_opt='alp')
                size = next(n for sig, n in sigmas.iteritems() if sig < .09)
                break
            except StopIteration:
                correlation.delete_events(len(correlation.TelRow[0]) / 2, max_ev)
        return round_up_to(size, 5)

    # =======================================================
    # region WRITE TREE
    def get_next_event(self):
        if self.AtEntry == self.NEntries:
            return False
        self.InTree.GetEntry(self.AtEntry)
        self.AtEntry += 1
        return True

    def set_branch_addresses(self):
        for name, branch in self.Branches.iteritems():
            self.NewTree.SetBranchAddress(name, branch)

    def save_tree(self):
        self.NewFile.cd()
        self.NewTree.Write()
        self.NewFile.Write()

    def clear_vectors(self):
        for vec in self.Branches.itervalues():
            vec.clear()

    def write_aligned_tree(self):
        offsets = self.find_offsets()
        self.NewFile = TFile(self.Converter.get_root_file_path(), 'RECREATE')
        self.NewTree = self.InTree.CloneTree(0)
        self.set_branch_addresses()
        self.start_pbar(self.NEntries)
        offset = 0
        while self.get_next_event():
            entry = self.AtEntry - 1
            self.ProgressBar.update(self.AtEntry)
            self.clear_vectors()
            if entry in offsets:
                offset += offsets[entry]
            if not offset:
                for name, lst in self.BranchLists.iteritems():
                    for value in lst[entry]:
                        self.Branches[name].push_back(value)
            else:
                ref = abs(offset) if offset < 0 else 0
                dut = offset if offset > 0 else 0
                # if we get out of range stop the loop
                try:
                    ref_plane = filter(lambda x: x < self.NDutPlanes, self.BranchLists['plane'][entry + ref])
                    dut_plane = filter(lambda x: x >= self.NDutPlanes, self.BranchLists['plane'][entry + dut])
                except IndexError:
                    break
                for i in xrange(len(ref_plane)):
                    for name in self.Branches.iterkeys():
                        self.Branches[name].push_back(self.BranchLists[name][entry + ref][i])
                for i in xrange(len(dut_plane)):
                    for name in self.Branches.iterkeys():
                        self.Branches[name].push_back(self.BranchLists[name][entry + dut][-len(dut_plane) + i])
            self.NewTree.Fill()
        self.ProgressBar.finish()
        self.save_tree()
    # endregion

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=self.Widgets, maxval=n)
        self.ProgressBar.start()
