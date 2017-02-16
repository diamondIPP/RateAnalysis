#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TFile, TH1F, vector
from collections import OrderedDict, Counter
from numpy import corrcoef, mean, sqrt
from Utils import set_root_output, log_message, log_critical, time, print_elapsed_time
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar


class PixAlignment:
    def __init__(self, converter):
        # main
        self.StartTime = time()
        self.Converter = converter
        self.Run = converter.Run
        self.NDutPlanes = 4
        # progress bar
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
        self.ProgressBar = None
        # files/trees
        self.InFile = TFile(converter.get_root_file_path())
        self.InTree = self.InFile.Get(self.Run.treename)
        self.NewFile = None
        self.NewTree = None
        # info
        self.Row1 = None
        self.Row2 = None
        self.load_rows()
        # alignment
        self.NEntries = int(self.InTree.GetEntries())
        self.AtEntry = 0
        # branches
        self.Branches = self.init_branches()

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
        dic['trigger_phase'] = vector('unsigned short')()
        return dic

    def load_rows(self):
        t = self.Run.log_info('Loading information from tree ... ', next_line=False)
        self.InTree.SetEstimate(self.InTree.Draw('plane', '', 'goff'))
        x, y = OrderedDict(), OrderedDict()
        p1, p2 = 2, 4
        n = self.InTree.Draw('plane:row:event_number', '', 'goff')
        planes = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
        rows = [int(self.InTree.GetV2()[i]) for i in xrange(n)]
        nrs = Counter([int(self.InTree.GetV3()[i]) for i in xrange(n)])
        n_ev = 0
        for ev, size in sorted(nrs.iteritems()):
            plane = planes[n_ev:size + n_ev]
            row = rows[n_ev:size + n_ev]
            if plane.count(p1) == 1:
                x[ev] = row[plane.index(p1)]
            if plane.count(p2) == 1:
                y[ev] = row[plane.index(p2)]
            n_ev += size
        self.Run.add_info(t)
        self.Row1 = x
        self.Row2 = y

    def check_alignment(self):
        t = self.Run.log_info('Checking aligment ... ', next_line=False)
        xt, yt = [], []
        n = 20
        for ev, row in self.Row1.iteritems():
            if ev in self.Row2:
                xt.append(row)
                yt.append(self.Row2[ev])
        correlations = [corrcoef(xt[j:(j + n)], yt[j:(j + n)])[0][1] for j in xrange(0, len(xt), n)]
        h = TH1F('h_ee', 'Event Alignment', int(sqrt(len(correlations))), 0, 1)
        for cor in correlations:
            h.Fill(cor)
        set_root_output(0)
        fit = h.Fit('gaus', 'qs', '', .6, 1)
        self.Run.format_histo(h, x_tit='Correlation Factor', y_tit='Number of Entries', y_off=1.4, stats=0)
        self.Run.save_histo(h, 'EventAlignmentControl', show=False, lm=.13, prnt=False)
        mean_, sigma = fit.Parameter(1), fit.Parameter(2)
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

    def find_offsets(self):
        t = self.Run.log_info('Scanning for precise offsets ... ', next_line=False)
        offsets = OrderedDict()
        n = 50
        thresh = .5
        # create n lists for events with -1, 0 and +1 offset
        n_offs = 2
        offs = [0] + [v for v in xrange(-n_offs, n_offs + 1) if v]
        xts = {off: OrderedDict() for off in offs}
        yts = {off: [] for off in offs}
        offset = 0
        checked = 3 * n
        for ev, row in self.Row1.iteritems():
            ev += offset
            for off in offs:
                this_ev = ev + off
                if this_ev in self.Row2:
                    xts[off][ev] = row
                    yts[off].append(self.Row2[this_ev])
            # check correlation of zero offset
            if len(yts[0]) % n == 0 and len(yts[0]) > checked:
                corr = corrcoef(xts[0].values()[(-2 * n):-n], yts[0][(-2 * n):-n])[0][1]
                # print len(yts[0]), xts[0].keys()[-2 * n], xts[0].keys()[-n], corr
                checked = len(yts[0])
                found_offset = False
                if corr < thresh:
                    # find exact event -> shift through three n buckets and take the last event if the correlation starts dropping
                    correlations = {xts[0].keys()[-n * 3 + k]: corrcoef(xts[0].values()[(-n * 4 + k):(-n * 3 + k)], yts[0][(-n * 4 + k):(-n * 3 + k)])[0][1] for k in xrange(2 * n)}
                    correlations = OrderedDict(sorted(correlations.iteritems()))
                    mean_ = mean(correlations.values()[:n])
                    mean_ = mean(correlations.values()[n:(3 * n / 2)]) if mean_ < thresh else mean_
                    off_event = correlations.keys()[correlations.values().index(filter(lambda x: x < mean_ - .1, correlations.itervalues())[0]) - 1]
                    # print off_event
                    # find offset
                    corrs = {corrcoef(xts[off].values()[(-3 * n / 2):], yts[off][(-3 * n / 2):])[0][1]: off for off in offs if off}
                    off, corr = corrs[max(corrs)], max(corrs)
                    if corr > thresh - .1:
                        offset += off
                        offsets[off_event] = off
                        # log_message('Found an offset of {v} at event {e}'.format(v=offsets[off_event], e=off_event))
                        found_offset = True
                        checked += n
                    if not found_offset:
                        off = self.find_detailed_offset(xts, yts, n, thresh)
                        if off is None:
                            log_critical('Something went wrong during finding the alignment offsets...')
                        else:
                            offset += off
                            checked += n
                            found_offset = True
                # reset old lists for speed improvements
                if len(yts[0]) >= 20 * n and not found_offset:
                    xts = {off: OrderedDict(sorted({ev: val for ev, val in xts[off].items()[(-3 * n):]}.iteritems())) for off in offs}
                    yts = {off: yts[off][(-3 * n):] for off in offs}
                    checked = 1
        self.Run.add_info(t)
        log_message('Found {n} offsets'.format(n=len(offsets)))
        return offsets

    @staticmethod
    def find_detailed_offset(x, y, n, thresh):
        # todo improve this... and maybe use it in general
        x_off = {off: lst.values()[(-3 * n / 2):] for off, lst in x.iteritems()}
        y_off = {off: lst[(-3 * n / 2):] for off, lst in y.iteritems()}
        corrs = {off: [corrcoef(x_off[off][i:(i + n / 5)], y_off[off][i:(i + n / 5)])[0][1] for i in xrange(0, n, n / 5)] for off in x_off.iterkeys()}
        offs = {lst.index(corr): off for off, lst in corrs.iteritems() for corr in lst if corr > thresh}
        if not len(offs):
            for off, corr in corrs.iteritems():
                print off, corr
            return
        return offs[max(offs)]

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
        print_banner('ALIGNING SHIT')
        # self.NewFile = TFile(self.Converter.get_root_file_path(), 'RECREATE')
        self.NewFile = TFile('test.root', 'RECREATE')
        self.NewTree = self.InTree.CloneTree(0)
        # self.NewTree = TTree(self.InTree.GetName, self.InTree.GetTitle())
        self.set_branch_addresses()
        self.start_pbar(self.NEntries)
        while self.get_next_event():
            self.ProgressBar.update(self.AtEntry)
            self.clear_vectors()
            for i in xrange(len(self.InTree.plane)):
                self.Branches['plane'].push_back(i)
                self.Branches['col'].push_back(i)
                self.Branches['row'].push_back(i)
                self.Branches['adc'].push_back(i)
                self.Branches['charge'].push_back(i)
            self.NewTree.Fill()
        self.save_tree()
    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=self.Widgets, maxval=n)
        self.ProgressBar.start()
