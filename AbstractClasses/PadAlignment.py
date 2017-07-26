#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TFile, vector, TProfile
from collections import OrderedDict
from numpy import mean
from Utils import log_message, time, print_elapsed_time
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar


class PadAlignment:
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
        self.IsAligned = self.check_alignment()
        if not self.IsAligned:
            # branches
            self.Branches = self.init_branches()
            self.BranchLists = {name: [] for name in self.Branches}
            # info
            self.ColSize = []
            self.PulserEvents = []
            self.load_variables()
            self.BucketSize = 30

    def __del__(self):
        self.InFile.Close()
        print_elapsed_time(self.StartTime, 'Pad Alignment')

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
        """ get all the telescope branches in vectors"""
        t = self.Run.log_info('Loading information from tree ... ', next_line=False)
        self.InTree.SetEstimate(self.InTree.Draw('plane', '', 'goff'))
        dic = {name: None for name in self.BranchLists}
        n = self.InTree.Draw('plane:row:col', '', 'goff')
        dic['plane'] = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
        dic['row'] = [int(self.InTree.GetV2()[i]) for i in xrange(n)]
        dic['col'] = [int(self.InTree.GetV3()[i]) for i in xrange(n)]
        n = self.InTree.Draw('adc:charge', '', 'goff')
        dic['adc'] = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
        dic['charge'] = [int(self.InTree.GetV2()[i]) for i in xrange(n)]
        n = self.InTree.Draw('@plane.size()', '', 'goff')
        self.ColSize = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
        n = self.InTree.Draw('Entry$', 'pulser', 'goff')
        self.PulserEvents = [int(self.InTree.GetV1()[i]) for i in xrange(n)]
        n_hits = 0
        for size in self.ColSize:
            for name, lst in dic.iteritems():
                self.BranchLists[name].append(lst[n_hits:size + n_hits])
            n_hits += size
        self.Run.add_info(t)

    def check_alignment(self, binning=1000):
        """ just check the zero correlation """
        nbins = self.NEntries / binning
        h = TProfile('h', 'Pulser Rate', nbins, 0, self.NEntries)
        self.InTree.Draw('(@col.size()>1)*100:Entry$>>h', 'pulser', 'goff')
        aligned = all(h.GetBinContent(bin_) < 40 for bin_ in xrange(h.GetNbinsX()))
        if not aligned:
            self.Run.log_info('Fast check found misalignment :-(')
        return aligned

    def find_offset(self, start, offset):
        means = {self.calc_mean_size(start, 1 + offset): 1, self.calc_mean_size(start, -1 + offset): -1}
        return next(means[key] for key in means.iterkeys() if key < self.Threshold)

    def calc_mean_size(self, start, off=0):
        return mean([self.ColSize[ev + off] > 1 for ev in self.PulserEvents[start:start + self.BucketSize]])

    def find_offsets(self):
        t = self.Run.log_info('Scanning for precise offsets ... ', next_line=False)
        n = self.BucketSize
        offsets = OrderedDict()
        offset = 0
        mean_sizes = []
        for i in xrange(0, len(self.PulserEvents), n):
            # check the number events during pulser events
            if i + n + offset > len(self.PulserEvents):
                break
            mean_sizes.append(self.calc_mean_size(i, offset))
            if mean_sizes[-1] > self.Threshold:
                # the aligned mean is two bunches before
                al_mean = mean_sizes[-3]
                # slide through the bucket before until this to find the off event
                for j in xrange(n):
                    sliding_mean = self.calc_mean_size((i - n + j), offset)
                    if sliding_mean > al_mean + .1:
                        off_event = self.PulserEvents[i + j - 1]
                        this_offset = self.find_offset(i + j - 1, offset)
                        offsets[off_event] = this_offset
                        # print off_event, offset
                        offset += this_offset
                        break

        self.Run.add_info(t)
        log_message('Found {n} offsets'.format(n=len(offsets)))
        return offsets

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
        macro = self.InFile.Get('region_information')
        if macro:
            macro.Write()
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
            if entry > self.NEntries - offset - 1:
                break
            for name, lst in self.BranchLists.iteritems():
                for value in lst[entry + offset]:
                    self.Branches[name].push_back(value)
            self.NewTree.Fill()
        self.ProgressBar.finish()
        self.save_tree()
    # endregion

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=self.Widgets, maxval=n)
        self.ProgressBar.start()
