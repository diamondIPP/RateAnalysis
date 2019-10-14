#!/usr/bin/env python
# --------------------------------------------------------
#       Class to align the DUT and REF events of the Rate Pixel Analysis
# created on February 13th 2017 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TFile, vector, TProfile
from Utils import *
from os import remove


class PadAlignment:
    def __init__(self, converter):

        # main
        self.StartTime = time()
        self.NDutPlanes = 4

        # converter
        self.Converter = converter
        self.Run = converter.Run

        print_banner('STARTING PAD EVENT ALIGNMENT OF RUN {r}'.format(r=self.Run.RunNumber))

        # files/trees
        self.InFile = TFile(converter.get_root_file_path())
        self.InTree = self.InFile.Get(self.Run.treename)
        self.NewFile = None
        self.NewTree = None

        # alignment
        self.NEntries = int(self.InTree.GetEntries())
        self.AtEntry = 0
        # branches
        self.Branches = self.init_branches()
        self.BranchLists = {name: [] for name in self.Branches}
        # info
        self.ColSize = []
        self.PulserEvents = []
        self.BucketSize = 30
        self.PulserRate = self.get_pulser_rate()
        self.Threshold = .4 if self.PulserRate > 5 else .2

    def __del__(self):
        self.InFile.Close()
        print_elapsed_time(self.StartTime, 'Pad Alignment')

    def run(self):
        if not self.is_aligned():
            self.load_variables()
            self.write_aligned_tree()
            if file_exists(self.Converter.ErrorFile):
                remove(self.Converter.ErrorFile)  # file is not needed anymore and just gets appended otherwise

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

    def is_aligned(self, binning=1000):
        """ just check the zero correlation """
        nbins = self.NEntries / binning
        h = TProfile('h', 'Pulser Rate', nbins, 0, self.NEntries)
        self.InTree.Draw('(@col.size()>1)*100:Entry$>>h', 'pulser', 'goff')
        aligned = all(h.GetBinContent(bin_) < 40 for bin_ in xrange(h.GetNbinsX()))
        self.Run.log_info('Fast check found misalignment :-(' if not aligned else 'Run {} is perfectly aligned :-)'.format(self.Run.RunNumber))
        return aligned

    def get_pulser_rate(self):
        h = TProfile('hprc', 'Pulser Rate', self.NEntries / 500, 0, self.NEntries)
        self.InTree.Draw('pulser:Entry$>>hprc', '', 'goff')
        values = [make_ufloat((h.GetBinContent(i), 1 / h.GetBinError(i))) for i in xrange(1, h.GetNbinsX() + 1) if h.GetBinContent(i) and h.GetBinError(i)]
        raw_mean = mean(values)
        values = [v for v in values if v < raw_mean + .1]
        m, s = mean_sigma(values)
        return m

    def find_offset(self, start, offset, n_events=None, n_offsets=5):
        offsets = sorted(range(-n_offsets, n_offsets + 1), key=lambda x: abs(x))
        offsets.pop(offsets.index(0))
        stop = self.BucketSize / 2 if n_events is None else n_events
        means = OrderedDict((i, self.calc_mean_size(start, i + offset, stop)) for i in offsets)
        # if we don't have beam in the beginning we cannot find out the offset
        # print ['{0:2.2f}'.format(i) for i in means.values()]
        if all(value < self.Threshold for value in means.itervalues()):
            return 0
        try:
            return next(key for key, value in means.iteritems() if value < self.Threshold)
        except StopIteration:
            return 0

    def find_zero_offset(self):
        """return the offset at the beginning of the run"""
        return self.find_offset(0, 0)

    def find_final_offset(self):
        off = 0
        while self.calc_mean_size(len(self.PulserEvents) - self.BucketSize - 1 - off, off) > self.Threshold:
            off += 1
        return off

    def get_error_offsets(self):
        self.Run.log_info('Using decoding errors to get the event alignment offsets')
        zero_offset = [(0, self.find_zero_offset())] if self.find_zero_offset() else []
        return OrderedDict(zero_offset + [(event_number, 1) for event_number in self.Converter.DecodingErrors])

    def calc_mean_size(self, start, off=0, n=None):
        n = n if n is not None else self.BucketSize
        return mean([self.ColSize[ev + off] > 3 for ev in self.PulserEvents[start:start + n]])

    def find_error_offset(self, i, offset):
        for ev in self.Converter.DecodingErrors:
            if self.PulserEvents[i - self.BucketSize / 2] < ev < self.PulserEvents[i + self.BucketSize / 2]:
                return ev, self.find_offset(self.PulserEvents.index(next(pev for pev in self.PulserEvents if pev > ev)), offset)
        return None, None

    def find_offsets(self):
        if self.find_final_offset() == len(self.Converter.DecodingErrors) + self.find_zero_offset() and self.Converter.DecodingErrors:
            return self.get_error_offsets()
        return self.find_shifting_offsets()

    def find_shifting_offsets(self, show=False):
        t = self.Run.log_info('Scanning for precise offsets ... ', next_line=False)
        n = self.BucketSize
        offset = self.find_offset(0, 0)
        # add first offset
        offsets = OrderedDict([(0, offset)] if offset else [])
        rates = [self.calc_mean_size(0)]
        i = 1
        g = self.Converter.Run.make_tgrapherrors('gse', 'Shifting Rate')
        g1 = self.Converter.Run.make_tgrapherrors('gse1', 'Shifting Rate')
        while i < len(self.PulserEvents) - abs(offset) - n:
            rate = self.calc_mean_size(i, offset)
            g.SetPoint(i - 1, i, rate)
            g1.SetPoint(i - 1, i, self.calc_mean_size(i, 1))
            # print i, '{0:1.2f}'.format(rate)
            if rate > self.Threshold:
                # first check if the event is in the decoding errors
                off_event, this_offset = self.find_error_offset(i, offset)
                if off_event is None:
                    # assume that the rate was good n/2 events before
                    good_rate = rates[-n / 2] if len(rates) > n / 2 else .1
                    for j, r in enumerate(rates[-n / 2:]):
                        if r > good_rate + .1:
                            # i + j - n/2 + n is the first misaligned event
                            off_event = self.PulserEvents[i + j - 1 + n / 2]
                            this_offset = self.find_offset(i + j - 1 + n / 2, offset)
                            if this_offset:
                                i += j - 1 + n / 2
                                break
                if this_offset:
                    # print 'Found offset:', off_event, this_offset
                    offsets[off_event] = this_offset
                    offset += this_offset
            rates.append(rate)
            i += 1
            if len(rates) > n:
                del rates[0]
        self.Run.format_histo(g, x_tit='Pulser Event', y_tit='Pulser Rate [%]', y_off=1.4)
        self.Run.format_histo(g1, color=self.Run.colors[0])
        self.Run.draw_histo(g, show=show)
        self.Run.draw_histo(g1, draw_opt='pl', canvas=get_last_canvas(), show=show)

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
        self.Run.log_info('successfully aligned the tree and saved it to {}'.format(self.NewFile.GetName()))

    def clear_vectors(self):
        for vec in self.Branches.itervalues():
            vec.clear()

    def write_aligned_tree(self):
        offsets = self.find_offsets()
        self.NewFile = TFile(self.Converter.get_root_file_path(), 'RECREATE')
        self.NewTree = self.InTree.CloneTree(0)
        self.set_branch_addresses()
        self.Run.start_pbar(self.NEntries)
        offset = 0
        while self.get_next_event():
            entry = self.AtEntry - 1
            self.Run.ProgressBar.update(self.AtEntry)
            self.clear_vectors()
            if entry in offsets:
                offset += offsets[entry]
            if entry > self.NEntries - offset - 1:
                break
            for name, lst in self.BranchLists.iteritems():
                for value in lst[entry + offset]:
                    self.Branches[name].push_back(value)
            self.NewTree.Fill()
        self.Run.ProgressBar.finish()
        self.save_tree()
    # endregion


if __name__ == '__main__':
    from Run import Run
    from Converter import Converter
    zrun = Run(259, test_campaign='201807', tree=False, verbose=True)
    zconv = Converter(zrun)
    z = PadAlignment(zconv)
    if not z.is_aligned():
        z.load_variables()
