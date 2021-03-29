#!/usr/bin/env python
# --------------------------------------------------------
#       general class for event alignment
# created on October 18th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TFile, vector, AddressOf
from helpers.utils import *


class EventAligment(object):
    def __init__(self, converter=None, verbose=True):

        # Time
        self.StartTime = time()

        # Converter
        self.Converter = converter
        self.Run = converter.Run

        # ROOT Files and Trees
        self.InFile = self.load_file()
        self.InTree = self.InFile.Get(self.Run.TreeName)
        self.NewFile = None
        self.NewTree = None

        # Info
        self.Verbose = verbose
        self.NEntries = int(self.InTree.GetEntries())
        self.MaxOffset = 5
        self.AtEntry = -1
        self.IsAligned = self.check_alignment_fast()

        # Branches
        self.Branches = self.init_branches()
        self.Variables = None
        self.NHits = None

        # Progress Bar
        self.PBar = PBar()

    def __del__(self):
        self.InFile.Close()
        print_elapsed_time(self.StartTime, 'Pad Alignment', show=self.Verbose)

    def __call__(self):
        self.run()

    def run(self):
        if not self.IsAligned:
            self.print_start()
            self.NHits = self.load_n_hits()
            self.Variables = self.load_variables()
            self.update_variables()
            if not self.check_alignment():  # make detailed check
                self.write_aligned_tree()
                if file_exists(self.Converter.ErrorFile):
                    remove(self.Converter.ErrorFile)  # file is not needed anymore and just gets appended otherwise

    # ----------------------------------------
    # region INIT
    def print_start(self):
        pass

    def load_file(self):
        f_name = self.Converter.get_eudaqfile_path()
        f_name = self.Converter.Run.RootFilePath if not file_exists(f_name) else f_name
        return TFile(f_name)

    def check_alignment(self):
        return self.IsAligned

    def check_alignment_fast(self):
        pass

    @staticmethod
    def init_branches():
        dic = OrderedDict()
        dic['plane'] = vector('unsigned short')()
        dic['col'] = vector('unsigned short')()
        dic['row'] = vector('unsigned short')()
        dic['adc'] = vector('short')()
        dic['charge'] = vector('unsigned int')()
        return dic

    def load_n_hits(self, n_entries=None, first_entry=0):
        self.InTree.SetEstimate(self.NEntries)
        return get_tree_vec(self.InTree, '@plane.size()', '', 'u1', choose(n_entries, self.NEntries), first_entry)

    def load_variables(self):
        """ get all the telescope branches in vectors"""
        t = self.Run.info('Loading information from tree ... ', endl=False)
        n_hits = self.InTree.Draw('plane', '', 'goff')
        self.InTree.SetEstimate(n_hits)
        plane, column, row, adc, charge = get_tree_vec(self.InTree, ['plane', 'col', 'row', 'adc', 'charge'], dtype=int)
        self.Run.add_to_info(t)
        return split(array([plane, column, row, adc, charge]), cumsum(self.NHits), axis=1)

    def reload(self):
        self.NHits = self.load_n_hits()
        self.Variables = self.load_variables()
        self.update_variables()

    def update_variables(self):
        """ get additional vectors"""
        pass
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region OFFSETS

    def find_start_offset(self):
        return 0

    def find_final_offset(self):
        return 0

    def load_error_offsets(self):
        self.Run.info('Using decoding errors to get the event alignment offsets')
        zero_offset = [(0, self.find_start_offset())] if self.find_start_offset() else []
        return OrderedDict(zero_offset + [(event_number, 1) for event_number in self.Converter.DecodingErrors])

    def find_manual_offsets(self):
        return {}

    def find_offsets(self):
        use_decoding_errors = (self.find_final_offset() == len(self.Converter.DecodingErrors) + self.find_start_offset()) and self.Converter.DecodingErrors.size
        errors = self.load_error_offsets() if use_decoding_errors else self.find_manual_offsets()
        info('Found {n} offsets'.format(n=len(errors)))
        return errors

    # endregion OFFSETS
    # ----------------------------------------

    # ----------------------------------------
    # region WRITE TREE
    def get_next_event(self):
        self.AtEntry += 1
        if self.AtEntry == self.NEntries:
            return False
        self.InTree.GetEntry(self.AtEntry)
        return True

    def set_branch_addresses(self):
        for name, branch in self.Branches.items():
            self.NewTree.SetBranchAddress(name, AddressOf(branch))

    def save_tree(self):
        self.NewFile.cd()
        self.NewTree.Write()
        macro = self.InFile.Get('region_information')
        if macro:
            macro.Write()
        self.NewFile.Write()
        self.Run.info('successfully aligned the tree and saved it to {}'.format(self.NewFile.GetName()))

    def clear_vectors(self):
        for vec in self.Branches.values():
            vec.clear()

    def show_branches(self):
        for i, j in self.Branches.items():
            print(i, list(j))

    def fill_branches(self, offset):
        pass

    def write_aligned_tree(self):
        offsets = self.find_offsets()
        self.NewFile = TFile(self.Converter.get_eudaqfile_path(), 'RECREATE')
        self.NewTree = self.InTree.CloneTree(0)
        self.set_branch_addresses()
        self.PBar.start(self.NEntries)
        offset = 0
        while self.get_next_event():
            self.PBar.update(self.AtEntry)
            self.clear_vectors()
            if self.AtEntry in offsets:
                offset += offsets[self.AtEntry]
            if self.AtEntry > self.NEntries - abs(offset) - 1:
                break
            self.fill_branches(offset)
            self.NewTree.Fill()
        self.PBar.finish()
        self.save_tree()
    # endregion WRITE TREE
    # ----------------------------------------


if __name__ == '__main__':

    from pad_run import PadRun
    from converter import Converter

    args = init_argparser(run=23, tc='201908')
    zrun = PadRun(args.run, testcampaign=args.testcampaign, tree=False, verbose=True)
    z = EventAligment(Converter(zrun))
