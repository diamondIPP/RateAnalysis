#!/usr/bin/env python
# --------------------------------------------------------
#       general class for event alignment
# created on October 18th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TFile
from helpers.utils import *
from numpy import sign, invert

MAX_SIZE = 255


class EventAligment(object):
    def __init__(self, converter=None):

        # Time
        self.StartTime = time()

        # Converter
        self.Converter = converter
        self.Run = converter.Run
        self.Draw = self.Run.Draw

        # ROOT Files and Trees
        self.InFile = self.load_file()
        self.InTree = self.InFile.Get(self.Run.TreeName)
        self.NewFile = None
        self.NewTree = None

        # Info
        self.NEntries = int(self.InTree.GetEntries())
        self.InTree.SetEstimate(self.NEntries)
        self.IsAligned = self.check_alignment()
        self.Offsets = {}
        self.FirstOffset = 0
        self.FinalOffset = 0

        # Branches
        self.Branches = self.init_branches()
        self.NHits = None
        self.Variables = None

        # Progress Bar
        self.PBar = PBar()

    def __call__(self):
        self.run()

    def run(self):
        if not self.IsAligned:
            self.print_start()
            self.NHits = self.load_n_hits()
            self.Variables = self.load_variables()
            if not self.check_alignment():  # make detailed check
                self.write_aligned_tree()
                if file_exists(self.Converter.ErrorFile):
                    remove(self.Converter.ErrorFile)  # file is not needed anymore and just gets appended otherwise
            eff = calc_eff(values=self.get_aligned(self.NewTree))[0]
            print_banner(f'{__class__.__name__} of run {self.Run.Number} finished in {get_elapsed_time(self.StartTime)} ({eff:.1f}% aligned)', color='green')

    def get_tree_vec(self, var, cut='', dtype=None, nentries=None, firstentry=0):
        return get_tree_vec(self.InTree, var, cut, dtype, nentries, firstentry)

    # ----------------------------------------
    # region INIT
    def print_start(self):
        pass

    def load_file(self):
        f_name = self.Converter.get_eudaqfile_path()
        return TFile(self.Converter.Run.RootFilePath if not file_exists(f_name) else f_name)

    def check_alignment(self, vec=None):
        self.Run.info(f'{calc_eff(values=invert(vec))[0]:.1f}% of the events are misaligned :-(' if not all(vec) else f'Run {self.Run.Number} is perfectly aligned :-)')
        return all(vec)

    def get_aligned(self, tree=None):
        pass

    @staticmethod
    def init_branches():
        return [('n_hits_tot', zeros(1, 'u1'), 'n_hits_tot/b'),
                ('plane', zeros(MAX_SIZE, 'u1'), 'plane[n_hits_tot]/b'),
                ('col', zeros(MAX_SIZE, 'u1'), 'col[n_hits_tot]/b'),
                ('row', zeros(MAX_SIZE, 'u1'), 'row[n_hits_tot]/b'),
                ('adc', zeros(MAX_SIZE, 'i2'), 'adc[n_hits_tot]/S'),
                ('charge', zeros(MAX_SIZE, 'float32'), 'charge[n_hits_tot]/F')]

    def get_hit_var(self):
        return 'n_hits_tot' if self.InTree.GetBranch('n_hits_tot') else '@col.size()'

    def load_n_hits(self, n_entries=None, first_entry=0):
        self.InTree.SetEstimate(self.NEntries)
        return get_tree_vec(self.InTree, self.get_hit_var(), '', 'u1', choose(n_entries, self.NEntries), first_entry)

    def load_variables(self):
        """ get all the telescope branches in vectors"""
        t = self.Run.info('Loading information from tree ... ', endl=False)
        self.InTree.SetEstimate(sum(self.NHits))
        data = get_tree_vec(self.InTree, ['plane', 'col', 'row', 'adc', 'charge'], dtype=['u1', 'u1', 'u1', 'i2', 'f2'])
        self.Run.add_to_info(t)
        return data

    def reload(self):
        self.NHits = self.load_n_hits()
        self.Variables = self.load_variables()
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region OFFSETS

    def find_first_offset(self):
        return 0

    def find_final_offset(self):
        return 0

    def load_error_offsets(self):
        self.Run.info('Using decoding errors to get the event alignment offsets')
        zero_offset = [(0, self.find_first_offset())] if self.find_first_offset() else []
        return OrderedDict(zero_offset + [(event_number, 1) for event_number in self.Converter.DecodingErrors])

    def find_manual_offsets(self):
        return {}

    def find_offsets(self, off, delta=1):
        use_decoding_errors = (self.find_final_offset() == len(self.Converter.DecodingErrors) + self.find_first_offset()) and self.Converter.DecodingErrors.size
        errors = self.load_error_offsets() if use_decoding_errors else self.find_manual_offsets()
        info('Found {n} offsets'.format(n=len(errors)))
        return errors

    # endregion OFFSETS
    # ----------------------------------------

    # ----------------------------------------
    # region WRITE TREE
    def save_tree(self):
        self.NewFile.cd()
        self.NewTree.Write()
        macro = self.InFile.Get('region_information')
        if macro:
            macro.Write()
        self.NewFile.Write()
        self.Run.info('successfully aligned the tree and saved it to {}'.format(self.NewFile.GetName()))

    def fill_branches(self, ev):
        pass

    def write_aligned_tree(self):
        set_root_output(False)
        self.find_offsets(self.FinalOffset + 5, max(sign(self.FinalOffset - self.FirstOffset), 1, key=abs))
        for name, o, leaf in self.Branches:  # remove old branches
            self.InTree.SetBranchStatus(name, 0)
        self.NewFile = TFile(self.Converter.get_eudaqfile_path(), 'RECREATE')
        self.NewTree = self.InTree.CloneTree(0)
        for name, o, leaf in self.Branches:  # set new branches
            self.NewTree.Branch(name, o, leaf)
        self.PBar.start(self.NEntries)
        offset = self.FirstOffset
        for ev, t in enumerate(self.InTree):
            self.PBar.update()
            if ev in self.Offsets:
                offset = self.Offsets[ev]
            if ev > self.NEntries - abs(offset) - 1:
                break
            self.fill_branches(ev + offset)
            self.NewTree.Fill()
        self.PBar.finish()
        self.save_tree()
    # endregion WRITE TREE
    # ----------------------------------------


if __name__ == '__main__':

    from pad_run import PadRun
    from converter import Converter

    # examples: (201508-442, ...)
    pargs = init_argparser(run=442)
    zrun = PadRun(pargs.run, testcampaign=pargs.testcampaign, load_tree=False, verbose=True)
    z = EventAligment(Converter(zrun))
