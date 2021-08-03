#!/usr/bin/env python
# --------------------------------------------------------
#       general class for event alignment
# created on October 18th 2019 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TFile
from helpers.utils import *
from helpers.draw import set_root_output
from numpy import invert, ones, append

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
        self.HitVar = 'n_hits_tot' if bool(self.InTree.GetBranch('n_hits_tot')) else '@col.size()'
        self.NEntries = int(self.InTree.GetEntries())
        self.InTree.SetEstimate(self.NEntries)
        self.IsAligned = self.check_alignment()
        self.Aligned = ones(self.NEntries, dtype='?')
        self.Offsets = {}
        self.FirstOffset = 0
        self.FinalOffset = 0

        # Branches
        self.Branches = self.init_branches()
        self.NHits = None
        self.NTot = None
        self.Variables = None

        # Progress Bar
        self.PBar = PBar(counter=True, t='s')

    def __call__(self):
        self.run()

    def run(self):
        if not self.IsAligned:
            self.print_start()
            self.NHits = self.load_n_hits()
            self.NTot = self.load_n_tot()
            self.Variables = self.load_variables()
            self.write_aligned_tree()
            if file_exists(self.Converter.ErrorFile):
                remove(self.Converter.ErrorFile)  # file is not needed anymore and just gets appended otherwise
            eff = calc_eff(values=self.get_aligned(self.NewTree))[0]
            print_banner(f'{__class__.__name__} of run {self.Run.Number} finished in {get_elapsed_time(self.StartTime)} ({eff:.1f}% aligned)', color='green')

    def get_tree_vec(self, var, cut='', dtype=None, nentries=None, firstentry=0):
        return get_tree_vec(self.InTree, var, cut, dtype, nentries, firstentry)

    def show_event(self, ev):
        print(f'NHits: {self.NHits[ev]}')
        for i, name in enumerate(self.get_tel_branches()):
            print(f'{name}: {self.Variables[i][self.NTot[ev]:self.NTot[ev + 1]]}')
        print(f'Trigger Phase: {self.Variables[-1][[2 * ev, 2 * ev + 1]]}')

    # ----------------------------------------
    # region INIT
    def print_start(self):
        print_banner(f'STARTING {self.__class__.__name__[:3].title()} EVENT ALIGNMENT OF RUN {self.Run.Number}')

    def load_file(self):
        f_name = self.Converter.get_eudaqfile_path()
        return TFile(self.Converter.Run.RootFilePath if not file_exists(f_name) else f_name)

    def check_alignment(self, *args, **kwargs):
        v = self.get_aligned(*args, **kwargs)
        self.Run.info(f'{calc_eff(values=invert(v))[0]:.1f}% of the events are misaligned :-(' if not all(v) else f'Run {self.Run.Number} is perfectly aligned :-)')
        return all(v)

    def get_aligned(self, *args, **kwargs):
        return array([])

    @staticmethod
    def init_branches():
        return {'n_hits_tot': (zeros(1, 'u1'), 'n_hits_tot/b'),
                'plane': (zeros(MAX_SIZE, 'u1'), 'plane[n_hits_tot]/b'),
                'col': (zeros(MAX_SIZE, 'u1'), 'col[n_hits_tot]/b'),
                'row': (zeros(MAX_SIZE, 'u1'), 'row[n_hits_tot]/b'),
                'adc': (zeros(MAX_SIZE, 'i2'), 'adc[n_hits_tot]/S'),
                'charge': (zeros(MAX_SIZE, 'float32'), 'charge[n_hits_tot]/F'),
                'trigger_phase': (zeros(2, 'u1'), 'trigger_phase[2]/b'),
                'aligned': (zeros(1, '?'), 'aligned/O')}

    @staticmethod
    def get_tel_branches():
        return ['plane', 'col', 'row', 'adc', 'charge']

    def load_n_hits(self, n_entries=None, first_entry=0):
        self.InTree.SetEstimate(self.NEntries)
        return get_tree_vec(self.InTree, self.HitVar, '', 'u1', choose(n_entries, self.NEntries), first_entry)

    def load_n_tot(self):
        return append(0, cumsum(self.NHits, dtype='i8'))

    def load_variables(self):
        """ get all the telescope branches in vectors"""
        t = self.Run.info('Loading information from tree ... ', endl=False)
        self.InTree.SetEstimate(self.NTot[-1])
        data = self.get_tree_vec(['plane', 'col', 'row', 'adc', 'charge'], dtype=['u1', 'u1', 'u1', 'i2', 'f2'])
        data += [self.get_tree_vec('trigger_phase', dtype='u1')]
        self.Run.add_to_info(t)
        return data

    def reload(self):
        self.NHits = self.load_n_hits()
        self.NTot = self.load_n_tot()
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

    def find_all_offsets(self, *args, **kwargs):
        return

    def set_aligned(self, bin_size=None):
        return
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

    def fill_branches(self, ev, offset):
        pass

    def write_aligned_tree(self):
        set_root_output(False)
        self.find_all_offsets()
        self.set_aligned()
        for name in self.Branches:  # remove old branches
            self.InTree.SetBranchStatus(name, 0)
        self.NewFile = TFile(self.Converter.get_eudaqfile_path(), 'RECREATE')
        self.NewTree = self.InTree.CloneTree(0)
        for name, (br, leaf) in self.Branches.items():  # set new branches
            self.NewTree.Branch(name, br, leaf)
        info('STEP 2: Writing the TTree ...')
        self.PBar.start(self.NEntries, counter=False)
        offset = self.FirstOffset
        for ev, t in enumerate(self.InTree):
            self.PBar.update()
            if ev in self.Offsets:
                offset = self.Offsets[ev]
            if ev > self.NEntries - abs(offset) - 1:
                break
            self.fill_branches(ev, offset)
            self.NewTree.Fill()
        self.PBar.finish()
        self.save_tree()
    # endregion WRITE TREE
    # ----------------------------------------


if __name__ == '__main__':

    from src.run import Run
    from converter import Converter

    # examples: (201508-442, ...)
    pargs = init_argparser(run=442)
    zrun = Run(pargs.run, testcampaign=pargs.testcampaign, load_tree=False, verbose=True)
    z = EventAligment(Converter(zrun))
