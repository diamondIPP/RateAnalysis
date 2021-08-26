# --------------------------------------------------------
#       sub analysis class
# created on Oct 4th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.analysis import Analysis, choose


class SubAnalysis(Analysis):
    """ small module to create all required fields for the subanalyses. """

    def __init__(self, analysis, sub_dir=None, pickle_dir='', dut=True):

        self.Ana = analysis
        self.Run = analysis.Run
        self.Tree = analysis.Tree
        self.NRocs = analysis.NRocs
        if dut and hasattr(analysis, 'DUT'):
            self.DUT = analysis.DUT
        try:
            self.Cut = analysis.Cut
            self.Bins = analysis.Bins
            self.StartEvent = analysis.StartEvent
        except AttributeError:
            pass
        super().__init__(analysis.TCString, sub_dir=choose(sub_dir, analysis.Draw.SubDir), pickle_dir=pickle_dir)
        self.Config = analysis.Config
        self.Draw.ServerDir = analysis.Draw.ServerDir

    def get_tree_vec(self, var, cut='', dtype=None, nentries=None, firstentry=0):
        return self.Run.get_tree_vec(var, cut, dtype, nentries, firstentry)

    def has_branch(self, branch):
        return self.Run.has_branch(branch)

    def get_t_var(self):
        return self.Ana.get_t_var()

    def get_t_args(self, rel_time=True):
        return self.Ana.get_t_args(rel_time)

    def get_event_at_time(self, seconds, rel=False):
        return self.Run.get_event_at_time(seconds, rel)

    def get_attenuator(self):
        return self.Ana.get_attenuator()

    def get_irradiation(self):
        return self.Ana.get_irradiation()

    def make_event_cut(self, events):
        return self.Ana.make_event_cut(events)


class PadSubAnalysis(SubAnalysis):

    def __init__(self, analysis, sub_dir=None, pickle_dir='', dut=True):

        self.Channel = analysis.Channel
        self.Polarity = analysis.Polarity
        self.DigitiserBinWidth = analysis.DigitiserBinWidth
        super().__init__(analysis, sub_dir, pickle_dir, dut)
