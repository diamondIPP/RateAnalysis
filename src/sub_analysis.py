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
        self.Bins = analysis.Bins
        self.Cut = analysis.Cut
        self.NRocs = analysis.NRocs
        self.StartEvent = analysis.StartEvent
        if dut and hasattr(analysis, 'DUT'):
            self.DUT = analysis.DUT
        super().__init__(analysis.TCString, sub_dir=choose(sub_dir, analysis.Draw.SubDir), pickle_dir=pickle_dir)

    def get_root_vec(self, n=0, ind=0, dtype=None, var=None, cut='', nentries=None, firstentry=0):
        return self.Run.get_root_vec(n, ind, dtype, var, cut, nentries, firstentry)

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


class PadSubAnalysis(SubAnalysis):

    def __init__(self, analysis, sub_dir=None, pickle_dir='', dut=True):

        self.Channel = analysis.Channel
        self.Polarity = analysis.Polarity
        super().__init__(analysis, sub_dir, pickle_dir, dut)
