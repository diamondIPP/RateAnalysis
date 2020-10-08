# --------------------------------------------------------
#       sub analysis class
# created on Oct 4th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.analysis import Analysis, choose


class SubAnanlysis(Analysis):
    """ small module to create all required fields for the subanalyses. """

    def __init__(self, analysis, results_dir=None, pickle_dir='', dut=True):

        self.Ana = analysis
        self.Run = analysis.Run
        self.Tree = analysis.Tree
        self.Bins = analysis.Bins
        self.Cut = analysis.Cut
        self.NRocs = analysis.NRocs
        if dut and hasattr(analysis, 'DUT'):
            self.DUT = analysis.DUT
        super().__init__(analysis.TCString, choose(results_dir, analysis.Draw.SubDir), pickle_dir)

    def get_root_vec(self, n=0, ind=0, dtype=None, var=None, cut='', nentries=None, firstentry=0):
        return self.Run.get_root_vec(n, ind, dtype, var, cut, nentries, firstentry)

    def has_branch(self, branch):
        return self.Run.has_branch(branch)

    def get_t_var(self):
        return self.Ana.get_t_var()

    def get_event_at_time(self, seconds, rel=False):
        return self.Run.get_event_at_time(seconds, rel)
