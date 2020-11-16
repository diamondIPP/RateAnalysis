# --------------------------------------------------------
#       sub analysis collection class
# created on Nov 14th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.analysis import Analysis, choose, do_nothing, array, file_exists


class SubCollection(Analysis):
    """ small module to create all required fields for the subcollections. """

    def __init__(self, ana_collection, sub_dir=None, pickle_dir=''):

        self.Ana = ana_collection
        self.Analyses = ana_collection.Analyses
        self.Runs = ana_collection.Runs
        self.DUT = ana_collection.DUT
        self.RunPlan = self.Ana.RunPlan
        super().__init__(ana_collection.TCString, sub_dir=choose(sub_dir, ana_collection.Draw.SubDir), pickle_dir=pickle_dir)
        self.Config = ana_collection.Config

    def get_times(self):
        return self.Ana.get_times()

    def get_fluxes(self):
        return self.Ana.get_fluxes()

    def get_x(self, vs_time=False, rel_time=False, rel_error=0., avrg=False):
        return self.Ana.get_x_var(vs_time, rel_time, rel_error, avrg)

    def get_x_args(self, vs_time, rel_time=False, x_range=None):
        return self.Ana.get_x_args(vs_time, rel_time, x_range)

    def get_anas(self, runs=None):
        return self.Ana.get_analyes(runs)

    def get_values(self, string, f, runs=None, pbar=None, avrg=False, picklepath=None, *args, **kwargs):
        runs = choose(runs, self.Ana.Runs)
        pbar = choose(pbar, 'redo' in kwargs and kwargs['redo'] or (True if picklepath is None else not all(file_exists(picklepath.format(run)) for run in runs)))
        values = []
        self.info('Generating {} ...'.format(string), prnt=pbar)
        self.PBar.start(len(runs)) if pbar else do_nothing()
        for ana in self.get_anas(runs):
            values.append(f(ana, *args, **kwargs))
            self.PBar.update() if pbar else do_nothing()
        return array(self.Ana.get_flux_average(array(values)) if avrg else values)
