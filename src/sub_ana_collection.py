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
        self.Analysis = self.Analyses[0]
        self.Runs = ana_collection.Runs
        self.DUT = ana_collection.DUT
        self.Ensemble = ana_collection.Ensemble
        self.RunPlan = ana_collection.RunPlan
        super().__init__(ana_collection.TCString, sub_dir=choose(sub_dir, ana_collection.Draw.SubDir), pickle_dir=pickle_dir)
        self.Config = ana_collection.Config

    def get_times(self):
        return self.Ana.get_times()

    def get_time_bins(self, bin_width=None, only_edges=False):
        return self.Ana.get_time_bins(bin_width, only_edges)

    def get_fluxes(self, pbar=False):
        return self.Ana.get_fluxes(pbar=pbar)

    def get_x(self, vs_time=False, avrg=False):
        return self.Ana.get_x_var(vs_time, avrg)

    def get_x_args(self, vs_time=False, rel_time=False, draw=False, **kwargs):
        return self.Ana.get_x_args(vs_time, rel_time, draw, **kwargs)

    def get_x_draw(self, vs_time=False):
        return self.Ana.get_x_draw(vs_time)

    def get_analyses(self, runs=None):
        return self.Ana.get_analyes(runs)

    def get_runs(self):
        return self.Ana.get_runs()

    def get_values(self, string, f, runs=None, pbar=None, avrg=False, picklepath=None, flux_sort=False, plots=False, *args, **kwargs):
        runs = choose(runs, self.Ana.Runs)
        pbar = choose(pbar, 'redo' in kwargs and kwargs['redo'] or (True if picklepath is None else not all(file_exists(picklepath.format(run)) for run in runs)))
        values = []
        self.info('Generating {} ...'.format(string), prnt=pbar)
        self.PBar.start(len(runs)) if pbar else do_nothing()
        for ana in self.get_analyses(runs):
            values.append(f(ana, *args, **kwargs))
            self.PBar.update() if pbar else do_nothing()
        return values if plots else array(self.Ana.get_flux_average(array(values))) if avrg else array(values, dtype=object)[self.get_fluxes().argsort() if flux_sort else ...]

    def get_plots(self, string, f, runs=None, pbar=None, avrg=False, picklepath=None, *args, **kwargs):
        return self.get_values(string, f, runs, pbar, avrg, picklepath, False, True, *args, **kwargs)

    def get_pickle_path(self, name='', suf='', sub_dir=None, dut=None, camp=None):
        return self.Analysis.make_simple_pickle_path(name, suf, sub_dir, run='{}', dut=dut, camp=camp)
