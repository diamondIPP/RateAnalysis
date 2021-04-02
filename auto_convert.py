#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert new files from the beam test
# created on August 15th 2016 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------


from os import stat
from os.path import basename
from helpers.utils import *
from src.run import Run
from src.run_selection import RunSelector
from glob import glob


class AutoConvert:

    def __init__(self, multi, first_run=None, end_run=None, test_campaign=None, verbose=False):

        self.Multi = multi

        self.Selection = RunSelector(testcampaign=test_campaign, verbose=verbose)
        self.StartAtRun = choose(first_run, self.find_last_converted())
        self.StopAtRun = 1e9 if not multi or end_run is None else int(end_run)
        self.Runs = self.load_runs()

    def find_last_converted(self):
        return max([int(remove_letters(basename(name))) for name in glob(join(self.Selection.Run.TCDir, 'root', '*', 'TrackedRun*.root'))])

    def load_runs(self):
        runs = array([run for run in self.Selection.get_runplan_runs() if not file_exists(self.Selection.get_final_file_path(run))], 'i2')
        return runs[(runs >= self.StartAtRun) & (runs <= self.StopAtRun)]

    def convert_run(self, run):

        self.Run.Converter.set_run(run, )

        # check if we have to convert the run
        if file_exists(self.Run.RootFilePath) or self.RunInfos[run]['runtype'] in ['test', 'crap', 'schrott']:
            return
        raw_file = self.Run.Converter.RawFilePath
        if not file_exists(raw_file, warn=True):
            return

        t = time()
        while file_is_beeing_written(raw_file):
            info('waiting until run {} is finished since {}'.format(run, get_running_time(t)), endl=False)
            sleep(1)
        r = Run(run, self.Run.TCString)
        return f'{run} --> {timedelta(seconds=round(time() - r.InitTime))}'

    def auto_convert(self):
        """Sequential conversion with check if the file is currently written. For usage during beam tests."""
        for run in self.Runs:
            if run >= self.FirstRun:
                self.convert_run(run)
            # wait until a new run was added to the run log
            t = time()
            while run == max(self.Runs):
                info('waiting for new run ... {} since {}'.format(run + 1, get_running_time(t)), endl=False)
                self.RunInfos = self.Run.load_run_info_file()
                sleep(1)

    def multi(self):
        """parallel conversion"""
        self.Selection.Run.info(f'Creating pool with {cpu_count()} processes')
        with Pool() as pool:
            runs = pool.starmap(Run, [(run, self.Selection.TCString, True, False) for run in self.Runs])
        print(runs)
        for run in runs:
            print(f'{run} --> {timedelta(seconds=round(run.TInit))}')

    def run(self):
        if not self.Runs.size:
            return info('There are no runs to convert :-)')
        self.multi() if self.Multi else self.auto_convert()

    def __call__(self, run):
        return self.convert_run(run)


def file_is_beeing_written(file_path):
    size1 = stat(file_path)
    sleep(4)
    size2 = stat(file_path)
    return size1 != size2


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-m', action='store_true', help='turn parallel processing ON')
    parser.add_argument('-tc', nargs='?', default=None)
    parser.add_argument('s', nargs='?', default=None, help='run number where to start, default [None], = stop if no end is provided', type=int)
    parser.add_argument('e', nargs='?', default=None, help='run number where to stop, default [None]')
    parser.add_argument('-v', action='store_false', help='turn verbose OFF')
    args = parser.parse_args()

    z = AutoConvert(args.m, args.s, args.e, args.tc, args.v)
    if z.Runs.size:
        print_banner(f'Starting {"multi" if z.Multi else "auto"} Conversion for runs {z.Runs[0]} - {z.Runs[-1]}', color='green')
    z.run()
    print_banner('Finished Conversion!', color='green')
