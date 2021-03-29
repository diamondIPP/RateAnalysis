#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert new files from the beam test
# created on August 15th 2016 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------


from os import stat
from os.path import basename
from helpers.utils import *
from src.run import Run
from glob import glob


class AutoConvert:

    def __init__(self, multi, first_run=None, end_run=None, test_campaign=None, verbose=False):

        self.Run = Run(None, testcampaign=test_campaign, tree=False, verbose=verbose)
        self.RunInfos = OrderedDict((int(key), value) for key, value in self.Run.load_run_info_file().items())
        self.Runs = array(list(self.RunInfos.keys()), 'i2')

        self.Multi = multi
        self.FirstRun = self.find_last_converted(first_run)
        self.EndRun = max(self.RunInfos) if end_run is None else int(end_run)
        self.Run.Converter.set_run(self.FirstRun)

    def find_last_converted(self, run=None):
        last = max([int(remove_letters(basename(name))) for name in glob(join(self.Run.TCDir, 'root', '*', 'TrackedRun*.root'))])
        return self.Runs[self.Runs >= int(run)][0] if run is not None else last

    def convert_run(self, run):

        self.Run.Converter.set_run(run)

        # check if we have to convert the run
        if file_exists(self.Run.RootFilePath) or self.RunInfos[run]['runtype'] in ['test', 'crap', 'schrott']:
            print('{}: final file exists'.format(run))
            return
        raw_file = self.Run.Converter.RawFilePath
        if not file_exists(raw_file, warn=True):
            return

        t = time()
        while file_is_beeing_written(raw_file):
            info('waiting until run {} is finished since {}'.format(run, get_running_time(t)), endl=False)
            sleep(1)
        print()
        Run(run, self.Run.TCString)

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
        n_cpus = cpu_count()
        self.Run.info('We got {0} CPUs\n'.format(n_cpus))
        self.Run.info('Creating pool with {0} processes\n'.format(n_cpus))
        self.Run.info('End conversion at run {}'.format(self.EndRun))

        pool = Pool(n_cpus)

        runs = [run for run in self.Runs if self.FirstRun <= run <= self.EndRun]
        results = [pool.apply_async(self, [run]) for run in runs]
        for res in results:
            print(res.get(timeout=2 * 24 * 60 * 60))

    def run(self):
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
    parser.add_argument('s', nargs='?', default=None, help='run number where to start, default [None], = stop if no end is provided')
    parser.add_argument('e', nargs='?', default=None, help='run number where to stop, default [None]')
    parser.add_argument('-v', action='store_false', help='turn verbose OFF')
    args = parser.parse_args()

    z = AutoConvert(args.m, args.s, args.e, args.tc, args.v)
    print_banner('Starting {m} Conversion at run {r}'.format(m='Multi' if z.Multi else 'Auto', r=z.FirstRun))
    z.run()
    print_banner('Finished Conversion!')
