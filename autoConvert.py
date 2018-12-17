#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert new files from the beam test
# created on August 15th 2016 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------


from os import stat, sys
from os.path import join, dirname, realpath
file_dir = dirname(realpath(__file__))
sys.path.append(join(file_dir, 'AbstractClasses'))
from Converter import Converter
from Run import Run
from multiprocessing import cpu_count
from argparse import ArgumentParser
from Utils import *


class AutoConvert:

    def __init__(self, multi, end_run=None, test_campaign=None, verbose=False):

        self.Converter = Converter(Run(None, test_campaign=test_campaign, tree=False, verbose=verbose))
        self.RunInfos = self.load_run_infos()

        self.Dir = file_dir
        self.LastConvertedFile = join(self.Dir, '.last_converted.txt')

        self.Multi = multi
        self.FirstRun = self.load_last_converted()
        self.EndRun = max(self.RunInfos) if end_run is None else end_run
        self.Converter.set_run(self.FirstRun)

    def load_last_converted(self):
        if not file_exists(self.LastConvertedFile):
            return 1
        with open(self.LastConvertedFile) as f:
            last_run = int(f.read())
            return next(run for run in self.RunInfos if run >= last_run)

    def save_last_converted(self, run, reset=False):
        with open(self.LastConvertedFile, 'w') as f:
            f.seek(0)
            f.write(str(run) if not reset else 0)

    def load_run_infos(self):
        return {int(key): value for key, value in self.Converter.Run.load_run_info_file().iteritems()}

    def convert_run(self, run):

        self.Converter.set_run(run)

        # check if we have to convert the run
        if file_exists(self.Converter.get_final_file_path()) or self.RunInfos[run]['runtype'] in ['test', 'crap', 'schrott']:
            return
        raw_file = self.Converter.get_raw_file_path()
        if not file_exists(raw_file):
            log_warning('The following raw file does not exist: {}'.format(raw_file))
            return

        t = time()
        while file_is_beeing_written(raw_file):
            log_message('waiting until run {} is finished since {}'.format(run, get_running_time(t)), overlay=True)
            sleep(1)
        print
        Run(run, self.Converter.Run.TCString)
        self.save_last_converted(run)

    def auto_convert(self):
        """Sequential conversion with check if the file is currently written. For usage during beam tests."""
        for run in self.RunInfos:
            if run >= self.FirstRun:
                self.convert_run(run)
            # wait until a new run was added to the run log
            t = time()
            while run == max(self.RunInfos):
                log_message('waiting for new run ... {} since {}'.format(run + 1, get_running_time(t)), overlay=True)
                self.RunInfos = self.load_run_infos()
                sleep(1)

    def multi(self):
        """parallel conversion"""
        n_cpus = cpu_count()
        self.Converter.Run.log_info('We got {0} CPUs\n'.format(n_cpus))
        self.Converter.Run.log_info('Creating pool with {0} processes\n'.format(n_cpus))

        pool = Pool(n_cpus)

        tasks = [(self, [run]) for run in [run for run in self.RunInfos if self.FirstRun <= run <= self.EndRun]]
        results = [pool.apply_async(execute, t) for t in tasks]
        for res in results:
            print res.get(timeout=10000)

    def run(self):
        self.multi() if self.Multi else self.auto_convert()

    def __call__(self, run):
        return self.convert_run(run)


def execute(func, args):
    func(*args)


def file_is_beeing_written(file_path):
    size1 = stat(file_path)
    sleep(4)
    size2 = stat(file_path)
    return size1 != size2


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-m', action='store_true')
    parser.add_argument('-tc', nargs='?', default='201807')
    parser.add_argument('e', nargs='?', default=None)
    parser.add_argument('-v', action='store_false')
    argms = parser.parse_args()

    z = AutoConvert(argms.m, argms.e, argms.tc, argms.v)
    print_banner('Starting {m} Conversion at run {r}'.format(m='Multi' if z.Multi else 'Auto', r=z.FirstRun))
    z.run()
    print_banner('Finished Conversion!')
