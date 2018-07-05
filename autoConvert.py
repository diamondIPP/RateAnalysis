#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert new files from the beam test
# created on August 15th 2016 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from AbstractClasses.Run import Run
from os import stat, system
from multiprocessing import cpu_count
from argparse import ArgumentParser
from AbstractClasses.Utils import *


class AutoConvert:

    def __init__(self, run, multi):
        self.Run = run
        self.Multi = multi
        self.NextRun = self.load_last_converted()
        self.MaxRun = 1000

    def load_last_converted(self):
        filename = join(self.Run.Dir, 'last_converted.txt')
        f_exists = file_exists(filename)
        f = open(filename, 'r+' if f_exists else 'w')
        next_run = int(f.read()) + 1 if f_exists else 1
        f.close()
        return next_run

    def save_last_converted(self, run, reset=False):
        f = open(join(self.Run.Dir, 'last_converted.txt'), 'w')
        f.seek(0)
        f.write(str(run) if not reset else 0)
        f.close()

    def convert_run(self, run):
        run_infos = {int(key): value for key, value in self.Run.load_run_info_file().iteritems()}
        is_converted = file_exists(self.Run.converter.get_final_file_path(run))
        is_crap = run_infos[run]['runtype'] in ['test', 'crap', 'schrott'] if run in run_infos else True
        raw_file = self.Run.converter.find_raw_file(run)

        # check if we have to convert the run
        next_raw_file = self.Run.converter.find_raw_file(run + 1, prnt=False)
        if run not in run_infos and next_raw_file and file_exists(next_raw_file):
            self.NextRun += 1
            return True
        if is_converted or is_crap or not raw_file:
            if int(sorted(run_infos.keys())[-1]) in xrange(self.NextRun, self.MaxRun):
                self.NextRun += 1
                return True

        # convert if the file is written atm
        if raw_file and not file_is_beeing_written(raw_file):
            run_infos = {int(key): value for key, value in self.Run.load_run_info_file().iteritems()}
            if run not in run_infos:
                return False
            system('cvlc --play-and-exit ~/Downloads/closing_time.mp3 --run-time=60')
            Run(run, self.Run.TCString)
            self.save_last_converted(run)
            self.NextRun += 1
            return True
        return False

    def auto_convert(self):
        """Sequential conversion with check if the file is currently written. For usage during beam tests."""
        n_loops = 0
        print
        while n_loops < 14400:
            log_message('checking for new run ... {0} {1}'.format(self.NextRun, '(' + str(n_loops) + ')' if n_loops else ''), overlay=True)
            for run in xrange(self.NextRun, self.MaxRun):
                if self.convert_run(run):
                    n_loops = 0
                    continue
                else:
                    n_loops += 1
                    break
            sleep(5)

    def multi(self):
        """parallel conversion"""
        n_cpus = cpu_count()
        print 'We got {0} CPUs\n'.format(n_cpus)

        print 'Creating pool with {0} processes\n'.format(n_cpus)
        pool = Pool(n_cpus)
        print

        tasks = [(self, [run]) for run in xrange(self.NextRun, 1000)]

        results = [pool.apply_async(execute, t) for t in tasks]
        for res in results:
            print res.get(timeout=5000)

    def run(self):
        self.multi() if self.Multi else self.auto_convert()

    def __call__(self, run):
        return self.convert_run(run)


def execute(func, args):
    func(*args)


def file_is_beeing_written(path):
    size1 = stat(path)
    sleep(4)
    size2 = stat(path)
    return size1 != size2


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-m', action='store_true')
    parser.add_argument('-tc', nargs='?', default='201807')
    arg = parser.parse_args()

    dummy = Run(1, arg.tc, tree=False)
    z = AutoConvert(dummy, arg.m)
    print_banner('Starting {m} Conversion at run {r}'.format(m='Multi' if z.Multi else 'Auto', r=z.NextRun))
    z.run()
