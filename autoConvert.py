#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert new files from the beam test
# created on August 15th 2016 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from time import sleep
from AbstractClasses.Utils import file_exists, print_banner, log_message
from os import stat, chdir, system
from json import load
from multiprocessing import cpu_count, Pool
from argparse import ArgumentParser
from os.path import realpath, dirname


tc = '201707-2'
prog_dir = dirname(realpath(__file__))
data_dir = '/data/psi_{0}_{1}'.format(tc[:4], tc[-2:])
raw_dir = '{0}/raw'.format(data_dir)
final_dir = '{0}/root/pads'.format(data_dir)
lc = '{dir}/last_converted.txt'.format(dir=prog_dir)
f_lc = open(lc, 'r+') if file_exists(lc) else open(lc, 'w')
next_run = int(f_lc.read()) + 1 if file_exists(lc) else 1
f_lc.close()
print_banner('Starting Pad Autoconversion at run {0}!'.format(next_run))


def load_runinfos():
    f = open('{0}/run_log.json'.format(data_dir))
    infos = load(f)
    f.close()
    return infos


run_infos = load_runinfos()


def save_last_converted(run, reset=False):
    f = open(lc, 'w')
    f.seek(0)
    f.write(str(run) if not reset else 0)
    f.close()


def file_is_beeing_written(path):
    size1 = stat(path)
    sleep(4)
    size2 = stat(path)
    return size1 != size2


def make_raw_run_str(run, old=False):
    return '{1}/run{2}{0}.raw'.format(str(run).zfill(5 if old else 6), raw_dir, tc[2:] if old else '')


def make_final_run_str(run):
    return '{1}/TrackedRun{0}.root'.format(run, final_dir)


def convert_run(run):
    global run_infos, next_run
    run_infos = load_runinfos()
    run_infos = {int(key): value for key, value in run_infos.iteritems()}
    if run not in run_infos or run_infos[run]['runtype'] in ['test', 'crap', 'schrott']:
        if int(sorted(run_infos.keys())[-1]) in xrange(next_run, 1000):
            next_run += 1
            return False
        else:
            return 3
    raw = make_raw_run_str(run)
    old_raw = make_raw_run_str(run, True)
    final = make_final_run_str(run)
    if not file_exists(final):

        if file_exists(raw) or file_exists(old_raw):
            raw = raw if file_exists(raw) else old_raw
            if not file_is_beeing_written(raw):
                chdir(prog_dir)
                cmd = 'python AbstractClasses/Run.py {0} -t -tc {1}'.format(run, tc)
                print cmd
                system(cmd)
            else:
                return 3
    save_last_converted(run)
    next_run += 1
    return 2


def auto_convert():
    n_loops = 0
    print
    while n_loops < 14400:
        log_message('checking for new run ... {0} {1}'.format(next_run, '(' + str(n_loops) + ')' if n_loops else ''), overlay=True)
        for run in xrange(next_run, 1000):
            typ = convert_run(run)
            if typ == 2:
                n_loops = 0
                continue
            elif typ == 3:
                n_loops += 1
                break

        sleep(5)


def execute(func, args):
    func(*args)


def multi():
    n_cpus = cpu_count()
    print 'We got {0} CPUs\n'.format(n_cpus)

    print 'Creating pool with {0} processes\n'.format(n_cpus)
    pool = Pool(n_cpus)
    print

    tasks = [(convert_run, [run]) for run in xrange(next_run, 1000)]

    results = [pool.apply_async(execute, t) for t in tasks]
    for res in results:
        print res.get(timeout=5000)

parser = ArgumentParser()
parser.add_argument('-m', action='store_true')
arg = parser.parse_args()

if arg.m:
    multi()
else:
    auto_convert()
