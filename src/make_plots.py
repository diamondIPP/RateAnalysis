#!/usr/bin/env python

from AnalysisCollection import AnalysisCollection
from Elementary import Elementary
from RunSelection import RunSelection
from glob import glob
from argparse import ArgumentParser
from Utils import *


def load_collection(plan, channel):
    sel = RunSelection()
    sel.select_runs_from_runplan(plan, ch=channel)
    threads = load_root_files(sel, plan)
    ana = AnalysisCollection(sel, threads=threads)
    ana.save_dir = '{info}_rp{rp}'.format(info=ana.make_info_string().strip('_'), rp=float(ana.RunPlan))
    ana.set_save_directory('PlotsFelix')
    return ana


def del_redundant_plots(res_dir, save_dir):
    for f in glob(join(res_dir, save_dir, '*', '*')):
        file_name = basename(f)
        for start in ['PulseHeightTime', 'PulseHeightZeroTime', 'PulserMean', 'PulseHeightFlu', 'PulseHeightZeroFlu']:
            if file_name.startswith(start):
                remove(f)


def create_plots():
    sel = RunSelection(verbose=False)
    sel.show_run_plans()
    for run_plan, info in sel.RunPlan.iteritems():
        for channel in xrange(1, sel.Run.get_n_diamonds(info['runs'][0])):
            print_banner('Starting AnalysisCollection for run plan {0} and ch {1}'.format(run_plan, channel), '-')
            z = load_collection(run_plan, channel)
            z.print_loaded()
            z.draw_ph_with_currents(show=False)
            z.draw_pulse_heights(show=False, redo=True)
            z.Pulser.draw_all_pulse_heights(show=False)
            z.draw_ph_distributions_below_flux(flux=80, show=False, save_plot=True)
            del_redundant_plots(z.ResultsDir, z.save_dir)
            z.close_files()
            z.__del__()


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('tc', nargs='?', default='201510')
    args = parser.parse_args()

    a = Elementary(args.tc, False)
    print_banner('STARTING RATE SCAN PLOT GENERATION FOR {}'.format(a.TCString))
    create_plots()
