#!/usr/bin/python

from AnalysisCollection import AnalysisCollection
from Elementary import Elementary
from RunSelection import RunSelection

tc = '201510'
a = Elementary(tc)
a.print_banner('STARTING RATE SCAN PLOT GENERATION')

only10 = True


def load_collection(plan, channel):
    sel = RunSelection(testcampaign=tc)
    sel.select_runs_from_runplan(plan)
    a.print_testcampaign()
    ana = AnalysisCollection(sel, channel)
    ana.save_dir = '{tc}_{dia}_{hv}_{rp}'.format(tc=ana.TESTCAMPAIGN, dia=ana.diamond_name, hv=ana.bias, rp=ana.run_plan)
    return ana


if only10:

    runplans10 = {3: [1], 5: [1], 8.1: [1, 2], 10.1: [1, 2]}
    upscans10 = {3: [1], 5: [1], 8.2: [1, 2], 10.2: [1, 2]}

    for rp, chs in runplans10.iteritems():

        for ch in chs:
            z = load_collection(rp, ch)

            z.draw_ph_with_currents(show=False)
            z.draw_pulse_heights(show=False)
            z.draw_pulser_info(show=False)
            z.close_files()

    for rp, chs in upscans10.iteritems():

        for ch in chs:
            z = load_collection(rp, ch)
            z.draw_signal_distributions(show=False, off=200)
            if rp == 8.2:
                z.draw_all_chi2s(show=False)
                z.draw_both_angles(show=False)
            z.close_files()

else:

    tc = '201508'
    a = Elementary(tc)

    runplans08 = {2: [2], 5.1: [1], 13: [2]}
    upscans08 = {2.2: [2], 5.2: [1], 13.1: [2]}

    for rp, chs in runplans08.iteritems():

        for ch in chs:
            z = load_collection(rp, ch)

            z.draw_ph_with_currents(show=False)
            z.draw_pulse_heights(show=False)
            z.draw_pulser_info(show=False)
            z.close_files()

    for rp, chs in upscans08.iteritems():

        for ch in chs:
            z = load_collection(rp, ch)
            z.draw_signal_distributions(show=False, off=200)
            z.close_files()
