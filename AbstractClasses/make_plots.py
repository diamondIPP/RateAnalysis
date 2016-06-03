#!/usr/bin/python

from AnalysisCollection import AnalysisCollection
from argparse import ArgumentParser
from Elementary import Elementary
from RunSelection import RunSelection
from Utils import *

main_parser = ArgumentParser()
main_parser.add_argument('runplan', nargs='?', default=3)
main_parser.add_argument('dia', nargs='?', default=1, type=int)
main_parser.add_argument('-tc', '--testcampaign', nargs='?', default='')
main_parser.add_argument('-t', '--type', nargs='?', default='full')
args = main_parser.parse_args()
tc = args.testcampaign if args.testcampaign.startswith('201') else None
run_plan = args.runplan
diamond = args.dia
a = Elementary(tc)
sel = RunSelection(testcampaign=tc)
sel.select_runs_from_runplan(run_plan)
message = 'STARTING PAD-ANALYSIS COLLECTION OF RUNPLAN {0}'.format(run_plan)
print '\n{delim}\n{msg}\n{delim}\n'.format(delim=len(str(message)) * '=', msg=message)
a.print_testcampaign()
z = AnalysisCollection(sel, diamond)

z.save_dir = '{tc}_{dia}_{hv}_{rp}'.format(tc=z.TESTCAMPAIGN, dia=z.diamond_name, hv=z.bias, rp=z.run_plan)

if args.type in ['full', 'f']:
    z.draw_pulse_heights(show=False)
    z.draw_ph_with_currents(show=False)

elif args.type in ['single', 's']:
    z.draw_signal_distributions(show=False, off=200)

elif args.type in ['special', 'x']:
    z.draw_all_chi2s(show=False)
    z.draw_both_angles(show=False)

else:
    log_warning('You entered an invalid type!')






