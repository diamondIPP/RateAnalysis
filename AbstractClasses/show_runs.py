#!/usr/bin/env python
from RunSelection import RunSelection
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('RunPlan', nargs='?', help='Displays the deltailed overview about a single scan', type=str, default=None)
parser.add_argument('-tc', '--testcampaign', default=None, help='TestCampaign')
parser.add_argument('-d', '--dia', type=str, default=None, help='diamond name')
args = parser.parse_args()

z = RunSelection(testcampaign=args.testcampaign)
if args.RunPlan is not None:
    z.print_testcampaign()
    z.select_runs_from_runplan(args.RunPlan)
    z.show_selected_runs()
else:
    z.show_run_plans(diamond=args.dia)
