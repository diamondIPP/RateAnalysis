#!/usr/bin/env python
from RunSelection import RunSelection
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('RunPlan', nargs='?', help='Displays the deltailed overview about a single scan', type=str, default=None)
parser.add_argument('-tc', '--testcampaign', type=str, default='201610', help='TestCampaign')
parser.add_argument('-d', '--dia', type=str, default=None, help='diamond name')
args = parser.parse_args()

testcampaign = str(args.testcampaign)
print 'TestCampaign: "{tc}"'.format(tc=testcampaign)

z = RunSelection(testcampaign=testcampaign)
if args.RunPlan is not None:
    z.select_runs_from_runplan(args.RunPlan)
    z.show_selected_runs()
else:
    z.show_run_plans(dia=args.dia)
