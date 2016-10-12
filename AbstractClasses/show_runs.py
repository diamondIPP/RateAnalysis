#!/usr/bin/env python
from AbstractClasses.RunSelection import RunSelection
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-R', '--RunPlan', help='Displays the deltailed overview about one Runplan', type=int, default=-1)
parser.add_argument('-tc', '--testcampaign', type=str, default='201610', help='TestCampaign')
args = parser.parse_args()

testcampaign = str(args.testcampaign)
print 'TestCampaingn: "{tc}"'.format(tc=testcampaign)
z = RunSelection(testcampaign=testcampaign)
if args.RunPlan != -1:
    print 'RunPlan:', args.RunPlan
    z.select_runs_from_runplan(args.RunPlan)
    z.show_selected_runs()
else:
    z.show_run_plans()
