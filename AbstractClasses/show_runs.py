#!/usr/bin/env python
from RunSelection import RunSelection
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-R","--RunPlan",help = "Displays the deltailed overview about one Runplan", type=int,default = -1)

parser.add_argument("-tc", "--testcampaign", type=str,default="201510",
                    help="TestCampaign: default=201510")
args = parser.parse_args()
print args

testcampaign=str(args.testcampaign)
print 'TestCampaingn: "{tc}"'.format(tc=testcampaign)
a = Elementary(testcampaign=testcampaign)
a.print_testcampaign()
z=RunSelection(testcampaign=testcampaign)
if args.RunPlan != -1:
    print 'RunPlan:',args.RunPlan
    z.select_runs_from_runplan(args.RunPlan)
    z.show_selected_runs()
else:
    z.show_run_plans()
