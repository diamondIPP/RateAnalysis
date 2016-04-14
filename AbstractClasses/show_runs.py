#!/usr/bin/env python
from RunSelection import RunSelection
from sys import argv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("RunPlan",help = "Displays the deltailed overview about one Runplan")

parser.add_argument("-tc", "--testcampaign", type=int,default=201510,
                    help="TestCampaign: default=201510")
args = parser.parse_args()
print args
for arg in args:
    print '*',arg

raw_input()
z=RunSelection()
if len(argv) > 1:
    z.set_test_campaign(str(args.testcampaign))
    z.select_runs_from_runplan(argv[1])
    z.show_selected_runs()
else:
    z.show_run_plans()
