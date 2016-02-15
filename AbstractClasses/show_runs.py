#!/usr/bin/env python
from RunSelection import RunSelection
from sys import argv

z=RunSelection()
if len(argv) > 1:
    z.select_runs_from_runplan(argv[1])
    z.show_selected_runs()
else:
    z.show_run_plans()
