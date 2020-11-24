RATEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
alias analyse='ipython -i $RATEDIR/analyse.py -- $@'
alias dutanalyse='ipython -i $RATEDIR/src/dut_analysis.py -- $@'
alias DrawCurrents='ipython -i $RATEDIR/src/currents.py -- $@'
alias ShowRunplans='python $RATEDIR/src/run_selection.py -s'
alias ShowSelection='python $RATEDIR/src/runplan_selection.py -s -v'
alias ShowSelections='python $RATEDIR/src/runplan_selection.py -sa -v'
alias RunSelection='ipython -i $RATEDIR/src/run_selection.py -- $@'
alias DiaSelection='ipython -i $RATEDIR/src/runplan_selection.py -- $@'
source $RATEDIR/../root6/rootinstall/bin/thisroot.sh
source $RATEDIR/venv/bin/activate
