DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
alias analyse='ipython -i $DIR/analyse.py -- $@'
alias dutanalyse='ipython -i $DIR/src/dut_analysis.py -- $@'
alias DrawCurrents='ipython -i ~/software/RateAnalysis/src/currents.py -- $@'
alias ShowRunplans='python ~/software/RateAnalysis/src/run_selection.py -s'
alias MasterSelection='python ~/software/RateAnalysis/src/run_selection.py -ms'
alias ShowSelection='python ~/software/RateAnalysis/src/runplan_selection.py -s -v'
alias ShowSelections='python ~/software/RateAnalysis/src/runplan_selection.py -sa -v'
alias RunSelection='ipython -i ~/software/RateAnalysis/src/run_selection.py -- $@'
alias DiaSelection='ipython -i ~/software/RateAnalysis/src/runplan_selection.py -- $@'
source $DIR/../root6/build/bin/thisroot.sh
source $DIR/venv/bin/activate
