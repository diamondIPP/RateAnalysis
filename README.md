# Analysis software for the diamond rate beam tests at PSI

## Requirements

- ROOT6.22+
- python3.6+

<ins>python packages:
- termcolor, numpy, scipy, screeninfo, progressbar, uncertainties, pytz, gtts, h5py

Install python virtual environment with all packages:
```shell
make venv
```

## Add aliases
- source .bash_aliases
- adds following aliases:
    - analyse (for all single runs and runplans)
    - DrawCurrents
    - ShowRunplans
    - MasterSelection
    - ShowSelection
    - ShowSelections
    - RunSelecion
    - DiaSelection
