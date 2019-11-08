#! /bin/bash

IFS='_'
cmds=($(src/DiamondRateScans.py -p -v $@))
for cmd in ${cmds[@]}
do
    echo ${cmd};
#    eval ${cmd};
    sleep 1;
done