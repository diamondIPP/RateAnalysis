#! /bin/bash

IFS='_'
cmds=($(src/runplan_selection.py -p -v $@))
for cmd in ${cmds[@]}
do
    echo ${cmd};
#    eval ${cmd};
    sleep 1;
done