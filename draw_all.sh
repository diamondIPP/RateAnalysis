#! /bin/bash

IFS='_'
cmds=($(AbstractClasses/DiamondRateScans.py -p -tc $1))
for cmd in ${cmds[@]}
do
    eval ${cmd};
    sleep 1;
done