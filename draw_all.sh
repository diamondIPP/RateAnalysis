#! /bin/bash

IFS='_'
cmds=($(AbstractClasses/DiamondRateScans.py -p -tc 201610))
for cmd in ${cmds[@]}
do
    eval $cmd;
done