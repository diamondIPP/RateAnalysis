#! /bin/bash

IFS='?'
cmds=($(python src/analysis_collection.py -p -v "$@"))

# parallel
parallel --timeout 3600 --bar "eval" ::: "${cmds[@]}"
wait

# sequential
#for cmd in "${cmds[@]}"
#do
#    echo "${cmd}";
#    eval "${cmd}";
#    sleep .1;
#done