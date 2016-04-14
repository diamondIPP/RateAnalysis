#!/bin/sh
echo $1 $2 $3
ipython -i AbstractClasses/PadAnalysis.py $1 -- $2 $3 $4
