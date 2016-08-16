#!/bin/bash
# alias.sh

shopt -s expand_aliases

## load the configuration parameters
source configParallel.ini
NUMFILES=`ls ${RAWFILESDIR} -1f | grep ${RAWFILESSTEM} | wc -l`
RAWFILES=`ls ${RAWFILESDIR}`
RAWFILES2=""
RUNNUMS=""
iii=0
for i in $RAWFILES; do
	temp=${i#${RAWFILESSTEM}}
	temp=${temp%${RAWEXT}}
	if (($temp>=${INI}&&$temp<=${FIN}))
	then
		RAWFILES2="${RAWFILES2} ${i}"
		RUNNUMS="${RUNNUMS} ${temp}"
		iii=$[$iii+1]
	fi
done
echo "Processing ${iii} files"

sequ1() {
    source configParallel.ini
    convertedfile0="${ROOTFILESDIR}/test${TESTCAMP}00${2}${ROOTEXT}"
    find ${convertedfile0} -size +1k -exec rm {} +
    convertedfile1="${ROOTFILESDIR}/${ROOTFILESSTEM}${2}${ROOTEXT}"
    find ${convertedfile1} -size -1000k -exec rm {} +
    if [ ! -f "${convertedfile1}" ]; then
	    echo "Converting run ${2}..."
	    echo IPY ${1} -- ${2}
	    IPY $1 -- $2 &> /dev/null
#	    echo $1 -- $2
	    echo "Finished converting run ${2}."
	fi
	sleep 1s
	find ${convertedfile0} -size +1k -exec rm {} +
	convertedfile1="${ROOTFILESDIR}/${ROOTFILESSTEM}${2}${ROOTEXT}"
	find ${convertedfile1} -size -1000k -exec rm {} +
	## sizefile=`wc -c "${convertedfile}" | cut -d' ' -f5`
	## if ((${sizefile}<=1000000))
	## then
}
## export sequ1 to make it usable by gnuparallel
export -f sequ1

sequ2() {
    source configParallel.ini
    convertedfile0="${ROOTFILESDIR}/test${TESTCAMP}00${2}${ROOTEXT}"
    find ${convertedfile0} -size +1k -exec rm {} +
    temp="${RESULTSFOLDER}/${TESTCAMP}_${2}"
    if [ ! -d "${temp}" ]; then
	    echo "Analysing run ${2}..."
	    echo IPY ${1} -- ${2}
	    IPY $1 -- $2 &> /dev/null
#	    echo $1 -- $2
	    echo "Finished analysing run ${2}."
	else
	    numplots=`find "${temp}/Plots" -maxdepth 1 -type f|wc -l`
        if ((${numplots}<${NUMELEMENTS}))
        then
            echo "Analysing run ${2}..."
            echo IPY ${1} -- ${2}
            IPY $1 -- $2 &> /dev/null
#           echo $1 -- $2
            echo "Finished analysing run ${2}."
        fi
	fi
	sleep 1s
	find ${convertedfile0} -size +1k -exec rm {} +
	convertedfile1="${ROOTFILESDIR}/${ROOTFILESSTEM}${2}${ROOTEXT}"
	find ${convertedfile1} -size -1000k -exec rm {} +
	## sizefile=`wc -c "${convertedfile}" | cut -d' ' -f5`
	## if ((${sizefile}<=1000000))
	## then
}
## export sequ2 to make it usable by gnuparallel
export -f sequ2
sleep 1s

echo "Starting parallel conversion..."

parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ1 ${ANAPATH} %%" ::: $RUNNUMS
wait
parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ1 ${ANAPATH} %%" ::: $RUNNUMS
wait
parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ1 ${ANAPATH} %%" ::: $RUNNUMS
wait
#parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ1 ${ANAPATH} %%" ::: $RUNNUMS
#wait
#parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ1 ${ANAPATH} %%" ::: $RUNNUMS
#wait
#parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ1 ${ANAPATH} %%" ::: $RUNNUMS
#wait

echo "Finished parallel conversion."
echo "Converting missing files..."

for i in $RUNNUMS; do
	temp="${ROOTFILESDIR}/${ROOTFILESSTEM}${i}${ROOTEXT}"
	if [ ! -f "${temp}" ]; then
		sequ1 ${ANAPATH} ${i}
	fi
done
sleep 1s

for i in $RUNNUMS; do
	temp="${ROOTFILESDIR}/${ROOTFILESSTEM}${i}${ROOTEXT}"
	if [ ! -f "${temp}" ]; then
		sequ1 ${ANAPATH} ${i}
	fi
done
sleep 1s

for i in $RUNNUMS; do
	temp="${ROOTFILESDIR}/${ROOTFILESSTEM}${i}${ROOTEXT}"
	if [ ! -f "${temp}" ]; then
		sequ1 ${ANAPATH} ${i}
	fi
done
sleep 1s

for i in $RUNNUMS; do
	temp="${ROOTFILESDIR}/${ROOTFILESSTEM}${i}${ROOTEXT}"
	if [ ! -f "${temp}" ]; then
		sequ1 ${ANAPATH} ${i}
	fi
done
sleep 1s

for i in $RUNNUMS; do
	temp="${ROOTFILESDIR}/${ROOTFILESSTEM}${i}${ROOTEXT}"
	if [ ! -f "${temp}" ]; then
		sequ1 ${ANAPATH} ${i}
	fi
done

echo "All files converted!"
echo "Starting files analysis..."

##parallel --timeout 200% --bar -I %% "sequ2 ${ANAPATH} %%" ::: $RUNNUMS
##parallel --timeout 200% --bar --dry-run -v -I %% --use-cpus-instead-of-cores "sequ2 ${ANAPATH} %%" ::: $RUNNUMS

parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ2 ${ANAPIXPATH} %%" ::: $RUNNUMS
wait
parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ2 ${ANAPIXPATH} %%" ::: $RUNNUMS
wait
parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ2 ${ANAPIXPATH} %%" ::: $RUNNUMS
wait
parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ2 ${ANAPIXPATH} %%" ::: $RUNNUMS
wait
parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ2 ${ANAPIXPATH} %%" ::: $RUNNUMS
wait
parallel --timeout 3600 --bar -I %% --use-cpus-instead-of-cores "sequ2 ${ANAPIXPATH} %%" ::: $RUNNUMS
wait

echo "Finished with parallel analysis"
echo "Analysing missing files..."

for i in $RUNNUMS; do
	temp="${RESULTSFOLDER}/${TESTCAMP}_${i}"
	if [ ! -d "${temp}" ]; then
		sequ2 ${ANAPIXPATH} ${i}
	else
	    numplots=`find "${temp}/Plots" -maxdepth 1 -type f|wc -l`
	    if ((${numplots}<${NUMELEMENTS}))
	    then
	        sequ2 ${ANAPIXPATH} ${i}
	    fi
	fi
done

echo "Finished all analysis!"
