#!/bin/bash

# This runs inside the container and builds the container
# We assume you did the first-time setup already

UBDLDIR=$1
ANADIR=$2
FILELIST=$3

let arrayid="$SLURM_ARRAY_TASK_ID+0"
echo $arrayid

cd $UBDLDIR
# setup container
source setenv_py3.sh
source configure.sh

#now move to working dir
cd $ANADIR
pwd

#split up lines into inputs
let line=${arrayid}+1
inline=`sed -n ${line}p ${FILELIST}`
IFS=' '
read -a strarr <<< "$inline"
kpsfile="${strarr[1]}"
echo ${kpsfile}
larcvfile="${strarr[3]}"
echo ${larcvfile}

./showerenergy ${kpsfile} ${larcvfile} outputs plots_${arrayid} > logs/log_${arrayid}.log

echo Finished
