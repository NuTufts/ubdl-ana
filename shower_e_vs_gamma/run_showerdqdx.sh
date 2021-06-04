#!/bin/bash

OFFSET=$1
STRIDE=$2
SAMPLE_NAME=$3
INPUTLIST=$4 # contains list with columns [fileid] [dlmerged] [dlana]
INPUTSTEM=$5

# we assume we are already in the container

WORKDIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/shower_e_vs_gamma/
UBDL_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl_py3/
PYSCRIPT_DIR=${UBDL_DIR}/larflow/larflow/Reco/test
OUTPUT_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/shower_e_vs_gamma/output/${SAMPLE_NAME}/showerdqdx/
OUTPUT_LOGDIR=${WORKDIR}/logdir/${SAMPLE_NAME}

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

echo "JOB ARRAYID: ${SLURM_ARRAY_TASK_ID}"

# LOCAL JOBDIR
local_jobdir=`printf /tmp/showerll_jobid%d_%04d_${SAMPLE_NAME} ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf showerdqdx_${SAMPLE_NAME}_jobid%d_%04d.log ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
echo "output logfile: "$local_logfile

#echo "SETUP CONTAINER/ENVIRONMENT"
cd ${UBDL_DIR}
alias python=python3
source /usr/local/root/root-6.22.06/bin/thisroot.sh
source setenv_py3.sh > /dev/null
source configure.sh > /dev/null
export PYTHONPATH=${PYSCRIPT_DIR}:${PYTHONPATH}
export LD_LIBRARY_PATH=${WORKDIR}:${LD_LIBRARY_PATH}
export PATH=${WORKDIR}:${PATH}
export OMP_NUM_THREADS=4

cd $local_jobdir

echo "STARTING TASK ARRAY ${SLURM_ARRAY_TASK_ID} for ${SAMPLE_NAME}" > ${local_logfile}


# run a loop
for ((i=0;i<${STRIDE};i++)); do

    # CALC JOB ID
    jobid=$(( ${start_jobid} + ${i} ))
    echo "JOBID ${jobid}" >> ${local_logfile}

    # GET FILE ID
    let lineno=${jobid}+1    
    let fileid=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $1 }'`
    fileidname=`printf fileid%04d ${fileid}`
    
    # subfolder dir
    let nsubdir=${fileid}/100
    subdir=`printf %03d ${nsubdir}`
  
    # GET INPUT FILENAMES
    dlmergedpath=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $2 }'` # dlmerged file
    kpsrecopath=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $3 }'`  # kpsana
    dlmergedbase=$(basename $dlmergedpath )
    kpsrecobase=$(basename $kpsrecopath )

    # OUTPUT FILES
    dqdxbase=$(echo $dlmergedbase | sed 's|'"${INPUTSTEM}"'|showerdqdx_'"${fileidname}"'|g')
    dqdxpath=${OUTPUT_DIR}/${subdir}/${dqdxbase}

    echo "fileidname: $fileidname" >> ${local_logfile}
    echo "dlmerged path: $dlmergedpath" >> ${local_logfile}
    echo "kpsreco path: $kpsrecopath" >> ${local_logfile}    
    echo "dlmerged base: $dlmergedbase" >> ${local_logfile}
    echo "kpsreco base: $kpsrecobase" >> ${local_logfile}
    echo "OUTPUT dqdx base: $dqdxbase" >> ${local_logfile}
    echo "OUTPUT dqdx path: $dqdxpath" >> ${local_logfile}    

    echo "JOBID ${jobid} running FILEID ${fileid} with file: ${dlmergedbase}"

    # copy input files to local folder
    scp $dlmergedpath .
    scp $kpsrecopath .

    #let larmatchok=`python ${WORKDIR}/check_larmatch.py ${baselm}`
    let larmatchok=1
    echo "RESULT OF LARMATCH CHECK: $larmatchok" >> ${local_logfile}

    if [ $larmatchok -eq 1 ]
    then
    
	CMD="python3 ${PYSCRIPT_DIR}/dev_showerbilineardqdx.py --input-kpsana ${kpsrecobase} --input-dlmerged ${dlmergedbase} --output ${dqdxbase} -mc"	
	echo $CMD
	echo $CMD >> ${local_logfile}
	$CMD >> ${local_logfile}

        # copy to subdir in order to keep number of files per folder less than 100. better for file system.
	echo "COPY output to "${dqdxpath} >> ${local_logfile}
	mkdir -p $OUTPUT_DIR/${subdir}/
	cp ${dqdxbase} ${OUTPUT_DIR}/${subdir}/
    else
	echo "LARMATCH put} NOT OK" >> ${local_logfile}
    fi

    rm ${PWD}/${dlmergedbase}
    rm ${PWD}/${kpsrecobase}    
    rm -f ${PWD}/${dqdxbase}
done

# copy log to logdir
cp $local_logfile ${OUTPUT_LOGDIR}/

# clean-up
cd /tmp
rm -r $local_jobdir
