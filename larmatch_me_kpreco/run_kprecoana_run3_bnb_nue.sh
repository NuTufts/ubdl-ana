#!/bin/bash

OFFSET=0
STRIDE=10
INPUTSTEM="merged_dlreco"

# we assume we are already in the container
SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge
WORKDIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_me_kpreco/
FILEIDLIST=${WORKDIR}/runlist_kprecoana_mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge.txt
UBDL_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ubdl/
OUTPUT_DIR=${WORKDIR}/output/${SAMPLE_NAME}
OUTPUT_LOGDIR=${WORKDIR}/logdir/${SAMPLE_NAME}

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

cudadev="cpu"
echo "JOB ARRAYID: ${SLURM_ARRAY_TASK_ID} : DEVICE = ${cudadev}"

# LOCAL JOBDIR
local_jobdir=`printf /tmp/kprecoana_jobid%d_%04d_${SAMPLE_NAME} ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf kprecoana_${SAMPLE_NAME}_jobid%d_%04d.log ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "output logfile: "$local_logfile

#echo "SETUP CONTAINER/ENVIRONMENT"
cd ${UBDL_DIR}
alias python=python3
source setenv_py3.sh > /dev/null
source configure.sh > /dev/null
export PATH=${WORKDIR}/build/src:${PATH}
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
    let fileid=`sed -n ${lineno}p ${FILEIDLIST} | awk '{ print $1 }'`
  
    # GET INPUT FILENAME
    inputfile=`sed -n ${lineno}p ${FILEIDLIST} | awk '{ print $2 }'`
    lminput=`sed -n ${lineno}p ${FILEIDLIST} | awk '{ print $3 }'`    
    baseinput=$(basename $inputfile )
    baselm=$(basename $lminput)
    echo "inputfile path: $inputfile" >> ${local_logfile}
    echo "larmatch path: $lminput" >> ${local_logfile}    
    echo "baseinput: $baseinput" >> ${local_logfile}
    echo "base larmatch: $baselm" >> ${local_logfile}    

    echo "JOBID ${jobid} running FILEID ${fileid} with file: ${baseinput}"

    # local outfile
    jobname=`printf jobid%04d ${jobid}`
    fileidname=`printf fileid%04d ${fileid}`
    local_outfile=$(echo $baseinput  | sed 's|'"${INPUTSTEM}"'|kprecoana_'"${fileidname}"'|g')
    local_basename=$(echo $baseinput | sed 's|'"${INPUTSTEM}"'|kprecoana_'"${fileidname}"'|g' | sed 's|.root||g')
    echo "outfile : "$local_outfile >> ${local_logfile}
    scp $inputfile $baseinput
    scp $lminput $baselm

    CMD="ana_kpreco ${baseinput} ${baselm} ${local_outfile}"
    echo $CMD
    echo $CMD >> ${local_logfile}
    $CMD >> ${local_logfile}

    # subfolder dir
    let nsubdir=${fileid}/100
    subdir=`printf %03d ${nsubdir}`

    # copy to subdir in order to keep number of files per folder less than 100. better for file system.
    echo "COPY output to "${OUTPUT_DIR}/${subdir}/ >> ${local_logfile}
    mkdir -p $OUTPUT_DIR/${subdir}/
    cp ${local_basename}*ana.root $OUTPUT_DIR/${subdir}/

    rm ${PWD}/*.root
done

# copy log to logdir
cp $local_logfile $OUTPUT_LOGDIR/

# clean-up
cd /tmp
rm -r $local_jobdir
