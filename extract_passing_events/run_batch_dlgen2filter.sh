#!/bin/bash

OFFSET=$1
STRIDE=$2
SAMPLE_NAME=$3
INPUTLIST=$4
OUTPUT_DIR=$5

INPUTSTEM=larflowreco

# we assume we are already in the container
WORKDIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/extract_passing_events/
UBDL_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl_py3/
RECO_TEST_DIR=${UBDL_DIR}/larflow/larflow/Reco/test/
OUTPUT_LOGDIR=${WORKDIR}/logdir/${SAMPLE_NAME}/

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

cudadev="cpu"
echo "JOB ARRAYID: ${SLURM_ARRAY_TASK_ID} : DEVICE = ${cudadev}"

# LOCAL JOBDIR
local_jobdir=`printf /tmp/dev_dlgen2filter_jobid%d_%04d_${SAMPLE_NAME} ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf dev_dlgen2filter_${SAMPLE_NAME}_jobid%d_%04d.log ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "output logfile: "$local_logfile

#echo "SETUP CONTAINER/ENVIRONMENT"
cd ${UBDL_DIR}
alias python=python3
source /usr/local/root/root-6.22.06/bin/thisroot.sh
source setenv_py3.sh > /dev/null
source configure.sh > /dev/null
export PYTHONPATH=${LARMATCH_DIR}:${PYTHONPATH}
export PATH=${WORKDIR}:${PATH}

cd $local_jobdir

echo "STARTING TASK ARRAY ${SLURM_ARRAY_TASK_ID} for ${SAMPLE_NAME}" > ${local_logfile}

# run a loop. 

for ((i=0;i<${STRIDE};i++)); do

    # CALC JOB ID
    jobid=$(( ${start_jobid} + ${i} ))
    echo "JOBID ${jobid}" >> ${local_logfile}

    # GET FILE PATH
    let lineno=${jobid}+1
    let fileid=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $1 }'`
    kpsana=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $2 }'`
    dlmerged=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $3 }'`
    kpslarlite=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $4 }'`
    larmatch=`sed -n ${lineno}p ${INPUTLIST} | awk '{ print $5 }'`
    
    base_kpsana=$(basename $kpsana )
    base_dlmerged=$(basename $dlmerged)
    base_kpslarlite=$(basename $kpslarlite)
    base_larmatch=$(basename $larmatch)

    base_joboutfile=$(echo $base_kpslarlite | sed 's|larflowreco_|dlgen2filtered_|g' | sed 's|_larlite.root|.root|g')
    
    echo "kpsana path: $kpsana" >> ${local_logfile}
    echo "dlmerged path: $dlmerged" >> ${local_logfile}
    echo "kpslarlite path: $kpslarlite" >> ${local_logfile}
    echo "larmatch path: $larmatch" >> ${local_logfile}
    echo "output file: $base_joboutfile" >> ${local_logfile}

    # COPY LOCAL
    cp $kpsana $base_kpsana
    cp $dlmerged $base_dlmerged
    cp $kpslarlite $base_kpslarlite
    cp $larmatch $base_larmatch

    cmd="extract_passing_events $base_kpsana $base_dlmerged $base_kpslarlite $base_larmatch $base_joboutfile"
    echo "RUN PROGRAM" >> ${local_logfile}
    echo $cmd >> ${local_logfile} 2>&1
    $cmd >> ${local_logfile}
    
    # COPY OUTPUT TO FINAL FOLDER
    let fileid_dir=${fileid}/100
    mkdir -p $OUTPUT_DIR/$fileid_dir/
    cp $base_joboutfile $OUTPUT_DIR/$fileid_dir/

    # CLEAN UP
    rm *.root
done

cp $local_logfile ${OUTPUT_LOGDIR}/

# clean-up
cd /tmp
rm -r $local_jobdir
