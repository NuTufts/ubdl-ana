#!/bin/bash

OFFSET=$1
STRIDE=$2
SAMPLE_NAME=$3

# we assume we are already in the container

WORKDIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/vertexing/
UBDL_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl/
INPUTLIST=${WORKDIR}/inputlists/${SAMPLE_NAME}.list
LARMATCH_DIR=${UBDL_DIR}/larflow/larmatchnet/
LARMATCH_RECO_DIR=${UBDL_DIR}/larflow/larflow/Reco/test
OUTPUT_DIR=${WORKDIR}/reco_outputdir/${SAMPLE_NAME}
LARMATCH_OUTDIR=${WORKDIR}/outputdir/${SAMPLE_NAME}
OUTPUT_LOGDIR=${WORKDIR}/reco_logdir/${SAMPLE_NAME}

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

# LOCAL JOBDIR
local_jobdir=`printf /tmp/kpsreco_jobid%04d_${SAMPLE_NAME} ${SLURM_ARRAY_TASK_ID}`
echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf kpsreco_${SAMPLE_NAME}_jobid%04d.log ${SLURM_ARRAY_TASK_ID}`
echo "output logfile: "$local_logfile

echo "SETUP CONTAINER/ENVIRONMENT"
cd ${UBDL_DIR}
source setenv.sh
source configure.sh
export PYTHONPATH=${LARMATCH_DIR}:${PYTHONPATH}

echo "GO TO JOBDIR"
cd $local_jobdir

echo "STARTING TASK ARRAY ${SLURM_ARRAY_TASK_ID} for ${SAMPLE_NAME}" > ${local_logfile}

# run a loop
for ((i=0;i<${STRIDE};i++)); do

    # jobid
    jobid=$(( ${start_jobid} + ${i} ))
    jobname=`printf jobid%04d ${jobid}`

    # subfolder dir
    let nsubdir=${jobid}/100
    subdir=`printf %03d ${nsubdir}`
  
    # GET DLMERGED INPUT FILENAME
    let lineno=${jobid}+1
    inputfile=`sed -n ${lineno}p ${INPUTLIST}`
    baseinput=$(basename $inputfile )
    echo "inputfile path: $inputfile"

    # LARMATCH RECO FILENAME
    larmatch_basename=$(echo $baseinput | sed 's|merged_dlreco|larmatch_kps|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}")
    larmatch_larlite_path=${LARMATCH_OUTDIR}/${subdir}/${larmatch_basename}_larlite.root
    larmatch_larcv_path=${LARMATCH_OUTDIR}/${subdir}/${larmatch_basename}_larcv.root # not used

    # KPSRECO FILENAME
    kpsreco_basename=$(echo $baseinput | sed 's|merged_dlreco|kpsreco|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}")
    kpsreco_rootname=${kpsreco_basename}.root

    echo "========================================"
    echo "JOB: ${jobid}"
    echo "[input] dlmerged: ${inputfile}"
    echo "[input] larmatch: ${larmatch_larlite_path}"
    echo "[output] kpsreco: ${kpsreco_rootname}"
    CMD="python ${LARMATCH_RECO_DIR}/run_kpsrecoman.py --input-dlmerged $inputfile --input-larflow ${larmatch_larlite_path} --output ${kpsreco_rootname} -tb"
    echo $CMD
    $CMD >> ${local_logfile} 2>&1

    # copy to subdir in order to keep number of files per folder less than 100. better for file system.
    echo "COPY output to "${OUTPUT_DIR}/${subdir}/
    mkdir -p $OUTPUT_DIR/${subdir}/
    cp ${kpsreco_basename}* $OUTPUT_DIR/${subdir}/
done

# copy log to logdir
cp $local_logfile $OUTPUT_LOGDIR/

# clean-up
cd /tmp
rm -r $local_jobdir
