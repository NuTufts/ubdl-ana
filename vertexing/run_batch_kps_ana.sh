#!/bin/bash

OFFSET=$1
STRIDE=$2
SAMPLE_NAME=$3
ISMC=$4

# we assume we are already in the container

WORKDIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/vertexing/
UBDL_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl/
INPUTLIST=${WORKDIR}/inputlists/${SAMPLE_NAME}.list
LARMATCHANA_DIR=${UBDL_DIR}/larflow/larflow/Ana/

LARMATCH_DIR=${UBDL_DIR}/larflow/larmatchnet/
GRID_SCRIPTS_DIR=${UBDL_DIR}/larflow/larmatchnet/grid_deploy_scripts/

RECO_OUTDIR=${WORKDIR}/reco_outputdir/${SAMPLE_NAME}
OUTPUT_DIR=${WORKDIR}/ana_outputdir/${SAMPLE_NAME}
OUTPUT_LOGDIR=${WORKDIR}/ana_logdir/${SAMPLE_NAME}

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

# LOCAL JOBDIR
local_jobdir=`printf /tmp/kpsana_jobid%04d_${SAMPLE_NAME} ${SLURM_ARRAY_TASK_ID}`
echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf kpsana_${SAMPLE_NAME}_jobid%04d.log ${SLURM_ARRAY_TASK_ID}`
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

    jobid=$(( ${start_jobid} + ${i} ))
    echo "JOBID ${jobid}"
  
    # GET INPUT FILENAME
    let lineno=${jobid}+1
    inputfile=`sed -n ${lineno}p ${INPUTLIST}`
    baseinput=$(basename $inputfile )
    echo "inputfile path: $inputfile"

    # jobname
    jobname=`printf jobid%04d ${jobid}`

    # subfolder dir
    let nsubdir=${jobid}/100
    subdir=`printf %03d ${nsubdir}`
    
    # KPSRECO FILENAME
    kpsreco_basename=$(echo $baseinput | sed 's|merged_dlreco|kpsreco|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}")
    kpsreco_rootname=${kpsreco_basename}_kpsrecomanagerana.root
    kpsreco_path=${RECO_OUTDIR}/${subdir}/${kpsreco_rootname}

    # get name of ana file
    local_outfile=$(echo $baseinput  | sed 's|merged_dlreco|kps_vertexana|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}.root")
    echo "outfile : "$local_outfile
    
    CMD="${LARMATCHANA_DIR}/./kpsreco_vertexana ${inputfile} ${kpsreco_path} ${local_outfile} ${ISMC}"
    echo $CMD
    $CMD >> ${local_logfile} 2>&1

    # copy to subdir in order to keep number of files per folder less than 100. better for file system.
    echo "COPY output to "${OUTPUT_DIR}/${subdir}/
    mkdir -p $OUTPUT_DIR/${subdir}/
    cp ${local_outfile} $OUTPUT_DIR/${subdir}/
done

# copy log to logdir
cp $local_logfile $OUTPUT_LOGDIR/

# clean-up
cd /tmp
rm -r $local_jobdir
