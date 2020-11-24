#!/bin/bash

OFFSET=$1
STRIDE=$2
SAMPLE_NAME=$3

# we assume we are already in the container

WORKDIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/makedata_spembednet/
UBDL_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl/
INPUTLIST=${WORKDIR}/inputlists/${SAMPLE_NAME}.txt
LARMATCH_DIR=${UBDL_DIR}/larflow/larmatchnet/
SPEMBED_DIR=${UBDL_DIR}/larflow/larflow/SpatialEmbed/test/
WEIGHTS_DIR=${UBDL_DIR}/larflow/larmatchnet/grid_deploy_scripts/larmatch_kps_weights/
WEIGHT_FILE=checkpoint.1974000th.tar
OUTPUT_DIR=${WORKDIR}/outputdir/${SAMPLE_NAME}
OUTPUT_LOGDIR=${WORKDIR}/logdir/${SAMPLE_NAME}

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

# LOCAL JOBDIR
local_jobdir=`printf /tmp/larmatch_kps_jobid%04d_${SAMPLE_NAME} ${SLURM_ARRAY_TASK_ID}`
echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf larmatch_kps_${SAMPLE_NAME}_jobid%04d.log ${SLURM_ARRAY_TASK_ID}`
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
    dirinput=$(dirname $inputfile)
    mcinfoname=$(echo $baseinput | sed 's|larcvtruth|mcinfo|g')
    mcinfofile=$dirinput/$mcinfoname
    echo "inputfile path: $inputfile"
    echo "mcinfo path: $mcinfofile"

    # local outfile
    jobname=`printf jobid%04d ${jobid}`
    local_outfile=$(echo $baseinput | sed 's|larcvtruth|larmatch_kps|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}.root")
    local_basename=$(echo $baseinput | sed 's|larcvtruth|larmatch_kps|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}")
    spembedout=$(echo $baseinput | sed 's|larcvtruth|spatialembed_voxels|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}.root")
    spembedout_base=$(echo $baseinput | sed 's|larcvtruth|spatialembed_voxels|g' | sed 's|.root||g' | xargs -I{} echo {}"-${jobname}")
    echo "outfile : "$local_outfile

    # copy input
    cp $inputfile $baseinput
    cp $mcinfofile $mcinfoname
    cp ${WEIGHTS_DIR}/${WEIGHT_FILE} ${WEIGHT_FILE}
    
    CMD="python $LARMATCH_DIR/deploy_kps_larmatch.py --supera $baseinput --weights ${WEIGHT_FILE} --output $local_outfile --min-score 0.3 --adc-name wiremc --chstatus-name wiremc --device-name cpu --use-unet -tb"
    echo $CMD
    $CMD >> ${local_logfile} 2>&1

    CMD2="python $SPEMBED_DIR/prep_spatialembed.py --input-larmatch ${local_basename}_larlite.root --input-larcv $baseinput --input-larlitetruth $mcinfoname -adc wiremc -o $spembedout"
    echo $CMD2
    $CMD2 >> ${local_logfile} 2>&1

    # subfolder dir
    let nsubdir=${jobid}/100
    subdir=`printf %03d ${nsubdir}`

    # copy to subdir in order to keep number of files per folder less than 100. better for file system.
    echo "COPY output to "${OUTPUT_DIR}/${subdir}/
    mkdir -p $OUTPUT_DIR/${subdir}/
    cp ${spembedout_base}_s3dembed.root $OUTPUT_DIR/${subdir}/
    rm ${local_basename}*
    rm ${spembedout_base}*
    rm ${baseinput}
    rm ${mcinfoname}
    rm ${WEIGHT_FILE}
done

# copy log to logdir
cp $local_logfile $OUTPUT_LOGDIR/

# clean-up
cd /tmp
rm -r $local_jobdir
