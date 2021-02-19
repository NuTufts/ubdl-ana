#!/bin/bash

OFFSET=$1
STRIDE=$2
SAMPLE_NAME=$3
INPUTLIST=$4

INPUTSTEM=larflowreco

# we assume we are already in the container
WORKDIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/dev_1e1p_selection/
UBDL_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl_py3/
RECO_TEST_DIR=${UBDL_DIR}/larflow/larflow/Reco/test/
OUTPUT_DIR=${WORKDIR}/output/${SAMPLE_NAME}/
OUTPUT_LOGDIR=${WORKDIR}/logdir/${SAMPLE_NAME}/

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

cudadev="cpu"
echo "JOB ARRAYID: ${SLURM_ARRAY_TASK_ID} : DEVICE = ${cudadev}"

# LOCAL JOBDIR
local_jobdir=`printf /tmp/dev_1e1p_sel_jobid%d_%04d_${SAMPLE_NAME} ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf dev_1e1p_sel_${SAMPLE_NAME}_jobid%d_%04d.log ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
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
# we copy files to tmp, then run the plotter on it?
touch input.list
outfile=`printf plots_1e1p_sel_${SAMPLE_NAME}_%03d.root ${SLURM_ARRAY_TASK_ID}`

for ((i=0;i<${STRIDE};i++)); do

    # CALC JOB ID
    jobid=$(( ${start_jobid} + ${i} ))
    echo "JOBID ${jobid}" >> ${local_logfile}

    # GET FILE PATH
    let lineno=${jobid}+1    
    inputfile=`sed -n ${lineno}p ${INPUTLIST}`
    baseinput=$(basename $inputfile )
    echo "inputfile path: $inputfile" >> ${local_logfile}
    echo "baseinput: $baseinput" >> ${local_logfile}

    echo "JOBID ${jobid} running FILEID ${fileid} with file: ${baseinput}"
    echo $inputfile >> input.list

done

make_1e1p_sel_plots input.list >> ${local_logfile} 2>&1
mv plots_1e1p_sel.root $outfile
outlist=`printf input_${SAMPLE_NAME}_%03d.list ${SLURM_ARRAY_TASK_ID}`
mv input.list $outlist

cp $outfile $OUTPUT_DIR/
cp $outlist $OUTPUT_DIR/

# clean-up
cd /tmp
rm -r $local_jobdir
