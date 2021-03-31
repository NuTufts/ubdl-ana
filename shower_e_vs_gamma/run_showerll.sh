#!/bin/bash

OFFSET=$1
STRIDE=$2
SAMPLE_NAME=$3
INPUTLIST=$4
INPUTSTEM=$5
FILEIDLIST=$6

# we assume we are already in the container

WORKDIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/shower_e_vs_gamma/
UBDL_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl_py3/
OUTPUT_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/shower_e_vs_gamma/output/${SAMPLE_NAME}/showerll/
PERFECT_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/shower_e_vs_gamma/output/${SAMPLE_NAME}/perfectreco/
OUTPUT_LOGDIR=${WORKDIR}/logdir/${SAMPLE_NAME}

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_ANA_DIR
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
local_logfile=`printf showerll_${SAMPLE_NAME}_jobid%d_%04d.log ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "output logfile: "$local_logfile

#echo "SETUP CONTAINER/ENVIRONMENT"
cd ${UBDL_DIR}
alias python=python3
source /usr/local/root/root-6.22.06/bin/thisroot.sh
source setenv_py3.sh > /dev/null
source configure.sh > /dev/null
export PYTHONPATH=${LARMATCH_DIR}:${PYTHONPATH}
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
    let fileid=`sed -n ${lineno}p ${FILEIDLIST} | awk '{ print $1 }'`
    fileidname=`printf fileid%04d ${fileid}`    

    # subfolder dir
    let nsubdir=${fileid}/100
    subdir=`printf %03d ${nsubdir}`
  
    # GET INPUT FILENAME
    inputfile=`sed -n ${lineno}p ${FILEIDLIST} | awk '{ print $2 }'`
    lminput=`sed -n ${lineno}p ${FILEIDLIST} | awk '{ print $3 }'`    
    baseinput=$(basename $inputfile )
    baselm=$(basename $lminput)
    perfectreco_base=$(echo $baseinput | sed 's|'"${INPUTSTEM}"'|perfectreco_'"${fileidname}"'|g' | sed 's|.root|_perfectreco.root|g')
    perfectreco_path=${PERFECT_DIR}/${subdir}/${perfectreco_base}

    echo "fileidname: $fileidname" >> ${local_logfile}
    echo "inputfile path: $inputfile" >> ${local_logfile}
    echo "larmatch path: $lminput" >> ${local_logfile}    
    echo "baseinput: $baseinput" >> ${local_logfile}
    echo "base larmatch: $baselm" >> ${local_logfile}
    echo "perfect reco: $perfectreco_path" >> ${local_logfile}

    echo "JOBID ${jobid} running FILEID ${fileid} with file: ${baseinput}"

    # local outfile
    jobname=`printf jobid%04d ${jobid}`
    local_outfile=$(echo $baseinput  | sed 's|'"${INPUTSTEM}"'|showerll_'"${fileidname}"'|g')
    echo "outfile : "$local_outfile >> ${local_logfile}
    scp $inputfile $baseinput
    scp $perfectreco_path .

    #let larmatchok=`python ${WORKDIR}/check_larmatch.py ${baselm}`
    let larmatchok=1
    echo "RESULT OF CHECK ${baselm}: $larmatchok" >> ${local_logfile}

    if [ $larmatchok -eq 1 ]
    then
    
	CMD="extract_variables $perfectreco_base $baseinput $local_outfile"
	echo $CMD
	echo $CMD >> ${local_logfile}
	$CMD >> ${local_logfile}

        # copy to subdir in order to keep number of files per folder less than 100. better for file system.
	echo "COPY output to "${OUTPUT_DIR}/${subdir}/ >> ${local_logfile}
	mkdir -p $OUTPUT_DIR/${subdir}/
	cp ${local_outfile} ${OUTPUT_DIR}/${subdir}/
    else
	echo "LARMATCH FILE ${lminput} NOT OK" >> ${local_logfile}
    fi
    #rm ${PWD}/${baseinput}
    #rm ${PWD}/${perfectreco_base}
done

# copy log to logdir
cp $local_logfile $OUTPUT_LOGDIR/

# clean-up
cd /tmp
rm -r $local_jobdir
