#!/bin/bash

OFFSET=0
STRIDE=50
SAMPLE_NAME=pi0select
INPUTSTEM=pi0select

# we assume we are already in the container

WORKDIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/
UBDL_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ubdl/
RECO_TEST_DIR=${UBDL_DIR}/larflow/larflow/Reco/test/
OUTPUT_DIR=${WORKDIR}/output_reco/larlite/
OUTPUT_ANA_DIR=${WORKDIR}/output_reco/ana/
OUTPUT_LOGDIR=${WORKDIR}/logdir/rec/
FILEIDLIST=${WORKDIR}/larflowreco_runlist.txt

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_ANA_DIR
mkdir -p $OUTPUT_LOGDIR

# WE WANT TO RUN MULTIPLE FILES PER JOB IN ORDER TO BE GRID EFFICIENT
start_jobid=$(( ${OFFSET} + ${SLURM_ARRAY_TASK_ID}*${STRIDE}  ))

cudadev="cpu"
echo "JOB ARRAYID: ${SLURM_ARRAY_TASK_ID} : DEVICE = ${cudadev}"

# LOCAL JOBDIR
local_jobdir=`printf /tmp/larflowreco_jobid%d_%04d_${SAMPLE_NAME} ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "local jobdir: $local_jobdir"
rm -rf $local_jobdir
mkdir -p $local_jobdir

# local log file
local_logfile=`printf larflowreco_${SAMPLE_NAME}_jobid%d_%04d.log ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#echo "output logfile: "$local_logfile

#echo "SETUP CONTAINER/ENVIRONMENT"
cd ${UBDL_DIR}
alias python=python3
source setenv_py3.sh > /dev/null
source configure.sh > /dev/null
export PYTHONPATH=${LARMATCH_DIR}:${PYTHONPATH}
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
    local_outfile=$(echo $baseinput  | sed 's|'"${INPUTSTEM}"'|larflowreco_'"${fileidname}"'|g')
    local_basename=$(echo $baseinput | sed 's|'"${INPUTSTEM}"'|larflowreco_'"${fileidname}"'|g' | sed 's|.root||g')
    echo "outfile : "$local_outfile >> ${local_logfile}
    scp $inputfile $baseinput
    scp $lminput $baselm

    let larmatchok=`python ${WORKDIR}/check_larmatch.py ${baselm}`
    echo "RESULT OF CHECK ${baselm}: $larmatchok" >> ${local_logfile}

    if [ $larmatchok -eq 1 ]
    then
    
	CMD="python $RECO_TEST_DIR/run_kpsrecoman.py --input-dlmerged ${baseinput}  --input-larflow ${baselm} --output ${local_outfile} -tb -mc --products rerun "
	echo $CMD >> ${local_logfile}
	$CMD >> ${local_logfile}

        # subfolder dir
	let nsubdir=${fileid}/100
	subdir=`printf %03d ${nsubdir}`

        # copy to subdir in order to keep number of files per folder less than 100. better for file system.
	echo "COPY output to "${OUTPUT_DIR}/${subdir}/ >> ${local_logfile}
	mkdir -p $OUTPUT_DIR/${subdir}/
	mkdir -p $OUTPUT_ANA_DIR/${subdir}/
	cp ${local_basename}*larlite.root $OUTPUT_DIR/${subdir}/
	cp ${local_basename}*ana.root $OUTPUT_ANA_DIR/${subdir}/
    else
	echo "LARMATCH FILE ${lminput} NOT OK" >> ${local_logfile}
    fi
    #echo "rm ${PWD}/${local_basename}*"
    #echo "rm ${PWD}/${baseinput}"
    #echo "rm ${PWD}/${baselm}"
    rm ${PWD}/${local_basename}*
    rm ${PWD}/${baseinput}
    rm ${PWD}/${baselm}
done

# copy log to logdir
cp $local_logfile $OUTPUT_LOGDIR/

# clean-up
cd /tmp
rm -r $local_jobdir
