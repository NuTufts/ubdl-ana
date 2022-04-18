#!/bin/bash

UBDL_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ubdl
WORKDIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/
OUTDIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/output_pi0_filtered
PYSCRIPT=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/filter_mc.py
INPUTLIST=${WORKDIR}/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list

stride=100
jobid=${SLURM_ARRAY_TASK_ID}
let startline=$(expr "${stride}*${jobid}")

mkdir -p $WORKDIR
jobworkdir=`printf "%s/working/pi0filtermc_jobid_%03d" $WORKDIR $jobid`
mkdir -p $jobworkdir
mkdir -p $OUTPUT_DIR

local_jobdir=`printf /tmp/pi0filtermc_jobid%03d $jobid`
rm -rf $local_jobdir
mkdir -p $local_jobdir

cd $local_jobdir
touch log_jobid${jobid}.txt
local_logfile=`echo ${local_jobdir}/log_jobid${jobid}.txt`

cd $UBDL_DIR
source setenv_py3.sh >> ${local_logfile} 2>&1
source configure.sh >>	${local_logfile} 2>&1
cd $local_jobdir

CMD="python3 ${PYSCRIPT}"
echo "SCRIPT: ${PYSCRIPT}" >> ${local_logfile} 2>&1
echo "startline: ${startline}" >> ${local_logfile} 2>&1

for i in {1..100}
do
    let lineno=$startline+$i
    dlmerged=`sed -n ${lineno}p $INPUTLIST`
    dlmerged_base=`basename ${dlmerged}`
    dlmerged_dir=`dirname ${dlmerged}`

    COMMAND="python3 ${PYSCRIPT} --dlmerged ${dlmerged} --filelist ${WORKDIR}/ListforAaroosh_pi0selection_bnb_overlay_run3.txt"
    echo $COMMAND
    echo $COMMAND >> ${local_logfile} 2>&1
    #$COMMAND >> ${local_logfile} 2>&1
    $COMMAND >> ${local_logfile}
    xmerged=`ls pi0select*_larlite.root | sed 's|larlite\.root|merged\.root|g'`
    hadd -f $xmerged pi0select_*.root
    cp pi0select_*merged.root ${OUTDIR}/
    rm pi0select_*.root
done

cp ${local_logfile} ${jobworkdir}/

cd /tmp
rm -r $local_jobdir

