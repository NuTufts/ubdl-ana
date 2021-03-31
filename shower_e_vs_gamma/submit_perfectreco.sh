#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=cheatreco
#SBATCH --output=cheatreco.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1-00:00:00
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
##SBATCH --partition=batch
##SBATCH --partition=preempt
#SBATCH --partition=wongjiradlab
##SBATCH --gres=gpu:p100:3
##SBATCH --partition ccgpu
##SBATCH --gres=gpu:t4:1
##SBATCH --nodelist=ccgpu01

container=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/shower_e_vs_gamma/
OFFSET=0
STRIDE=40

SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge
INPUTFILE=/cluster/tufts/wongjiradlab/nutufts/dlgen2prod/run3inputlists/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list
INPUTSTEM=merged_dlreco
FILEIDLIST=${RUN_DLANA_DIR}/runlist_reco_mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.txt

module load singularity

# CPU MODE
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_perfectreco.sh $OFFSET $STRIDE $SAMPLE_NAME ${INPUTFILE} ${INPUTSTEM} ${FILEIDLIST}"

