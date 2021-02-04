#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=ubdlvareco
#SBATCH --output=ubdl_vareco.log
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
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/va_study/
OFFSET=0
STRIDE=200

SAMPLE_NAME=mcc9_v29e_dl_run3_G1_extbnb_dlana
INPUTFILE=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/va_study/inputlists/mcc9_v29e_dl_run3_G1_extbnb_dlana_MRCNN_INPUTS_LIST.txt
INPUTSTEM=merged_dlana
FILEIDLIST=${RUN_DLANA_DIR}/runlist_vareco_${SAMPLE_NAME}.txt

module load singularity

# CPU MODE
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_vareco_data_cpu.sh $OFFSET $STRIDE $SAMPLE_NAME ${INPUTFILE} ${INPUTSTEM} ${FILEIDLIST}"

