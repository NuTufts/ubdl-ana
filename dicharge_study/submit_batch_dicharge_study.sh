#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=dicharge
#SBATCH --output=dicharge_batch.log
#SBATCH --mem-per-cpu=2000
#SBATCH --time=1:00:00
#SBATCH --array=0-9

container=/cluster/tufts/wongjiradlab/twongj01/ubdl/larflow/larmatchnet/grid_deploy_scripts/singularity_pytorch1.3cpu.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/dicharge_study/
OFFSET=0
STRIDE=10

SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge # 15519
DLMERGED_STEM=merged_dlreco

module load singularity
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_dicharge_study.sh $OFFSET $STRIDE $SAMPLE_NAME ${DLMERGED_STEM}"

