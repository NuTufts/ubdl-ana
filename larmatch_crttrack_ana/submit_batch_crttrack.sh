#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=crttrack
#SBATCH --output=crttrack_tmw.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1:00:00
#SBATCH --array=0-99

container=/cluster/tufts/wongjiradlab/twongj01/ubdl/larflow/larmatchnet/grid_deploy_scripts/singularity_pytorch1.3cpu.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/larmatch_crttrack_ana/
OFFSET=0
STRIDE=1

SAMPLE_NAME=mcc9_v29e_dl_run3_G1_extbnb # 34254

module load singularity
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_crttrack.sh $OFFSET $STRIDE $SAMPLE_NAME ${RUN_DLANA_DIR}"

