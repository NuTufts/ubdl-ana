#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=dedxana
#SBATCH --output=dedxana_tmw.log
#SBATCH --mem-per-cpu=2000
#SBATCH --time=30:00
#SBATCH --array=0-25

container=/cluster/tufts/wongjiradlab/twongj01/ubdl/larflow/larmatchnet/grid_deploy_scripts/singularity_pytorch1.3cpu.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/larmatch_1m1p_ana/
OFFSET=0
STRIDE=100

SAMPLE_NAME=mcc9_v29e_dl_run3_G1_bnb_dlfilter_1m1p # 2503

module load singularity
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_dedxana.sh $OFFSET $STRIDE $SAMPLE_NAME ${RUN_DLANA_DIR}"

