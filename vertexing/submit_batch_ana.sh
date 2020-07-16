#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=kpsana_batch_mc
#SBATCH --output=kps_ana_mc_tmw.log
#SBATCH --mem-per-cpu=2000
#SBATCH --time=1:00:00
#SBATCH --array=50-199

container=/cluster/tufts/wongjiradlab/twongj01/ubdl/larflow/larmatchnet/grid_deploy_scripts/singularity_pytorch1.3cpu.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/vertexing/
OFFSET=0
STRIDE=1

SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_intrinsic_nue_LowE # 580


module load singularity
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_kps_ana.sh $OFFSET $STRIDE $SAMPLE_NAME"

