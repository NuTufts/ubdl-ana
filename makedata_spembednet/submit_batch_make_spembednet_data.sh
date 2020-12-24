#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=sp3demb
#SBATCH --output=kps_batch_bnbnu_rerun.log
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --array=0

container=/cluster/tufts/wongjiradlab/twongj01/ubdl/larflow/larmatchnet/grid_deploy_scripts/singularity_pytorch1.3cpu.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/makedata_spembednet/
OFFSET=0
STRIDE=10

#SAMPLE_NAME=mcc9_v13_bnbnue_corsika_training
SAMPLE_NAME=mcc9_v13_bnbnue_corsika_valid


module load singularity
#srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_kps_larmatch.sh $OFFSET $STRIDE $SAMPLE_NAME"
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_truehit_training_voxels.sh $OFFSET $STRIDE $SAMPLE_NAME"

