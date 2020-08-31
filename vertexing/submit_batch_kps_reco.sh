#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=kpsreco
#SBATCH --output=kpsreco_tmw.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1-00:00:00
#SBATCH --array=0-9

container=/cluster/tufts/wongjiradlab/twongj01/ubdl/larflow/larmatchnet/grid_deploy_scripts/singularity_pytorch1.3cpu.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/vertexing/
OFFSET=0
STRIDE=10

#SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_intrinsic_nue_LowE # 580 files
#SAMPLE_NAME=mcc9_v28_wctagger_extbnb # 19864 files
#SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge # 15519 # missing 29,84
SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge # 3331

module load singularity
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_kps_reco.sh $OFFSET $STRIDE $SAMPLE_NAME"

