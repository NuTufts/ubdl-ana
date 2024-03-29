#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=dedxana
#SBATCH --output=dedxana_tmw.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1:00:00
#SBATCH --array=0-99

container=/cluster/tufts/wongjiradlab/twongj01/ubdl/larflow/larmatchnet/grid_deploy_scripts/singularity_pytorch1.3cpu.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/larmatch_dedx_ana/
OFFSET=0
STRIDE=10

#SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge # 1551
#SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge
SAMPLE_NAME=mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge_goodlist

module load singularity
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_dedxana_showers.sh $OFFSET $STRIDE $SAMPLE_NAME ${RUN_DLANA_DIR}"

