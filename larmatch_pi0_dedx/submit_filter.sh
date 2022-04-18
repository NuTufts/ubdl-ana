#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=filterpi0
#SBATCH --output=filterpi0_bnbnu.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=2:00:00
#SBATCH --array=2-155
##SBATCH --partition=preempt
#SBATCH --partition=batch
##SBATCH --partition=wongjiradlab
#SBATCH --error=gridlog_filterpi0.%j.%N.err

container=/cluster/tufts/wongjiradlabnu/larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
DATA_PREP_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/

module load singularity/3.5.3
cd /cluster/tufts/

# mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge.list: nfiles=15519, njobs=155
singularity exec ${container} bash -c "cd ${DATA_PREP_DIR} && source run_filter_mc.sh"


