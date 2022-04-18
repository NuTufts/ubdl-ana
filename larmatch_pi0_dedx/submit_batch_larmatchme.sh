#!/bin/bash

# slurm submission script for making larmatch training data
#SBATCH --job-name=larmatchme
#SBATCH --output=larmatchme_mcc9_run3_bnbnu_sub0.log
#SBATCH --mem-per-cpu=6000
#SBATCH --time=8:00:00
#SBATCH --array=0-3
#SBATCH --cpus-per-task=4
##SBATCH --partition=batch
#SBATCH --partition=wongjiradlab
##SBATCH --partition=preempt
##SBATCH --gres=gpu:p100:3
##SBATCH --partition ccgpu
##SBATCH --gres=gpu:a100:1
##SBATCH --nodelist=ccgpu01
#SBATCH --error=gridlog/griderr_larmatcheme_deploy_mcc9_run3_bnbnu_sub0.%j.%N.err

container=/cluster/tufts/wongjiradlabnu//larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
SCRIPTS_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/

module load singularity/3.5.3
# GPU MODE
#srun singularity exec --nv ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_kps_larmatch.sh $OFFSET $STRIDE $SAMPLE_NAME ${INPUTFILE} ${INPUTSTEM} ${FILEIDLIST}"
# CPU MODE
cd /cluster/tufts/
singularity exec ${container} bash -c "cd ${SCRIPTS_DIR} && source run_batch_larmatchminkowski.sh"

