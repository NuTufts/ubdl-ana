#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=kprecoana
#SBATCH --output=kprecoana.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1:00:00
#SBATCH --array=0
#SBATCH --cpus-per-task=1
##SBATCH --partition=batch
#SBATCH --partition=preempt
##SBATCH --partition=wongjiradlab
##SBATCH --gres=gpu:p100:3
##SBATCH --partition ccgpu
##SBATCH --gres=gpu:t4:1
##SBATCH --nodelist=ccgpu01
#SBATCH --error=gridlog/gridlog_kprecoana_extbnb.%j.%N.err

container=/cluster/tufts/wongjiradlabnu//larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
RUN_DIR="/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_me_kpreco/"
module load singularity/3.5.3

# CPU MODE
#srun singularity exec --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab ${container} bash -c "cd ${RUN_DIR} && source run_kprecoana.sh"
#srun singularity exec --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab ${container} bash -c "cd ${RUN_DIR} && source run_kprecoana_run3_bnb_nue.sh"
srun singularity exec --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab ${container} bash -c "cd ${RUN_DIR} && source run_kprecoana_run3_G1_extbnb.sh"

