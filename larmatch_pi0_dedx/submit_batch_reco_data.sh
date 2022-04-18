#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=lfreco
#SBATCH --output=larflow_reco_data_bnbnu.log
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1-00:00:00
#SBATCH --array=0-3
#SBATCH --cpus-per-task=1
##SBATCH --partition=batch
##SBATCH --partition=preempt
#SBATCH --partition=wongjiradlab
##SBATCH --gres=gpu:p100:3
##SBATCH --partition ccgpu
##SBATCH --nodelist=ccgpu01
#SBATCH --error=gridlog/gridlog_larflowreco_pi0select.%j.%N.err

container=/cluster/tufts/wongjiradlabnu//larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif
SCRIPTS_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/ana/larmatch_pi0_dedx/

module load singularity/3.5.3

# CPU MODE
singularity exec --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab ${container} bash -c "cd ${SCRIPTS_DIR} && source run_batch_larflowreco_data_cpu.sh"

