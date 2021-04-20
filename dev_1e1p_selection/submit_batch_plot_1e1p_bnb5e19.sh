#!/bin/bash

# slurm submission script for making larmatch training data

#SBATCH --job-name=plotLEEdev
#SBATCH --output=plot_LEE_dev_sub0.txt
#SBATCH --mem-per-cpu=2000
#SBATCH --time=1:00:00
#SBATCH --array=0-99
#SBATCH --cpus-per-task=1
##SBATCH --partition=batch
#SBATCH --partition=preempt
##SBATCH --partition=wongjiradlab
##SBATCH --gres=gpu:p100:3
##SBATCH --partition ccgpu
##SBATCH --gres=gpu:t4:1
##SBATCH --nodelist=ccgpu01

container=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
RUN_DLANA_DIR=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/dev_1e1p_selection/
OFFSET=0
STRIDE=117

SAMPLE_NAME=mcc9_v28_wctagger_bnb5e19
INPUTFILE=/cluster/tufts/wongjiradlab/twongj01/ubdl-ana/dev_1e1p_selection/bnb5e19.txt
DLMERGED_LIST=/cluster/tufts/wongjiradlab//nutufts/dlgen2prod/run1inputlists/mcc9_v28_wctagger_bnb5e19_filelist.txt
module load singularity

# CPU MODE
srun singularity exec ${container} bash -c "cd ${RUN_DLANA_DIR} && source run_batch_1e1p_plots.sh $OFFSET $STRIDE $SAMPLE_NAME ${INPUTFILE} 1 0 ${DLMERGED_LIST}"

