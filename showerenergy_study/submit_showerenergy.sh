#!/bin/bash

#SBATCH --job-name=showerenergy
#SBATCH --output=log-showerenergy
##SBATCH --partition batch
#SBATCH --partition=wongjiradlab
#SBATCH --time=0-2:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-1999


ANADIR=/cluster/tufts/wongjiradlab/kmason03/ubdl-ana/showerenergy_study
UBDLDIR=/cluster/tufts/wongjiradlab/kmason03/ubdl
CONTAINER=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
FILELIST=/cluster/tufts/wongjiradlab/kmason03/ubdl-ana/showerenergy_study/runlist_dlgen2filter_mcc9_v29e_dl_run3b_bnb_intrinsic_nue_overlay_nocrtremerge.txt

module load singularity
singularity exec ${CONTAINER} bash -c "cd ${ANADIR} && source run_showerenergy.sh ${UBDLDIR} ${ANADIR} ${FILELIST} "
