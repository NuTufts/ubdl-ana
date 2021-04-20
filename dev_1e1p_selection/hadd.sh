#!/bin/bash

out=$1
input=$2

thisdir=$PWD

container=/cluster/tufts/wongjiradlab/larbys/larbys-containers/ubdl_depsonly_py3.6.11_u16.04_cu11_pytorch1.7.1.simg
ubdldir=/cluster/tufts/wongjiradlab/twongj01/ubdl_py3
module load singularity
singularity exec $container bash -c "cd $ubdldir && source setenv_py3.sh && source configure.sh && cd $thisdir && hadd $out $input"
