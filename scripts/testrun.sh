#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe impi 2
#$ -q c1normal

LIBRARY_PATH=$LIBRARY_PATH:/home/hal9000/library/lib
export LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hal9000/library/lib
export LD_LIBRARY_PATH


cd ..
. /fefs/opt/x86_64/intel/parallel_studio_xe_2017/impi/2017.2.174/bin64/mpivars.sh

mpirun -bootstrap sge /home/hal9000/workspace/mp-lamp/lamp \
    --item ./samples/cont_data/synth_100_10.data --pos ./samples/cont_data/synth_100_10.class \
    --show_progress true --log true -discretize true -ratio 0.95
#date
