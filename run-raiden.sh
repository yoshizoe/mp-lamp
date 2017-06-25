#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe impi 8
#$ -q c1normal

# Get variables
# This wacky line is required for whatever reasons.
. /fefs/opt/x86_64/intel/parallel_studio_xe_2017/impi/2017.2.174/bin64/mpivars.sh
#source /home/hal9000/.bashrc
LIBRARY_PATH=$LIBRARY_PATH:/home/hal9000/library/lib
export LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hal9000/library/lib
export LD_LIBRARY_PATH

# TODO: choose which data to test
cd ..
INSTDIR=./samples/cont_data
WDIR=/home/hal9000/workspace/mp-lamp

if [ $method = "mahito" ]
then
#    mpirun -bootstrap sge $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize false
#    mpirun -bootstrap sge $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize false
else
    mpiexec $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize true -r 0.95
#    mpirun -bootstrap sge $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize true -r 0.95
fi


