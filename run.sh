#!/bin/bash

#cd ..
INSTDIR=./samples/cont_data
WDIR=/home/yuu/workspace/mp-lamp

if [ $method == "mahito" ]
then
#    mpirun -bootstrap sge $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize false
    mpiexec $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize false
else
    mpiexec $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize true -ratio 0.95
#    mpirun -bootstrap sge $WDIR/cont-lamp --item ${INSTDIR}/${instance}.data --pos ${INSTDIR}/${instance}.class  --show_progress true --log true -discretize true -r 0.95
fi


