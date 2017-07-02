#!/bin/bash

#scons -j 4 debug=0
#cp mp_build/opt/mp-main/cont-lamp .
scons -j 4 debug=10
cp mp_build/dbg/mp-main/cont-lamp .
