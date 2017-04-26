#!/bin/bash

scons -j 4 debug=1
cp mp_build/dbg/mp-main/cont-lamp .
cp mp_build/dbg/mp-main/bin-lamp .
