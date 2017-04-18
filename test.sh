#!/bin/bash

./make.sh
./run.sh > results.txt

haspattern=`cat results.txt | grep "# number of significant patterns=1" | wc -l`

if [ $haspattern = "1" ]
then
    echo "TEST PASSED"
    git-cola
else
    cat results.txt
    echo "!!TEST FAILED!!"
fi
