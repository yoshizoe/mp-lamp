#!/bin/bash

./make.sh
mpiexec -np 2 ./bin-lamp --item samples/sample_data/sample_item.csv --pos samples/sample_data/sample_expression_over1.csv --show_progress true --log true > results.txt

haspattern=`cat results.txt | grep "# number of significant patterns=1" | wc -l`

if [ $haspattern = "1" ]
then
    echo "TEST PASSED"
else
    cat results.txt
    echo "!!TEST FAILED!!"
fi
