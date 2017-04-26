#!/bin/bash

if [ -f "made.txt" ]
then
    rm made.txt
fi
./make.sh 2>&1 | tee -a made.txt

successfullymade=`cat made.txt | grep "error" | wc -l`

if [ $successfullymade != "0" ]
then
    echo "build unsuccessful"
    exit 0
fi

if [ -f "results.txt" ]
then
    rm results.txt
fi

./run.sh 2>&1 | tee -a results.txt

haspattern=`cat results.txt | grep "# number of significant patterns=1" | wc -l`

if [ $haspattern = "1" ]
then
    echo "TEST PASSED"
    git-cola
else
    cat results.txt
    echo "!!TEST FAILED!!"
fi
