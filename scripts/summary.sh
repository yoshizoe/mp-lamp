#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -m e
#$ -M yuu.jinnai@riken.jp
PATH=$PATH:/home/hal9000/library/bin
export PATH
LIBRARY_PATH=$LIBRARY_PATH:/home/hal9000/library/lib
export LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hal9000/library/lib
export LD_LIBRARY_PATH

#fresults=`echo $results | awk -F '__' '{print}'`
echo "results=$results"
#echo "fresults=$fresults"
# TODO: Add a method to generate tags.
# TODO: Add a method to add new tags when necessary.

if [ ! -f $summary ]
then
    for result in $results
    do
	cat $result | awk -f ./parse-header.awk > $summary   
	break
    done
fi

for result in $results
do
    echo $result
#    cat $result
    cat $result | awk -f ./parse.awk >> $summary
done


Rscript --vanilla plot_summary.R ${summary}
#cd ../../scripts

d=`date`
#git stage ${summary} ${summary}.pdf ../src/R/.RData
git add -u
git commit -m "$d pdf autocommitted: $summary"
#git push
