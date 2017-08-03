#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

DATA=$(sed -n -e "$SGE_TASK_ID p" split_files.txt)
HMM1=$1
HMM2=$2
EVAL=$5

echo "input file: " $DATA
echo "HHM1: " $HMM1
echo "HMM2: " $HMM2
echo "eval: " $EVAL

echo

echo "writing HMM1 to: " $DATA.$3
nhmmscan --noali --cpu 8 --incT 6 --tblout $DATA.$3 -E $EVAL $HMM1 $DATA >/dev/null
echo "writing HMM2 to: " $DATA.$4
nhmmscan --noali --cpu 8 --incT 6 --tblout $DATA.$4 -E $EVAL $HMM2 $DATA >/dev/null
echo "complete"