#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G


DIR=$1
FASTA=${1}.fa
ID=$2
RFILE=$3
REG1=$4
REG2=$5
LOWQUAL=$6

cd $DIR

Rscript $RFILE $DIR $REG1 $REG2 $FASTA $ID $LOWQUAL
