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
FPL=$6
RPL=$7

cd $DIR

Rscript $RFILE $DIR $REG1 $REG2 $FASTA $ID $FPL $RPL
