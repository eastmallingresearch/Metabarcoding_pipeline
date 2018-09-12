#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

INFILE=$1
OUTDIR=$2
PREFIX=$3
OTU=$4
EP=${5:-plus}

usearch -otutab $INFILE -db $OUTDIR/$PREFIX.${OTU}.fa -sample_delim . -strand $EP -id 0.97 -biomout $OUTDIR/$PREFIX.${OTU}_table.biom -otutabout $OUTDIR/$PREFIX.${OTU}_table.txt -notmatched $OUTDIR/$PREFIX.${OTU}_nomatch.fa -userout $OUTDIR/$PREFIX.${OTU}_hits.out -userfields query+target

