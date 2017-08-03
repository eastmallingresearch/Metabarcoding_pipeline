#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

dir=$1
OUTDIR=$2
PREFIX=$3

for SCRIPT_DIR; do true; done

Rscript $SCRIPT_DIR/merge_biom.R $OUTDIR $dir ".*otu_table.*"  $PREFIX

