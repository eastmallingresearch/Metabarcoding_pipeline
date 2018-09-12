#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

dir=$1
OUTDIR=$2
PREFIX=$3
VER=$4

for SCRIPT_DIR; do true; done

if  ! (( $VER )) ; then VER=""; fi

Rscript $SCRIPT_DIR/merge_otu_tables.R $dir ".*otu_table.*" $OUTDIR  ${PREFIX}${VER}

