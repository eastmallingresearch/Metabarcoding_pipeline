#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G
#$ -pe smp 8

# SCRIPT_DIR=$1;shift
OUTDIR=$1/data/$2
PREFIX=$3

cd $TMP

usearch -calc_distmx ${OUTDIR}/${PREFIX}.otus.fa -distmxout ${PREFIX}.phy -distmo fractdiff -format phylip_square -threads 8
#usearch -calc_distmx ${OUTDIR}/${PREFIX}.zotus.fa -distmxout ${PREFIX}.z.phy -distmo fractdiff -format phylip_square -threads 8

cp * $OUTDIR/.
