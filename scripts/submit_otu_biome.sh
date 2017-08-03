#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

dir=$1
OUTDIR=$2
PREFIX=$3

for SCRIPT_DIR; do true; done

$SCRIPT_DIR/otu_to_biom.pl $dir/row_biom $dir/col_biom $dir/data_biom >${OUTDIR}/${PREFIX}.otu_table.biom
$SCRIPT_DIR/biom_maker.pl ${OUTDIR}/${PREFIX}.taxa ${PREFIX}.otu_table.biom >${OUTDIR}/${PREFIX}.taxa.biom

