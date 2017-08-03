#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

SCRIPT_DIR=$1;shift
OUTDIR=$1/data/$2
PREFIX=$3

shift;shift;shift

cd $OUTDIR

#### Assign Taxonomy
usearch9 -utax ${PREFIX}.otus.fa -db $SCRIPT_DIR/../taxonomies/utax/${PREFIX}_ref.udb -strand both -utaxout ${PREFIX}.reads.utax -rdpout ${PREFIX}.rdp -alnout ${PREFIX}.aln.txt $@
usearch9 -utax ${PREFIX}.zotus.fa -db $SCRIPT_DIR/../taxonomies/utax/${PREFIX}_ref.udb -strand both -utaxout z${PREFIX}.reads.utax -rdpout z${PREFIX}.rdp -alnout z${PREFIX}.aln.txt $@

cat ${PREFIX}.rdp|$SCRIPT_DIR/mod_taxa.pl > ${PREFIX}.taxa 
cat z${PREFIX}.rdp|$SCRIPT_DIR/mod_taxa.pl > z${PREFIX}.taxa 
