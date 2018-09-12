#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1;shift
OUTFILE=$1;shift
OUTDIR=$1;shift
PRIMERS=$1;shift
MINL=$1;shift
QUAL=$1;shift
FPL=$1;shift
RPL=$1;shift
SCRIPT_DIR=$1;shift

mkdir -p $OUTDIR/filtered 
mkdir -p $OUTDIR/unfiltered 

cd $TMP

usearch -search_oligodb $F -db $PRIMERS -strand both -userout ${F}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 

grep primer1 -v ${F}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl $F > ${F}.t2.fastq


mv ${F}.t2.fastq $OUTDIR/unfiltered/${OUTFILE}.r1.unfiltered.fastq



