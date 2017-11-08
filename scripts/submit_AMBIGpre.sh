#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4
PRIMERS=$5
ADAPTERS=$6
MINL=$7
MAXDIFF=$8
SCRIPT_DIR=$9

mkdir -p $OUTDIR/joined 
mkdir -p $OUTDIR/unjoined 

cd $TMP


usearch -fastq_mergepairs $F -reverse $R -fastqout ${OUTFILE}.t1  -fastq_pctid 0 -fastq_maxdiffs $(($MINL*${MAXDIFF}/100)) -fastq_minlen $MINL -fastq_minovlen 0 #-fastq_trunctail 25
usearch -search_oligodb ${OUTFILE}.t1 -db $ADAPTERS -strand both -userout ${OUTFILE}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
awk -F"\t" '{print $1}' ${OUTFILE}.t1.txt|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${OUTFILE}.t1 > ${OUTFILE}.t2
mv ${OUTFILE}.t2 $OUTDIR/joined/${OUTFILE}.unfiltered.fastq


(( MINL= MINL - 100 ))
usearch -search_oligodb $F -db $PRIMERS -strand both -userout ${OUTFILE}.F.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
usearch -search_oligodb $R -db $PRIMERS -strand both -userout ${OUTFILE}.R.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 

grep primer1 -v ${OUTFILE}.F.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl $F > ${OUTFILE}.F.t2.fastq
grep primer2 -v ${OUTFILE}.R.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl $R > ${OUTFILE}.R.t2.fastq

mv ${OUTFILE}.F.t2.fastq $OUTDIR/unjoined/${OUTFILE}.r1.unfiltered.fastq
mv ${OUTFILE}.R.t2.fastq $OUTDIR/unjoined/${OUTFILE}.r2.unfiltered.fastq

