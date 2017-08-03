#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1
R=$2
OUTFILE=$3
OUTDIR=$4
PRIMERS=$5
MINL=$6
MAXR2=$7
QUAL=$8
SCRIPT_DIR=$9

mkdir -p $OUTDIR/filtered 
mkdir -p $OUTDIR/unfiltered 

cd $TMP

usearch -search_oligodb $F -db $PRIMERS -strand both -userout ${F}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 
usearch -search_oligodb $R -db $PRIMERS -strand both -userout ${R}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 

grep primer1 -v ${F}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl $F > ${F}.t2.fastq
grep primer2 -v ${R}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl $R > ${R}.t2.fastq

usearch -fastq_filter ${F}.t2.fastq -fastq_trunclen $MAXR2 -fastq_minlen $MINL -fastq_maxee $QUAL -fastaout ${OUTFILE}_t1.fa  #-fastq_truncqual 26
usearch -fastq_filter ${R}.t2.fastq -fastq_minlen $MINL -fastq_trunclen $MAXR2 -fastq_maxee $QUAL -fastaout ${OUTFILE}_t2.fa #-fastq_truncqual 26

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <${OUTFILE}_t1.fa > ${OUTFILE}_R1.fa
sed -i -e '1d' ${OUTFILE}_R1.fa
sed -i -e 's/ .*//' ${OUTFILE}_R1.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <${OUTFILE}_t2.fa > ${OUTFILE}_R2.fa
sed -i -e '1d' ${OUTFILE}_R2.fa
sed -i -e 's/ .*//' ${OUTFILE}_R2.fa

mv ${OUTFILE}_R1.fa $OUTDIR/fasta/.
mv ${OUTFILE}_R2.fa $OUTDIR/fasta/.
mv ${F}.t2.fastq $OUTDIR/unfiltered/${OUTFILE}.r1.unfiltered.fastq
mv ${R}.t2.fastq $OUTDIR/unfiltered/${OUTFILE}.r2.unfiltered.fastq

rm ${F}.t1.txt ${R}.t1.txt ${OUTFILE}_t1.fa ${OUTFILE}_t2.fa


