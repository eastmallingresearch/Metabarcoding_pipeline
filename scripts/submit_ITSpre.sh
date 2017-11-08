#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1;shift
R=$1;shift
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
usearch -search_oligodb $R -db $PRIMERS -strand both -userout ${R}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 

grep primer1 -v ${F}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl $F > ${F}.t2.fastq
grep primer2 -v ${R}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl $R > ${R}.t2.fastq

usearch -fastq_filter ${F}.t2.fastq -fastq_trunclen $(( MINL - $FPL ))  -fastq_stripleft $FPL -fastq_maxee $QUAL -fastaout ${OUTFILE}_t1.fa  
usearch -fastq_filter ${R}.t2.fastq -fastq_trunclen $(( MINL - $RPL )) -fastq_stripleft $RPL -fastq_maxee $QUAL -fastaout ${OUTFILE}_t2.fa  

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


