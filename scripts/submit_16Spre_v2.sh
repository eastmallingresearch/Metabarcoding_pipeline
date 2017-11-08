#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1;shift
R=$1;shift
OUTFILE=$1;shift
OUTDIR=$1;shift
ADAPTERS=$1;shift
MINL=$1;shift
MAXDIFF=$1;shift
QUAL=$1;shift
FPL=$1;shift
RPL=$1;shift
SCRIPT_DIR=$1;shift

LEN=$(( MINL - RPL ))

LABEL=${OUTFILE}.

mkdir -p $OUTDIR/filtered 
mkdir -p $OUTDIR/unfiltered 

cd $TMP

#cd $OUTDIR 

usearch -fastq_mergepairs $F -reverse $R -fastqout ${OUTFILE}.t1  -fastq_pctid 0 -fastq_maxdiffs $(($MINL*${MAXDIFF}/100)) -fastq_minlen $MINL -fastq_minovlen 0 #-fastq_trunctail 25

#usearch9 -fastq_mergepairs $F -reverse $R -fastqout ${OUTFILE}.t1  -fastq_maxdiffpct $MAXDIFF -fastq_maxdiffs $(($MINL*${MAXDIFF}/100)) -fastq_minlen $MINL 
usearch -search_oligodb ${OUTFILE}.t1 -db $ADAPTERS -strand both -userout ${OUTFILE}.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 

#cat ${OUTFILE}.t1.txt|awk -F"\t" '{print $1}'|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${OUTFILE}.t1 > ${OUTFILE}.t2
awk -F"\t" '{print $1}' ${OUTFILE}.t1.txt|sort|uniq|$SCRIPT_DIR/adapt_delete.pl ${OUTFILE}.t1 > ${OUTFILE}.t2

awk  -v SL="$FPL" -v SR="$RPL" -F" " '{if(NR % 2 == 0){$1=substr($1,(SL+1),(length($1)-SL-SR))};print $1}' ${OUTFILE}.t2 > ${OUTFILE}.t3

usearch -fastq_filter ${OUTFILE}.t3 -fastq_maxee $QUAL -relabel $LABEL -fastaout ${OUTFILE}.t3.fa

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <${OUTFILE}.t3.fa > ${OUTFILE}.filtered.fa

sed -i -e '1d' ${OUTFILE}.filtered.fa

mv ${OUTFILE}.filtered.fa $OUTDIR/filtered/.
mv ${OUTFILE}.t3 $OUTDIR/unfiltered/${OUTFILE}.unfiltered.fastq

rm ${OUTFILE}.t1.txt ${OUTFILE}.t1 ${OUTFILE}.t3.fa
