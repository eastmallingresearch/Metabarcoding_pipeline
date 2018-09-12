#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

F=$1;shift
R=$1;shift
OUTFILE=$1;shift
OUTDIR=$1;shift
MINL=$1;shift
MAXDIFF=$1;shift


cd $TMP

usearch -fastq_mergepairs $F -reverse $R -fastqout ${OUTFILE}.t1  -fastq_pctid 0 -fastq_maxdiffs $(($MINL*${MAXDIFF}/100)) -fastq_minlen $MINL -fastq_minovlen 0 #-fastq_trunctail 25

mv ${OUTFILE}.t1 $OUTDIR/${OUTFILE}.unfiltered.fastq