#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

INFILE=$1
OUTDIR=$2
PREFIX=$3
OTU=$4
EP=${5:-plus}


#if [ -z $EP ]; then
usearch -otutab $INFILE -db $OUTDIR/$PREFIX.${OTU}.fa -strand $EP -id 0.97 -biomout $OUTDIR/$PREFIX.${OTU}_table.biom -otutabout $OUTDIR/$PREFIX.${OTU}_table.txt -notmatched $OUTDIR/$PREFIX.${OTU}_nomatch.fa -userout $OUTDIR/$PREFIX.${OTU}_hits.out -userfields query+target
#usearch -otutab $INFILE -db $OUTDIR/$PREFIX.zotus.fa -strand $EP -id 0.97 -biomout $OUTDIR/$PREFIX.zotu_table.biom -otutabout $OUTDIR/$PREFIX.zotu_table.txt -notmatched $OUTDIR/$PREFIX.znomatch.fa -userout $OUTDIR/$PREFIX.zhits.out -userfields query+target
#	usearch -usearch_global $INFILE -db $OUTDIR/$PREFIX.otus.fa -strand plus -id 0.97 -biomout $OUTDIR/$PREFIX.otu_table.biom -otutabout $OUTDIR/$PREFIX.otu_table.txt -output_no_hits -userout $OUTDIR/$PREFIX.hits.out -userfields query+target
#else
#	usearch -otutab $INFILE -db $OUTDIR/$PREFIX.otus.fa -strand both -id 0.97 -biomout $OUTDIR/$PREFIX$EP.otu_table.biom -otutabout $OUTDIR/$PREFIX$EP.otu_table.txt -notmatched $OUTDIR/$PREFIX$EP.nomatch.fa -userout $OUTDIR/$PREFIX$EP.hits.out -userfields query+target
#	usearch -otutab $INFILE -db $OUTDIR/$PREFIX.zotus.fa -strand both -id 0.97 -biomout $OUTDIR/$PREFIX$EP.zotu_table.biom -otutabout $OUTDIR/$PREFIX$EP.zotu_table.txt -notmatched $OUTDIR/$PREFIX$EP.znomatch.fa -userout $OUTDIR/$PREFIX$EP.zhits.out -userfields query+target
#	usearch -usearch_global $INFILE -db $OUTDIR/$PREFIX.otus.fa -strand both -id 0.97 -biomout $OUTDIR/$PREFIX$EP.otu_table.biom -otutabout $OUTDIR/$PREFIX$EP.otu_table.txt -output_no_hits -userout $OUTDIR/$PREFIX$EP.hits.out -userfields query+target
#fi