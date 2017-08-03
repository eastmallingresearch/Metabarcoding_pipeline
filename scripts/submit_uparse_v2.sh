#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

OUTDIR=$1/data/$2
PREFIX=$3
FILTDIR=$OUTDIR/$PREFIX/filtered
SL=$4
SR=$5

for SCRIPT_DIR; do true; done

#cd $OUTDIR

cd $TMP

#### Concatenate files
cat ${FILTDIR}/*.fa > ${PREFIX}.temp.fa
echo $(ls -gh ${PREFIX}.temp.fa)

#### Remove multiplex primers and pad reads to same length
X=`awk '{if ($1!~/>/){mylen=mylen+length($0)}else{print mylen;mylen=0};}' ${PREFIX}.temp.fa|awk '$0>x{x=$0};END{print x}'`
usearch -fastx_truncate ${PREFIX}.temp.fa -stripleft $SL -stripright $SR -trunclen $X -padlen $X -fastaout ${PREFIX}.fa
#rm ${PREFIX}.temp.fa

#### Dereplication
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${PREFIX}.fa|$SCRIPT_DIR/get_uniq.pl > ${PREFIX}.sorted.fasta 
#rm ${PREFIX}.fa

#### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))

usearch -cluster_otus ${PREFIX}.sorted.fasta -otus ${PREFIX}.otus.fa -relabel OTU -minsize 4

mv ${PREFIX}.otus.fa $OUTDIR/.

usearch -unoise3 ${PREFIX}.sorted.fasta -zotus ${PREFIX}.zotus.fa #-relabel OTU #-minampsize 8

#usearch -unoise ${PREFIX}.sorted.fasta -tabbedout ${PREFIX}.txt -fastaout ${PREFIX}.otus.fa -relabel OTU #-minampsize 8

#perl -pi -e 's/uniq.*/OTU . ++$n/ge' ${PREFIX}.otus.fa

#rm ${PREFIX}.sorted.fasta


mv ${PREFIX}.zotus.fa $OUTDIR/.
