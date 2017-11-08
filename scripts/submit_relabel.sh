#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

F=$1
OUTFILE=$2
OUTDIR=$3
SCRIPT_DIR=$4

# LABEL=${OUTFILE}.

#for F in *R1.fa;do  
 OUTFILE=$(echo $F|awk -F"_" '{print $1,$2}' OFS="_")
 usearch -fastx_relabel $F -prefix ${OUTFILE}. -fastaout $OUTFILE.t1.fa
 awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  $OUTFILE.t1.fa > ${OUTFILE}.filtered.fa
 sed -i -e '1d' ${OUTFILE}.filtered.fa
 mv ${OUTFILE}.filtered.fa $OUTDIR/.
 rm $OUTFILE.t1.fa
#done
