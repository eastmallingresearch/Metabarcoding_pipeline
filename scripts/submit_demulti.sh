#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

R1=$1
shift
R2=$1
shift

for SCRIPT_DIR; do true; done

cd $TMP

F=$(echo $R1|awk -F"/" '{print $NF}')
R=$(echo $R2|awk -F"/" '{print $NF}')

OUTDIR=$(echo $R1|sed "s/$F//")


mkfifo $F.de
mkfifo $R.de

zcat -f -- $R1 > $F.de &
zcat -f -- $R2 > $R.de &

${SCRIPT_DIR}/demulti_v2.pl $F.de $R.de $@

rm $F.de $R.de

cp *.fastq $OUTDIR/.
