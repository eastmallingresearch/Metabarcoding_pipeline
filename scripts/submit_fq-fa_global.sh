#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

INFILE=$(sed -n -e "$SGE_TASK_ID p" $1)
dir=$2
SL=$3
SR=$4
DB=$5
FILE=$(echo $INFILE|awk -F"/" '{print $NF}')
PREFIX=$(echo $FILE|awk -F"." '{print $1}')

for SCRIPT_DIR; do true; done

if [[ "$FILE" =~ "r1" ]]; then
	awk -v S="$PREFIX" -v SL="$SL" -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){$1=substr($1,(SL+1));print $1}}' $INFILE >$dir/$PREFIX.t1.fa
elif [[ "$FILE" =~ "r2" ]]; then 
	awk -v S="$PREFIX" -v SR="$SR" -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){$1=substr($1,(SR+1));print $1}}' $INFILE >$dir/$PREFIX.t2.fa
else
	awk -v S="$PREFIX" -v SL="$SL" -v SR="$SR" -F" " '{if(NR % 4 == 1){print ">" S "." count+1 ";"$1;count=count+1} if(NR % 4 == 2){$1=substr($1,(SL+1),(length($1)-SL-SR));print $1}}' $INFILE >$dir/$PREFIX.t1.fa
fi

usearch -usearch_global $dir/$PREFIX.t1.fa -db $DB -strand plus -id 0.97 -otutabout $dir/${PREFIX}.otu_table.txt