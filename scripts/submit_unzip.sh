#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FILE=$1
shift;
TYPE=${1:-1}

for SCRIPT_DIR; do true; done

if [ $TYPE == 1 ]; then
	pigz -d $FILE
else
	OUTFILE=$(echo $FILE|awk -F"." '{NF--;print}' OFS=".")
	zcat $FILE > $OUTFILE
fi

