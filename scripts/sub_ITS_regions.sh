#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G

MYOTUS=$1;shift
HMM1=$1;shift
HMM2=$1;shift
EVAL=$1;shift
SCRIPT_DIR=$1;shift

cd $TMP

nhmmscan --noali --cpu 8 --incT 6 --tblout ssu  -E $EVAL $HMM1 $MYOTUS >/dev/null
nhmmscan --noali --cpu 8 --incT 6 --tblout 58ss -E $EVAL $HMM2 $MYOTUS >/dev/null

Rscript $SCRIPT_DIR/ITS_regions.R $MYOTUS $TMP/ssu $TMP/58ss

