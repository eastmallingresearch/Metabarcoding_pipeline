#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

DATA=$1
OUT=$2

usearch -unoise $DATA -tabbedout ${OUT}.txt -fastaout ${OUT}.fa
perl -pi -e 's/uniq.*/OTU . ++$n/ge' ${OUT}.fa


