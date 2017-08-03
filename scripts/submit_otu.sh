#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

OPTIND=1 

OUTDIR=$1/data/$2
PREFIX=$3
UNFILTDIR=$OUTDIR/$PREFIX/unfiltered
SL=$4
SR=$5

for SCRIPT_DIR; do true; done

while getopts ":hc:" options; do
	case "$options" in
	c)  
 	    program=$OPTARG
	    ;;
	h)  
	    exit 1
 	    ;;
	\?) 
	    echo "Invalid option: -$OPTARG" >&2
	    exit 1
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    exit 1
      	    ;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

case $program in
submit_fastq_fasta)
	INFILE=$(sed -n -e "$SGE_TASK_ID p" $1)
	dir=$2
	SL=$3
	SR=$4
	FILE=$(echo $INFILE|awk -F"/" '{print $NF}')
	PREFIX=$(echo $FILE|awk -F"." '{print $1}')

	for SCRIPT_DIR; do true; done

	if [[ "$FILE" =~ "r1" ]]; then
		$SCRIPT_DIR/fq2fa_v2.pl $INFILE $dir/$PREFIX.t1.fa $PREFIX $SL 0
	elif [[ "$FILE" =~ "r2" ]]; then 
		$SCRIPT_DIR/fq2fa_v2.pl $INFILE $dir/$PREFIX.t2.fa $PREFIX 0 $SR
	else
		$SCRIPT_DIR/fq2fa_v2.pl $INFILE $dir/$PREFIX.t1.fa $PREFIX $SL $SR
	fi
	
	exit 1
	;;
submit_cat_files)
	dir=$1
	cat $dir/*.t1.fa > $dir/t1.fa
	cat $dir/*.t2.fa > $dir/t2.fa
	exit 1
	;;
submit_global_search)
	INFILE=$1
	OUTDIR=$2
	PREFIX=$3
	EP=$4

	if [ -z $EP ]; then
		usearch -usearch_global $INFILE -db $OUTDIR/$PREFIX.otus.fa -strand plus -id 0.97 -biomout $OUTDIR/$PREFIX.otu_table.biom -otutabout $OUTDIR/$PREFIX.otu_table.txt -output_no_hits -userout $OUTDIR/$PREFIX.hits.out -userfields query+target
	else
		usearch -usearch_global $INFILE -db $OUTDIR/$PREFIX.otus.fa -strand both -id 0.97 -biomout $OUTDIR/$PREFIX$EP.otu_table.biom -otutabout $OUTDIR/$PREFIX$EP.otu_table.txt -output_no_hits -userout $OUTDIR/$PREFIX$EP.hits.out -userfields query+target
	fi
	
	exit 1
	;;
submit_search_hits)
	INFILE=$(sed -n -e "$SGE_TASK_ID p" $1)
	OUTDIR=$2
	HITS=$3
	S=$(echo $INFILE|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')

	for SCRIPT_DIR; do true; done

	grep "$S.*\*" $HITS|awk -F";" '{print $2}'|awk -F" " '{print $1}'|$SCRIPT_DIR/seq_select_v2.pl $OUTDIR/t2.fa >> $OUTDIR/t3.fa

	exit 1
	;;
submit_tidy)
	for thing in "$@";do
		echo removing $thing
		rm -r $thing
	done
	exit 1
	;;
*)
	echo "Invalid options: $program" >&2
	exit 1
esac
