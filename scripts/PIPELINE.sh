#!/bin/bash

OPTIND=1 

SCRIPT_DIR=$(readlink -f ${0%/*})
S=$SCRIPT_DIR


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
demulti|demultiplex)
	LOC=$1
	shift
 	for f in $LOC
	do
		R1=$f     
		R2=$(echo $R1|sed 's/_R1_/_R2_/')    
		qsub $SCRIPT_DIR/submit_demulti.sh $R1 $R2 $@ $SCRIPT_DIR
	done
	exit 1
	;;
16Spre|OOpre)
	LOC=$1; shift
	for f in $LOC
	do
		R1=$f
   		R2=$(echo $R1|sed 's/_R1_/_R2_/')
		#S=$(echo $f|awk -F"_" -v D=$(echo $LOC|awk -F"/" '{print($(NF-3))}') '{print $2"D"D}')	
		S=$(echo $f|awk -F"/" '{print $NF'}|awk -F"_" {'print $1,$2'} OFS="_")	
		qsub $SCRIPT_DIR/submit_16Spre_v2.sh $R1 $R2 $S $@ $SCRIPT_DIR
	done
	exit 1
	;;
ITSpre)
	LOC=$1
	shift
	for f in $LOC
	do     
		R1=$f;     
		R2=$(echo $R1|sed 's/_R1_/_R2_/');
		#S=$(echo $f|awk -F"_" -v D=$(echo $LOC|awk -F"/" '{print($(NF-3))}') '{print $2"D"D}');
		S=$(echo $f|awk -F"/" '{print $NF'}|awk -F"_" {'print $1,$2'} OFS="_")	
		qsub $SCRIPT_DIR/submit_ITSpre.sh $R1 $R2 $S $@ $SCRIPT_DIR
	done
	exit 1
	;;
NEMpre)
	LOC=$1;	shift
	for f in $LOC
	do     
		R1=$f;     
		R2=$(echo $R1|sed 's/_R1_/_R2_/');
		S=$(echo $f|awk -F"/" '{print $NF'}|awk -F"_" {'print $1,$2'} OFS="_")	
		qsub $SCRIPT_DIR/submit_NEMpre.sh $R1 $R2 $S $@ $SCRIPT_DIR
	done
	exit 1
	;;
AMBIGpre)
	LOC=$1
	shift
	for f in $LOC
	do     
		R1=$f;     
		R2=$(echo $R1|sed 's/_R1/_R2/');
		#S=$(echo $f|awk -F"_" -v D=$(echo $LOC|awk -F"/" '{print($(NF-3))}') '{print $2"D"D}');
		S=$(echo $f|awk -F"/" '{print $NF'}|awk -F"_" {'print $1,$2'} OFS="_")	
		qsub $SCRIPT_DIR/submit_AMBIGpre.sh $R1 $R2 $S $@ $SCRIPT_DIR
	done
	exit 1
	;;

procends)
	LOC=$1
	shift	
	READ=$1
	shift

	cd $LOC
	for f in *${READ}.fa
	do
		d=$(echo $f|awk -F"." '{print $1}')
		mkdir $d
		split -l 2000 $f -a 3 -d ${d}/$f.
		cd $d
		find $PWD -name '*.fa.*' >split_files.txt
		TASKS=$(wc -l split_files.txt|awk -F" " '{print $1}')
        		qsub -t 1-$TASKS:1 $SCRIPT_DIR/submit_nscan.sh $@
		cd ..
	done
	exit 1
	;;
ITS)
	LOC=$1
	shift
	for d in $LOC
	do
		S=$(echo $d|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}')
		qsub $SCRIPT_DIR/submit_ITS.sh $d $S $@
	done
	exit 1
	;;
ITS_regions)
	qsub $SCRIPT_DIR/sub_ITS_regions.sh $@ $SCRIPT_DIR
	exit 1
	;;

merge_hits)
	qsub $SCRIPT_DIR/submit_merge_hits.sh $@ $SCRIPT_DIR
	exit 1
	;;
UPARSE|uparse)
	qsub $SCRIPT_DIR/submit_uparse_v2.sh $@ $SCRIPT_DIR
	exit 1
	;;
UCLUS|uclus)
	qsub $SCRIPT_DIR/submit_uparse.sh $@ $SCRIPT_DIR
	exit 1
	;;
OTU2)
	OUTDIR=$1/data/$2
	PREFIX=$3
	UNFILTDIR=$OUTDIR/$PREFIX/unfiltered
	SL=$4
	SR=$5
	R2=${6:-false}
	JOBNAME=OTU_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

	cd $OUTDIR
	dir=`mktemp -d -p $OUTDIR`
	find $UNFILTDIR -name '*.fastq' >$dir/files.txt
	TASKS=$(wc -l $dir/files.txt|awk -F" " '{print $1}')
	qsub -N ${JOBNAME}_1 -t 1-$TASKS:1 $SCRIPT_DIR/submit_otu.sh -c submit_fastq_fasta $dir/files.txt $dir $SL $SR $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_1 -N ${JOBNAME}_2 $SCRIPT_DIR/submit_otu.sh -c submit_cat_files $dir
	qsub -hold_jid ${JOBNAME}_2 -N ${JOBNAME}_3 $SCRIPT_DIR/submit_otu.sh -c submit_global_search $dir/t1.fa $OUTDIR $PREFIX
	
	if $R2; then
		find $UNFILTDIR -name '*.r2.*' >$dir/R2.files.txt			
		TASKS=$(wc -l $dir/R2.files.txt|awk -F" " '{print $1}')
		qsub -hold_jid ${JOBNAME}_3 -N ${JOBNAME}_4 -t 1-$TASKS:1 $SCRIPT_DIR/submit_otu.sh -c submit_search_hits $dir/files.txt $dir $OUTDIR/$PREFIX.hits.out $SCRIPT_DIR
		qsub -hold_jid ${JOBNAME}_4 -N ${JOBNAME}_5 $SCRIPT_DIR/submit_otu.sh -c submit_global_search $dir/t3.fa $OUTDIR $PREFIX 2
		qsub -hold_jid ${JOBNAME}_5 $SCRIPT_DIR/submit_otu.sh -c submit_tidy $dir $PREFIX.hits.out ${PREFIX}2.hits.out OTU_*_1.* OTU_*_4.*
	else
		qsub -hold_jid ${JOBNAME}_3 $SCRIPT_DIR/submit_otu.sh -c submit_tidy $dir $PREFIX.hits.out OTU_*_1.* 
	fi  

	exit 1
	;;
OTU|otu)
	OUTDIR=$1/data/$2
	PREFIX=$3
	UNFILTDIR=$OUTDIR/$PREFIX/unfiltered
	SL=$4
	SR=$5
	R2=${6:-false}

	if $R2; then
		EP=both
	else
		EP=plus
	fi

	JOBNAME=OTU_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

	cd $OUTDIR
	dir=`mktemp -d -p $OUTDIR`
	find $UNFILTDIR -name '*.fastq' >$dir/files.txt
	TASKS=$(wc -l $dir/files.txt|awk -F" " '{print $1}')
	qsub -N ${JOBNAME}_1 -t 1-$TASKS:1 $SCRIPT_DIR/submit_fastq_fasta.sh $dir/files.txt $dir $SL $SR $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_1 -N ${JOBNAME}_2 $SCRIPT_DIR/submit_cat_files.sh $dir $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_2 -N ${JOBNAME}_3 $SCRIPT_DIR/submit_global_search.sh $dir/t1.fa $OUTDIR $PREFIX otus plus
	qsub -hold_jid ${JOBNAME}_2 -N ${JOBNAME}_4 $SCRIPT_DIR/submit_global_search.sh $dir/t1.fa $OUTDIR $PREFIX zotus plus
	
	if $R2; then
		find $UNFILTDIR -name '*.r2.*' >$dir/R2.files.txt			
		TASKS=$(wc -l $dir/R2.files.txt|awk -F" " '{print $1}')
		qsub -hold_jid ${JOBNAME}_3 -N ${JOBNAME}_4 -t 1-$TASKS:1 $SCRIPT_DIR/submit_search_hits.sh $dir/files.txt $dir $OUTDIR/$PREFIX.hits.out $SCRIPT_DIR
		qsub -hold_jid ${JOBNAME}_4 -N ${JOBNAME}_5 $SCRIPT_DIR/submit_global_search.sh $dir/t3.fa $OUTDIR $PREFIX otus both
		qsub -hold_jid ${JOBNAME}_4 -N ${JOBNAME}_6 $SCRIPT_DIR/submit_global_search.sh $dir/t3.fa $OUTDIR $PREFIX zotus both
		qsub -hold_jid ${JOBNAME}_5,${JOBNAME}_6 $SCRIPT_DIR/submit_tidy.sh $dir $PREFIX.hits.out ${PREFIX}2.hits.out OTU_*_1.* OTU_*_4.*
	else
		qsub -hold_jid ${JOBNAME}_3,${JOBNAME}_4 $SCRIPT_DIR/submit_tidy.sh $dir $PREFIX.hits.out OTU_*_1.*
	fi  

	exit 1
	;;
OTUS)
	OUTDIR=$1/data/$2
	PREFIX=$3
	SL=$4
	SR=$5
	OTU=$6
    VER=${7:-0}
	UNFILTDIR=$OUTDIR/$PREFIX/unfiltered	
	
	JOBNAME=OTU_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

	cd $OUTDIR
	dir=`mktemp -d -p $OUTDIR`
	find $UNFILTDIR -name '*.fastq' >$dir/files.txt
	TASKS=$(wc -l $dir/files.txt|awk -F" " '{print $1}')
	qsub -N ${JOBNAME}_1 -t 1-$TASKS:1 -tc 10 $SCRIPT_DIR/submit_fq-fa_global.sh $dir/files.txt $dir $SL $SR ${OUTDIR}/${OTU} $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_1 -N ${JOBNAME}_2 $SCRIPT_DIR/submit_cat_global.sh $dir $OUTDIR $PREFIX $VER $SCRIPT_DIR
	# qsub -hold_jid ${JOBNAME}_2 -N ${JOBNAME}_3 $SCRIPT_DIR/submit_otu_biome.sh $dir $OUTDIR $PREFIX $SCRIPT_DIR
	  qsub -hold_jid ${JOBNAME}_2 $SCRIPT_DIR/submit_tidy.sh $dir OTU_*_1.*
	exit 1
	;;
tax_assign)
	qsub $SCRIPT_DIR/submit_taxonomy.sh $SCRIPT_DIR $@ 
	exit 1
	;;
dist)
	qsub $SCRIPT_DIR/submit_dist.sh $SCRIPT_DIR $@
	exit 1
	;;
denoise)
	qsub $SCRIPT_DIR/submit_denoise.sh $@ $SCRIPT_DIR
	exit 1
	;;
unzip)
	qsub $SCRIPT_DIR/submit_unzip.sh $@ $SCRIPT_DIR
	exit 1
	;;
qcheck)
	qsub $SCRIPT_DIR/submit_qcheck.sh $@ $SCRIPT_DIR
	exit 1
	;;
TEST)
	OUTDIR=$1/data/$2
	PREFIX=$3
	UNFILTDIR=$OUTDIR/$PREFIX/unfiltered
	SL=$4
	SR=$5
	R2=${6:-false}
	JOBNAME=OTU_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

	cd $OUTDIR
	dir=`mktemp -d -p $OUTDIR`
	find $UNFILTDIR -name '*.fastq' >$dir/files.txt
	TASKS=$(wc -l $dir/files.txt|awk -F" " '{print $1}')
	qsub -N ${JOBNAME}_1 -t 1-$TASKS:1 $SCRIPT_DIR/submit_fq-fa_global.sh $dir/files.txt $dir $SL $SR $OUTDIR/${PREFIX}.otus.fa $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_1 -N ${JOBNAME}_2 $SCRIPT_DIR/submit_cat_global.sh $dir $OUTDIR $PREFIX $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_2 -N ${JOBNAME}_3 $SCRIPT_DIR/submit_otu_biome.sh $dir $OUTDIR $PREFIX $SCRIPT_DIR
	qsub -hold_jid ${JOBNAME}_3 $SCRIPT_DIR/submit_tidy.sh $dir OTU_*_1.*
	exit 1
	;;
*)
	echo "Invalid program: $program" >&2
	exit 1
esac
	
