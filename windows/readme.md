## Introduction
Because I've been asked...

## Demultiplexing
```
P1F=CCTACGGGNGGCWGCAG # Bacteria
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA # Fungi
P2R=ATATGCTTAAGTTCAGCGGG

demulti_v2.pl forward_read reverse_read mismatches $P1F $P1R $P2F $P2R
# set mismatches to number of allowed mismatches in index
# Outputs three sets of files (F and R), primer pair1, primer pair2 and unassigned 
```
usearch version:
```
usearch -fastx_demux reads.fq -index index.fq -barcodes bar.fa -filename_suffix .fq
```
I don't like this at all  - requires an index and barcode file and is not read pair aware.  
My version is far simpler (it's written in Perl, so can be made to run on Windows )

Note to self - still need to work on the rest of this. The pipeline uses a lot of shell/awk/sed scripts which will be a bit of a pain to convert to Windows. I have donea  bit of cdm/power shell scripting - 10 years ago... 

## Preprocessing
``` #16S
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
"$ARDERI/data/$RUN/$SSU/fastq/*R1*.fastq" \
$ARDERI/data/$RUN/$SSU \
$ARDERI/metabarcoding_pipeline/primers/adapters.db \
$MINL $MINOVER $QUAL
```

``` #ITS
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
"$ARDERI/data/$RUN/$SSU/fastq/*R1*.fastq" \
$ARDERI/data/$RUN/$SSU/fasta \
$ARDERI/metabarcoding_pipeline/primers/primers.db \
$MINL $MAXR2 $QUAL; 
```

## ITS specific (if wanted - R1 only)
Find ssu/5.8s/lsu
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
$ARDERI/data/$RUN/$SSU/fasta \
R1 \
$ARDERI/metabarcoding_pipeline/hmm/ssu_end.hmm 	\
$ARDERI/metabarcoding_pipeline/hmm/58s_start.hmm \
ssu 58ss 20
``` 
Remove ssu/5.8s/lsu
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$ARDERI/data/$RUN/$SSU/fasta/*R1" \
 $ARDERI/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
 "*.\\.ssu" \
 "*.\\.58"
 ```
 
 correct names 
 (obviously this is not going to work under windows as it is using awk 
 
 maybe could do the same with powershell - been a long time since I used it)
  ```
 for f in /home/deakig/projects/Oak_decline/data/run/ITS/unfiltered/*r1*; do
F=$(echo $f|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}')
L=$(echo $f|awk -F"/" '{print $NF}'|awk -F"." '{print $1}' OFS=".") 
awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $F.tmp && mv $F.tmp $F
done
```
 
## UPARSE pipeline

### Cluster
```
#denoise
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $ARDERI $RUN $SSU $FPL $RPL

# or clustering with cluser_otu
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UCLUS $ARDERI $RUN $SSU $FPL $RP
```

### Assign taxonomy
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL
```

### Create OTU table
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL
```
