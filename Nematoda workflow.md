# Nematode workflow

## Conditions
SSU determines the file location
FPL is forward primer length
RPL is reverse primer length

```shell
# nematodes
SSU=NEM 
FPL=23
RPL=18

# all
MINL=100
MAXL=300
QUAL=1
```
## Pre-processing
Script will:<br>
1. Remove reads with both forward and reverse primers<br>
2. Remove reads with adapter contamination<br>
3. Filter for quality and minimum length (with UTRIM)<br>
4. Convert FASTQ to single line FASTA

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c NEMpre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/nematode.db \
 $MINL $MAXL $QUAL
```


### SSU/58S/LSU removal 
Prbably best to skip this part - doesn't work too well (or at all) with the primer set we use.  
I've removed it form the pipeline for now.

If you really want to run it, it's the same as the fungal version, but will require additional HMM files.  

### Move (forward) fasta and rename headers
```shell
for f in $PROJECT_FOLDER/data/$RUN/$SSU/fasta/*.fa; do 
 F=$(echo $f|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $f|awk -F"/" '{print $NF}'|awk -F"." '{print $1}' OFS=".") ;
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $f > $F.tmp && mv $F.tmp $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$F;
done
```

## UPARSE

### Cluster 

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE \ $PROJECT_FOLDER $RUN $SSU 0 0
```
### Assign taxonomy

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign \ $PROJECT_FOLDER $RUN $SSU 
```

### Create OTU tables

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU \ $PROJECT_FOLDER $RUN $SSU $FPL $RPL true
```

### [Bacteria workflow](../master/16S%20%20workflow.md)  
### [Fungi workflow](../master//ITS%20workflow.md)  
### [Oomycete workflow](../master/Oomycota%20workflow.md)
