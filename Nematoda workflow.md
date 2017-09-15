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
MINL=200
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
 $PROJECT_FOLDER/data/$RUN/$SSU/fasta \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/nematode.db \
 $MINL $MAXL $QUAL
```

#### Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2

```shell
mkdir -p $PROJECT_FOLDER/data/$RUN/$SSU/filtered
find $PROJECT_FOLDER/data/$RUN/$SSU/fasta -type f -name *_R*|xargs -I myfile mv myfile $PROJECT_FOLDER/data/$RUN/$SSU/filtered/.

cd $PROJECT_FOLDER/data/$RUN/$SSU/filtered
for f in $PROJECT_FOLDER/data/$RUN/$SSU/filtered/*R1.fa
do
    R1=$f
    R2=$(echo $R1|sed 's/_R1\.fa/_R2\.fa/')
    S=$(echo $f|awk -F"_" '{print $1}'|awk -F"/" '{print $NF}')
    $PROJECT_FOLDER/metabarcoding_pipeline/scripts/catfiles_v2.pl $R1 $R2 $S;
done

mkdir R1
mkdir R2
mv *_R1* R1/.
mv *_R2* R2/
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


###[16S workflow](../master/16S%20%20workflow.md)
###[Statistical analysis](../master/statistical%20analysis.md)
