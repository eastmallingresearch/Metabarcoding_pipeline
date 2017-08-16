# Oomycota workflow

## Conditions
SSU determines the file location
FPL is forward primer length
RPL is reverse primer length

```shell
# oomycota
SSU=OO 
FPL=21
RPL=20

# all
MINL=300
MINOVER=5
QUAL=0.5
```

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size min_join_overlap max_errrors 

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OOpre \
  "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
  $PROJECT_FOLDER/data/$RUN/$SSU \
  $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
  $MINL $MINOVER $QUAL
```
### SSU and 5.8S removal 

# move files to keep consistent with Fungal ITS workflow
```shell
mv $PROJECT_FOLDER/data/$RUN/$SSU/filtered/* $PROJECT_FOLDER/data/$RUN/$SSU/fasta.
rename 's/filtered\.//' $PROJECT_FOLDER/data/$RUN/OO/fasta/*.fa
```

#### Identify SSU and 5.8S regions

This will create a large number of array jobs on the cluster

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
 $PROJECT_FOLDER/data/$RUN/$SSU/fasta \
 "" \
 $PROJECT_FOLDER/metabarcoding_pipeline/hmm/others/Oomycota/ssu_end.hmm \
 $PROJECT_FOLDER/metabarcoding_pipeline/hmm/others/Oomycota/58s_start.hmm \
 ssu 58ss 20
```

#### Remove identified SSU and 5.8S regions
```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
  "$PROJECT_FOLDER/data/$RUN/$SSU/fasta/*[^fa]" \
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
  "*.\\.ssu" "*.\\.58"
```

There's a slight problem with one of the scripts and the fasta names...
```shell
find $PROJECT_FOLDER/data/$RUN/$SSU/fasta -type f -name *r1.fa|xargs -I myfile mv myfile $PROJECT_FOLDER/data/$RUN/$SSU/filtered/.

#mkdir $PROJECT_FOLDER/data/$RUN/$SSU/filtered/intermediate
#mv *filtered* $PROJECT_FOLDER/data/$RUN/$SSU/filtered/intermediate/.
#find . -maxdepth 1 -type d -name "S*" -exec mv '{}' intermediate \;

rename 's/\.r1//' *.fa

#for f in *.fa; do
#	sed -i -e 's/ .*//' $f
#done
```

## OTU assignment 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting 

### Cluster 
```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN $SSU 0 0
```
### Assign taxonomy
```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU 
cp $SSU.otus.fa $SSU_v2.otus.fa
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU_v2
rm $SSU_v2.otus.fa
```

### Create OTU tables
```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN $SSU $FPL $RPL
```

###[16S workflow](../master/16S%20%20workflow.md)
###[ITS workflow](../master//ITS%20workflow.md)
###[Nematode workflow](../master/Nematoda%20workflow.md)
###[Statistical analysis](../master/statistical%20analysis.md)



