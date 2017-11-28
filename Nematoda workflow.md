# Nematode workflow

Pipeline updated to use PE workflow

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
MINL=150
MINOVER=10
QUAL=0.5
```

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size min_join_overlap max_errrors 

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c NEMpre \
  "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
  $PROJECT_FOLDER/data/$RUN/$SSU \
  $PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
  $MINL $MINOVER $QUAL $FPL $RPL
```

## UPARSE

### Cluster 

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE \ $PROJECT_FOLDER $RUN $SSU 23 18
```
### Assign taxonomy

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign \ $PROJECT_FOLDER $RUN $SSU 
```

### Create OTU tables

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU \ $PROJECT_FOLDER $RUN $SSU $FPL $RPL
```

### [Bacteria workflow](../master/BAC%20%20workflow.md)  
### [Fungi workflow](../master//FUN%20workflow.md)  
### [Oomycete workflow](../master/Oomycota%20workflow.md)
