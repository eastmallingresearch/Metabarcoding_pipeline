# Bacteria workflow

## Conditions
SSU determines the file location
FPL is forward primer length
RPL is reverse primer length

```shell
# all
MINL=300
MINOVER=5
QUAL=0.5

#bacteria
SSU=BAC
FPL=17
RPL=21
```

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size percent_diff max_errrors 

```shell
$PROJECT_FILE/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
"$PROJECT_FILE/data/$RUN/$SSU/fastq/*R1*.fastq" \
$PROJECT_FILE/data/$RUN/$SSU \
$PROJECT_FILE/metabarcoding_pipeline/primers/adapters.db \
$MINL $MINOVER $QUAL
```
## UPARSE

### Cluster 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting 

```shell
#denoise
$PROJECT_FILE/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FILE $RUN $SSU $FPL $RPL

#clustering with cluser_otu
#$PROJECT_FILE/metabarcoding_pipeline/scripts/PIPELINE.sh -c UCLUS $PROJECT_FILE $RUN $SSU $FPL $RPL
```

#### Work around for usearch bug 10.1
```shell
sed -i -e 's/Zotu/OTU/' 16S.zotus.fa
```

### Assign taxonomy
NOTE:- I still need to build nematode utax taxonomy database from Silva_SSU.

```shell
$PROJECT_FILE/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FILE $RUN $SSU 
```

### OTU evolutionary distance

Output a phylogentic tree in phylip format (both upper and lower triangles)
(usearch9 doesn't work)
```shell
usearch8.1 -calc_distmx 16S.otus.fa -distmxout 16S.phy -distmo fractdiff -format phylip_square
```

### Create OTU table 

```shell
$PROJECT_FILE/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FILE $RUN $SSU $FPL $RPL
```

If unfiltered data is too much for usearch(32) to handle :

```shell
$PROJECT_FILE/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FILE $RUN $SSU $FPL $RPL
```
(this may actually be quicker than OTU, need to check)

### Combine biome and taxa

biom_maker.pl will take a hacked rdp taxonomy file (mod_taxa.pl) and UPARSE biome and output a combined taxa and biome file to standard out.

```shell
$PROJECT_FILE/metabarcoding_pipeline/scripts/biom_maker.pl ITS.taxa ITS.otu_table.biom >ITS.taxa.biom
```

### Poor quality data work round 
Occasionally, due to v.poor reverse read quality, joining of f+r reads fails for the vast majority. The following will cluster f+r reads separately and then merge read counts which align to the same OTU. I've dropped the clustering down to 0.95 similarity - both reads aligning to the same OTU at this similarity, I'd suggest is pretty good evidence they're the same. 
I've also added a rev compliment routine to fq2fa_v2.pl, means the reverse reads can be called as plus strand by usearch_global.

```shell
for f in $PROJECT_FILE/data/$RUN/16S/fastq/*R1*; do
 R1=$f
 R2=$(echo $R1|sed 's/\_R1_/\_R2_/')
 S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
 $PROJECT_FILE/metabarcoding_pipeline/scripts/fq2fa_v2.pl $R1 $PROJECT_FILE/data/$RUN/16S.r1.unfiltered.fa $S $fpl 0
 $PROJECT_FILE/metabarcoding_pipeline/scripts/fq2fa_v2.pl $R2 $PROJECT_FILE/data/$RUN/16S.r2.unfiltered.fa $S $rpl 30 rev
done
usearch9 -usearch_global 16S.r1.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r1.txt -userfields query+target+id
usearch9 -usearch_global 16S.r2.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r2.txt -userfields query+target+id
$PROJECT_FILE/metabarcoding_pipeline/scripts/PIPELINE.sh -c merge_hits $PROJECT_FILE/metabarcoding_pipeline/scripts/merge_hits.R hits.r1.txt hits.r2.txt 16S.otu_table.txt
$PROJECT_FILE/metabarcoding_pipeline/scripts/otu_to_biom.pl row_biom col_biom data_biom >16S.otu_table.biom
rm row_biom col_biom data_biom
```

### [ITS workflow](../master//ITS%20workflow.md)
### [Oomycete workflow](../master/Oomycota%20workflow.md)
### [Nematode workflow](../master/Nematoda%20workflow.md)
### [Statistical analysis](../master/statistical%20analysis.md)
