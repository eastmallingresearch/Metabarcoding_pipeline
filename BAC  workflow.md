# Bacteria workflow

## Conditions
SSU determines the file location  
FPL is forward primer length  
RPL is reverse primer length  

```shell
#bacteria
SSU=BAC
FPL=17
RPL=21

MINL=300
MINOVER=5
QUAL=0.5

```

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size percent_diff max_errrors 

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
"$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
$PROJECT_FOLDER/data/$RUN/$SSU \
$PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
$MINL $MINOVER $QUAL $FPL $RPL 
```
## UPARSE

### Cluster 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory certain steps. Alternative scripts for dereplication and sorting are supplied. The final two command arguments will strip bases, left and right respectively.

```shell
#denoise
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN $SSU 0 0
#clustering with cluser_otu
#$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UCLUS $PROJECT_FOLDER $RUN $SSU $FPL $RPL
```

#### Work around for usearch bug 10.1
```shell
sed -i -e 's/Zotu/OTU/' 16S.zotus.fa
```

### Assign taxonomy

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU 
```

### OTU evolutionary distance

Output a phylogentic tree in phylip format (both upper and lower triangles)
(usearch9 doesn't work - 10 does work, but the syntax is slightly diffrent, need to update this...)
```shell
usearch -calc_distmx BAC.zotus.fa -distmxout BACZ.phy -format phylip_square
usearch -calc_distmx BAC.otus.fa -distmxout BAC.phy -format phylip_square
```

### Create OTU table 

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN $SSU $FPL $RPL
```

If unfiltered data is too much for usearch(32) to handle:

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $PROJECT_FOLDER $RUN $SSU $FPL $RPL
```


### Combine biome and taxa

biom_maker.pl will take a hacked rdp taxonomy file (mod_taxa.pl) and UPARSE biome, and output a combined taxa and biome file to standard out.

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/biom_maker.pl 16S.taxa 16S.otu_table.biom > 16S.taxa.biom
```

### Poor quality data work round 
Occasionally, due to v.poor reverse read quality, joining of f+r reads fails for the vast majority. The following will cluster f+r reads separately and then merge read counts which align to the same OTU. I've dropped the clustering down to 0.95 similarity - both reads aligning to the same OTU at this similarity, I'd suggest is pretty good evidence they're the same. 
I've also added a rev compliment routine to fq2fa_v2.pl, means the reverse reads can be called as plus strand by usearch_global.
(this needs testing with usearch10)

```shell
for f in $PROJECT_FOLDER/data/$RUN/16S/fastq/*R1*; do
 R1=$f
 R2=$(echo $R1|sed 's/\_R1_/\_R2_/')
 S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/fq2fa_v2.pl $R1 $PROJECT_FOLDER/data/$RUN/16S.r1.unfiltered.fa $S $fpl 0
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/fq2fa_v2.pl $R2 $PROJECT_FOLDER/data/$RUN/16S.r2.unfiltered.fa $S $rpl 30 rev
done
usearch9 -usearch_global 16S.r1.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r1.txt -userfields query+target+id
usearch9 -usearch_global 16S.r2.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r2.txt -userfields query+target+id
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c merge_hits $PROJECT_FOLDER/metabarcoding_pipeline/scripts/merge_hits.R hits.r1.txt hits.r2.txt 16S.otu_table.txt
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/otu_to_biom.pl row_biom col_biom data_biom >16S.otu_table.biom
rm row_biom col_biom data_biom
```

### [Fungi workflow](../master//FUN%20workflow.md)  
### [Nematode workflow](../master/Nematoda%20workflow.md)
### [Oomycete workflow](../master/Oomycota%20workflow.md)
