# 16/18S workflow

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
SSU=16S
FPL=17
RPL=21
```

## Pre-processing
Script will join PE reads (with a maximum % difference in overlap) remove adapter contamination and filter on minimum size and quality threshold.
Unfiltered joined reads are saved to unfiltered folder, filtered reads are saved to filtered folder.

16Spre.sh forward_read reverse_read output_file_name output_directory adapters min_size percent_diff max_errrors 

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c 16Spre \
"$ARDERI/data/$RUN/$SSU/fastq/*R1*.fastq" \
$ARDERI/data/$RUN/$SSU \
$ARDERI/metabarcoding_pipeline/primers/adapters.db \
$MINL $MINOVER $QUAL
```
## UPARSE

### Cluster 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting 

```shell
#denoise
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $ARDERI $RUN $SSU $FPL $RPL

#clustering with cluser_otu
#$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UCLUS $ARDERI $RUN $SSU $FPL $RPL
```

#### Work around for usearch bug 10.1
```shell
sed -i -e 's/Zotu/OTU/' 16S.zotus.fa
```

### Assign taxonomy
NOTE:- I still need to build nematode utax taxonomy database from Silva_SSU.

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $ARDERI $RUN $SSU 
```

### OTU evolutionary distance

Output a phylogentic tree in phylip format (both upper and lower triangles)
(usearch9 doesn't work)
```shell
usearch8.1 -calc_distmx 16S.otus.fa -distmxout 16S.phy -distmo fractdiff -format phylip_square
```

### Create OTU table 

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL
```

If unfiltered data is too much for usearch(32) to handle :

```shell
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTUS $ARDERI $RUN $SSU $FPL $RPL
```
(this may actually be quicker than OTU, need to check)

### Combine biome and taxa

biom_maker.pl will take a hacked rdp taxonomy file (mod_taxa.pl) and UPARSE biome and output a combined taxa and biome file to standard out.

```shell
$ARDERI/metabarcoding_pipeline/scripts/biom_maker.pl ITS.taxa ITS.otu_table.biom >ITS.taxa.biom
```

### Poor quality work round (oops, make that poor quality <i> data</i> work round)

Occasionally, due to v.poor reverse read quality, joining of f+r reads fails for the vast majority. The following will cluster f+r reads separately and then merge read counts which align to the same OTU. I've dropped the clustering down to 0.95 similarity - both reads aligning to the same OTU at this similarity, I'd suggest is pretty good evidence they're the same. 
I've also added a rev compliment routine to fq2fa_v2.pl, means the reverse reads can be called as plus strand by usearch_global.

```shell
for f in $ARDERI/data/$RUN/16S/fastq/*R1*; do
 R1=$f
 R2=$(echo $R1|sed 's/\_R1_/\_R2_/')
 S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
 $ARDERI/metabarcoding_pipeline/scripts/fq2fa_v2.pl $R1 $ARDERI/data/$RUN/16S.r1.unfiltered.fa $S $fpl 0
 $ARDERI/metabarcoding_pipeline/scripts/fq2fa_v2.pl $R2 $ARDERI/data/$RUN/16S.r2.unfiltered.fa $S $rpl 30 rev
done
usearch9 -usearch_global 16S.r1.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r1.txt -userfields query+target+id
usearch9 -usearch_global 16S.r2.unfiltered.fa -db 16S.otus.fa -strand plus -id 0.95 -userout hits.r2.txt -userfields query+target+id
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c merge_hits $ARDERI/metabarcoding_pipeline/scripts/merge_hits.R hits.r1.txt hits.r2.txt 16S.otu_table.txt
$ARDERI/metabarcoding_pipeline/scripts/otu_to_biom.pl row_biom col_biom data_biom >16S.otu_table.biom
rm row_biom col_biom data_biom
```

### [ITS workflow](../master//ITS%20workflow.md)
### [Oomycete workflow](../master/Oomycota%20workflow.md)
### [Nematode workflow](../master/Nematoda%20workflow.md)
### [Statistical analysis](../master/statistical%20analysis.md)




### OLD
```
### Not implemented (but keep incase something breaks)
#### Concatenate files
#cat $METAGENOMICS/data/$RUN/16S/filtered/*filtered* > $METAGENOMICS/data/$RUN/16S.t.fa
#### Truncate and pad (Remove multiplex primers and pad reads to same length.)
#X=`cat 16S.t.fa|awk '{if ($1!~/>/){mylen=mylen+length($0)}else{print mylen;mylen=0};}'|awk '$0>x{x=$0};END{print x}'`
#usearch8.1 -fastx_truncate 16S.t.fa -stripleft 17 -stripright 21 -trunclen $X -padlen $X -fastaout 16S.fa
#rm 16S.t.fa
#### Dereplication
#cat 16S.fa|awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'|$METAGENOMICS/scripts/get_uniq.pl > #16S.sorted.fasta 
#rm 16S.fa
#### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))
#usearch8.1 -cluster_otus 16S.sorted.fasta -otus 16S.otus.fa -uparseout 16S.out.up -relabel OTU -minsize 2

#usearch8.1 -derep_fulllength 16S.fa -fastaout 16S.uniques.fasta -sizeout 
#usearch8.1 -sortbysize 16S.uniques.fasta -fastaout 16S.sorted.fasta -minsize 2
#rm 16S.fa 16S.uniques.fasta

 	#$METAGENOMICS/data/$RUN/$SSU/filtered \
	#$METAGENOMICS/data/$RUN \
	#$SSU $FPL $RPL

#### Assign Taxonomy
usearch9 -utax 16S.otus.fa -db $METAGENOMICS/taxonomies/utax/16s_ref.udb -strand both -utaxout 16S.reads.utax -rdpout 16S.rdp -alnout 16S.aln.txt
cat 16S.rdp|$METAGENOMICS/scripts/mod_taxa.pl > 16S.taxa

#### Concatenate unfiltered reads (Unfiltered fastq will need to be converted to fasta first )
for f in $METAGENOMICS/data/$RUN/$SSU/unfiltered/*.fastq
do
	S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
	$METAGENOMICS/scripts/fq2fa_v2.pl $f $METAGENOMICS/data/$RUN/$SSU.unfiltered.fa $S $FPL $RPL
done
#### Make table (Creates an OTU table of read counts per OTU per sample)
usearch9 -usearch_global $SSU.unfiltered.fa -db $SSU.otus.fa -strand plus -id 0.97 -biomout $SSU.otu_table.biom -otutabout $SSU.otu_table.txt
```
