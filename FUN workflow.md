# Fungi workflow
This workflow is for none overlapping ITS primers, i.e. covering both ITS regions.
The Oomycete workflow is more appopriate if the primers do overlap.

```shell
SSU=FUN
FPL=23 
RPL=21

MINL=200
QUAL=1
```

## Pre-processing
Script will:<br>
1. Remove reads with both forward and reverse primers<br>
2. Remove reads with adapter contamination<br>
3. Remove primers
4. Filter for quality and minimum length<br>
5. Convert FASTQ to single line FASTA

If not running SSU/58S/LSU removal, set FPL to 68  and RPL to 66 for the below. 

```shell

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITSpre \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fastq/*R1*.fastq" \
 $PROJECT_FOLDER/data/$RUN/$SSU \
 $PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
 $MINL $QUAL $FPL $RPL
```

### SSU/58S/LSU removal 

It is debatable whether this is necessary - and it can take a while to run (on a buzy cluster). The SSU/LSU shared region should be about 45 bases after the forward or reverse primers are removed. This can be removed in the pre-processing step. It's more tricky to get  the start postion for the 58S shared region due to variable length of the ITS region. Just trimming off the final 60 odd bases gives good results. It's debatable whether this is of any use either - especially if using denoising (potentially with 97% clustering it could be an issue).

If not using any further preprocessing the  below should be run to get forward reads in the correct format for the UPARSE stages
```
for F in $PROJECT_FOLDER/data/$RUN/$SSU/fasta/*_R1.fa; do 
 FO=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $F|awk -F"/" '{print $NF}'|awk -F"_" '{print $1}') ;
 echo $F
 echo $FO
 echo $L
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $F > $FO.tmp && mv $FO.tmp $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$FO;
done
```

I've split this into a forward only and a forward and reverse pipeline.  
The forward pipeline will need to be run for both (except where stated)

#### Forward pipeline

##### Identify SSU, 5.8S  and LSU regions

This will create a large number of array jobs on the cluster

```shell
# forward reads
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
$PROJECT_FOLDER/data/$RUN/$SSU/fasta \
R1 \
$PROJECT_FOLDER/metabarcoding_pipeline/hmm/ssu_end.hmm \
$PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_start.hmm \
ssu 58ss 20
```

##### Remove SSU, 5.8S  and LSU regions
```shell
# forward reads
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fasta/*R1" \
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_SSU_58Ss.R \
 "*.\\.ssu" \
 "*.\\.58"
 $FPL $RPL
```

##### Run if using forward only
Move fasta and rename fasta header
```shell
for D in $PROJECT_FOLDER/data/$RUN/$SSU/fasta/*1; do 
 F=$(echo $D|awk -F"/" '{print $NF}'|awk -F"_" '{print $1".r1.fa"}'); 
 L=$(echo $D|awk -F"/" '{print $NF}'|awk -F"." '{print $1}' OFS=".") ;
 awk -v L=$L '/>/{sub(".*",">"L"."(++i))}1' $D/$F > $F.tmp && mv $F.tmp $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$F;
# mv $PROJECT_FOLDER/data/$RUN/$SSU/fasta/${L}_R1/$F $PROJECT_FOLDER/data/$RUN/$SSU/filtered/$L; 
done
```

#### Forward and reverse pipeline

##### Identify SSU, 5.8S  and LSU regions

This will create a large number of array jobs on the cluster

```shell
# reverse reads
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c procends \
 $PROJECT_FOLDER/data/$RUN/$SSU/fasta \
 R2 \
 $PROJECT_FOLDER/metabarcoding_pipeline/hmm/lsu_start.hmm \
 $PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_end.hmm \
 lsu 
```

##### Remove SSU, 5.8S  and LSU regions and merge output
If reverse read quality was poor and it was necessary to truncate reads to get more than a couple of reads past set LOWQUAL to TRUE

LOWQUAL keeps reads which lack 5.8S homology - this is necessary as trimming will in most instances have removed the homologous region
```shell
# reverse reads
LOWQUAL=FALSE   
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS \
 "$PROJECT_FOLDER/data/$RUN/$SSU/fasta/*R2" \
 $PROJECT_FOLDER/metabarcoding_pipeline/scripts/rm_58Se_LSU_v2.R \
 "*.\\.58" \
 "*.\\.lsu" \
 $LOWQUAL
```

##### Return ITS1 where fasta header matches ITS2, unique ITS1 and unique ITS2

```shell
find $PROJECT_FOLDER/data/$RUN/$SSU/fasta -type f -name *.r*.fa|xargs -I myfile mv myfile $PROJECT_FOLDER/data/$RUN/$SSU/filtered/.

cd $PROJECT_FOLDER/data/$RUN/$SSU/filtered
for f in $PROJECT_FOLDER/data/$RUN/$SSU/filtered/*r1.fa
do
    R1=$f
    R2=$(echo $R1|sed 's/\.r1\.fa/\.r2\.fa/')
    S=$(echo $f|awk -F"." '{print $1}'|awk -F"/" '{print $NF}')
    $PROJECT_FOLDER/metabarcoding_pipeline/scripts/catfiles_v2.pl $R1 $R2 $S;
done

mkdir R1
mkdir R2
mv *r1* R1/.
mv *r2* R2/.
```

## UPARSE
FPL=23 
RPL=21

### Cluster 
This is mostly a UPARSE pipeline, but usearch (free version) runs out of memory for dereplication and subsequent steps. I've written my own scripts to do the dereplication and sorting 

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $PROJECT_FOLDER $RUN $SSU 0 0
```

Work around for usearch bug 10.1
```shell
sed -i -e 's/Zotu/OTU/' FUN.zotus.fa
```

### Assign taxonomy
```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c tax_assign $PROJECT_FOLDER $RUN $SSU 
```

### OTU evolutionary distance

Output a phylogentic tree in phylip format (both upper and lower triangles)
(usearch9 doesn't work - haven't tested v10)
```shell
usearch -calc_distmx FUN.otus.fa -distmxout FUN.phy -distmo fractdiff -format phylip_square
```

### Create OTU tables

Concatenates unfiltered reads, then assigns forward reads to OTUs. For any non-hits, attemps to assign reverse read (ITS2) to an OTU. 

For forward only remove the true (or set to false) at the end of the command

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $PROJECT_FOLDER $RUN $SSU $FPL $RPL true
```


### Combine biome and taxa

biom_maker.pl will take a hacked rdp taxonomy file (mod_taxa.pl) and UPARSE biome and output a combined taxa and biome file to standard out.

```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/biom_maker.pl ITS.taxa ITS.otu_table.biom >ITS.taxa.biom
```

### [Bacteria workflow](../master/BAC%20%20workflow.md)  
### [Nematode workflow](../master/Nematoda%20workflow.md)
### [Oomycete workflow](../master/Oomycota%20workflow.md)
