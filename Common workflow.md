# Common workflow

Set up project folder.

If the project has multiple sequencing runs, RUN should be set to location where files are to be stored.

```shell
PROJECT_FOLDER=~/projects/my_project_folder
RUN=.
ln -s $PROJECT_FOLDER/metabarcoding_pipeline $MBPL

mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
mkdir $PROJECT_FOLDER/data/$RUN/quality
mkdir $PROJECT_FOLDER/data/$RUN/ambiguous
mkdir -p $PROJECT_FOLDER/data/$RUN/16S/fastq
mkdir $PROJECT_FOLDER/data/$RUN/16S/filtered
mkdir $PROJECT_FOLDER/data/$RUN/16S/unfiltered
mkdir -p $PROJECT_FOLDER/data/$RUN/ITS/fastq
mkdir $PROJECT_FOLDER/data/$RUN/ITS/filtered
mkdir $PROJECT_FOLDER/data/$RUN/ITS/unfiltered
```

## Decompress files

The demultiplexing step will accept gz compressed files - so this step may not be necessary

```shell
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*.gz; do 
	$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c unzip $FILE
done
```

## QC
Qualtiy checking with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```shell
for FILE in $$PROJECT_FOLDER/data/$RUN/fastq/*; do 
	$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $$PROJECT_FOLDER/data/$RUN/quality
done
```

## Demultiplexing

We have multiplexed 16S and ITS PCR reactions in same sequencing run which can be seperated by the index
Run demulti.pl to demultiplex these into fungal and bacterial fastq files. Ambiguous reads are written to two (f & r) seperate files.

Running something like the below should give a good indication of what index_1 and index_2 should be - this is useful if you don't knwo what the primer sequences are and to get a feel of how many mismatches (if necesary) to use. 
```shell
sed -n '2~4p' $(ls|head -n1)|grep -x "[ATCG]\+"|cut -c-16|sort|uniq| \
tee zzexpressions.txt|xargs -I%  grep -c "^%" $(ls|head -n1) >zzcounts.txt
```

Any sequence which has too many mismatches, or none mathching primers is removed to a file x.ambigous.fq

demultiplex can accept any number of primer pairs (though for this project only 2 primer pairs are multiplexed)

<table>
Primers:
<tr><td><td>Forward<td>Reverse</tr>
<tr><td>16S<td>CCTACGGGNGGCWGCAG<td>GACTACHVGGGTATCTAATCC</tr>
<tr><td>ITS<td>CTTGGTCATTTAGAGGAAGTAA<td>ATATGCTTAAGTTCAGCGGG</tr>
<tr><td>OO<td>GAAGGTGAAGTCGTAACAAGG<td>AGCGTTCTTCATCGATGTGC</tr>
<tr><td>Nem<td>CGCGAATRGCTCATTACAACAGC<td>GGCGGTATCTGATCGCC</tr>
</table>


```shell
#bacteria and fungi
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

$PROJECT_FOLDER/mbpl/scripts/PIPELINE.sh -c demultiplex \
	"$PROJECT_FOLDER/data/$RUN/fastq/*_R1_*" 0 \
	$P1F $P1R $P2F $P2R


mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/16S/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/ITS/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.


#nematode and oomycete
P1F=CGCGAATRGCTCATTACAACAGC
P1R=GGCGGTATCTGATCGCC
P2F=GAAGGTGAAGTCGTAACAAGG
P2R=AGCGTTCTTCATCGATGTGC

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
	"$PROJECT_FOLDER/data/$RUN/fastq/*Nem*_R1_*" 0 \
	$P1F $P1R $P2F $P2R


mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/NEM/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/OO/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.
```
### Ambiguous data
Ambiguous data should not be used for OTU clustering/denoising, but it can be counted in the OTU tables.
Requires converting to FASTA with approprite labels

### [16S workflow](../master/16S%20%20workflow.md)
### [ITS workflow](../master//ITS%20workflow.md)
### [Oomycete workflow](../master/Oomycota%20workflow.md)
### [Nematode workflow](../master/Nematoda%20workflow.md)
### [Statistical analysis](../master/statistical%20analysis.md)
