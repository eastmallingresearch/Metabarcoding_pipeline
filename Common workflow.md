# Common workflow

## Set up project folders

If the project has multiple sequencing runs, RUN should be set to location where files are to be stored.

```shell
PROJECT_FOLDER=~/projects/my_project_folder
mkdir -p $PROJECT_FOLDER
ln -s $PROJECT_FOLDER/metabarcoding_pipeline $MBPL

RUN=.
mkdir -p $PROJECT_FOLDER/data/$RUN/fastq
mkdir $PROJECT_FOLDER/data/$RUN/quality
mkdir $PROJECT_FOLDER/data/$RUN/ambiguous

RIB="BAC FUN OO NEM"
for s in $RIB; do
mkdir -p $PROJECT_FOLDER/data/$RUN/$s/fastq
mkdir $PROJECT_FOLDER/data/$RUN/$s/filtered
mkdir $PROJECT_FOLDER/data/$RUN/$s/unfiltered
mkdir $PROJECT_FOLDER/data/$RUN/$s/fasta
done
```

Copy raw Fastq files (or link to) to $PROJECT_FOLDER/data/$RUN/fastq

## Decompress files (not required, unless using none multiplexed data)

Append 2 to decompress sym links

```shell
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*.gz; do 
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c unzip $FILE 2
done
```

## QC
Qualtiy checking with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```shell
for FILE in $PROJECT_FOLDER/data/$RUN/fastq/*; do 
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/data/$RUN/quality
done
```

## Demultiplexing

Script demulti.pl demultiplexs mixed (e.g. ITS and 16S) libraries based on the primer sequence. Number of acceptable mismatches in the primer sequence can be specified (0 by default). Any sequence which has too many mismatches, or none mathching primers is written to ambiguous.fq (f & r seperately). The script accepts multiple primer pairs.

<table>
Possible primers:
<tr><td><td>Forward<td>Reverse</tr>
<tr><td>Bacteria<td>CCTACGGGNGGCWGCAG<td>GACTACHVGGGTATCTAATCC</tr>
<tr><td>Fungi<td>CTTGGTCATTTAGAGGAAGTAA<td>ATATGCTTAAGTTCAGCGGG</tr>
<tr><td>Oomycete<td>GAAGGTGAAGTCGTAACAAGG<td>AGCGTTCTTCATCGATGTGC</tr>
<tr><td>Nematode<td>CGCGAATRGCTCATTACAACAGC<td>GGCGGTATCTGATCGCC</tr>
</table>

If the primers are unknown, running something like the below should give a good indication of what they are (run from folder of decompressed fastq files). It will also give a good indication of how many mismatches (if any) to use for demulti.pl. 
```shell
sed -n '2~4p' $(ls|head -n1)|grep -x "[ATCG]\+"|cut -c-16|sort|uniq| \
tee zzexpressions.txt|xargs -I%  grep -c "^%" $(ls|head -n1) >zzcounts.txt
```

Run below to demultiplex
```shell
# Demultiplex bacterial and fungal amplicons
P1F=CCTACGGGNGGCWGCAG
P1R=GACTACHVGGGTATCTAATCC
P2F=CTTGGTCATTTAGAGGAAGTAA
P2R=ATATGCTTAAGTTCAGCGGG

$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c demultiplex \
"$PROJECT_FOLDER/data/$RUN/fastq/*_R1_*" 0 \
$P1F $P1R $P2F $P2R

mv $PROJECT_FOLDER/data/$RUN/fastq/*ps1* $PROJECT_FOLDER/data/$RUN/BAC/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ps2* $PROJECT_FOLDER/data/$RUN/FUN/fastq/.
mv $PROJECT_FOLDER/data/$RUN/fastq/*ambig* $PROJECT_FOLDER/data/$RUN/ambiguous/.

# Demultiplex nematode and oomycete amplicons
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
Would require converting to FASTA with approprite labels  - the below should do this will join PE and remove adapters/primers
```shell
$PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c AMBIGpre \
"$PROJECT_FOLDER/data/$RUN/ambiguous/*R1*.fastq" \
$PROJECT_FOLDER/data/$RUN/ambiguous \
$PROJECT_FOLDER/metabarcoding_pipeline/primers/primers.db \
$PROJECT_FOLDER/metabarcoding_pipeline/primers/adapters.db \
300 5
```

### [16S workflow](../master/16S%20%20workflow.md)
### [ITS workflow](../master//ITS%20workflow.md)
### [Oomycete workflow](../master/Oomycota%20workflow.md)
### [Nematode workflow](../master/Nematoda%20workflow.md)
### [Statistical analysis](../master/statistical%20analysis.md)
