# Index 
 1. Description
 3. HMM Preperation for ITS analysis  
 5. [Common workflow](../master/Common%20workflow.md)
 7. [Bacterial 16S workflow](../master/16S%20%20workflow.md)
 18. [Fungal ITS workflow](../master/ITS%20workflow.md)
 20. [Oomycete ITS workflow](../master/Oomycota%20workflow.md)
 22. [Nematode 18S workflow](../master/Nematoda%20workflow.md)
 24. [Statistical analysis](../master/statistical%20analysis.md)
 

## Description
Metabarcoding pipeline for Illumina MiSeq data. This pipeline is base in the main part on usearch (v10 currently) 32 bit (http://www.drive5.com/usearch/) and is designed to run on a Sun Grid Engine cluster. 

Tested for bacterial 16S, fungal ITS, Oomycete ITS and Nematode amplicons

## Setup pipeline
```shell
# set MBPL variable to pipeline folder
MBPL=~/metabarcoding_pipeline
# to set permanetly for future (bash) shell sessions (be careful with this, if you have settings in ~/.profile they will no longer load)
echo export MBPL=~/metabarcoding_pipeline >>~/.bash_profile
```

## HMM Preperation for ITS analysis

The ITS pipelines are more involved and include scripts for removing common regions of 18S (SSU, 5.8S and LSU). The current implementation uses hidden markov models provided with ITSx (http://microbiology.se/software/itsx/) of these regions and HHMMER v 3.1b2 (http://hmmer.janelia.org/) to find them within the seqeunces. Scripts are provided to then remove the regions. 

The HMM files need to be pepared before using the pipeline

The script cut_hmm.pl splits the combined hmm into individual files, which allows the four parts (ssu end, 58S start, 58 end and lsu start) to be in their own files. This allows a roughly three to four fold speed increase for the whole pipeline.  

NOTE: This is dependent on the hmm files downloaded from ITSx - if they change the format it may no longer work.


```shell
perl $MBPL/scripts/cut_hmm.pl $MBPL/hmm/F.hmm $MBPL/hmm/chopped_hmm Fungi

cd $MBPL/hmm/chopped_hmm

cat *SSU*> t1
cat *58S_start* > t2
cat *58S_end* > t3
cat *LSU* > t4

hmmconvert t1 > $MBPL/hmm/Fun/ssu_end.hmm
hmmconvert t2 > $MBPL/hmm/Fun/58s_end.hmm
hmmconvert t3 > $MBPL/hmm/Fun/58s_start.hmm
hmmconvert t4 > $MBPL/hmm/Fun/lsu_start.hmm

cd $MBPL/hmm/Fun
rm -r $MBPL/hmm/chopped_hmm

for f in *.hmm
do
	sed -i -e'/^LENG/a MAXL  90' $f
done

hmmpress ssu_end.hmm
hmmpress 58s_end.hmm
hmmpress 58s_start.hmm
hmmpress lsu_start.hmm
```
#### NOTES
Files copied to $MBPL/hmm  
Repeat for O.hmm for oomycetes (or for any of the other HMMs you want to include) and set Fungi to Oomycota in the call to cut_hmm.pl . The output files from hmmpress will need to be copied to another location e.g. $MBL/hmm/OO

## Taxonomy reference databases
Assigning taxonomy to OTUs requires a reference database(s) and these will need to be configured for use with the pipeline.

The Unite V7 (fungi) and RDP trainset 15 (bacteria) reference database were downloaded from  
http://drive5.com/usearch/manual/utax_downloads.html  
Configuration has been tested with usearch8.1 (probably works the same for 10)
```shell
usearch8.1 -makeudb_utax refdb.fa -output 16s_ref.udb -report 16s_report.txt
usearch8.1 -makeudb_utax refdb.fa -utax_trainlevels kpcofgs ‑utax_splitlevels NVpcofgs -output ITS_ref.udb -report ITS_report.txt
```

### Oomycota
The oomycota database combines three sets of data; 1) a subset (stamenopiles) of the silva_ssu database https://www.arb-silva.de/browser/, 2) a supplied 3rd party database 3) and a usearch specific Unite database (https://unite.ut.ee/sh_files/utax_reference_dataset_28.06.2017.zip)

NOTE: It is now possible to download just the stramenopiles subset from silva - the below code may need modifying 

The silva_ssu and 3rd party databases required slight modification before use with usearch

```shell
# convert silva_ssu multiline to single line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < SILVA_123_SSURef_Nr99_tax_silva.fasta | tail -n +2 >silva.fa

grep Peronosporomycetes -A1 --no-group-separator silva.fa | sed -e 's/ Eukaryota;SAR;Stramenopiles;Peronosporomycetes/;tax=tax=k:SAR;p:Heterokontophyta;c:Oomycota;/' > oomycota.silva.fa


# combine and replace fasta headers with headers including full taxonomy
awk -F";" 'NR==FNR{a[$1]=$0;next;}a[$1]{$0=a[$1]}1' Oomycota.txt Oomycota.fasta > Oomycota_new.fasta


usearch9 -makeudb_sintax Oomycota_new.fasta -output Oomycota.udp

# convert multiline to single line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < oomycetes.ITS1.fa | tail -n +2 > out.fasta

```
### Nematodes
Nematode database is also a subset of Silva_ssu  
Again note that as the required section of the database can be downloaded directly the below (prior to the usearch command) may not be required

```shell
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < SILVA_123_SSURef_Nr99_tax_silva.fasta | tail -n +2 >silva.fa

grep Nematoda -A1 --no-group-separator silva.fa | sed -e 's/ Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;/;tax=k:Metazoa;p:/'| \
awk -F";" '{
	if(NF>1){
		if(NF==5) {print $1";"$2";" $3";c:"$4";s:"$5}
		if(NF==6) {print $1";"$2";" $3";c:"$4";o:"$5";s:"$6}
		if(NF==7) {print $1";"$2";" $3";c:"$4";o:"$6";s:"$7}
		
	} else {print $1}
}'  > nem_tax.fasta

grep Eumetazoa -A1 --no-group-separator silva.fa > Eumetazoa.fa

grep -n -A1 Nematoda Eumetazoa.fa | \
sed -n 's/^\([0-9]\{1,\}\).*/\1d/p' | \
sed -f - Eumetazoa.fa|awk -F";" '{if(NF>1){print $1";tax=k:Metazoa;p:"$6}else {print $1}}' > nonem_tax.fasta

cat nem_tax.fasta nonem_tax.fasta > Eumetazoa_tax.fasta

usearch -makeudb_sintax nem_tax.fasta -output NEM_ref.udp
usearch -makeudb_sintax Eumetazoa_tax.fasta -output xNEM_ref.udp
```

___
### [Common workflow](../master/Common%20workflow.md)
