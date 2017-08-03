# Index 
 1. HMM Preperation for ITS analysis  
 2. [Common workflow](../master/Common%20workflow.md)
 6. [16S workflow](../master/16S%20%20workflow.md)
 18. [ITS workflow](../master/ITS%20workflow.md)
 20. [Oomycete workflow](../master/Oomycota%20workflow.md)
 22. [Nematode workflow](../master/Nematoda%20workflow.md)
 24. [Statistical analysis](../master/statistical%20analysis.md)
 
## HMM Preperation for ITS analysis
Using HHMMER v 3.1b2 (http://hmmer.janelia.org/)

Used HMM files from ITSx (http://microbiology.se/software/itsx/)

```shell
perl $ARDERI/metabarcoding_pipeline/scripts/cut_hmm v.3.1 $ARDERI/metabarcoding_pipeline/hmm/chopped_hmm fungi
cd $ARDERI/metabarcoding_pipeline/hmm/chopped_hmm
cat *SSU*> t1
cat *58S_start* > t2
cat *58S_end* > t3
cat *LSU* > t4
hmmconvert t1 > ssu_end.hmm
hmmconvert t2 > 58s_end.hmm
hmmconvert t3 > 58s_start.hmm
hmmconvert t4 > lsu_start.hmm

rm t1 t2 t3 t4

for f in *.hmm
do
	sed -i -e'/^LENG/a MAXL  90' $f
done

hmmpress ssu_end.hmm
hmmpress 58s_end.hmm
hmmpress 58s_start.hmm
hmmpress lsu_start.hmm
```
Ouptut files were copied to $METAGENOMICS/hmm. Hacked the HMM files to include a MAXL satement (required) and manually split out SSU,58S and LSU into seperate files (only fungal hmms are implemented in this pipeline)

## Taxonomy reference databases
Reference databases were downloaded from:
http://drive5.com/usearch/manual/utax_downloads.html
(Unite V7 and RDP trainset 15)
```shell
usearch8.1 -makeudb_utax refdb.fa -output 16s_ref.udb -report 16s_report.txt
usearch8.1 -makeudb_utax refdb.fa -utax_trainlevels kpcofgs ‑utax_splitlevels NVpcofgs -output ITS_ref.udb -report ITS_report.txt
```

Oomycota database was created from a subset of the silva_ssu (stamenopiles) database
```shell
#combine and replace fasta headers with headers including full taxonomy
awk -F";" 'NR==FNR{a[$1]=$0;next;}a[$1]{$0=a[$1]}1' Oomycota.txt Oomycota.fasta > Oomycota_new.fasta

usearch9 -makeudb_sintax Oomycota_new.fasta -output Oomycota.udp

# convert multiline to single line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < oomycetes.ITS1.fa | tail -n +2 > out.fasta

```
Nematode database is also a subset of Silva_ssu - I've made a second taxonomy file containing other, non nematode eukaryotes.
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

usearch9 -makeudb_sintax nem_tax.fasta -output nematode.udp
usearch9 -makeudb_sintax Eumetazoa_tax.fasta -output nematode2.udp
```

___
### [Common workflow](../master/Common%20workflow.md)
### [16S workflow](../master/16S%20%20workflow.md)
### [ITS workflow](../master//ITS%20workflow.md)
### [Statistical analysis](../master/statistical%20analysis.md)

## Combine samples from multiple NGS runs

Biom table for samples from multiple NGS runs are required.

This will mean the names of each fasta will need to be made unique and the sequence lengths will need to be set to the same.

Something like the below will copy samples with the wholename string to a new location. Uses original fastq file name and reconsructing the sample ID for each sample used in the workflow.
```shell
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/16S/filtered/"$3"D"$1".filtered.fa 16S/filtered/."}' > runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/16S/unfiltered/"$3"D"$1".unfiltered.fastq 16S/unfiltered/."}' >> runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/ITS/filtered/"$3"D"$1".fa ITS/filtered/."}' >> runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R1*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/ITS/unfiltered/"$3"D"$1".r1.unfiltered.fastq ITS/unfiltered/."}' >> runme.sh
find .. -type f -wholename "*[0-9]/fastq/G[0|O]*R2*"|awk -F"/" '{print $2"_"$4}'|awk -F"_" '{print "cp ../"$1"/ITS/unfiltered/"$3"D"$1".r2.unfiltered.fastq ITS/unfiltered/."}' >> runme.sh
./runme.sh
```


