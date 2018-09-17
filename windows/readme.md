## Introduction
Metabarcoding pipeline for Windows (because I've been asked...)

The pipeline should be able to run on Cygwin without any modification (maybe).

## Pure windows version
Scripts are run in the following order:  

decompress??  
demulti_v3.pl (this will need decompressed input - via zcat currently. Will need a windows tool to handle gz files??)  
submit_16Spre_v2.sh - this uses usearch, awk, sed and perl  
adapt_delete.pl stdin (sorted/unique usearch -search_oligodb -userout via awk scripts)  
submit_uparse_v2.sh   
  dereplication and sorting via get_uniq.pl (this takes single line fasta as input - easy to make it take multiline)  
  usearch  
submit_taxonomy.sh  
  usearch   
  mod_taxa.pl/mod_taxa_sintax.pl (uses stdin - can be adjusted)  
submit_cat_files.sh  
  cat only  
submit_global_search.sh  
  usearch  
  
The quickest method (outside a Cygwin implementation) would be to use a Perl + usearch implementation.
Shouldn't be too hard to implement - I'll stay away from any powershell/vbscript for now.

## Perl for Windows installation

https://learn.perl.org/installing/windows.html

Set path to include perl\bin and perl\site\bin
```cmd
REM set perl paths
setx path "c:\perl\bin;c:\perl\site\bin"

REM install perl module installer
cpan App::cpanminus

REM install modules
cpanm Scalar::Util
cpanm List::Util
cpanm List::UtilsBy
cpanm IO::Uncompress::Gunzip
```

## Decompression
Perl includes modules for decompression...
 

## Demultiplexing (and primer matching)
This should always be run-even if not demuliplexing-to ensure both forawrd and reverse primers are present in a read pair 

demulti_v3.pl can handle both gz and uncompressed input fastq

```cmd
demulti_v3.pl FORWARD_READ REVERSE_READ [MAX_ALLOWED_DIFF] FORWARD_PRIMER_1 REVERSE_PRIMER_1 [FORWARD_PRIMER_n] [REVERSE_PRIMER_n]
```
MAX_ALLOWED_DIFF has a max of 9 (default 0 if not set). Greater than 9, units digit indicates no. allowed mismtches and common adapter sequence will be identified in F/R reads as well as primer.

No limit to number of primer pairs 

Output files will have the same name as the input, but appended with ps1.fastq/ps2.fastq etc. or ambig.fastq

```cmd
set P1F=CCTACGGGNGGCWGCAG
set P1R=GACTACHVGGGTATCTAATCC
set P2F=CTTGGTCATTTAGAGGAAGTAA
set P2R=ATATGCTTAAGTTCAGCGGG

perl demulti_v3.pl sample1_forward.fq.gz sample1_reverse.fq.gz 0 %P1F% %P1R% %P2F% %P2R%

```

## Preprocessing

16S
```cmd

powershell.exe -noprofile -executionpolicy bypass -file 16S_pre.ps1

REM merge f + r
usearch -fastq_mergepairs FORWARD_READ -reverse REVERSE_READ -fastqout OUTFILE.t1.fq  -fastq_pctid 0 -fastq_maxdiffs (MINL * MAXDIFF)/100 -fastq_minlen MINL -fastq_minovlen 0 -minhsp 0

REM find adapter contamination
usearch -search_oligodb OUTFILE.t1.fq -db ADAPTERS_DB -strand both -userout OUTFILE.t1.txt -userfields query+target+qstrand+diffs+tlo+thi+trowdots 

REM remove adapter contaminated reads and f + r primers (removes primers by length rather than by match - guaranteed exact if demultiplex set to 0)
perl adapt_delete.pl OUTFILE.t1.txt OUTFILE.t1.fq FORWARD_PRIMER_LENGTH REVERSE_PRIMER_LENGTH > SAMPLE_NAME.unfiltered.fa

REM filter on quality
usearch -fastq_filter SAMPLE_NAME.unfiltered.fa -fastq_maxee QUALITY_SCORE -relabel SAMPLE_NAME -fastaout SAMPLE_NAME.filtered.fa

REM convert to single line fasta - not required yet (or at all?)
REM awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <${OUTFILE}.t3.fa > ${OUTFILE}.filtered.fa
REM sed -i -e '1d' ${OUTFILE}.filtered.fa
```

ITS
```cmd
REM ITS to do
```

 
## UPARSE pipeline

### Cluster
```
type %FILTDIR%/*.fa > %PREFIX%.fa

#### Dereplication
REM min_count would usually be set to a number > 0; default is 1
get_uniq.pl min_count < %PREFIX%.fa > %PREFIX%.sorted.fasta 

#### Clustering (Cluster dereplicated seqeunces and produce OTU fasta (also filters for chimeras))

usearch -cluster_otus ${PREFIX}.sorted.fasta -otus ${PREFIX}.otus.fa -relabel OTU -minsize 4

mv ${PREFIX}.otus.fa $OUTDIR/.

usearch -unoise3 ${PREFIX}.sorted.fasta -zotus ${PREFIX}.zotus.fa #-relabel OTU #-minampsize 8

#usearch -unoise ${PREFIX}.sorted.fasta -tabbedout ${PREFIX}.txt -fastaout ${PREFIX}.otus.fa -relabel OTU #-minampsize 8

#perl -pi -e 's/uniq.*/OTU . ++$n/ge' ${PREFIX}.otus.fa

#rm ${PREFIX}.sorted.fasta


mv ${PREFIX}.zotus.fa $OUTDIR/.


#denoise
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UPARSE $ARDERI $RUN $SSU $FPL $RPL

# or clustering with cluser_otu
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c UCLUS $ARDERI $RUN $SSU $FPL $RP



```

### Assign taxonomy
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL
```

### Create OTU table
```
$ARDERI/metabarcoding_pipeline/scripts/PIPELINE.sh -c OTU $ARDERI $RUN $SSU $FPL $RPL
```
