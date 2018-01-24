# requires a qiime formatted taxonomy file (with % identity)
# probably best using a ssu/5.8s filtered otu file
for F in  $PROJECT_FOLDER/data/$RUN/FUN.otus.fa; do   
  $PROJECT_FOLDER/metabarcoding_pipeline/scripts/PIPELINE.sh -c ITS_regions \
  $F \
  $PROJECT_FOLDER/metabarcoding_pipeline/hmm/ssu_end.hmm \
  $PROJECT_FOLDER/metabarcoding_pipeline/hmm/58s_start.hmm \
  20
done

# assign taxonomy with qiime/uclust
assign_taxonomy.py \
-r $PROJECT_FOLDER/metabarcoding_pipeline/taxonomies/its/qiime/sh_refs_qiime_ver7_dynamic_01.12.2017.fasta \
-t $PROJECT_FOLDER/metabarcoding_pipeline/taxonomies/its/qiime/sh_taxonomy_qiime_ver7_dynamic_01.12.2017.txt \
-i V2_FUN.otus.fa \
-m uclust \ 
-o . \
--similarity 0.1

# FUNguild requires a specific input file which is a combination of the qiime txt and log output files (and OTU table)
grep -E "^H" V2_FUN.otus_tax_assignments.log|awk -F"[\t_]" '{split($10,a,"\.");print $9,$4"%",a[1],a[2]"_"$11}' OFS="\t" > log.txt

# then from R
library(dplyr)

tax <- read.table("V2_FUN.otus_tax_assignments.txt",sep="\t",header=F)
log <- read.table("log.txt",sep="\t",header=F)
countData <- read.table("../FUN.otus_table.txt",sep="\t",header=T,row.names=1, comment.char = "")
