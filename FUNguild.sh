# FUNGuild has been updated and simplfied, and no longer requires a qiime format taxonmy file
# The following will convert a the standard taxonomy file from the metabarcoding pipeline to FUNGuild compatible
awk -F"," '{print $1,$2,$4,$6,$8,$10,$12,$14}' OFS="\t"  < FUN.sintax.taxa > temp.taxa

# Then to run:
FUNGuild.py guild -taxa temp.taxa
# Output will be a file names "temp.guilds.txt" which can be analysed in R




##### OLD STUFF - NO LONGER REQUIRED #####
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
cat V2_FUN.otus_tax_assignments.txt|mod_taxa_qiime.pl > tax.txt

# then from R
library(dplyr)

tax <- read.table("tax.txt",sep=",",header=F)
log <- read.table("log.txt",sep="\t",header=F)
countData <- read.table("../FUN.otus_table.txt",sep="\t",header=T, comment.char = "")

log$perc <- as.numeric(sub("%","",log$V2))
test <- log  %>% group_by(V1) %>% summarise(perc=max(perc))
test2 <- left_join(test,log)
test2<-test2[!duplicated(test2$V1),]
log_tax <- left_join(test2,tax,by=c("V1"="V1"))
counts_log_tax <- left_join(countData,log_tax,by=c("X.OTU.ID","V1))
counts_log_tax$taxonomy <- paste(counts_log_tax$V2.x,counts_log_tax$V8,counts_log_tax$V3.x,counts_log_tax$V4.x,"reps",paste(counts_log_tax$V2.y,counts_log_tax$V3.y,counts_log_tax$V4.y,counts_log_tax$V5,counts_log_tax$V6,counts_log_tax$V7,counts_log_tax$V8,sep=";") ,sep="|")

