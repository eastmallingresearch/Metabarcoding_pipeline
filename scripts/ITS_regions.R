library(Biostrings)
library(data.table)
args <- commandArgs(TRUE) 

RHB <- args[1]
ssu <- args[2]
ss <-  args[3]

print(getwd())

# get ssu data
X <- read.table(ssu,skip=2,header=F) # annoyingly fread can't handle multiple whitespace as a seperator
X <- X[with(X, ave(V14,V3, FUN=max)==V14),]
X <- X[!duplicated(X[,c(3,8)]),c(3,8)]
X <- X[with(X, ave(V8,V3,FUN=max)==V8),]
colnames(X) <- c("seq","end")
r_start <- data.table(X,key="seq")

# get 58ss data
X <- read.table(ss,skip=2,header=F)
X <- X[with(X, ave(V14,V3, FUN=max)==V14),]
X <- X[!duplicated(X[,c(3,7)]),c(3,7)]
X <- X[with(X, ave(V7,V3,FUN=min)==V7),]
colnames(X) <- c("seq","start")
r_end <- data.table(X,key="seq")

# read OTUs into DNAStringSet
myOTUs <- readDNAStringSet(RHB)

# rearrange to be consistent with the data tables
myOTUs <- myOTUs[order(myOTUs@ranges@NAMES)]

# get the OTU names
mytable <- data.table(seq=names(myOTUs),key="seq")

# add the 58ss start position
mytable <- setkey(r_end[mytable],seq)

# add the ssu end position
mytable <- setkey(r_start[mytable],seq)

# for any OTUs lacking a 58ss start position set ITS region to end of string
mytable$start[is.na(mytable$start)] <- width(myOTUs[mytable$seq[is.na(mytable$start)]])

# for any OTUs lacking an ssu end set ITS position to start of string
mytable$end[is.na(mytable$end)] <- 0

# output the ITS region postions
write.table(mytable,paste0(RHB,".regions.txt"),sep="\t",quote=F,row.names=F)

# Find OTUs with odd looking ssu and 58s regions and set to full length sequence
mytable[!(mytable$start-mytable$end)>10,2:3] <- data.frame(0,width(myOTUs[mytable[!(mytable$start-mytable$end)>10]$seq]))

# create a range object from mytable
ITS_IR <- IRanges(start=mytable$end+1,end=mytable$start-1,names=mytable$seq)

# extract ranges from OTU string set
ITS <- DNAStringSet(myOTUs,start=ITS_IR@start,width=(ITS_IR@width))

# write fasta output
writeXStringSet(ITS,sub("(\\/)(?!.*\\/)","\\/V2_",RHB,perl=T))
