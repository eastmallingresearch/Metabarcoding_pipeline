#### Load libraries
library(Biostrings)

args <- commandArgs(TRUE) #get command line variables
#args <- c(".","*.\\.ssu","*.\\.58ss","../S77_R1.fa", "ITS1")
# 1: folder
# 2: regex start file names
# 3: regex end file names
# 4: sample fasta
# 5: sample ID
print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])


if(exists(args[7])) {
	args[4] <- sub("\\.fa","\\.filtered\\.fa",args[4])
}

load_tables <- function (files,func,End=F) {
  tables <- lapply(files, load_table,func=func,End=End)
  do.call(rbind, tables)
}

load_table <- function(file,func,End) {
 tryCatch({
    X <- read.table(file,skip=2,header=F)
    X <- X[with(X, ave(V14,V3, FUN=max)==V14),]
    if (isTRUE(End)) {
      X <- X[!duplicated(X[,c(3,8)]),c(3,8)]
      X <- X[with(X, ave(V8,V3,FUN=func)==V8),]
      colnames(X) <- c("seq","end")
    } else {
      X <- X[!duplicated(X[,c(3,7)]),c(3,7)]
      X <- X[with(X, ave(V7,V3,FUN=func)==V7),]
      colnames(X) <- c("seq","start")
    }
      return(X) 
    }, error = function(err) {
  } )
} 

print("loading ssu tables")
r_start <- load_tables(list.files(args[1],args[2],full.names=T),max,T)
print("loading 58S tables")
r_end   <- load_tables(list.files(args[1],args[3],full.names=T),min)

print("reading DNA table")
myfasta <- readDNAStringSet(args[4])

library(data.table)
library(dplyr)
mytable <- data.table(seq=names(myfasta))
mytable <- as.data.table(left_join(mytable,r_start))
mytable <- as.data.table(left_join(mytable,r_end))
mytable$start[is.na(mytable$start)] <- (width(myfasta[mytable$seq[is.na(mytable$start)]])-as.numeric(args[6]))
mytable$end[is.na(mytable$end)] <- as.numeric(args[7])

print("merging tables")
#mytable <- merge(r_start,r_end,by.all=T,all.x=T)
#if(exists(args[6])){mytable$start<-myfasta[mytable$seq]@ranges@width}

print("removing NAs")
#mytable <- na.omit(mytable)
mytable <- mytable[((mytable$start-mytable$end)>40),]
#myfasta <- myfasta[mytable$seq]

print("set ITS_IR")
ITS_IR <- IRanges(start=mytable$end+1,end=mytable$start-1,names=mytable$seq)

print("set ITS")
ITS <- DNAStringSet(myfasta,start=ITS_IR@start,width=(ITS_IR@width-1))
#ITS <- ITS[ITS@ranges@width>=140]


print("write ITS")
writeXStringSet(ITS,paste(args[5],".r1.fa",sep=""))

