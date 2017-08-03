#### Load libraries
library(data.table)
####

args <- commandArgs(TRUE)

tmpdir <- paste0(args[2],"/")
qq <- lapply(list.files(tmpdir,args[3],full.names=T),function(x) fread(x))
qm <- as.data.table(melt(qq,id.vars="#OTU ID",))

colnames(qm)[1] <- c("OTU")
qx <- dcast(qm,OTU~variable,sum)
write.table(qx,paste0(args[1],"/",args[4],".otu_table.txt"),sep="\t",row.names=F,quote=F,na="0") 

qm2 <- melt(qx,id.vars="OTU")
qm2 <- qm2[value!=0]
qm2<-qm2[order(variable)]
qm2<-qm2[,variable:=as.factor(variable)]
qm2<-qm2[,sam_pos:=(as.numeric(variable)-1)]
qm2 <- qm2[order(OTU)]
qm2<-qm2[,OTU:=as.factor(OTU)]
qm2<-qm2[,otu_pos:=(as.numeric(OTU)-1)]
write.table(qm2[,5:3,with=F],paste0(tmpdir,"data_biom"),sep=",",row.names=F,col.names=F,quote=F,na="0")
write.table(levels(qm2$OTU),paste0(tmpdir,"row_biom"),sep=",",row.names=F,col.names=F,quote=F)
write.table(levels(qm2$variable),paste0(tmpdir,"col_biom"),sep=",",row.names=F,col.names=F,quote=F)
