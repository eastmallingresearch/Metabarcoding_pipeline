#### Load libraries
library(data.table)
####

args <- commandArgs(TRUE)

#print(args)

# directory containing otu tables to merge
tmpdir <- paste0(args[1],"/")

# read otu tables into a list of data tables
qq <- lapply(list.files(tmpdir,args[2],full.names=T),function(x) fread(x))

# full outer join the otu tables
merged_table <- Reduce(function(...) merge(..., all=T,by = "#OTU ID"), qq)

# write combined OTU table
write.table(merged_table,paste0(args[3],"/",args[4],".otu_table.txt"),sep="\t",row.names=F,col.names=T,quote=F,na="0")
