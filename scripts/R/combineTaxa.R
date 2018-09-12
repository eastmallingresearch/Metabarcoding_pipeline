combineTaxa2 <-
function(
	taxData, 
	rank="species", 
	# cut-off for matching ranks
	confidence=0.95,
	returnFull=F
) {
	
	# data tables and dplyr are going to be used 
	require(plyr)
	require(data.table)
	
	# new filter 
	levelled <- taxaConfVec(taxData,confidence,which(colnames(taxData)==rank))
	taxData <- data.table(taxData[grep("[^\\(].[^\\)]$",levelled),],keep.rownames="OTU")

	# group any remaining OTUs by the rank value, will return the OTUs as a character list
	taxData <- ddply(taxData,rank,summarize,OTUS=list(as.character(OTU)))
	
	# return taxa associated with more than one OTU
	if(!returnFull) {
		return(taxData[lapply(taxData[,2],function(x) length(unlist(x)))>1,])
	}

	taxData
}

combineTaxa <-
function(
	# path to the taxonomy file
	path, 
	# regex to specify filter for species (or other rank) column ( (s),(g),(f) and etc.)
	rank="species", 
	# cut-off for matching ranks
	confidence=0.95, 
	# column order MUST match k,p,c,o,f,g,s,k_conf,p_conf,c,conf,o_conf,f_conf,g_conf,s_conf
	# use column_order to reorder the columns to the required use -99 to keep the same order(or some other "large" number)  
	column_order=c(1,3,5,7,9,11,13,2,4,6,8,10,12,14),
	# whether the taxonomy file contains a header
	header=F, 
	# taxonomy file seperator
	sep=",",
	 # column with row names, set to NULL if no row names 
	row.names=1,
	# pass any further arguments to read.table
	...
	
) {
	
	# data tables and dplyr are going to be used 
	require(plyr)
	require(data.table)
	
	# read taxonomy file into a data frame
	taxData <- read.table(path,header=header,sep=sep,row.names=row.names,...)

	# reorder columns
	taxData<-taxData[,column_order]

	# add best "rank" at confidence and tidy-up the table
	taxData<-phyloTaxaTidy(taxData,confidence)

	# create a regex string out of the rank value
	rank_reg <- paste0("\\(",tolower(substr(rank,1,1)),"\\)")

	# filter the data for ranks at at least the given confidence
	taxData <- data.table(taxData[grep(rank_reg,taxData$rank),],keep.rownames="OTU")

	# group any remaining OTUs by the rank value, will return the OTUs as a character list
	taxData <- ddply(taxData,~rank,summarize,OTUS=list(as.character(OTU)))
	
	# return taxa associated with more than one OTU
	taxData[lapply(taxData[,2],function(x) length(unlist(x)))>1,]

}

combCounts <- 
function(
	# results from combineTaxa function
	combinedTaxa,
	# associated count data table
	countData
) {
	
	start<-nrow(countData)+1
	countData<-rbind(countData,t(sapply(combinedTaxa[,2],function(x) colSums(rbind(countData[rownames(countData)%in%unlist(x),],0)))))
	end <- nrow(countData)
	rownames(countData)[start:end] <- lapply(combinedTaxa[,2],function(x) paste(x[[1]][1],length(unlist(x)),sep="_"))
	countData <- countData[!row.names(countData)%in%unlist(combinedTaxa[,2]),]
  	countData[complete.cases(countData),]
}


combTaxa <- 
function (
	# results from combineTaxa function
	combinedTaxa,
	# the taxonomy table which was used for combineTaxa
	taxData
) { 
	keep <- taxData[unlist(lapply(combinedTaxa[,2],"[[",1)),]
	rownames(keep)<-lapply(combinedTaxa[,2],function(x) paste(x[[1]][1],length(unlist(x)),sep="_"))
	taxData<- taxData[!row.names(taxData)%in%unlist(combinedTaxa[,2]),]
	rbind(taxData,keep)
}
