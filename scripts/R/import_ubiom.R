import_ubiom <- function (
	locX,
	locY,
	locZ
){
	options(stringsAsFactors = FALSE)
	countData <- read.table(locX,header=T,sep="\t", comment.char="")
	rownames(countData ) <- countData [,1]
	countData <- countData [,-1]
	taxonomy <- read.csv(locY,header=F)
	taxonomy <- taxonomy [,c(1,2,4,6,8,10,12,14)]
	rownames(taxonomy) <- taxonomy[,1]
	taxonomy <- taxonomy[,-1]
	colnames(taxonomy) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
	colData <- read.table(locZ,sep="\t",header=T)
	rownames(colData) <- colData [,1]
	colData <- colData[,-1,drop=FALSE]
	countData <- countData[,rownames(colData)]
	ls.biom <- list(countData,colData, taxonomy)
	names(ls.biom) <- c("countData","colData","taxonomy")
	return(ls.biom)
}