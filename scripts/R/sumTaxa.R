sumTaxa <- function(
	obj,
	taxon="phylum",
	design="condition",
	proportional=F
){
# sums by sample data 
	suppressPackageStartupMessages(require(plyr))
	suppressPackageStartupMessages(require(reshape2))
	tx <- obj[[2]][,taxon]
	dtx <- cbind(obj[[1]],tx)
	md <- melt(dtx,id="tx")
	obj[[3]]$all <- "all"
	md$variable <- mapvalues(md$variable,from=rownames(obj[[3]]), to=as.character(obj[[3]][,design]))
	nd <- dcast(md,...~variable,sum)
	colnames(nd)[1] <- taxon
	if(proportional) {
		nd[-1] <-  prop.table(as.matrix(nd[,-1]),2)*100
	}
	return(nd)
}