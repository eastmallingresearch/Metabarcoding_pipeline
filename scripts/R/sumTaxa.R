sumTaxa <- function(
	obj,
	taxon="phylum",
	design="all",
	proportional=F,
	conf=0.65,
	meanDiff=F,
	weighted=T
){
# sums by sample data 
	suppressPackageStartupMessages(require(plyr))
	suppressPackageStartupMessages(require(dplyr))
	suppressPackageStartupMessages(require(tibble))
	suppressPackageStartupMessages(require(reshape2))

	obj[[2]][,taxon] <- taxaConfVec(obj[[2]][,-8],conf=conf,level=which(colnames(obj[[2]])==taxon))
	dtx <- left_join(rownames_to_column(obj[[1]]),rownames_to_column(obj[[2]][,taxon,drop=F]))
	rownames(dtx) <- dtx[,1]
   	md <- melt(dtx,id=c(taxon,"rowname"))

	obj[[3]]$all <- "all"
	if(length(design)>1) {
		obj[[3]][[paste(design,collapse=":")]] <- factor(apply(as.data.frame(obj[[3]][, design,drop = FALSE]), 1, paste, collapse = " : "))
		design <- paste(design,collapse=":")		
	}

	if(meanDiff&length(levels(as.factor(obj[[3]][[design[1]]])))!=2) {
		warning("Design has more than 2 levels mean difference may be difficult to interpret")	
	}

	if (!meanDiff) {
		md <- md[,setdiff(names(md), "rowname")]
	}

	md$variable <- mapvalues(md$variable,from=rownames(obj[[3]]), to=as.character(obj[[3]][,design]))
	nd <- dcast(md,...~variable,sum)

	if(meanDiff) {
		temp <- nd[,c(-1,-2)]
		cols <- ncol(temp)
		nnd <- as.data.frame(lapply(seq(1,(cols-1)),function(i) lapply(seq((i+1),cols),function(j) 
		{x=as.data.frame(abs(temp[,i]-temp[,j]));
		 colnames(x)[1] <-paste(colnames(temp)[i],colnames(temp)[j],sep="|");
		 if(weighted){x<-x*((temp[,i]+temp[,j])/sum(temp[,c(i,j)]))}
		 return(x)}))) 
		nd <- aggregate(nnd,by=list(nd[[taxon]]),sum)
	}

	colnames(nd)[1] <- taxon
	if(proportional) {
		nd[-1] <-  prop.table(as.matrix(nd[,-1]),2)*100
	}
	return(nd)
}

sumTaxaAdvanced <-
function (
	obj,# list(countdata,taxData,colData)
	taxon="phylum", 	# taxon (str) is the taxonomic level of interest
	conf=0.65,
	design="all", 		# condition (str) describes how the samples should be grouped (must be column of sample data)
	proportional=T,	# proportional (bool) whether the graph should use proportional or absolute values
	cutoff=1, 	# cutoff (double) for proportional graphs. 
	topn=0, 		# topn (int)taxons to display (by total reads) for non-prortional graphs. T
	others=T, 	# combine values less than cutoff/topn into group "other"
	ordered=F 	# order by value (max to min)
) {
	
	obj$taxData <- 
	taxa_sum <- sumTaxa(obj,taxon=taxon,design=design,conf=conf)
	taxa_sum$taxon[grep("\\(",taxa_sum$taxon)] <- taxa_sum$taxon[sub("\\(.*"," incertae sedis",taxa_sum$taxon)]

	if(!topn) {
		obj[[3]]$MLUflop <- 1 #assigns the MLU flop digit
		tx <- sumTaxa(obj,taxon=taxon,"MLUflop")
		tx[,-1] <- prop.table(as.matrix(tx[,-1]),2)*100
		txk <- tx[tx[,2]>=cutoff,1]
	} else {
		taxa_sum[,ncol(taxa_sum)+1]<- 0
		taxa_sum <- taxa_sum[order(rowSums(taxa_sum[,-1]),decreasing=T),]
		taxa_sum <- taxa_sum[,-ncol(taxa_sum)]	
		txk <- taxa_sum[1:topn,1]
	}
	
	if(proportional) {
		taxa_sum[,-1] <- prop.table(as.matrix(taxa_sum[,-1]),2)*100
	}

	taxa_cut <- taxa_sum[taxa_sum[,1]%in%txk,]
	taxa_cut <- taxa_cut[order(taxa_cut[,1],decreasing=F),]
	if(others) {
		taxa_cut <- rbind(taxa_cut,setNames(data.frame(x="others" ,t(colSums(as.data.frame(taxa_sum[!taxa_sum[,1]%in%txk,-1])))),names(taxa_cut)))
	}
	taxa_cut <- na.omit(taxa_cut)
	taxa_cut[,1] <- as.factor(taxa_cut[,1])
	if(ordered) {
		taxa_cut <- taxa_cut[order(rowSums(as.data.frame(taxa_cut[,-1])),decreasing=T),]
	}
	
	return(taxa_cut)

}
