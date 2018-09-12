plotTaxa <- function(
	obj=mybiom, 	# obj (phloseq) a phyloseq object which must include taxonomy and sample data (or alternatively an S3 list)
	taxon="phylum", 	# taxon (str) is the taxonomic level of interest
	design="all", 		# condition (str) describes how the samples should be grouped (must be column of sample data)
	proportional=T,	# proportional (bool) whether the graph should use proportional or absolute values
	cutoff=1, 	# cutoff (double) for proportional graphs. 
	topn=0, 		# topn (int)taxons to display (by total reads) for non-prortional graphs. T
	others=T, 	# combine values less than cutoff/topn into group "other"
	ordered=F, 	# order by value (max to min)
	stack="stack",	# stacked graph, set to "dodge" for side by side
	type=2, 		# type: (1) by sample (2) by taxonomy 
	fixed=F, 		# fixed is a ggplot parameter to apply coord_fixed(ratio = 0.1)
	ncol=1, 		# ncol is a ggplot paramter to use n columns for the legend
	trans=T,		# set to False if using original/prior transformed counts (useful for larger OTU tables)
	bw=F,		# set graph to black and white 
	ylab="",
	ylims=NULL,	# 
	margins=c(0.2,0.2,0.2,1.5),
	textSize=14,
	returnData=F,
	legend=T,
	NEW=T,
    conf=0.65,
	transform=function(o,design,...) 
	{	
		suppressPackageStartupMessages(require(DESeq2))
		dots <- list(...)
		
		if(!is.null(dots$calcFactors)) {
			calcFactors <- dots$calcFactors
			dots$calcFactors<-NULL
			if(length(dots)>=1) {
				 assay(varianceStabilizingTransformation(ubiom_to_des(o,design=design,calcFactors=calcFactors),unlist(dots)))
			} else {
				 assay(varianceStabilizingTransformation(ubiom_to_des(o,desing=design,calcFactors=calcFactors)))
			}
		} else {
			if(NEW) {
				assay(varianceStabilizingTransformation(o$dds))
			} else {
				assay(varianceStabilizingTransformation(ubiom_to_des(o,design=as.formula(design)),...))
			}
		}
	}, # data transformation function 
	... # arguments to pass to transform function (obviously they could just be set in the function, but this looks neater)
) {
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(scales))

	if(isS4(obj)) {
		obj <- phylo_to_ubiom(obj)
	} 


	if(trans) {
		temp <- design
		idx <- grep(design,colnames(obj[[3]]))
		if(length(unique(obj[[3]][idx]))<=1) {
			design<-1
		}
		obj[[1]] <- as.data.frame(transform(obj,as.formula(paste0("~",design)),...))
		#obj[[1]] <- as.data.frame(transform(obj,as.formula("~1"),...))
		design<-temp
	}
	#obj[[1]][obj[[1]]] <- obj[[1]][obj[[1]]]+abs(min(obj[[1]][obj[[1]]]))

	obj[[1]][obj[[1]]<0] <- 0


	if(NEW) {
		obj <- list(countData=obj$countData,taxData=obj$taxData,colData=obj$colData)
	}

	

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
	taxa_cut <- taxa_cut[order(taxa_cut[,1],decreasing=T),]
	if(others) {
		taxa_cut <- rbind(taxa_cut,setNames(data.frame(x="others" ,t(colSums(taxa_sum[!taxa_sum[,1]%in%txk,-1]))),names(taxa_cut)))
	}
	taxa_cut <- na.omit(taxa_cut)
	taxa_cut[,1] <- as.factor(taxa_cut[,1])
	if(ordered) {
		taxa_cut[,ncol(taxa_cut)+1] <- 0
		taxa_cut[,1] <- reorder(taxa_cut[,1],-rowSums(taxa_cut[,-1]))
		taxa_cut <- taxa_cut[,-ncol(taxa_cut)]
	}

	if(returnData) {
		return(taxa_cut)
	}

	md2 <- melt(taxa_cut,id=colnames(taxa_cut)[1])

	md2$variable <- factor(md2$variable, levels=levels(md2$variable)[order(levels(md2$variable))]  )

	md2$value <- as.numeric(md2$value)

	if (type==1) {
		g <- ggplot(md2,aes_string(x=md2[,2],y=md2[,3],fill=taxon))
	} else {
		colnames(md2) <- c("taxa",design,"value")
		g <- ggplot(md2,aes_string(x=as.factor(md2[,1]),y=md2[,3],fill=design))
	}
	

	if(bw) {
		g<-g+geom_bar(stat="identity",colour="white",fill="black",position="stack")
	} else {
		g <- g + geom_bar(stat="identity",colour="white",position="stack")		
	}

	g <- g  + xlab("")

	if (fixed) {
		g <- g  + coord_fixed(ratio = 0.1)
	} 

	scaleFUN<-function(x) sprintf("%.0f", x)

	g <- g + scale_y_continuous(expand = c(0, 0),labels = scaleFUN,limits=ylims)
	g <- g + ylab(ylab)
	g <- g + guides(fill=guide_legend(ncol=ncol))
	g <- g + theme_blank()
	g <- g + theme(
		axis.text.x = element_text(angle = 45, hjust = 1,size=textSize),
		plot.margin=unit(margins,"cm"), 
	    axis.line.y = element_line(colour = "black",size=1),
		axis.ticks.x=element_blank(),
		text=element_text(size=textSize),
		axis.title.y=element_text(size=(textSize-2))

	)
	if (type==1) {
		g <- g + theme(legend.text=element_text(face="italic"))
	} else if (type==2) {
		#g <- g + theme(axis.text.x=element_text(face="italic"))
	}
	if(!legend) {
		g<- g+guides(fill=FALSE)
	}
	
	return(g)
}


