plotPCA <- function (	
	obj, 
	design = "condition",
	labelby,
	ntop = 500,
	pcx = 1,
	pcy = 2, 
	returnData = FALSE,
	cofix=F,
	trans=T,
	addLabels=F,
	filter=1,
	filterFun=function(o,f){return(o)},
	transform=function(o,design,...) 
	{	
		suppressPackageStartupMessages(require(DESeq2))
		dots <- list(...)
		if(!is.null(dots$calcFactors)) {
			calcFactors <- dots$calcFactors
			dots$calcFactors<-NULL
			dots$object <- phylo_to_des(o,as.formula(design),calcFactors=calcFactors)
			assay(do.call(varianceStabilizingTransformation, dots))
		} else {
			assay(varianceStabilizingTransformation(phylo_to_des(o,as.formula(design)),...))
		}
	},...
) {
	suppressPackageStartupMessages(require(genefilter))
	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(DESeq2))

	if(trans) {
		obj@otu_table@.Data <- transform(obj,as.formula(paste0("~",design)),...)
	}
	
	obj <- filterFun(obj,filter)    	

 	if (returnData) {
 		#d <- pca$x
 		#attr(d, "percentVar") <- percentVar
 		pca <- prcomp(t(otu_table(obj)))
 		pca$percentVar <- pca$sdev^2/sum(pca$sdev^2)
 		return(pca)
	}
 
 	rv <- rowVars(otu_table(obj))
 
	colData <- sample_data(obj)
	#suppressWarnings(as.data.frame(as.matrix(obj@sam_data)))
	select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	pca <- prcomp(t(otu_table(obj)[select, ]))
	percentVar <- pca$sdev^2/sum(pca$sdev^2)
	
	if (!all(design %in% names(colData))) {
		stop("the argument 'design' should specify columns of colData")
	}
	design.df <- as.data.frame(colData[, design,drop = FALSE])
	group <- if (length(design) > 1) {
		factor(apply(design.df, 1, paste, collapse = " : "))
	}
	else {
		as.factor(sample_data(obj)[[design]])
	}
	
	if (!missing(labelby)) {
		shape <- as.factor(sample_data(obj)[[labelby]])
		d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,design.df,shape = shape)
		colnames(d)[grep("shape", colnames(d))] <- labelby
	} else {
		d <- data.frame(PC1 = pca$x[, pcx], PC2 = pca$x[, pcy], group = group,design.df)
	}

	colnames(d)[grep("group", colnames(d))] <- design

	if(cofix) {
		d[,1] <- d[,1] * percentVar[pcx]
		d[,2] <- d[,2] * percentVar[pcy]
	}

	g <- ggplot()
	g <- g + coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
	if (!missing(labelby)) {
	g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group, shape=shape),size=3)
	g <- g + scale_shape_discrete(name=labelby)
	} else {
	g <- g + geom_point(data=d, mapping=aes(x=PC1, y=PC2, colour=group),size=3)
	}
	g <- g + scale_colour_discrete(name=design)
	g <- g + xlab(paste0("PC",pcx,": ", round(percentVar[pcx] * 100), "% variance"))
	g <- g + ylab(paste0("PC",pcy,": ", round(percentVar[pcy] * 100), "% variance"))
	return(g)
}