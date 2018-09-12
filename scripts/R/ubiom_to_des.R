ubiom_to_des <- function(
	obj, 
	design=~1,
	fit=F,
	filter,
	calcFactors=function(d)
	{
		sizeFactors(estimateSizeFactors(d))
	},
	...
){
	suppressPackageStartupMessages(require(DESeq2))

	invisible(mapply(assign, names(obj), obj,MoreArgs=list(envir = environment())))
	
	colData <- colData[colnames(countData),,drop = FALSE]

	if(!missing(filter)) {
		filter <- eval(filter)
		colData <- droplevels(colData[filter,])
		countData <- countData[,filter]
	}

	dds <- 	suppressWarnings(DESeqDataSetFromMatrix(countData,colData,design))

	sizeFactors(dds) <- calcFactors(dds)

    if (fit) {
     	return(DESeq(dds,...))
    } else {
    	return(dds)
    }
} 
