ubiom_to_des <- function(
	obj, 
	colData="colData",
	countData="countData",
	design=~1,
	fit=F,
	calcFactors=function(d)
	{
		sizeFactors(estimateSizeFactors(d))
	},
	...
){
	suppressPackageStartupMessages(require(DESeq2))
	
	obj[[colData]] <- obj[[colData]][colnames(obj[[countData]]),,drop = FALSE]

	dds <- 	suppressWarnings(DESeqDataSetFromMatrix(obj[[countData]],obj[[colData]],design))

	sizeFactors(dds) <- calcFactors(dds)

    	if (fit) {
    	 	return(DESeq(dds,...))
    	} else {
    		return(dds)
    	}
} 
