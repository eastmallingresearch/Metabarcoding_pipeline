phylo_to_des <- function(
	obj,
	design=~1,
	fit=F,
	obj2=NA,	
	calcFactors=function(d,o)
	{
		sizeFactors(estimateSizeFactors(d))
	},
	...
){
	suppressPackageStartupMessages(require(phyloseq))
	suppressPackageStartupMessages(require(DESeq2))
	dds <-  phyloseq_to_deseq2(obj,design)
	sizeFactors(dds) <- calcFactors(dds,obj2)
    	if (fit) {
    	 	return(DESeq(dds,...))
    	} else {
    		return(dds)
    	}
} 