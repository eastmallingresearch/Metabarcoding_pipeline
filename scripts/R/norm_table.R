norm_table <- 
function(phylo,calcFactors=function(o,...)estimateSizeFactorsForMatrix(o,...),...)
{

	sizeFactors <- calcFactors(otu_table(phylo),...)
	if(missing(sizeFactors)) {
		sizeFactors <- sample_data(phylo)$sizeFactors
	}
	
	 t(t(otu_table(phylo))/sizeFactors)
}

