phylo_to_ubiom <- function(
	obj
){
	suppressPackageStartupMessages(require(phyloseq))
	list(
		countData=as.data.frame(obj@otu_table@.Data),
		taxonomy=as.data.frame(obj@tax_table@.Data),
		colData=as.data.frame(suppressWarnings(as.matrix(sample_data(obj)))) # suppresses a warning from the matrix call about loss of S4 state
	)
}