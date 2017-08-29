ubiom_to_phylo <- function(
	obj
){
	suppressPackageStartupMessages(require(phyloseq))
	phyloseq(
		otu_table(obj[[1]],taxa_are_rows=T),
	 	tax_table(as.matrix(obj[[2]])),
	 	sample_data(obj[[3]])
	 )
}