des_to_phylo <- function(
	obj
){
	suppressPackageStartupMessages(require(phyloseq))
	suppressPackageStartupMessages(require(DESeq2))
	X<-
	phyloseq(
		otu_table(assay(obj),taxa_are_rows=T),
	 	sample_data(as.data.frame(colData(obj)))
	 )
}