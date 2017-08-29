fltTaxon <- function(
	obj,
	taxon="phylum",
	out="phylo"
){
# same as phyloseq tax_glom (drops NA columns returned by tax_glom), but works on S3 biom data (i.e. ubiom)
# perhaps tax_glom is good with big datasets as fltTaxon is miles faster - aggregate will get pretty slow for large datasets
	if(class(obj)[[1]]=="phyloseq") {
		obj <- phylo_to_ubiom(obj)
	}
	n <- which(colnames(obj[[2]])==taxon)
	x <- aggregate(obj[[1]],by=obj[[2]][,1:n],sum)
	ls.biom <- list(x[,(n+1):ncol(x)],x[,1:n],obj[[3]])
	names(ls.biom) <- c("countData","taxonomy","colData")
	if(out=="phylo") {
		return(ubiom_to_phylo(ls.biom))
	}
	return(ls.biom)	
}