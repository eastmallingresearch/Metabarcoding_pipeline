ubiome_to_taxmap <- function(
	obj
){
  suppressPackageStartupMessages(require(metacoder))
  counts   <- as.data.frame(obj[[1]], stringsAsFactors = FALSE)
  meta     <- obj[[2]]
  tax_data <- as.data.frame(obj[[3]], stringsAsFactors = FALSE)
  
	phyloseq(
		otu_table(obj[[1]],taxa_are_rows=T),
	 	tax_table(as.matrix(obj[[3]])),
	 	sample_data(obj[[2]])
	 )
}
function (obj, class_regex = "(.*)", class_key = "taxon_name")
{   
    datasets <- list()
    mappings <- c()
    possible_ranks <- unique(unlist(strsplit(taxa::ranks_ref$ranks,
        split = ",")))
    if (!is.null(obj@otu_table)) {
        otu_table <- obj@otu_table
        if (!otu_table@taxa_are_rows) {
            otu_table <- t(otu_table)
        }
        otu_table <- as.data.frame(otu_table, stringsAsFactors = FALSE)
        otu_table <- cbind(data.frame(otu_id = rownames(otu_table),
            stringsAsFactors = FALSE), otu_table)
        datasets <- c(datasets, list(otu_table = otu_table))
        mappings <- c(mappings, c(`{{name}}` = "{{name}}"))
    }
    if (!is.null(obj@sam_data)) {
        sam_data <- as.data.frame(as.list(obj@sam_data), stringsAsFactors = FALSE)
        if (!is.null(rownames(obj@sam_data))) {
            sam_data <- cbind(sample_id = rownames(obj@sam_data),
                sam_data)
        }
        sam_data[] <- lapply(sam_data, as.character)
        datasets <- c(datasets, list(sample_data = sam_data))
        mappings <- c(mappings, NA)
    }
    if (!is.null(obj@phy_tree)) {
        datasets <- c(datasets, list(phy_tree = obj@phy_tree))
        mappings <- c(mappings, NA)
    }
    if (!is.null(obj@refseq)) {
        refseq <- as.character(obj@refseq)
        datasets <- c(datasets, list(ref_seq = refseq))
        mappings <- c(mappings, c(`{{name}}` = "{{name}}"))
    }
    tax_cols <- colnames(tax_data)
    output <- taxa::parse_tax_data(tax_data = tax_data, datasets = datasets,
        class_cols = tax_cols, mappings = mappings, named_by_rank = TRUE,
        class_regex = class_regex, class_key = class_key)
    output$filter_taxa(output$taxon_names() != "NA")
    if ("otu_table" %in% names(output$data)) {
        otu_tab_index <- which(names(output$data) == "otu_table")
        output$data <- c(output$data[otu_tab_index], output$data[-otu_tab_index])
    }
    return(output)
}
