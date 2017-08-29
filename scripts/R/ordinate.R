ordinate <- 
function (physeq, method = "DCA", distance = "bray", formula = NULL,
    ...)
{
    library(phyloseq)
    if (inherits(physeq, "formula")) {
        .Deprecated(msg = paste0("First argument, `physeq`, as formula is deprecated.\n",
            "There is now an explicit `formula` argument.\n",
            "Please revise method call accordingly."))
        formchar = as.character(physeq)
        if (length(formchar) < 3) {
            stop("Need both sides of formula in this deprecated syntax... Revisit ordinate() documentation / examples.")
        }
        physeq <- get(as.character(physeq)[2])
        newFormula = as.formula(paste0("~", formchar[length(formchar)]))
        return(ordinate(physeq, method = method, distance = distance,
            formula = newFormula, ...))
    }
    method_table <- c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS",
        "MDS", "PCoA")
    if (inherits(physeq, "character")) {
        if (physeq == "help") {
            cat("Available arguments to methods:\n")
            print(c(method_table))
            cat("Please be exact, partial-matching not supported.\n")
            cat("Can alternatively provide a custom distance.\n")
            cat("See:\n help(\"distance\") \n")
            return()
        }
        else if (physeq == "list") {
            return(c(method_table))
        }
        else {
            cat("physeq needs to be a phyloseq-class object, \n")
            cat("or a character string matching \"help\" or \"list\". \n")
        }
    }
    if (!inherits(physeq, "phyloseq") & !inherits(physeq, "otu_table")) {
        stop("Expected a phyloseq object or otu_table object.")
    }
    if (method == "DCA") {
        return(decorana(veganifyOTU(physeq), ...))
    }
    if (method %in% c("CCA", "RDA")) {
        return(cca.phyloseq(physeq, formula, method, ...))
    }
    if (method == "CAP") {
        return(capscale.phyloseq(physeq, formula, distance, ...))
    }
    if (method == "DPCoA") {
        return(DPCoA(physeq, ...))
    }
    if (inherits(distance, "dist")) {
        ps.dist <- distance
    }
    else if (class(distance) == "character") {
        vegdist_methods <- c("manhattan", "euclidean", "canberra",
            "bray", "kulczynski", "jaccard", "gower", "altGower",
            "morisita", "horn", "mountford", "raup", "binomial",
            "chao")
        if (method == "NMDS" & distance %in% vegdist_methods) {
            return(metaMDS(veganifyOTU(physeq), distance, ...))
        }
        ps.dist <- distance(physeq, distance, ...)
    }
    if (method %in% c("PCoA", "MDS")) {
        return(pcoa(ps.dist))
    }
    if (method == "NMDS") {
        return(metaMDS(ps.dist))
    }
}
