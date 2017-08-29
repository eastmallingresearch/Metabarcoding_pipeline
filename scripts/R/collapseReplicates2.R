# collapses samples with unequal libraries to mean and recalculates size factors 

collapseReplicates2 <- 
function (object, groupby, run, renameCols = TRUE)
{
    if (!is.factor(groupby))
        groupby <- factor(groupby)
    groupby <- droplevels(groupby)
    stopifnot(length(groupby) == ncol(object))
    sp <- split(seq(along = groupby), groupby)
    sizefactors <- sapply(sp, function(i) prod(sizeFactors(object)[i,drop=F]))
    countdata <- sapply(sp, function(i) 
	rowMeans(counts(object,normalize=T)[,i, drop = FALSE]*prod(sizeFactors(object)[i,drop=F])))
    mode(countdata) <- "integer"
    colsToKeep <- sapply(sp, `[`, 1)
    collapsed <- object[, colsToKeep]
    dimnames(countdata) <- dimnames(collapsed)
    assay(collapsed) <- countdata
    if (!missing(run)) {
        stopifnot(length(groupby) == length(run))
        colData(collapsed)$runsCollapsed <- sapply(sp, function(i) paste(run[i],
            collapse = ","))
    }
    if (renameCols) {
        colnames(collapsed) <- levels(groupby)
    }
    sizeFactors(collapsed) <- sizefactors

    #stopifnot(sum(as.numeric(assay(object))) == sum(as.numeric(assay(collapsed))))
    collapsed
}
