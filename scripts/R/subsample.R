# subsamples OTU data
# set depth to rarefy data, or use some a vector for a proportional scheme (e.g., if you want to combine different numbers of biological replicate per condition)
# sampling can be done with or without replacement

subsample <- function(
  obj,
  depth=1000,
  replace=T,
  sdr=.Random.seed
){
  if(length(depth)==1) {
    depth <- sapply(1:ncol(obj),function(i) min(depth,sum(obj[,i])))
  }

  l <- lapply(1:ncol(obj),function(i) {
    set.seed(sdr)
    if(replace) {
      otu <- sample(rownames(obj),depth[i],replace=T,prob=(obj[,i]/sum(obj[,i])))
	  } else {
      otu <- sample(unlist(sapply(seq_along(obj[,i]),function(x) {rep(rownames(obj)[x],obj[x,i])})))
      otu <- sample(otu,depth[i])
    }
    table(otu)
  })
  cts <- Reduce(function(...) merge(..., all=T,by = "otu"), l)
  rownames(cts) <- cts[,1]
  cts <- cts[-1]
  cts[is.na(cts)] <- 0
  colnames(cts) <- colnames(obj)
  cts
}