fread.phylip <- 
function(path,...){
	dt <- fread(path,...)
	m <- as.matrix(dt[,-1])
	rownames(m) <- colnames(m) <- dt[[1]]
	return(m)
}