batchEffectDes <-
function(obj,batch) {

	rld <- log2(counts(obj,normalize=T)+1)
	mypca <- prcomp(t(rld))

	myformula <- as.formula(paste("mypca$x~",batch,sep=""))

	pc.res <- resid(aov(myformula,obj@colData))

	mu <- colMeans(t(rld))

	Xhat <- pc.res %*% t(mypca$rotation)
	Xhat <- t(scale(Xhat, center = -mu, scale = FALSE))
	Xhat <- (2^Xhat)-1

	Xhat[Xhat<0] <- 0
	Xhat <- round(Xhat,6)

	return(Xhat)
}
