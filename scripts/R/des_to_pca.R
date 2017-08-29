des_to_pca <- function(obj) {
	X<-prcomp(t(assay(varianceStabilizingTransformation(obj))))
	X$percentVar<-X$sdev^2/sum(X$sdev^2)
	X
}