dt_to_df <- function(
	DT,
	row_names=1
){
	DF <- as.data.frame(DT)
	row.names(DF) <- DF[,row_names]
	DF <- DF[,-row_names]
	DF
} 
