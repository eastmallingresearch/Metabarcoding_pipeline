# functions for use with data table

#sets NA to 0 by reference
unsetNA = function(DT) {
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}

#sets 0 to NA by ref
setNA = function(DT) {
  for (j in seq_len(ncol(DT)))
    set(DT,which(DT[[j]]==0),j,NA)
}
  
