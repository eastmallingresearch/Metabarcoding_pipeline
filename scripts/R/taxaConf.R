function (obj, conf = 0.65, level = 7) 
{
  rank <- apply(obj, 1, function(x) {
    l <- (level + 7)
    y <- abs(as.numeric(x[8:l]))
    i <- last(which(y >= conf))
    i[is.na(i)] <- 1
    if (i == 1) {
      s = "(k)"
    }
    else if (i == 2) {
      s = "(p)"
    }
    else if (i == 3) {
      s = "(c)"
    }
    else if (i == 4) {
      s = "(o)"
    }
    else if (i == 5) {
      s = "(f)"
    }
    else if (i == 6) {
      s = "(g)"
    }
    else {
      s = "(s)"
    }
    ret <- paste(x[i], s, sep = "")
    return(ret)
  })
  X <- suppressMessages(as.data.frame(as.matrix(obj)))
  X$rank <- as.character(rank)
  return(X)
}

taxaConfVec <-
function (obj,conf=0.65,level=7){
	rank <- apply(obj,1,function(x){
		l <- (level+7)
		y <-abs(as.numeric(x[8:l]))
		i <- last(which(y>=conf))
		i[is.na(i)] <-1 
		##edit
  	    s=""
	    if(i!=level) {
		 if(i==1){s="(k)"
		 } else if(i==2) {s="(p)"
		 } else if(i==3) {s="(c)"
		 } else if(i==4) {s="(o)"	
		 } else if(i==5) {s="(f)"
		 } else if(i==6) {s="(g)"
		 } else {s="(s)"}
        }
		ret <- paste(x[i],s,sep="")
 		return(ret)
	})
    return(as.character(rank))
}
