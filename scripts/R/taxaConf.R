taxaConf <- 
function (obj,conf=0.65,level=7){

  isEmpty <- function(x) {
    return(length(x)==0)
  }
  
  	rank <- apply(obj,1,function(x){
		l <- (level+7)
		y <-abs(as.numeric(x[8:l]))
		i <-  last(which(y>=conf))
		#i <- which.max(y<conf)-1 # which.max can't distinguish between all or zero passing
		#if(i<(level-2)&(y[(i+2)]>=conf)) {
		#	i <- max.col((t(y)<conf),"last")-1
		#}
		i[isEmpty(i)] <-1 
		#print(i)
		# i[i<1] <- sum(y>conf);
		##edit
		if(i==1){s="(k)"
		}else if(i==2) {s="(p)"
		}else if(i==3) {s="(c)"
		}else if(i==4) {s="(o)"	
		}else if(i==5) {s="(f)"
		}else if(i==6) {s="(g)"
		}else {s="(s)"}
		ret <- paste(x[i],s,sep="")
		##end edit
 		return(ret)
	})
	X<-suppressMessages(as.data.frame(as.matrix(obj)))
	X$rank <- as.character(rank)
	return(X)
	#return(tax_table(as.matrix(X)))
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
