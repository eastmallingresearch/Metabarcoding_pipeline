plotOrd <- function (
	obj,
	colData,
	design=NULL, # column(s) of colData
	shapes=NULL, # column(s) of colData
	label=NULL,  # column(s) of colData
	facet=NULL,  # column(s) of colData. This doesn't add a layer to the graph just adds facet as a data column, to use call g + facet_wrap(~facet) 
	plot="point", # or "label"
	labelSize=4, # for text label to point
	labelPosition=c(-1.5,0.5), # for text label to point
	sublabels=F,
	cluster=NULL,
	continuous=F,
	colourScale=c(low="red", high="yellow"), #greyscale low="#000000", high="#DCDCDC"
	cbPalette=F,
	pointSize=2,
	textSize=11,
	textFont="Helvetica",
	xlims=NULL,
	ylims=NULL,
	legend="right",
	title=NULL,
	xlabel,
	ylabel,
	axes=c(1,2),
	dimx=1,
	dimy=2,
	alpha=NULL,

	exclude=T, # sometimes it can be useful to include a vector of points to exclude
	noPlot=F, 
	...
) {

	suppressPackageStartupMessages(require(ggplot2))
	suppressPackageStartupMessages(require(viridis))
	
	if(missing(obj)|missing(colData)) return(print("Error : please specify both obj and colData"))

	ef <- function(X) {
		cat("WARNING: Incorrect columns specified \"",X,"\"",sep="")
		return(NULL)
	}

	design <- tryCatch(colnames(colData[,design,drop=F]),error=function(e)ef(design))
	shapes <- tryCatch(colnames(colData[,shapes,drop=F]),error=function(e)ef(shapes))
	label  <- tryCatch(colnames(colData[,label,drop=F]), error=function(e)ef(label))
	facet  <- tryCatch(colnames(colData[,facet,drop=F]), error=function(e)ef(facet))

	# check if label is set if using label as a plot type
	ll<-T
	if(tolower(plot)=="label") {
		ll<-F
		shapes<-NULL
		alpha<-NULL
		if (!length(label)) {
			print("No label column specified defaulting to first column of colData")
			label=colnames(colData)[1]
		}
 
	}

	obj     <- obj[!rownames(obj)%in%exclude,]
	colData <- colData[!rownames(colData)%in%exclude,]
	
	if(!is.null(title)){if(title=="debugging"){invisible(mapply(assign, ls(),mget(ls()), MoreArgs=list(envir = globalenv())));return(title)}}

	d <- as.data.frame(obj[,axes])

	if(!is.null(cluster)) {
	#	km <- kmeans(obj,...)
	#	colData$Cluster<-as.factor(km$cluster)
	} 

	x = colnames(d)[1]
	y = colnames(d)[2]

	aes_map <-aes_string(x=x,y=y)

	if (length(design)) {
		colour <- if (length(design) > 1) {
			factor(apply(as.data.frame(colData[, design,drop = FALSE]), 1, paste, collapse = " : "))
		}
		else {
			as.factor(colData[[design]])
		}
		d <- cbind(d,colour=colour)
		if(continuous) {
			aes_map <- modifyList(aes_map,aes_string(colour=as.number(colour)))
		} else {
			aes_map <- modifyList(aes_map,aes_string(colour=colour))
		}
	}

	if (length(shapes)) {
		shape <- if (length(shapes) > 1) {
			factor(apply(as.data.frame(colData[, shapes,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[shapes]])
		}
		d <- cbind(d,shapes = shape)
		aes_map <- modifyList(aes_map,aes_string(shape="shapes"))
	}

	if(length(label)) {
		label <- if (length(label) > 1) {
			factor(apply(as.data.frame(colData[, label,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[label]])
		}
		d <- cbind(d,label = label)
		d$label[sublabels] <- NA
		aes_map <- modifyList(aes_map,aes_string(label="label"))
	}

	if(length(facet)) {
		facet <- if (length(facet) > 1) {
			factor(apply(as.data.frame(colData[, facet,drop = FALSE]), 1, paste, collapse = " : "))
		} else {
			as.factor(colData[[facet]])
		}
		d <- cbind(d,facet = facet)
		aes_map <- modifyList(aes_map,aes_string(facet=facet))
	}


	cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

	if(noPlot) return(list(data=d,aes=aes_map))

	### ggplot ###
	g <- ggplot(data=d,aes_map) 
	
	g <- g + if(tolower(plot)=="label"){
		geom_label(size=(labelSize))
	} else {
		geom_point(size=pointSize,na.rm = TRUE)
	}
	g <- g + coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, expand = TRUE)

	g <- g + theme_classic_thin(textSize,textFont) %+replace% theme(legend.position=legend)

	if(!is.null(design)) {
		if(continuous) {
			g <- g + scale_colour_gradient(low=colourScale[1], high=colourScale[2],name=design)
		} else {
			if(cbPalette) {
				g<-g+scale_colour_manual(values=cbbPalette)	+ guides(colour=guide_legend(title=design))
			} else {
				g<-g+scale_colour_viridis(discrete=TRUE)+ guides(colour=guide_legend(title=design))
			}
		}
	}

	if (!is.null(shapes)) {
		g <- g + scale_shape_discrete(name=shapes)
	}

	if(!is.null(alpha)) g <-g+ geom_point(size=(pointSize+(pointSize*1.5)),alpha=alpha)

	if(!is.null(cluster)) {
			g<-g+stat_ellipse(geom="polygon", level=cluster, alpha=0.2)
	}

	if (!missing(xlabel)) {g <- g + xlab(xlabel)}
	if (!missing(ylabel)) {g <- g + ylab(ylabel)}

	if(length(label)&ll) { 
		g <- g + geom_text(size=(labelSize), vjust=labelPosition[1], hjust=labelPosition[2],check_overlap = TRUE)
	}

	return(g)
}
