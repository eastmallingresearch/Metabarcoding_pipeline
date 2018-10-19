#' Ordination plot
#'
#' plotOrd can help in making pretty PCA/Ordination plots.
#' 
#' This is a function for making simple ggplot ordination plots from two column data with 
#' additional metadata. Based on \code{DESeq2::plotPCA}, but with extra functionality for labeling
#' points, and etc. All the additional functionality can be applied post return of ggplot object.
#' 
#'
#' @param obj a dataframe containing xy coordinates, e.g. principal component
#'   scores.
#' @param colData dataframe containing sample metadata.
#' @param design column(s) of colData to discriminate sample type (colour).
#' @param shapes column(s) of colData to discriminate sample type (shape).
#' @param label column(s) of colData to use for sample labeling.
#' @param facet column(s) of colData. Adds a faceting column to the 
#'  returned ggplot object. To use call g + facet_wrap(~facet). 
#' @param plot Plot either [point] or [label] (Default = "point").
#' @param labelSize Text size for labels(Default = 4).
#' @param labelPosition Label position relative to point(Default = c(-1.5,0.5)).
#' @param sublabels a numeric vector of labels to remove (Default = F).
#' @param cluster Set to turn on clustering, value is stat_sllipse confidence.
#' @param continuous T/F whether design is a continuos variable (default FALSE).
#' @param colourScale Vector used for continuous colour scales (Low to High)
#'   (Default = c(low="red", high="yellow")) #greyscale low="#000000", high="#DCDCDC".
#' @param cbPalette Use a predefined colour blind palette.
#' 	Max eight factors allowable in design (Default = F).
#' @param pointSize The size of plot points (Default = 2).
#' @param textSize The base text size for the graph (Default = 11).
#' @param textFont The base text font for the graph (Default = "Helvetica").
#' @param xlims, ylims Numeric vectors of axis limits, 
#'  e.g. c(-8,8) (Default unlimited).
#' @param legend Position of legend. 
#'   Set to "none" to remove legend (Default = "right").
#' @param legendDesign Display legend for design (Default=True).
#' @param legendShape Display legend for shapes (Default=True).
#' @param title Title (Default is to not use a title).
#' @param xlabel, ylabel Set axis labels.
#' @param axes Columns of obj to use for plotting (Default = c(1,2)).
#' @param alpha Numeric value, "prettifies" points by adding an extra outer circle 
#'   with given alpha value.
#' @param exclude vector of points to exclude (Default show all points).
#' @param noPlot T/F if set do not plot return list of: 
#'   [1] selected obj axes and [2] aesthetics (Default FALSE).
#' @param ... additional parameters (unused).
#' @return A a ggplot scatter plot of the axes taken from obj, 
#' colours as per design and shapes as per shapes (unless noPlot set to TRUE).
#' @examples
#' d <- data.frame(PCA1=runif(10,-8,8),PCA2=runif(10,-4,6))
#' m <- data.frame(Sample=seq(1,10),
#'	Condition=rep(c("H","D"),5),
#'	Site=c(rep("A",5),rep("B",5))) 
#' 
#' plotOrd(obj=d,colData=m,design="Condition")
#'
#' plotOrd(obj=d,colData=m,design="Condition",shapes="Site")
#'
#' plotOrd(obj=d,colData=m,design="Condition",xlims=c(-2,2), label="Sample")

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
	legendDesign=T,
	legendShape=T,
	title=NULL,
	xlabel,
	ylabel,
	axes=c(1,2),
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
			aes_map <- modifyList(aes_map,aes(colour=as.number(colour)))
		} else {
			aes_map <- modifyList(aes_map,aes_string(colour="colour"))
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
		aes_map <- modifyList(aes_map,aes_string(facet="facet"))
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
			g <- g + scale_colour_gradient(low=colourScale[1], high=colourScale[2],name=design,guide=legendDesign)
		} else {
			if(cbPalette) {
				g<-g+scale_colour_manual(values=cbbPalette,guide=legendDesign) + guides(colour=guide_legend(title=design))
			} else {
				g<-g+scale_colour_viridis(discrete=TRUE,guide=legendDesign) + guides(colour=guide_legend(title=design))
			}
		}
	}

	if (!is.null(shapes)) {
		g <- g + scale_shape_discrete(name=shapes,guide=legendShape)
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
