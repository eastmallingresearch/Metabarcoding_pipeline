theme_blank <- function(base_size = 16, base_family = "")
{
	theme_bw(base_size = base_size, base_family = base_family) %+replace% 
	theme(
		panel.border = element_blank(), 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank()
	)
}

theme_classic_thin <- function(base_size = 16, base_family = "") 
{
	theme_blank(base_size = base_size, base_family = base_family) %+replace% 
	theme(	
		axis.line.x = element_line(size=0.3,colour = "black"),
		axis.line.y = element_line(size=0.3,colour = "black"),
		axis.text = element_text(colour = "black")
	)
}


theme_facet_blank <- function(base_size = 16, base_family = "",angle=-90,t=2,r=0,b=0,l=0,hjust=0,vjust=0.5,text_angle="vertical")
{
	# angle=-90,hjust=0,vjust=0.5
	# angle=90,hjust=1,vjust=0.5
	# angle=45, hjust=1, vjust=1
	# angle=45, hjust=0.5, vjust=0.5 (centered)
	# angle=-45, hjust=0.5, vjust=0.5 (centered)
	# angle=-45, hjust=0, vjust=1
	# angle =0,hjust=0.5
	
	theme_classic_thin(base_size = base_size, base_family = base_family) %+replace% 
	theme(
		panel.border = element_rect(colour = "black", fill=NA, size=0.5),
		axis.text.x = element_text(angle = angle, margin = margin(t=t,r=r,b=b,l=l),hjust=hjust,vjust=vjust)
	)
}



#  theme(
#	legend.title = element_blank(),
#	legend.position="bottom", 
#   	legend.direction="horizontal",
#   	legend.key = element_rect(colour = NA),
#	text=element_text(size=18)
#  )
