theme_blank <- function(base_size = 11, base_family = "")
{
	theme_bw(base_size = base_size, base_family = base_family) %+replace% 
	theme(
		panel.border = element_blank(), 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank()
	)
}

theme_classic_thin <- function(base_size = 11, base_family = "") 
{
	theme_blank(base_size = base_size, base_family = base_family) %+replace% 
	theme(	
		axis.line.x = element_line(size=0.3,colour = "black"),
		axis.line.y = element_line(size=0.3,colour = "black"),
		axis.text = element_text(colour = "black")
	)
}




#  theme(
#	legend.title = element_blank(),
#	legend.position="bottom", 
#   	legend.direction="horizontal",
#   	legend.key = element_rect(colour = NA),
#	text=element_text(size=18)
#  )