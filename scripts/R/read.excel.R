library(tidyverse)
library(data.table)
library(readxl)    
read.excel <- function(filename, filetype = as.data.table) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    names(x) <- sheets
    x <- Filter(Negate(function(x) is.null(unlist(x))), x)
    x <- lapply(x, filetype)
    x
}
