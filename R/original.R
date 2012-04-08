## to reduce burden to import too much packages
## I redfined/move functions in some other packages to here
.all_aesthetics <- c("adj", "alpha", "angle", "bg", "cex", "col",  "color", "colour", "fg",
                     "fill", "group", "hjust",  "label", "linetype", "lower", "lty",
                     "lwd",  "max", "middle", "min", "order", "pch", "radius", "sample",
                     "shape" , "size",  "srt",  "upper", "vjust", "weight", "width",
                     "x" ,  "xend", "xmax", "xmin", "xintercept",  "y",
                     "yend", "ymax", "ymin", "yintercept",  "z")

ggplot2_is.constant <- function(x){
    sapply(x, function(x) "I" %in% all.names(asOneSidedFormula(x)))
}

## for ggplot2
ggplot2_is_calculated_aes <- function(aesthetics){
  match <- "\\.\\.([a-zA-z._]+)\\.\\."
  stats <- rep(FALSE, length(aesthetics))
  grepl(match, sapply(aesthetics, deparse))
}


ggplot2_strip_dots <- function(aesthetics){
  match <- "\\.\\.([a-zA-z._]+)\\.\\."
  strings <- lapply(aesthetics, deparse)
  strings <- lapply(strings, gsub, pattern = match, replacement = "\\1")
  lapply(strings, function(x) parse(text = x)[[1]])
}
