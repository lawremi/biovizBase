## to reduce burden to import too much packages
## I redfined/move functions in some other packages to here

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
