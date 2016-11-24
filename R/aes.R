parseArgsForAes <- function(args){
    res <- if ("mapping" %in% names(args)) {
               args[["mapping"]]
           } else {
               Find(function(x) {
                   class(x) == "uneval"
               }, args)
           }

    if (is.null(res))
        res <- structure(list(), class="uneval")
    
  idx <- ggplot2_is_calculated_aes(res)
  res[idx] <- ggplot2_strip_dots(res[idx])
  res
}

parseArgsForNonAes <- function(args){
  lst <- unlist(lapply(args, function(x) class(eval(x, parent.frame())) != "uneval"))
  args[lst]
}

