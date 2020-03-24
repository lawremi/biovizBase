parseArgsForAes <- function(args){
    res <- if ("mapping" %in% names(args)) {
               args[["mapping"]]
           } else {
               Find(function(x) {
                   is(x, "uneval")
               }, args)
           }

    if (is.null(res))
        res <- structure(list(), class="uneval")
    
  idx <- ggplot2_is_calculated_aes(res)
  res[idx] <- ggplot2_strip_dots(res[idx])
  res
}

parseArgsForNonAes <- function(args){
    Filter(function(x) !is(eval(x, parent.frame()), "uneval"), args)
}

