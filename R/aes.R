parseArgsForAes <- function(args){
  aes.lst <- unlist(lapply(args, function(x){
    class(eval(x, parent.frame())) == "uneval"
  }))
  if(length(aes.lst)){
    idx <- base::which(aes.lst)
    if(length(idx))
      res <- eval(args[[idx]])
    else
      res <- list()
  }else{
    res <- list()
  }
  idx <- ggplot2_is_calculated_aes(res)
  res[idx] <- ggplot2_strip_dots(res[idx])
  res
}

parseArgsForNonAes <- function(args){
  lst <- unlist(lapply(args, function(x) class(eval(x, parent.frame())) != "uneval"))
  args[lst]
}

