parseArgsForAes <- function(args){
  ## aesthetics <- plyr:::compact(args[.all_aesthetics])
  ## aesthetics <- aesthetics[!ggplot2_is.constant(aesthetics)]
  ## aes_names <- names(aesthetics)
  ## aesthetics <- rename_aes(aesthetics)
  ## class(aesthetics) <- "uneval"
  
  aes.lst <- unlist(lapply(args, function(x){
    class(eval(x, parent.frame())) == "uneval" 
  }))
  
  if(length(aes.lst)){
    idx <- base::which(aes.lst)
    if(length(idx)){
      res <- args[[idx]]
      ## res <- eval(args[[idx]])
    }else{
      res <- list()
    }
  }else{
    res <- list()
  }
  ## res <- c(res, aesthetics)
  idx <- ggplot2_is_calculated_aes(res)
  res[idx] <- ggplot2_strip_dots(res[idx])
  res
}

parseArgsForNonAes <- function(args){
  lst <- unlist(lapply(args, function(x) class(eval(x, parent.frame())) != "uneval"))
  args[lst]
}

