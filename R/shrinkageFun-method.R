setMethod("shrinkageFun", "GenomicRanges", function(obj, max.gap = 1L){
  obj <- obj[order(as.numeric(end(obj)))]
  shrinkage <- pmax(0, width(obj) - max.gap)
  shrinksum <- c(0L, cumsum(shrinkage))
  function(x) {
    res <- shift(x, -shrinksum[findInterval(start(x), end(obj)) + 1L])
    metadata(res)$coord <- "truncate_gaps"
    metadata(res)$max.gap <- max.gap
    metadata(res)$gaps <- obj
    values(res)$.ori <- x
    res
  }  
})

setMethod("shrinkageFun", "IRanges", function(obj, max.gap = 1L){
  obj <- obj[order(as.numeric(end(obj)))]
  shrinkage <- pmax(0, width(obj) - max.gap)
  shrinksum <- c(0L, cumsum(shrinkage))
  function(x) {
    res <- shift(x, -shrinksum[findInterval(start(x), end(obj)) + 1L])
    metadata(res)$coord <- "truncate_gaps"    
    metadata(res)$max.gap <- max.gap
    metadata(res)$gaps <- obj    
    values(res)$.ori <- x    
    res
  }  
})

is_coord_truncate_gaps<- function(obj){
  if(!is.null(metadata(obj)$coord))
    res <- metadata(obj)$coord == "truncate_gaps"
  else
    res <- FALSE
  res
}
