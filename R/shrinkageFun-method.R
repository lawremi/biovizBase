setMethod("shrinkageFun", "GenomicRanges", function(obj, max.gap = 1L){
  obj <- obj[order(as.numeric(end(obj)))]
  shrinkage <- pmax(0, width(obj) - max.gap)
  shrinksum <- c(0L, cumsum(shrinkage))
  function(x) {
    shift(x, -shrinksum[findInterval(start(x), end(obj)) + 1L])
  }  
})

setMethod("shrinkageFun", "IRanges", function(obj, max.gap = 1L){
  obj <- obj[order(as.numeric(end(obj)))]
  shrinkage <- pmax(0, width(obj) - max.gap)
  shrinksum <- c(0L, cumsum(shrinkage))
  function(x) {
    shift(x, -shrinksum[findInterval(start(x), end(obj)) + 1L])
  }  
})

