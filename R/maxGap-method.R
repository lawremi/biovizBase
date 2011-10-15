setMethod("maxGap", "GenomicRanges", function(obj, ratio = 0.0025){
  ratio * width(range(ranges(obj)))
})
