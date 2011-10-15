setMethod("getFragLength", "character", function(obj, model, which,
                                                 type = "remove.gaps",
                                                 maxGap = 0L){
  require(Rsamtools)
  if(missing(which))
    which <- range(reduce(model))
  if(tools::file_ext(obj) == "bam")
    res <- scanBam(obj, param = ScanBamParam(which = which, 
                          flag = scanBamFlag(isProperPair = TRUE,
                              isFirstMateRead = TRUE)))

  if(!length(res[[1]]$pos))
    stop("no paired end reads found in the the data")
  res <- res[[1]]
  m1 <- GRanges(res$rname, IRanges(res$pos, width = res$qwidth),
                mapid = 1:length(res$pos), isize = res$isize)
  ## suppose they have the same width
  m2 <- GRanges(res$rname, IRanges(res$mpos, width = res$qwidth),
                mapid = 1:length(res$pos), isize = res$isize)
  m <- c(m1, m2)
  res <- getFragLength(m, model, type = type, maxGap = 0L)
  res
})

setMethod("getFragLength", "GenomicRanges", function(obj, model, maxGap = 0L,
                                                     type = "remove.gaps"){
  model <- reduce(model)
  if(type == "remove.gaps"){
    grl <- split(obj, values(obj)$mapid)
    gr.res <- unlist(range(grl))
    names(gr.res) <- names(grl)
    mapid <- names(gr.res)
    idx <- match(mapid, values(obj)$mapid)
    values(gr.res)$isize <- values(obj)$isize[idx]
    g.gap <- gaps(c(ranges(obj), ranges(model)))
    cut.fun <- shrinkageFun(g.gap, maxGap)
    obj.cut <- cut.fun(obj)
    grl <- split(obj.cut, values(obj.cut)$mapid)
    grl.r <- range(grl)
    .fragLength <- unlist(width(grl.r))
    values(gr.res)$.fragLength <- .fragLength[names(gr.res)]
    gr.res
  }
  list(gr = gr.res,
       cut.fun = cut.fun,
       gr.cut = unlist(grl.r),
       ## model.cut = cut.fun(model),
       ## model = model,
       reads = obj,
       ## gaps = g.gap,
       fragLength = .fragLength)
})

