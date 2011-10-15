setMethod("addSteppings",c("GenomicRanges"),function(obj, group.name, extend.size = 0){
  if(!missing(group.name)){
    if(! (group.name %in% colnames(values(obj))))
      stop("group.name must be one of the elementMetadta column")
  }else{
    group.name <- NULL
  }
  ## levels only make sense on different seqnames
  ## so split them first
  obj.lst <- split(obj,as.character(seqnames(obj)))
  lv <- endoapply(obj.lst,function(x){
    if(!is.null(group.name)){
      x.n <- split(x, values(x)[,group.name])
      irs <- unlist(range(ranges(x.n, ignore.strand = TRUE)))
      irs.new <- resize(irs, fix = "center", width = width(irs) + extend.size)
      irs.new <- sort(irs.new)
      .lvs <- disjointBins(irs.new)
      values(x)$.levels <- .lvs[as.character(values(x)[,group.name])]
      x
    }else{
      irs <- ranges(x)
      values(x)$.levels <- as.numeric(disjointBins(resize(irs, fix = "center",
                                  width = width(irs) + extend.size)))
      x
  }})
  res <- unlist(lv)
  res
})

