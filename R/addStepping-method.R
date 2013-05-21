setMethod("addStepping",c("GenomicRanges"),function(obj, group.name, extend.size = 0,
                                                    fix = "center",
                                                     group.selfish = TRUE){
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
      if(!group.selfish){
        x.n <- split(x, values(x)[,group.name])
        irs <- unlist(range(ranges(x.n, ignore.strand = TRUE)))
        irs.new <- resize(irs, fix = fix, width = width(irs) + extend.size)
        irs.new <- sort(irs.new)
        .lvs <- disjointBins(irs.new)
        values(x)$stepping <- .lvs[as.character(values(x)[,group.name])]
        x
      }else{
      values(x)$stepping <- as.numeric(as.factor(values(x)[,group.name]))
      x
    }
    }else{
      irs <- ranges(x)
      values(x)$stepping <- as.numeric(disjointBins(resize(irs, fix = "center",
                                  width = width(irs) + extend.size)))
      x
  }})
  res <- unlist(lv)
  res
})

