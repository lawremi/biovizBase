setMethod("spliceSummary", c("GenomicRanges", "GRangesList"), function(obj,
                                                             model,
                                                             weighted = TRUE){
  lst <- lapply(model, function(md){
    logi <- isMatchedWithModel(md, obj)
    if(!length(logi))
      logi <- rep(FALSE, length(obj))
    logi
  })
  res <- do.call(rbind, lst)
  if(weighted){
    res.s <- apply(res, 2, sum)
    idx.all <- res.s > 1
    res.n <- as.numeric(res)
    res.m <- matrix(res.n, nrow = nrow(res))
    res.m[,idx.all] <- res.m[,idx.all]/res.s[idx.all]
    freq <- apply(res.m, 1, sum)
  }else{
    freq <- apply(res, 1, sum)
  }
  names(freq) <- names(model)
  freq
})

setMethod("spliceSummary", c("character", "GRangesList"), function(obj,
                                                                   model,
                                                                   weighted = TRUE){

  if(tools::file_ext(obj) == "bam"){
    ## bf <- BamFile(obj)
    require(Rsamtools)    
  }else{
    stop("Only support file with .bam extention name now") 
  }
  
  ## if(missing(which))
  which <- range(reduce(unlist(model)))
  
  ga <- readBamGappedAlignments(obj, param = ScanBamParam(which = which),
                                use.name = TRUE)
  ## reduce model
  ## assumption here is that: contains unoverlaped exons.
  ## must contain exon id            
  exons <- reduce(unlist(model), ignore.strand = TRUE);
  ## first leave only junction read
  idx <- sapply(cigar(ga), isJunctionRead)
  ga.junction <- ga[idx]
  gr.junction <- as(ga.junction, "GRanges")
  spliceSummary(gr.junction, model = model, weighted = weighted)
})


isJunctionRead <- function(cigar){
  grepl("N", cigar)
}

## model is a GRanges containing exons, sgr is a single gr object
isMatchedWithModel <- function(model, gr){
  ## require on the same chromosome now
  if(length(unique(as.character(seqnames(model)))) > 1)
    stop("only support single chromosome summary now")
  if(unique(as.character(seqnames(model))) !=
     unique(as.character(seqnames(gr))))
    stop("Cannot be mapped due to different space")
  model <- sort(model)
  ## do we consider strand information?
  countOverlaps(ranges(gr), ranges(model)) == 2
}

setMethod("spliceSummary", c("GenomicRanges", "GenomicRanges"), function(obj,
                                                                         model,
                                                                       model_id = NULL){
  ov <- findOverlaps(ranges(obj), ranges(model))
  mm <- matchMatrix(ov)
  res <- split(mm[,2], mm[,1])
  idx <- unlist(lapply(split(mm[,2], mm[,1]),length)) > 1
  res <- lapply(res[idx], function(x){
    if(!is.null(model_id))
      x <- values(model)[,model_id][x]
    paste(head(x, 1), tail(x, 1), sep = "-")
  })
  nms <- names(table(unlist(res)))
  res <- as.numeric(table(unlist(res)))
  names(res) <- nms
  res
})


setMethod("spliceSummary", c("character", "GenomicRanges"), function(obj,
                                                                   model,
                                                                   model_id = NULL){

  if(tools::file_ext(obj) == "bam"){
    require(Rsamtools)    
  }else{
    stop("Only support file with .bam extention name now") 
  }
  ## if(missing(which))
  which <- range(model)
  
  ga <- readBamGappedAlignments(obj, param = ScanBamParam(which = which),
                                use.name = TRUE)
  ## first leave only junction read
  grl <- grglist(ga)
  idx <- elementLengths(grl) > 1
  ga.junction <- ga[idx]
  gr.junction <- as(ga.junction, "GRanges")
  spliceSummary(gr.junction, model = model, model_id = model_id)
})

