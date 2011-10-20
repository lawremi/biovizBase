## setGeneric("fetch", function(obj, ...) standardGeneric("fetch"))

## ======================================================================
##  For "TranscriptDb" object
## ======================================================================
## FIXME: to remove the dependency, I change S4 to simple check
## So don't have to import class
## setMethod("fetch", "TranscriptDb", function(obj, which,...,gene.id,
##                                             type = c("all", "single",
##                                               "exons.reduce",
##                                               "exons.all",
##                                               "exons.geneid")){
  
fetch <- function(obj, which, ..., gene.id,
                  resize.extra = 10,
                  include.level = TRUE,
                  use.name = TRUE,
                  type = c("all", "single", "exons.reduce", "exons.all", "exons.geneid",
                            "gapped.pair", "raw", "all")){

  .txdb.type <- c("all", "single", "exons.reduce", "exons.all", "exons.geneid")
  .bamfile.type <- c("gapped.pair", "raw", "all")

  ## type <- match.arg(type)
  if(length(type) >1){
  if(is(obj,"TranscriptDb"))
    type <- "all"
  if(is(obj,"BamFile"))
    type <- "gapped.pair"
  }
  
  if(is(obj, "TranscriptDb")){
    if(!(type %in% .txdb.type))
      stop("type for TranscriptDb must be ", .txdb.type)
    require(GenomicFeatures)
    if(type == "exons.geneid"){
      message("return exons within one gene...")
      ## this require gene id, only exons
      exons.gene <- exonsBy(obj, by = "gene")
      if(missing(gene.id))
        stop("Please specify the gene id")
      res <- exons.gene[[as.character(gene.id)]]
    }
    if(type == "exons.reduce"){
      message("Parsing exons...")
      exons <- exonsByOverlaps(obj, which, ...)
      exons <- reduce(exons)
      values(exons)$type <- "exon"
      res <- exons
    }
    if(type == "exons.all"){
      message("Parsing exons...")
      exons <- exonsByOverlaps(obj, which, ...)
      values(exons)$type <- "exon"
      res <- exons
    }
    if(type %in% c("all","single")){
      message("Parsing exons...")
      exons <- exonsBy(obj, "tx")
      message("Parsing cds...")
      cdss <- cdsBy(obj, "tx")
      ## message("Parsing introns...")
      ## introns <- intronsByTranscript(obj)
      message("Parsing Transcirpts...")
      ## tx <- transcriptsByOverlaps(obj, which)
      tx <- transcripts(obj)
      tx <- subsetByOverlaps(tx, which)
      ## based on tx_id
      txids <- values(tx)$tx_id
      lst <- lapply(txids, function(id){
        id <- as.character(id)
        exons.cur <- exons[[id]]
        values(exons.cur) <- data.frame(tx_id = id, type = "exon")
        ## need to check cds
        cds.cur <- cdss[[id]]
        seqs <- unique(as.character(seqnames(exons.cur)))
        exon_union <- range(exons.cur)
        introns.cur <- GRanges(seqs,IRanges(gaps(ranges(exons.cur))))
        values(introns.cur) <- data.frame(tx_id = id, type = "intron")
        if(!is.null(cds.cur)){
          values(cds.cur) <- data.frame(tx_id = id, type = "cds")
          gaps.cur <- GRanges(seqs,IRanges(gaps(ranges(cds.cur))))
          values(gaps.cur) <- data.frame(tx_id = id, type = "gap")
          cds_union <- range(cds.cur)
          utrs <- setdiff(exon_union, cds_union)
          values(utrs) <- data.frame(tx_id = id, type = "utr")
          gr <- c(exons.cur, introns.cur, cds.cur, gaps.cur, utrs)
        }else{
          utrs <- exons.cur
          values(utrs)$type <- factor("utr")
          gaps.cur <- introns.cur
          values(gaps.cur)$type <- factor("gap")
          ## not sure about gaps
          gr <- c(exons.cur, introns.cur, utrs, gaps.cur)
        }
        gr
      })
      res <- do.call("c", lst)
      res <- keepSeqlevels(res, unique(as.character(seqnames(res))))
    }
    if(type == "single"){
      cds.s <- reduce(res[values(res)$type == "cds"])
      values(cds.s)$type <- factor("cds")
      exon.s <- reduce(res[values(res)$type == "exon"])
      values(exon.s)$type <- factor("exon")
      utr.s <- setdiff(range(exon.s), range(cds.s))
      values(utr.s)$type <- factor("utr")
      gap.s <- gaps(cds.s, start = min(start(cds.s)),
                    end = max(end(cds.s)))
      values(gap.s)$type <- factor("gap")
      res <- c(cds.s, utr.s, gap.s)
    }
    message("Done")
  }
  if(is(obj, "GappedAlignments")){
    require(Rsamtools)
    ## require(Rsamtools)
    if(!missing(which))
      obj <- subsetByOverlaps(obj, which)
    if(!length(obj))
      stop("No entry in the GappedAlignments data")
    if(is.null(names(obj)))
      stop("Missing qname, please make use.name = TRUE, when reading GappedAlignments")

    grg <- grglist(obj)
    grg.u <- unlist(grg)
    ## right way..
    values(grg.u)$qname <- rep(names(grg), times = elementLengths(grg))
    ## values(grg.u)$qname <- names(grg)
    ## message("extracting information...")
    ## values(grg.u)$read.group <- rep(1:length(grg), times = elementLengths(grg))
    values(grg.u)$junction <- rep(ifelse(elementLengths(grg)>1, TRUE,
                                         FALSE),
                                  times = elementLengths(grg))
    ## grl <- split(grg.u, values(grg.u)$qname)
    ## irs <- unlist(range(ranges(grl)))
    ## length(obj)
    ## sort based on junction first
    if(include.level)
      grg.u <- addSteppings(grg.u, group.name = "qname")
      ## values(grg.u)$.levels <- rep(disjointBins(irs), times = elementLengths(grl))
    ## seqlevels(grg.u) <- unique(as.character(seqnames(grg.u)))
    grg.u <- keepSeqlevels(grg.u, unique(as.character(seqnames(grg.u))))
    names(grg.u) <- NULL
    res <- grg.u
  }
  if(is(obj, "BamFile")){
    if(!(type %in% .bamfile.type))
      stop("type for TranscriptDb must be ", .bamfile.type)
    require(Rsamtools)
    if(type == "gapped.pair"){
      message("Read GappedAlignments from BamFile...")
      ga <- readBamGappedAlignments(obj,
                                    param = ScanBamParam(which = which),
                                    use.name = use.name, ...)
      res <- fetch(ga)
    }

    if(type == "raw"){
      message("Read Raw from BamFile by calling scanBam...")
      res <- scanBamGRanges(obj, which, ...)
    }

    if(type == "all"){
      message("Read Raw from BamFile by calling scanBam...")
      res.mb <- scanBamGRanges(obj, which,
                               flag = scanBamFlag(isFirstMateRead = TRUE))
      message("Read GappedAlignments from BamFile...")
      ga <- readBamGappedAlignments(obj,
                                    param = ScanBamParam(which = which),
                                    use.name = use.name, ...)
      res.gp <- fetch(ga)
      message("Combine...")
      nms <- values(res.mb)$qname
      ## fow now just, using isize
      isize <- values(res.mb)$issize
      names(isize) <- nms
      values(res.gp)$isize <- isize[values(res.gp)$qname]
      res <- res.gp
    }
  }
  res
}

scanBamGRanges <- function(file, which, ...){
  require(Rsamtools)
  param <- ScanBamParam(which = which, ...)
  message("scanBam...")
  bam <- scanBam(file, param = param)
  message("Coerce list to GRanges")
  .names <- c("rname", "strand", "pos", "qwidth")
  ##  check
  if(!all(.names %in% Rsamtools:::bamWhat(param)))
    stop("Param must contains rname, strand, pos, qwidth")
  dfnm <- setdiff(Rsamtools:::bamWhat(param), .names)
  ## bam <- bam[[1]]  
  res <- lapply(1:length(bam), function(i){
    bm <- bam[[i]]
    df <- as(bm, "DataFrame")[,dfnm]
    if(length(bm$rname))
      GRanges(bm$rname, IRanges(start = bm$pos, width = bm$qwidth),
              strand = bm$strand, df)
    else
      GRanges()
  })
  bam.gr <- do.call("c", res)
}



