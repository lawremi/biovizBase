## setGeneric("fetch", function(obj, ...) standardGeneric("fetch"))
fetch <- function(obj, which, ..., gene.id,
                  truncate.gaps = FALSE,
                  truncate.fun = NULL,
                  ratio = 0.0025,
                  resize.extra = 10,
                  include.level = TRUE,
                  use.name = TRUE,
                  columns = c("tx_id", "tx_name","gene_id"),
                  type = c("all", "single", "exons.reduce", "exons.all",
                    "exons.geneid",  "gapped.pair", "raw", "all")){
  .txdb.type <- c("all", "single", "exons.reduce",
                  "exons.all", "exons.geneid")
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
    if(is.list(which)){
      message("Parsing exons based on which(list) arguments")
      temp <- exons(obj, vals = which, columns = columns)
      which <- range(temp)
    }    
    ## 1st set all the sequences to be inactive:
    isActiveSeq(obj)[seqlevels(obj)] <- FALSE
    ## Then set only "chr9" to be active:
    seqnms <- as.character(unique(seqnames(which)))
    isActiveSeq(obj)[seqnms] <- TRUE
    ## You can see which are active like this:
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
      message("Parsing transcripts...")
      tx <- transcripts(obj, columns = columns)
      tx <- subsetByOverlaps(tx, which)
      message("Aggregating...")      
      txids <- values(tx)$tx_id
      lst <- lapply(txids, function(id){
        id <- as.character(id)
        tx_nm <- values(tx)[values(tx)$tx_id == id, "tx_name"]
        gene_id <- values(tx)[values(tx)$tx_id == id, "gene_id"]
        gene_id <- paste(unlist(gene_id), sep = ",")
        if(!length(gene_id))
          gene_id <- NA
        exons.cur <- exons[[id]]
        values(exons.cur) <- data.frame(tx_id = id,
                                        tx_name = tx_nm,
                                        gene_id = gene_id,
                                        type = "exon")
        
        cds.cur <- cdss[[id]]
        seqs <- unique(as.character(seqnames(exons.cur)))
        exon_union <- reduce(exons.cur)
        ir.g <- IRanges(gaps(ranges(exons.cur)))
        if(length(ir.g)){
          introns.cur <- GRanges(seqs,ir.g)
          values(introns.cur) <- data.frame(tx_id = id,
                                            tx_name = tx_nm,
                                            gene_id = gene_id,
                                            type = "intron")

        }else{
          introns.cur <- GRanges()
        }
        if(!is.null(cds.cur)){
          values(cds.cur) <- data.frame(tx_id = id,
                                        tx_name = tx_nm,
                                        gene_id = gene_id,
                                        type = "cds")

          cds_union <- reduce(cds.cur)
          utrs <- setdiff(exon_union, cds_union)
          if(length(utrs))
            values(utrs) <- data.frame(tx_id = id,
                                       tx_name = tx_nm,
                                       gene_id = gene_id,
                                       type = "utr")
          else
            utrs <- GRanges()
          ir.g2 <- gaps(reduce(c(ranges(cds_union), ranges(utrs))))
          if(length(ir.g2)){
            gaps.cur <- GRanges(seqs,ir.g2)
          values(gaps.cur) <- data.frame(tx_id = id,
                                         tx_name = tx_nm,
                                        gene_id = gene_id,
                                         type = "gap")
          }else{
            gaps.cur <- GRanges()
          }
          gr <- c(exons.cur, introns.cur, cds.cur, gaps.cur, utrs)
        }else{
          utrs <- exons.cur
          values(utrs)$type <- factor("utr")
          gaps.cur <- introns.cur
          if(length(gaps.cur)){
            values(gaps.cur)$type <- factor("gap")
          }
          gr <- c(exons.cur, introns.cur, utrs, gaps.cur)
        }
        gr
      })
      res <- do.call("c", lst)
      isActiveSeq(obj)[seqlevels(obj)] <- TRUE
      if(length(res))
        res <- keepSeqlevels(res, unique(as.character(seqnames(res))))
      else
        res <- GRanges()
    }
    if(type == "single"){
      if(length(res)){
      cds.s <- reduce(res[values(res)$type == "cds"])
      values(cds.s)$type <- factor("cds")
      exon.s <- reduce(res[values(res)$type == "exon"])
      values(exon.s)$type <- factor("exon")
      utr.s <- setdiff(exon.s, cds.s)
      values(utr.s)$type <- factor("utr")
      gap.s <- gaps(cds.s, start = min(start(cds.s)),
                    end = max(end(cds.s)))
      values(gap.s)$type <- factor("gap")
      res <- c(cds.s, utr.s, gap.s)
    }
    }
    message("Done")
  }
  if(is(obj, "GappedAlignments")){
    ## require(Rsamtools)
    ## require(Rsamtools)
    if(!missing(which))
      obj <- subsetByOverlaps(obj, which)
    if(!length(obj))
      stop("No entry in the GappedAlignments data")
    if(is.null(names(obj)))
      stop("Missing qname, please make use.name = TRUE, when reading GappedAlignments")
    grg <- grglist(obj)
    grg.u <- stack(grg, ".grl.name")
    message("extracting information...")
    values(grg.u)$junction <- rep(ifelse(elementLengths(grg)>1, TRUE, FALSE),
                                  times = elementLengths(grg))
    names(grg.u) <- NULL
    res <- grg.u
  }
  if(is(obj, "BamFile")){
    if(!(type %in% .bamfile.type))
      stop("type for TranscriptDb must be ", .bamfile.type)
    ## require(Rsamtools)
    if(type == "gapped.pair"){
      message("Read GappedAlignments from BamFile...")
      ga <- readBamGappedAlignments(obj,
                                    param = ScanBamParam(which = which),
                                    use.names = use.name, ...)
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
                                    use.names = use.name, ...)
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
  if(truncate.gaps){
    if(is.null(truncate.fun)){
      if("gap" %in% unique(values(res)$type))
        idx <- values(res)$type %in% c("utr", "cds")
      res.s <- reduce(res[idx], ignore.strand = TRUE)
      truncate.fun <- shrinkageFun(gaps(res.s, min(start(res.s)), max(end(res.s))),
                                   maxGap(gaps(res.s, min(start(res.s)), max(end(res.s))),
                                          ratio = ratio))
    }
    res <- truncate.fun(res)
  }
  res
}

scanBamGRanges <- function(file, which, ...){
  ## require(Rsamtools)
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


