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

## functions downbelow from workshop by Michael Lawrence
## will be implmented by into his package, so I will replace this function later
splicefun <- function(files, txdb, which, id,  xlim, txdb.chr.pre = character(), weight = 1){
  elementGaps <- function(x) {
    x_flat <- unlist(x, use.names = FALSE)
    egaps <- gaps(ranges(x))
    first_segment <- start(PartitioningByWidth(x))
    sn <- seqnames(x_flat)[first_segment][togroup(egaps)]
    strand <- strand(x_flat)[first_segment][togroup(egaps)]
    relist(GRanges(sn, unlist(egaps, use.names = FALSE), 
                   strand, seqlengths = seqlengths(x)), 
           egaps)
  }

  pairReadRanges <- function(reads) {
    pairs <- split(unlist(reads, use.names=FALSE),
                   factor(names(reads)[togroup(reads)], 
                          unique(names(reads))))
    metadata(pairs) <- metadata(reads)
    xs <- values(reads)$XS
    has_xs <- !is.na(xs)
    pair_xs <- setNames(rep.int(NA, length(pairs)), 
                        names(pairs))
    pair_xs[names(reads)[has_xs]] <- xs[has_xs]
    values(pairs)$XS <- unname(pair_xs)
    pairs
  }  

  strandFromXS <- function(pairs) {
    xs <- values(pairs)$XS
    strand <- ifelse(!is.na(xs) & xs != "?", xs, "*")
    strand(pairs) <- relist(Rle(strand, elementLengths(pairs)), 
                            pairs)
    pairs
  }

  gr2key <- function(x) {
    paste(seqnames(x), start(x), end(x), strand(x), 
          sep = ":")
  }

  key2gr <- function(x, ...) {
    key_mat <- matrix(unlist(strsplit(x, ":", fixed=TRUE)), 
                      nrow = 4)
    GRanges(key_mat[1,],
            IRanges(as.integer(key_mat[2,]), 
                    as.integer(key_mat[3,])),
            key_mat[4,], ...)
  }
  findIsoformOverlaps <- function(pairs) {
    splices <- values(pairs)$splices
    hits <- findOverlaps(pairs, tx)
    hit_pairs <- ranges(pairs)[queryHits(hits)]
    hit_splices <- ranges(splices)[queryHits(hits)]
    hit_tx <- ranges(tx)[subjectHits(hits)]
    read_within <- 
      elementLengths(setdiff(hit_pairs, hit_tx)) == 0L
    tx_within <- 
      elementLengths(intersect(hit_tx, hit_splices)) == 0L
    compatible <- read_within & tx_within
    compat_hits <- hits[compatible]
    reads_unique <- tabulate(queryHits(compat_hits), 
                             queryLength(compat_hits)) == 1L
    unique <- logical(length(hits))
    unique[compatible] <- reads_unique[queryHits(compat_hits)]
    strand_specific <- 
      all(strand(pairs) != "*")[queryHits(hits)]
    values(hits) <- DataFrame(strand_specific,
                              compatible,
                              unique)
    list(hits = hits,
         compatible = compatible,
         strand_specific = strand_specific)
  }
  countIsoformHits <- function(hits, compatible, strand_specific) {
    countByTx <- function(x) {
      tabulate(subjectHits(hits)[x], subjectLength(hits))
    }
    compatible_strand <- 
      countByTx(with(values(hits), 
                     compatible & strand_specific))
    counts <- DataFrame(compatible_strand,
                        lapply(values(hits)[-1], countByTx))
    counts
  }

  summarizeSplices <- function(reads) {
    splices <- values(reads)$splices
    splices_flat <- unlist(splices, use.names = FALSE)
    if(length(splices_flat)){    
    splice_table <- table(gr2key(splices_flat))
    splice_summary <- 
      key2gr(names(splice_table), 
             score = as.integer(splice_table),
             novel = !names(splice_table) %in% tx_keys,
             seqlengths = seqlengths(splices))
  }else{
    splice_summary <- GRanges()
  }
    splice_summary
  }

  getUniqueReads <- function(reads, hits) {
    sel <- values(hits)$unique & 
    subjectHits(hits) %in% c(1, 4)
    reads[unique(queryHits(hits)[sel])]
  }

  message("Parsing gene structure from txdb...")
  if(!missing(id)){
  aldoa_gr <- exons(txdb, vals = list(gene_id = id),
                    columns = c("tx_id", "gene_id"))
  aldoa_gr <- keepSeqlevels(aldoa_gr, unique(as.character(seqnames(aldoa_gr))))
  ## FIXME later
  nms <- as.character(names(seqlengths(aldoa_gr)))
  nms.new <- paste(txdb.chr.pre, nms, sep = "")
  names(nms.new) <- nms
  aldoa_gr <- renameSeqlevels(aldoa_gr, nms.new)
  aldoa_range <- range(aldoa_gr)
  wh <- aldoa_range  
  aldoa_vals <- values(aldoa_gr)
  tx <- multisplit(aldoa_gr, aldoa_vals$tx_id)
  tx_to_val <- match(names(tx), unlist(aldoa_vals$tx_id))
  values(tx)$gene_id <- 
    rep(unlist(aldoa_vals$gene_id), 
        elementLengths(aldoa_vals$tx_id))[tx_to_val]
  values(tx)$tx_id <- names(tx)
}else if(!missing(which) & is(which, "GRanges")){
    isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
    seqnms <- as.character(unique(seqnames(which)))
    seqnms <- gsub(txdb.chr.pre, "", seqnms)
    isActiveSeq(txdb)[seqnms] <- TRUE
    ## aldoa_gr <- exons(txdb, columns = c("tx_id", "gene_id"))
    ## aldoa_gr <- keepSeqlevels(aldoa_gr, unique(as.character(seqnames(aldoa_gr))))
    ## nms <- as.character(names(seqlengths(aldoa_gr)))
    ## nms.new <- paste(txdb.chr.pre, nms, sep = "")
    ## names(nms.new) <- nms
    ## aldoa_gr <- renameSeqlevels(aldoa_gr, nms.new)
    ## aldoa_gr <- subsetByOverlaps(aldoa_gr, which)
    ## aldoa_range <- range(aldoa_gr)
    aldoa_range <- which
    exons_grl <- exonsBy(txdb)
    gr.l <- stack(exons_grl, ".sample")
    gr.l <- keepSeqlevels(gr.l, unique(as.character(seqnames(gr.l))))
    nms <- as.character(names(seqlengths(gr.l)))
    nms.new <- paste(txdb.chr.pre, nms, sep = "")
    names(nms.new) <- nms
    gr.l <- renameSeqlevels(gr.l, nms.new)
      exons_grl <- split(gr.l, values(gr.l)$.sample)
     ans <- subsetByOverlaps(exons_grl, which)
    values(ans)$tx_id <- names(ans)
    tx_gr <- transcripts(txdb, columns = c("tx_id", "gene_id"))
    values(ans)$gene_id <- 
      drop(values(tx_gr)$gene_id)[match(names(ans), 
                                        values(tx_gr)$tx_id)]
    tx <- ans
    ## wh <- aldoa_range
    ## aldoa_vals <- values(aldoa_gr)
    ## tx <- multisplit(aldoa_gr, aldoa_vals$tx_id)
    ## tx_to_val <- match(names(tx), unlist(aldoa_vals$tx_id))
    ## values(tx)$gene_id <- 
    ##   rep(unlist(aldoa_vals$gene_id), 
    ##       elementLengths(aldoa_vals$tx_id))[tx_to_val]
    ## values(tx)$tx_id <- names(tx)
}

  message("Parsing bam files")
  bamFiles <- Rsamtools::BamFileList(files)

  message("Analysing...")
  introns <- elementGaps(tx)
  introns_flat <- unlist(introns, use.names = FALSE)
  tx_keys <- gr2key(introns_flat)

  readReadRanges <- function(bam) {
    param <- ScanBamParam(tag = "XS", which = aldoa_range)
    ga <- readGappedAlignments(path(bam), 
                               use.names = TRUE, 
                               param = param)
    reads <- grglist(ga)
    metadata(reads)$bamfile <- bam
    splices <- elementGaps(reads)
    values(splices)$XS <- values(reads)$XS
    pairs <- pairReadRanges(reads)
    pairs <- strandFromXS(pairs)
    splices <- pairReadRanges(splices)
    splices <- strandFromXS(splices)
    values(pairs)$splices <- splices
    pairs
  }


  N <- length(bamFiles)
  nms <- names(bamFiles)
  lst_hits <- lst_counts <- lst_splices <- lst_rr <- list()
  message("parsing hits/counts/splcies")
  for(i in 1:N){
    bf <- bamFiles[[i]]
    rr <- readReadRanges(bf)
    lst_rr <- c(lst_rr, rr)
    hits.lst <- findIsoformOverlaps(rr)
    lst_hits <- c(lst_hits, hits.lst$hits)
    lst_counts <- c(lst_counts, countIsoformHits(hits.lst$hits, hits.lst$compatible, hits.lst$strand_specific))
    lst_splices <- c(lst_splices, summarizeSplices(rr))
  }
  names(lst_hits) <- nms
  names(lst_counts) <- nms
  names(lst_splices) <- nms
  names(lst_rr) <- nms
  ## normal <- readReadRanges(bamFiles)
  ## tumor <- readReadRanges(bamFiles$tumor)
  
  ## normal_hits <- findIsoformOverlaps(normal)
  ## normal_counts <- countIsoformHits(normal_hits)
  ## normal_splices <- summarizeSplices(normal)

  ## tumor_hits <- findIsoformOverlaps(tumor)
  ## tumor_counts <- countIsoformHits(tumor_hits)
  ## tumor_splices <- summarizeSplices(tumor)


###################################################
### code chunk number 39: combine-samples
###################################################
  assays <- do.call(mapply, c(list(cbind), lst_counts, list(SIMPLIFY = FALSE)))
  ## assays <- mapply(cbind, normal_counts, SIMPLIFY = FALSE)
  ## assays <- mapply(cbind, normal_counts, tumor_counts, 
  ##                  SIMPLIFY = FALSE)
  colData <- DataFrame(tumorStatus =  names(bamFiles))
  rownames(colData) <- colData$tumorStatus
  se <- SummarizedExperiment(assays, tx, colData)

###################################################
### code chunk number 40: order-se
###################################################
  uc <- assay(se, "unique")
  uc_ord <- order(rowSums(uc), decreasing = TRUE)
  uc_top <- uc[head(uc_ord, 2),]
  ## fisher.test(uc_top)$estimate

###################################################
### code chunk number 41: get-unique-reads
###################################################
  lst_uniq <- lapply(1:N, function(i){
    getUniqueReads(lst_rr[[i]], lst_hits[[i]])
  })
  names(lst_uniq) <- nms
  ## normal_uniq <- getUniqueReads(normal, normal_hits)
  ## tumor_uniq <- getUniqueReads(tumor, tumor_hits)
  both_uniq <- do.call(mstack, lapply(lst_uniq, unlist))

###################################################
### code chunk number 42: combine-splices
###################################################
  lst_uniq_splices <- lapply(lst_uniq, summarizeSplices)
  names(lst_uniq_splices) <- nms
  uniq_splices <- do.call(mstack, lst_uniq_splices)
  all_splices <- do.call(mstack, lst_splices)
  novel_splices <- all_splices[values(all_splices)$novel]
  ## novel_splices <-  all_splices[values(all_splices)$novel &
  ##                               values(all_splices)$score == 9]
  ## novel_splices <-  all_splices[values(all_splices)$novel &
  ##                               values(all_splices)$score == 9]
  ## uniq_novel_splices <- c(uniq_splices, novel_splices)
  ## uniq_novel_splices <- c(uniq_splices, all_splices)
  uniq_novel_splices <- c(uniq_splices, novel_splices)  
  both_uniq <- keepSeqlevels(both_uniq, unique(as.character(seqnames(both_uniq))))
  .wt <- max(width(uniq_novel_splices))/max(coverage(both_uniq)) * weight
  .wt <- as.numeric(.wt)
  aes.res <- do.call(aes,list(size = substitute(score), 
                              height = substitute(width/.wt, list(.wt = .wt)),
                              color = substitute(novel)))

  message("Constructing splicing graphics....")
  p.novel <- do.call(geom_arch, c(list(data = uniq_novel_splices,
                                    ylab = "coverage",
                                       rect.height = 0), list(aes.res)))
  p.s <- ggplot() + p.novel + do.call(stat_coverage, list(data = both_uniq, facets = name~.))
  message("Constructing gene model....")
  if(length(txdb.chr.pre))
    wh <- GRanges(gsub(txdb.chr.pre, "", as.character(seqnames(wh))),
                  ranges(wh))
  ## browser()
  tx_un <- stack(tx, ".sample")
  tx_cur <- keepSeqlevels(tx_un, unique(as.character(seqnames(tx_un))))
  tx_cur <- tx_cur[, setdiff(colnames(values(tx_cur)), c("tx_id", "gene_id"))]
  tx_cur <- split(tx_cur, values(tx_cur)$.sample)
  tx_track <- do.call(autoplot, list(object = tx_cur, geom = "alignment", ylab = ""))
  ## tx_track <- autoplot(tx_16, geom = "alignment", ylab = "")
  ## tx_track <- do.call(autoplot, list(object = txdb, which = wh))
  if(missing(xlim))
    xlim <- c(start(wh), end(wh))
  tracks(p.s, tx_track, xlim  = xlim, heights = c(3, 1))
}


