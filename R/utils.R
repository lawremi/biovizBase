containLetters <- function(obj, all = FALSE){
  obj <- as.character(obj)
  obj <- tolower(obj)
  obj <- unlist(strsplit(obj,""))
  if(!all){
    res <- any(obj%in%letters)
  }else{
    res <- all(obj%in%letters)
  }
  res
}

GCcontent <- function(obj, ..., view.width, as.prob = TRUE){
  require(BSgenome)
  seqs <- getSeq(obj, ..., as.character = FALSE)
  if(missing(view.width)){
    res <- letterFrequency(seqs, letters="CG", as.prob = as.prob)
  }else{
    res <- lapply(seqs, function(x) letterFrequencyInSlidingView(x, view.width,
                                              letters="CG", as.prob = as.prob))
    if(length(res) == 1)
      res <- unlist(res)
  }
  res
}


centroidPos <- function(obj){
  if(!isIdeogram(obj))
    stop("Please pass an ideogram, please check function getIdeogram")
  obj.s <- split(obj, as.character(seqnames(obj)))
  lst <- lapply(obj.s,function(x){
    arms <- substr(values(x)$name, 1,1)
    idx <- which(arms == "q")[1]
    data.frame(x = start(x)[idx], y = 5, seqnames = unique(as.character(seqnames(x))))
  })
  df <- do.call(rbind, lst)
  df
}

## ## utils to generate pair-end
pspanGR <- function(file, region, sameChr = TRUE, isize.cutoff = 170){
  ## FIXME: move unmated?
  bam <- scanBam(file, param=ScanBamParam(which = region),
                       flag = scanBamFlag(hasUnmappedMate = FALSE))
  bam <- bam[[1]]
  bamrd <- GRanges(bam$rname, IRanges(bam$pos, width = bam$qwidth),
                      strand = bam$strand,
                      mseqname = bam$mrnm,
                      mstart = bam$mpos,
                      isize = bam$isize)
  ## why negative?sometime
  bamrd <- bamrd[abs(bam$isize) >= isize.cutoff]
  if(sameChr){
    idx <- as.character(seqnames(bamrd)) == values(bamrd)$mseqname 
    bamrd <- bamrd[idx]
  }
  if(length(bamrd)){
  p1 <- GRanges(seqnames(bamrd),
                ranges(bamrd))
  p2 <- GRanges(values(bamrd)$mseqname,
                IRanges(values(bamrd)$mstart, width = 75))
  pspan <- punion(p1, p2, fill.gap = TRUE)
  pgaps <- pgap(ranges(p1), ranges(p2))
  return(list(pspan = pspan, pgaps = pgaps, p1 = p1, p2 = p2))
}else{
  return(NULL)
}
}


## retrun a new GRanges, with new coord
transformGRangesForEvenSpace <- function(gr){
  ## make sure it's reduced one
  ## TODO: release a error
  ## if(anyDuplicated(matchMatrix(findOverlaps(ranges(gr)))[,1]))
  ##   stop("Your GRanges object should not have overlaped intervals")
  ## get whole coordinate
  st <- start(range(gr))
  ed <- end(range(gr))
  wd <- width(range(gr))
  ## we need to find midpoint first
  N <- length(gr)
  wid <- wd/N
  x.new <- st + wid/2 + (1:N - 1) * wid
  values(gr)$x.new <- x.new
  gr  
}

transformGRangesToDfWithTicks <- function(gr, fixed.length = NULL){
  if(is.null(fixed.length)){
    ## fixed.length <- getSpaceByData(gr)
    fixed.length <- getSpaceBySeqlengths(gr)
  }else{
    if(is.null(names(fixed.length)))
      stop("Please name your fixed.length with seqnames")
  }
  ## chr.order <- unique(as.character(seqnames(gr.snp)))
  chr.order <- names(seqlengths(gr))
  sft <- getShiftSpace(fixed.length[chr.order])
  df <- as.data.frame(gr)
  df$midpoint <- (df$start + df$end)/2
  ticks.pos <- fixed.length/2
  ## ticks.pos <- data.frame(pos = ticks.pos, seqnames = )
  df <- shiftDfBySpace(df, sft)
  df$seqnames <- factor(as.character(df$seqnames),
                        levels = chr.order)
  ticks <- sft[names(ticks.pos)] + ticks.pos
  list(df = df, ticks = ticks)
}

## compute default length, data-wide
getSpaceByData <- function(gr){
  grl <- split(gr, seqnames(gr))
  ## ifseqlengths(gr)
  res <- unlist(width(range(grl)))
  res
}

getSpaceBySeqlengths <- function(gr){
  res <- seqlengths(gr)
  if(any(is.na(res))){
    idx <- is.na(res)
    res <- getSpaceByData(gr)[names(res[idx])]
  }
  res
}

getShiftSpace <- function(val){
  N <- length(val)
  ## avoid integer overflow, covert to numeric
  res <- c(0, cumsum(as.numeric(val[-N])))
  res <- res + 1
  names(res) <- names(val)
  res
}

shiftDfBySpace <- function(df, space){
  ## GRanges need to be coerced to df
  ## overcome the integer limits of ranges
  ## make it vectorized
  chrs <- as.character(df$seqnames)
  shifts <- space[chrs]
  df$start <- df$start + shifts
  df$end <- df$end + shifts
  if("midpoint" %in% colnames(df))
    df$midpoint <- df$midpoint + shifts
  df
}


genSymbols <- function(org){
  ## TODO, add gene id
  if(!is(org, "OrgDb")){
    stop("org must be a OrgDb object")
  }
  pkg.name <- paste("package", org$packageName, sep = ":")
  sym.name <- grep("SYMBOL$", ls(pkg.name), value = TRUE)
  chr.name <- grep("CHR$", ls(pkg.name), value = TRUE)
  chrloc.name <- grep("CHRLOC$", ls(pkg.name), value = TRUE)
  chrlocend.name <- grep("CHRLOCEND$", ls(pkg.name), value = TRUE)
  en.name <- grep("ENSEMBL$", ls(pkg.name), value = TRUE)
  ## a little hack
  sym.obj <- eval(as.name(sym.name))
  chr.obj <- eval(as.name(chr.name))
  chrloc.obj <- eval(as.name(chrloc.name))
  chrlocend.obj <- eval(as.name(chrlocend.name))
  en.obj <- eval(as.name(en.name))
  symbol_id <- mappedRkeys(sym.obj)
  message("get genen id...")
  gene_id <- mget(symbol_id, revmap(sym.obj), ifnotfound=NA)
  gene_id <- unlist2(gene_id)
  gene_id <- gene_id[!is.na(gene_id)]
  message("get chromosome information...")
  chr <- toTable(chr.obj[gene_id])
  chrloc <- toTable(chrloc.obj[gene_id])
  chrlocend <- toTable(chrlocend.obj[gene_id])
  chrs <- chrloc
  chrs$end_location <- chrlocend$end_location
  res <- chrs
  if(sum(duplicated(gene_id)))
    stop("Duplicated gene_id")# no duplicate
  idx <- match(res$gene_id, gene_id)
  res$symbol <- names(gene_id[idx])
  if(sum(with(res, start_location > 0 & end_location < 0)))
    stop("Wierd pattern found")
  if(sum(with(res, start_location < 0 & end_location > 0))) 
   stop("Wierd pattern found")
  message("get strand information...")
  res$sense <- with(res, start_location>=0&end_location>=0)
  res.sense <- res[res$sense,]
  res.antisense <- res[!res$sense,]
  res$strand <- strand(!res$sense)
  message("get EMSEMBL ID...")
  ensembl_id <- unlist(mget(gene_id, en.obj))
  idx <- match(res$gene_id, names(ensembl_id))
  res$ensembl_id <- ensembl_id[idx]
  message("generating a GRanges object for genen symbols...")
  genesymbol <- with(res,GRanges(seqnames = paste("chr",Chromosome,sep = ""),
                                 IRanges(abs(start_location), abs(end_location)),
                                 strand = strand,
                                 symbol = symbol,
                                 ensembl_id = ensembl_id))
  names(genesymbol) <- values(genesymbol)$symbol
  message("Done")
  genesymbol
}
