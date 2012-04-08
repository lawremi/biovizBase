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

setGeneric("getGaps", function(obj, ...) standardGeneric("getGaps"))
setMethod("getGaps", "GRanges", function(obj, group.name = NULL, facets = NULL){
if(!length(facets))
    facets <- as.formula(~seqnames)
  allvars <- all.vars(as.formula(facets))
  allvars.extra <- allvars[!allvars %in% c(".", "seqnames")]
if(!group.name %in% colnames(values(obj)))
  stop(group.name, " is not in obj")
  grl <- splitByFacets(obj, facets)  
  grl <- lapply(grl, function(dt){
      res <- split(dt, values(dt)[,group.name])
      gps.lst <- lapply(res, function(x){
      if(length(x) > 1){
        gps <- gaps(ranges(x))
        if(length(gps)){
          seqs <- unique(as.character(seqnames(x)))
          ir <- gps
          gr <- GRanges(seqs, ir)
          values(gr)[,"stepping"] <- unique(values(x)[,"stepping"])
          values(gr)[,allvars.extra] <- rep(unique(values(x)[, allvars.extra]),
                                            length(gr))
          
          gr
        }else{
          NULL
        }}else{
          NULL
        }
    })
      idx <- which(!unlist(lapply(gps.lst, is.null)))
      gps <- do.call(c, gps.lst[idx])
  })
  grl <- grl[!unlist(lapply(grl, is.null))]
  ## grl
  if(length(grl)){
    ## res <- unlist(do.call(GRangesList, grl))a
    res <- unlist(do.call(GRangesList, do.call(c, grl)))
    values(res)$type <- "gaps"
    res <- resize(res, width = width(res) + 2, fix = "center")
  }else{
    res <- GRanges()
  }
  res  
})



getIdeoGR <- function(gr){
  if(!is(gr, "GenomicRanges"))
    stop("require GenomicRanges")
  if(all(is.na(seqlengths(gr)))){
    warning("geom(ideogram) need valid seqlengths information for accurate mapping,
                 now use reduced information as ideogram... ")
    res <- reduce(gr, ignore = TRUE)
    start(res) <- 1
    res
  }else{
    res <- as(seqinfo(gr), "GenomicRanges")
  }
  res
}

## get major/mionr scale, y value, major text?
getScale <- function(gr, unit = NULL, n = 100, type = c("M", "B", "sci")){
  type <- match.arg(type)
  if(is.null(unit)){
    unit <- sum(as.numeric(width(gr)))/n
    unit <- 10^floor(log10(unit))
  }
  ## not like normal scale
  res <- split(gr, seqnames(gr))
  grl <- endoapply(res, function(gr){
    st <- 1
    ed <- end(gr)
    major.pos <- seq(st, ed, by = 5*unit)
    minor.pos <- seq(st, ed, by = unit)
    minor.pos <- setdiff(minor.pos, major.pos)
    if(type == "M"){
      texts <- c(paste(as.character(round(major.pos/1e6, digits = 1)), "M",
                       sep = ""),
               rep("", length(minor.pos)))
    }
    if(type == "B"){
      texts <- c(paste(as.character(round(major.pos/1e9, digits = 1)), "B",
                       sep = ""),
               rep("", length(minor.pos)))
    }
    if(type == "sci"){
    texts <- c(as.character(format(major.pos, scientific = TRUE, digit = 1)),
               rep("", length(minor.pos)))
    }
    GRanges(seqnames(gr),
            IRanges(start = c(major.pos, minor.pos),
                    width = 1),
            type = c(rep("major", length(major.pos)), rep("minor", length(minor.pos))),
            scale.y = c(rep(3, length(major.pos)),
                         rep(1, length(minor.pos))),
            text.major = texts)
  })
  res <- unlist(grl)
  seqlengths(res) <- seqlengths(gr)
  res
}

getNR <- function(x, type = c("NUSE", "RLE"),range = 0, ...){
  require(affyPLM)
  compute.nuse <- function(which) {
    nuse <- apply(x@weights[[1]][which, ], 2, sum)
    1/sqrt(nuse)
  }

  type <- match.arg(type)
  model <- x@model.description$modelsettings$model

  if (type == "NUSE") {
    if (x@model.description$R.model$which.parameter.types[3] == 
        1 & x@model.description$R.model$which.parameter.types[1] == 
        0) {
      grp.rma.se1.median <- apply(se(x), 1, median, 
                                  na.rm = TRUE)
      res <- grp.rma.rel.se1.mtx <- sweep(se(x), 1, grp.rma.se1.median, 
                                          FUN = "/")

    }
    else {
      which <- indexProbesProcessed(x)
      ses <- matrix(0, length(which), 4)
      if (x@model.description$R.model$response.variable == 
          1) {
        for (i in 1:length(which)) ses[i, ] <- compute.nuse(which[[i]])
      }
      else {
        stop("Sorry I can't currently impute NUSE values for this PLMset object")
      }
      grp.rma.se1.median <- apply(ses, 1, median)
      res <- grp.rma.rel.se1.mtx <- sweep(ses, 1, grp.rma.se1.median, 
                                          FUN = "/")
      
    }
  }
  if(type == "RLE"){
    if (x@model.description$R.model$which.parameter.types[3] == 
        1) {
      medianchip <- apply(coefs(x), 1, median)
      res <- sweep(coefs(x), 1, medianchip, FUN = "-")
    }
    else {
      stop("It doesn't appear that a model with sample effects was used.")
    }  
  }
  res
}



getFormalNames <- function(..., remove.dots = TRUE){
  res <- lapply(list(...), function(fun){
    if(is.function(fun)){
      res <- names(formals(fun))
      if(remove.dots)
        res <- res[ res != "..."]
      res
    }else{
      stop("arguments passed to getFormalNames must be functions")
    }
  })
  res <- unlist(res)
  res <- res[!duplicated(res)]
  res
}


subsetArgsByFormals <- function(args, ..., remove.dots = TRUE){
  .formals <- getFormalNames(..., remove.dots = remove.dots)
  res <- args[names(args) %in% .formals]
  res
}

flatGrl <- function(object, indName = "grl_name"){
  idx <- togroup(object)
  gr <- stack(object, indName)
  values(gr) <-   cbind(values(gr), values(object)[idx,,drop = FALSE])
  gr
}
