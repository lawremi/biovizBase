## rules
## 1. add coord into metadata
## 2. never revise oringal
## 3. treak with labels
setGeneric("transformToGenome", function(data, ...) standardGeneric("transformToGenome"))
setMethod("transformToGenome", "GRanges", function(data, space.skip = 0.1, chr.weight = NULL){
  seqs.l <- seqlengths(data)
  if(all(!is.na(seqs.l))){
    chr.l <- seqs.l
    seqs.suml <- sum(as.numeric(seqs.l))
  }else{
    ## if no seqlengths are found, use data range
    chr.l <- max(end(split(data, seqnames(data))))
    seqs.suml <- sum(as.numeric(chr.l))
  }
  chr.l0 <- chr.l
  N <- length(chr.l)
  .start <- start(data)
  .end <- end(data)
  ## proportion specified for certain chromosomes
  if(!is.null(chr.weight)){
    if(any(!is.numeric(chr.weight)))
      stop("chr.weight must be a vector of numeric values.")
    if(sum(chr.weight) > 1)
      stop("sum of chr.weight must be equal or less than 1.")
    if(length(chr.weight) <= N){
      if(!all(names(chr.weight) %in% names(chr.l)))
        stop("chr.weight names must be a subset of seqinfo chromosome names")
      if(length(chr.weight) < N){
        .chrs <- setdiff(names(chr.l), names(chr.weight))
        .w <- chr.l[.chrs]/sum(as.numeric(chr.l[.chrs])) * (1 - sum(chr.weight))
        names(.w) <- .chrs
        chr.weight <- c(chr.weight, .w)
      }
    }else{
      stop("length of chr.weight cannot exceed length of seqinfo")
    }
    chr.l <- seqs.suml * chr.weight
    .start <- .start/chr.l0[as.character(seqnames(data))] * chr.l[as.character(seqnames(data))]
    .end <- .end/chr.l0[as.character(seqnames(data))] * chr.l[as.character(seqnames(data))]
  }  
  chr.l.back <- chr.l
  space.skip <- space.skip * seqs.suml
  skps <- space.skip * ((1:length(seqs.l))-1)
  names(skps) <- names(seqlengths(data))
  nms <- names(chr.l)
  chr.l <- cumsum(as.numeric(chr.l))
  chr.l2 <- c(0, chr.l[-length(chr.l)])
  .breaks <- chr.l2 + chr.l.back/2  + skps
  names(.breaks) <- nms
  names(chr.l2) <- nms
  sts.new <- .start + skps[as.character(seqnames(data))] +
    chr.l2[as.character(seqnames(data))] 
  ed.new <- .end + skps[as.character(seqnames(data))] +
    chr.l2[as.character(seqnames(data))] 
  max.chr <- rev(names(chr.l2))[1]
  x.max <- chr.l.back[max.chr] + skps[max.chr] + chr.l2[max.chr] + space.skip
  names(x.max) <- "genome"
  metadata(data)$breaks <- .breaks
  metadata(data)$labels <- names(.breaks)
  values(data)$.start <- round(sts.new)
  values(data)$.end <- round(ed.new)
  metadata(data)$space.skip <- space.skip
  metadata(data)$coord <- "genome"
  metadata(data)$x.max <- x.max
  data
})


setMethod("transformToGenome", "GRangesList", function(data, space.skip = 0.1, chr.weight = NULL){
  obj <- stack(data, ".group")
  obj <- transformToGenome(obj, space.skip = space.skip, chr.weight = NULL)
  obj
})
  
is_coord_genome <- function(data){
  if(!is.null(metadata(data)$coord)){  
   res <- metadata(data)$coord == "genome"
  }else{
    res <- FALSE
  }
  res
}

rescaleGr <- function(gr, maxSize = 1e8){
  if(!is_coord_genome(gr))
    stop("must pass a GRanges of coord genome")
  x.max <- metadata(gr)$x.max
  if(x.max > maxSize){
    values(gr)$.rescale.idx <- 1:length(gr)    
    gr <- gr[order(values(gr)$.start),]    
    .breaks <- metadata(gr)$breaks
    metadata(gr)$rescale <- TRUE
    .v <- rescale(c(1, .breaks, values(gr)$.start, values(gr)$.end), to = c(1, maxSize))
    .breaks <- .v[2:(1 + length(.breaks))]
    .v <- .v[-(1:(1+length(.breaks)))]
    N <- length(.v)
    res <- GRanges(seqnames(gr), IRanges(start  = .v[1:(N/2)],
                                     end = .v[(N/2+1):N]))
    values(res) <- values(gr)
    ## res <- sort(res)
    metadata(res) <- metadata(gr)
    metadata(res)$breaks <- .breaks
    res <- res[order(values(res)$.rescale.idx)]
  }else{
    res <- GRanges(seqnames(gr), IRanges(start  = values(gr)$.start,
                                     end = values(gr)$.end))
    values(res) <- values(gr)
    metadata(res) <- metadata(gr)
  }
  res
}

## then need a transformation to circlular view
## data is a GRanges dataect
transformToCircle <- function(data, x = NULL, y = NULL, ylim = NULL,
                      radius = 10, trackWidth =10, direction = c("clockwise", "anticlockwise"),
                      mul = 0.05){

  if(!".midpoint" %in% colnames(values(data)))
    values(data)$.midpoint <- (values(data)$.start + values(data)$.end)/2
  if(!length(y)){
    values(data)$.circle.equal.y <- 1
    y <- as.character(".circle.equal.y")
    ## stop("y is missing")
  }else if(!is.character(y) && !is.numeric(y)){
    stop("y must be character or numeric")
  }
  if(!length(x)){
    x <- ".midpoint"
  }else{
    if(x == "midpoint" & (!"midpoint" %in% colnames(values(data))))
      x <- ".midpoint"
  }

  if(length(ylim)){
    idx <- values(data)[,y] <= max(ylim) & values(data)[,y] >= min(ylim)
    data <- data[idx]
  }
  
  temp.x <- values(data)[,x] 
  
  if(is.character(y))
    temp.y <- values(data)[,y]
  if(is.numeric(y)){
    temp.y <- y
  }
  

  if(length(ylim))
    temp.y.r <- expand_range(ylim, mul = mul)
  else
    temp.y.r <- expand_range(range(temp.y), mul = mul)
  temp.y <- (temp.y - min(temp.y.r))/diff(temp.y.r) * trackWidth + radius
  
  ## always from 1 to max
  direction <- match.arg(direction)
  if(direction == "clockwise")
    dir <- 1
  else
    dir <- -1
  ## need to use max of "ideo"
  if(is_coord_genome(data)){
    x.max <- metadata(data)$x.max
  }else{
    x.max <- max(end(data))
  }
  angle.unit <- pi* 2/ x.max
  x.new <- temp.y * sin(angle.unit * (temp.x - 1) * dir)
  y.new <- temp.y * cos(angle.unit * (temp.x - 1)* dir)
  values(data)$.circle.x <- x.new
  values(data)$.circle.y <- y.new
  values(data)$.circle.angle <- angle.unit * (temp.x-1) * dir
  data
}

## interpolate segment or rectangle
## rectangle only use stepping as y

transformToRectInCircle <- function(data, y = NULL, space.skip = 0.1,
                                    trackWidth = 10, radius = 10,
                                    direction = c("clockwise", "anticlockwise"),
                                    n = 100, mul = 0.05,
                                    chr.weight = NULL){
  if(is_coord_genome(data))
    x.max <- values(data)$x.max
  else
    x.max <- max(end(data))
  x.unit <- x.max/n
  direction <- match.arg(direction)
  data.back <- data

    ## need to consider the space
  if(length(y)){
    temp.y <- as.character(values(data)[,as.character(y)])
    if(is.character(temp.y)){
      temp.y <- as.numeric(as.factor((temp.y)))
    }
    if(!all(check.integer(temp.y)))
      stop("geom(bar) require specified 'y', must be integer which indicates the levels")
    values(data)$stepping <- temp.y
  }else{
      values(data)$stepping <- disjointBins(data, ignore.strand = TRUE)
  }
  names(data) <- NULL
  df <- as.data.frame(data)

  lst <- lapply(1:nrow(df), function(i){
    if(df[i,"width"] > 1){
      .n <- round(df[i, "width"]/x.unit)
      if(.n< 5) .n <- 5
      inter.fun <- function(x, y) approx(x, y, n = .n)
      res.x <- as.integer(inter.fun(c(df[i,"start"], df[i,"end"]), c(0, 0))$x)
    }else{
      res.x <- c(df[i, "start"], df[i, "start"])
    }
    res <- data.frame(.biovizBase.new.x = c(res.x, rev(res.x)))
    N <- nrow(res)
    res$.biovizBase.group <- i
    res$.biovizBase.level <- c(rep(df[i, "stepping"] - 0.4, N/2),
                               rep(df[i, "stepping"] + 0.4, N/2))
    res <- suppressWarnings(cbind(res, df[i,]))
    res$.int.id <- seq_len(nrow(res))
    res
  })
  res <- do.call("rbind", lst)
  .y <- ".biovizBase.level"
  res.gr <- GRanges(res$seqnames, IRanges(start = res$.biovizBase.new.x,
                                          width = 1),
                    strand = res$strand)
  values(res.gr) <- subset(res, select = -c(start, end, width, strand,seqnames))
  seqlengths(res.gr) <- seqlengths(data)
  res <- transformToGenome(res.gr, space.skip = space.skip, chr.weight = chr.weight)
  res <- transformToCircle(res, y = .y, radius = radius, trackWidth = trackWidth,
                   direction = direction, mul = mul)
  res
}

transformToBarInCircle <- function(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                     direction = c("clockwise", "anticlockwise"),
                                   n = 100, mul = 0.05, chr.weight = NULL){
  
  if(is_coord_genome(data))
    x.max <- values(data)$x.max
  else
    x.max <- max(end(data))
  x.unit <- x.max/n
  direction <- match.arg(direction)
  data.back <- data
    ## need to consider the space
  if(length(y)){
    values(data)$stepping <- values(data)[, as.character(y)]
  }else{
    if("score" %in% colnames(values(data))){
      message("use 'score' for default y, or you can specify it in aes()")
     values(data)$stepping <- values(data)$score
    }else{
      stop("please provide aes(y = ), or use 'score' as column")
    }
  }
  names(data) <- NULL
  df <- as.data.frame(data)
  lst <- lapply(1:nrow(df), function(i){
    if(df[i,"width"] > 1){
      .n <- round(df[i, "width"]/x.unit)
      if(.n< 5) .n <- 5
      inter.fun <- function(x, y) approx(x, y, n = .n)
      res.x <- as.integer(inter.fun(c(df[i,"start"], df[i,"end"]), c(0, 0))$x)
    }else{
      res.x <- c(df[i, "start"], df[i, "start"])
    }
    res <- data.frame(.biovizBase.new.x = c(res.x, rev(res.x)))
    N <- nrow(res)
    res$.biovizBase.group <- i
    res$.biovizBase.new <- c(rep(0, N/2),
                             rep(df[i, "stepping"], N/2))
    res <- suppressWarnings(cbind(res, df[i,]))
    res$.int.id <- seq_len(nrow(res))
    res
  })
  res <- do.call("rbind", lst)
  .y <- ".biovizBase.new"
  res.gr <- GRanges(res$seqnames, IRanges(start = res$.biovizBase.new.x,
                                          width = 1),
                    strand = res$strand)
  values(res.gr) <- subset(res, select = -c(start, end, width, strand,seqnames))
  seqlengths(res.gr) <- seqlengths(data)
  res <- transformToGenome(res.gr, space.skip = space.skip, chr.weight = chr.weight)
  res <- transformToCircle(res, y = .y, radius = radius, trackWidth = trackWidth,
                   direction = direction, mul = mul)
  res
}

## ok, segment allow user to use a flexible y
transformToSegInCircle <- function(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                      direction = c("clockwise", "anticlockwise"),
                      n = 100, chr.weight = NULL){
  if(is_coord_genome(data))
    x.max <- values(data)$x.max
  else
    x.max <- max(end(data))
  x.unit <- x.max/n

  direction <- match.arg(direction)
  ## do the linear interpolatoin first

  if(!length(y)){
    values(data)$stepping <- disjointBins(ranges(data))
  }
  names(data) <- NULL
  df <- as.data.frame(data)
  lst <- lapply(1:nrow(df), function(i){
    if(df[i,"width"] > 1){
      .n <- round(df[i, "width"]/x.unit)
      if(.n< 5) .n <- 5
      inter.fun <- function(x, y) approx(x, y, n = .n)
      res.x <- as.integer(inter.fun(c(df[i,"start"], df[i,"end"]), c(0, 0))$x)
    }else{
      res.x <- c(df[i, "start"], df[i, "start"])
    }
    res <- data.frame(.biovizBase.new.x = res.x)
    res$.biovizBase.group <- i
    if(!length(y)){          
      res$.biovizBase.level <- df[i, "stepping"]
    }
    N <- nrow(res)
    df.extra <- do.call("rbind", lapply(1:N, function(k) df[i,]))
    res <- cbind(res, df.extra)
  })
  res <- do.call("rbind", lst)
  if(!length(y))
    .y <- ".biovizBase.level"
  else
    .y <- as.character(y)

  res.gr <- GRanges(res$seqnames, IRanges(start = res$.biovizBase.new.x,
                                          width = 1),
                    strand = res$strand)
  values(res.gr) <- subset(res, select = -c(start, end, width, strand,seqnames))
  seqlengths(res.gr) <- seqlengths(data)
  res <- transformToGenome(res.gr, space.skip, chr.weight = chr.weight)
  res <- transformToCircle(res, y = .y, radius = radius, trackWidth =trackWidth,
                   direction = direction)
  res
}


transformToSegInCircle2 <- function(data, y = NULL, space.skip = 0.1,
                                    trackWidth = 10, radius = 10,
                                    direction = c("clockwise", "anticlockwise"),
                                    n = 100, mul = 0.05, chr.weight = NULL){

  if(is_coord_genome(data))
    x.max <- values(data)$x.max
  else
    x.max <- max(end(data))
  x.unit <- x.max/n
  direction <- match.arg(direction)
  data.back <- data

    ## need to consider the space
  if(length(y)){
    temp.y <- as.character(values(data)[,as.character(y)])
    if(is.character(temp.y)){
      temp.y <- as.numeric(as.factor((temp.y)))
    }
    if(!all(check.integer(temp.y)))
      stop("geom(bar) require specified 'y', must be integer which indicates the levels")
    values(data)$stepping <- temp.y
  }else{
      values(data)$stepping <- disjointBins(data, ignore.strand = TRUE)
  }
  names(data) <- NULL
  df <- as.data.frame(data)

  lst <- lapply(1:nrow(df), function(i){
    if(df[i,"width"] > 1){
      res.x <- c((df[i, "start"] + df[i, "end"])/2, (df[i, "start"] + df[i, "end"])/2)      
    }else{
      res.x <- c(df[i, "start"], df[i, "start"])
    }
    res <- data.frame(.biovizBase.new.x = res.x)
    N <- nrow(res)
    res$.biovizBase.group <- i
    res$.biovizBase.level <- c(df[i, "stepping"] - 0.4,
                               df[i, "stepping"] + 0.4)
    res <- suppressWarnings(cbind(res, df[i,]))
    res$.int.id <- seq_len(nrow(res))
    res
  })
  res <- do.call("rbind", lst)
  .y <- ".biovizBase.level"
  res.gr <- GRanges(res$seqnames, IRanges(start = res$.biovizBase.new.x,
                                          width = 1),
                    strand = res$strand)
  values(res.gr) <- subset(res, select = -c(start, end, width, strand,seqnames))
  seqlengths(res.gr) <- seqlengths(data)
  res <- transformToGenome(res.gr, space.skip = space.skip, chr.weight = chr.weight)
  res <- transformToCircle(res, y = .y, radius = radius, trackWidth = trackWidth,
                   direction = direction, mul = mul)
  res  
}

## for a special GRanges
transformToLinkInCircle <- function(data, linked.to, space.skip = 0.1, trackWidth = 10,
                                    radius = 10, link.fun = function(x, y, n = 100)
                                    bezier(x, y, evaluation = n),
                                    direction = c("clockwise", "anticlockwise"),
                                    chr.weight = NULL){

  direction <- match.arg(direction)
  if(missing(linked.to))
    stop("linked.to must be provided and be a GRanges")
  ## trace id
  N <- length(data)
  values(data)$.biovizBase.idx <- 1:N
  values(values(data)[, linked.to])$.biovizBase.idx <- 1:N
  obj <- transformToGenome(data, space.skip, chr.weight = chr.weight)
  obj <- transformToCircle(obj, y = 0, radius = radius, trackWidth =trackWidth,
                   direction = direction)
  obj <- obj[order(values(obj)$.biovizBase.idx)]
  names(obj) <- NULL
  df <- as.data.frame(obj)
  linktodata <- values(data)[,linked.to]
  ## missing y
  linktodata <- transformToGenome(linktodata, space.skip, chr.weight = chr.weight)
  ## keep order
  linktodata <- transformToCircle(linktodata, radius = radius,
                          y = 0,
                          trackWidth =trackWidth,
                          direction = direction)
  linktodata <- linktodata[order(values(linktodata)$.biovizBase.idx)]  
  linkdf <- as.data.frame(linktodata)
  extra.df <- subset(df, select = -c(.circle.x, .circle.y))
  linkdf2 <- data.frame(from.x = df$.circle.x,
                        from.y = df$.circle.y,
                        to.x = linkdf$.circle.x,
                        to.y = linkdf$.circle.y)
  linkdf2 <- cbind(linkdf2, extra.df)
  lst <- lapply(1:nrow(linkdf2), function(i){
    ndf <- link.fun(c(linkdf2[i, "from.x"],0,linkdf2[i,"to.x"]),
                    c(linkdf2[i, "from.y"],0,linkdf2[i,"to.y"]))
    ndf <- as.data.frame(do.call("cbind", ndf))
    N <- nrow(ndf)
    ndf$.biovizBase.group <- i
    ## need to past other meta data
    extra.df <- subset(linkdf2[i,], select = -c(from.x, from.y, to.x, to.y))
    extra.df2 <- do.call("rbind", lapply(1:N, function(i) extra.df))
    ndf <- cbind(ndf, extra.df2)
  })
  res <- do.call("rbind", lst)
  res
}


## ddplyFun <- function(data, .fun, ..., window = 1e5){
##   grl <- split(data, seqnames(data))
##   idx <- elementNROWS(grl) > 0
##   grl <- endoapply(grl[idx], function(gr){
##     ## let's make it more works for multiple seqnames
##     bks <- seq(1, max(end(gr)), by = window)
##     N <- length(bks)
##     seqn <- unique(seqnames(gr))
##     if(length(seq) > 1)
##       stop("Multiple seqnames found")
##     if(!N%%2){
##       qgr <- GRanges(seqn, IRanges(start = bks[1:(N-1)],
##                                    end = bks[2:N]))
##     }else{
##       bks <- c(bks, max(end(gr)))
##       N <- length(bks)
##       qgr <- GRanges(seqn, IRanges(start = bks[1:(N-1)],
##                                    end = bks[2:N]))
##     }
##     if(max(end(qgr)) < max(end(gr))){
##       qgr <- c(qgr, GRanges(seqn, IRanges(max(end(qgr))+ 1,
##                                           max(end(gr)))))
##     }
##     of <- findOverlaps(gr, qgr, select = "first")
##     idx <- !is.na(of)
##     df <- as.data.frame(gr[idx])
##     df$of <- of[idx]
##     resdf <- plyr::ddply(df, .(of), .fun, ...)
##     idx <- as.numeric(resdf$of)
##     ## idx <- idx[!is.na(idx)]
##     res <- qgr[idx]
##     values(res) <- as.data.frame(resdf[,2])
##     res
##   })
##   res <- unlist(grl[!is.null(grl)])
##   names(res) <- NULL
##   sort(res)
## }

setGeneric("transformToDf", function(data, ...) standardGeneric("transformToDf"))

## setMethod("transformToDf", c("eSet"), function(model, data){
##   object <- model
##   pdata <- phenoData(object)
##   df <- as(object, "data.frame")
##   df$sampleNames <- row.names(df)
##   df.m <- melt(df, id.vars =  c(varLabels(pdata), "sampleNames"))
##   df.m
## })

setMethod("transformToDf", c("GRanges"), function(data){
  vals <- values(data)
  if(length(vals)){
    idx <- !unlist(lapply(vals@listData, function(x) is(x, "List") & !is(x, "DNAStringSet")))
    if(sum(!idx))
      warning(colnames(vals)[!idx], " column has been dropped for the reason that the corecion of class List is not supported")
    df <- as.data.frame(data[,idx])
  }else{
    df <- as.data.frame(data)
  }
  df$midpoint <- (df$start+df$end)/2
  df
})



transformDfToGr <- function(data, seqnames = NULL, start = NULL, end = NULL,
                            width = NULL, strand = NULL,
                            to.seqnames = NULL, to.start = NULL, to.end = NULL,
                            to.width = NULL, to.strand = NULL,
                            linked.to = "to.gr"){
  
  ## this losing seqinfo....
  colnms <- colnames(data)
  if(is.null(seqnames)){
    if("seqnames" %in% colnms)
      seqnames <- "seqnames"
    else
      stop("Please sepicify which column represent the seqnames")
  }

  seqnames <- data[, seqnames]
  if(length(c(start, end, width)) <2){
    stop("Pleaase sepeicfy two of start/end/width")
  }

  if(is.null(end)){
    if("end" %in% colnms)
      end <- data[,"end"]
  }else{
    end <- data[, end]
  }

  if(is.null(start)){
    if("start" %in% colnms)
      start <- data[,"start"]
  }else{
    start <- data[, start]
  }

  if(is.null(width)){
    if("width" %in% colnms)
      width <- data[,"width"]
  }else{
    width <- data[, width]
  }

  gr <- GRanges(as.character(seqnames), IRanges(start = start, width = width, end = end))
  if(is.null(strand)){
    if("strand" %in% colnms){
      strand <- data[, "strand"]
      strand(gr) <- strand
    }}else{
      strand(gr) <- data[, strand]
    }

  values(gr) <- data[,!colnames(data) %in% c(start, end, width, seqnames, strand)]
  ## linked ranges
  if(!is.null(to.seqnames)){
    message("try to creating a linked GRanges as column to.gr")
    if(length(c(to.start, to.end, to.width)) > 1){
      if(is.null(to.end)){
        if("to.end" %in% colnms)
          to.end <- data[,"to.end"]
      }else{
        to.end <- data[, to.end]
      }

      if(is.null(to.start)){
        if("to.start" %in% colnms)
          to.start <- data[,"to.start"]
      }else{
        to.start <- data[, to.start]
      }

      if(is.null(width)){
        if("to.width" %in% colnms)
          to.width <- data[,"to.width"]
      }else{
        to.width <- data[, to.width]
      }
      to.gr <- GRanges(seqnames = as.character(data[,to.seqnames]),
                       IRanges(start = to.start,
                               end = to.end,
                               width = to.width))
      if(!is.null(to.strand))
        strand(to.gr) <- data[,to.strand]
      values(gr)[,linked.to] <- to.gr
    }
  }
  gr
  
}

check.integer <- function(x){
  sapply(x, function(x)
    !length(grep("[^[:digit:]]", as.character(x))))
}


## retrun a new GRanges, with new coord
transformGRangesForEvenSpace <- function(gr){
  st <- start(range(gr))
  ed <- end(range(gr))
  wd <- width(range(gr))
  N <- length(gr)
  wid <- wd/N
  x.new <- st + wid/2 + (1:N - 1) * wid
  values(gr)$x.new <- x.new
  metadata(gr)$coord <- "even"
  gr
}

transformToEven <- function(data){
  st <- start(range(data))
  ed <- end(range(data))
  wd <- width(range(data))
  N <- length(data)
  wid <- wd/N
  x.new <- st + wid/2 + (1:N - 1) * wid
  gr <- GRanges(start = x.new, width = wid)
  values(gr) <- values(data)
  ## values(gr)$x.new <- x.new
  values(gr)$.ori <- data
  metadata(gr)$coord <- "even"
  data
}

## setGeneric("transform")
## setMethod("transform", "GRanges", function(_data, space.skip = 0.1,
##                                            maxSize = 1e9,
##                                            max.gap = 1L,
##                                            method = c("genome",
##                                              "truncate_gap",
##                                              "even")){
##   method <- match.arg(method)
##   obj <- switch(method,
##                 genome = {
##                   transformToGenome(_data, space.skip = space.skip)
##                 },
##                 truncate_gaps = {
                  
##                 },
##                 even = {
##                   transformToEven(_data, space.skip = space.skip)
##                 })
## }

## autoplot(gr, layout = "circle", geom = "ideo", rect.inter.n = 10, space.skip = 0.2,
##          aes(fill = .ori.seqnames))
## autoplot(gr, coord = "genome")
## autoplot(gr, layout = "circle",space.skip = 0.2, geom = "bar", grid = TRUE)
## gr2 <- transformToGenome(gr, space.skip = 0)
## gr2
## metadata(gr2)$coord <- NULL
## autoplot(gr2)

setGeneric("transformToArch", function(data, ...) standardGeneric("transformToArch"))
setMethod("transformToArch", "GRanges", function(data, width = 1){
  grl <- split(data, 1:length(data))
  grl <- endoapply(grl, function(gr){
    sts <- c(start(gr)-width, end(gr))
    res <- GRanges(seqnames(gr), IRanges(sts, width = width))
    values(res) <- values(gr)[c(1, 1),]
    seqlengths(res) <- seqlengths(gr)
    res
  })
  grl
})

setMethod("transformToArch", "GRangesList", function(data, width = 1){
  coord <- metadata(data)$coord
  res <- unlist(endoapply(data, function(gr){
    gps <- gaps(gr, start = min(start(gr)), end = max(end(gr)))
    gps <- gps[strand(gps) == "*"]
    values(gps) <- values(gr[1,])
    gps
  }))
  metadata(res)$coord <- coord
  res
})

## ## coord
## library(GenomicRanges)
## gr1 <- GRanges("chr1", IRanges(start = 1:10, width = 5))
## gr2 <- GRanges("chr1", IRanges(start = 100:110, width = 5))
## gr3 <- GRanges("chr1", IRanges(start = 200:210, width = 5))
## grl <- GenomicRangesList(gr1, gr2, gr3)


## setGeneric("transformToOrigin", function(data, ...) standardGeneric("transformToOrigin"))
## setMethod("transformToOrigin", "GRanges", function(data, ...){
##   ori <- metadata(gr)$origin
##   if(is.null(ori))
##     ori <- 1
##   values(gr)$.ori <- data
##   metadata(gr)$coord <- "even"
## })

## setMethod("transformToOrigin", "GenomicRangesList", function(data, ...){
  
##   values(gr)$.ori <- data
##   metadata(gr)$coord <- "even"
## })


