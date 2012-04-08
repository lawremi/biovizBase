transformToGenome <- function(data, space.skip = 0.1){
  data <- sort(data)
  seqs.l <- seqlengths(data)
  if(all(!is.na(seqs.l))){
    chr.l <- seqs.l
    seqs.suml <- sum(as.numeric(seqs.l))
  }else{
    ## if no seqlengths are found, use data range
    chr.l <- max(end(split(data, as.character(seqnames(data)))))
    seqs.suml <- sum(as.numeric(chr.l))
  }
  chr.l.back <- chr.l
  space.skip <- space.skip * seqs.suml
  skps <- space.skip * ((1:length(seqs.l)))
  names(skps) <- names(seqlengths(data))
  nms <- names(chr.l)
  chr.l <- cumsum(as.numeric(chr.l))
  chr.l2 <- c(0, chr.l[-length(chr.l)])
  names(chr.l2) <- nms
  
  sts.new <- start(data) + skps[as.character(seqnames(data))] +
    chr.l2[as.character(seqnames(data))]
  ed.new <- end(data) + skps[as.character(seqnames(data))] +
    chr.l2[as.character(seqnames(data))]

  gr.new <- GRanges("genome", IRanges(sts.new, ed.new),
                    strand = strand(data))
  values(gr.new) <- values(data)
  values(data) <- NULL
  values(gr.new)$.ori <- data
  max.chr <- rev(names(chr.l2))[1]
  x.max <- chr.l.back[max.chr] + skps[max.chr] + chr.l2[max.chr]
  names(x.max) <- "genome"
  seqlengths(gr.new) <- x.max
  metadata(gr.new)$space.skip <- space.skip
  metadata(gr.new)$coord <- "genome"
  gr.new
}

is_coord_genome <- function(data){
  if(!is.null(metadata(data)$coord)){  
   res <- metadata(data)$coord == "genome"
  }else{
    res <- FALSE
  }
  res
}


## then need a transformation to circlular view
## data is a GRanges dataect
transformToCircle <- function(data, x = NULL, y = NULL,
                      radius = 10, trackWidth =10, direction = c("clockwise", "anticlockwise"),
                      mul = 0.05){
  if(!".midpoint" %in% colnames(values(data)))
    values(data)$.midpoint <- start(data) + width(data)/2
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
    
  temp.x <- values(data)[,x]
  
  if(is.character(y))
    temp.y <- values(data)[,y]
  if(is.numeric(y)){
    temp.y <- y
  }
  temp.y.r <- expand_range(range(temp.y), mul = mul)
  temp.y <- (temp.y - min(temp.y.r))/diff(temp.y.r) * trackWidth + radius
  
  ## always from 1 to max
  direction <- match.arg(direction)
  if(direction == "clockwise")
    dir <- 1
  else
    dir <- -1
  ## need to use max of "ideo"
  if(!is.null(metadata(data)$transformed) && metadata(data)$transformed){
    x.max <- seqlengths(data)
  }else{
    x.max <- max(end(data))
  }
  angle.unit <- pi* 2/ x.max
  x.new <- temp.y * sin(angle.unit * temp.x * dir)
  y.new <- temp.y * cos(angle.unit * temp.x * dir)
  values(data)$.circle.x <- x.new
  values(data)$.circle.y <- y.new
  values(data)$.circle.angle <- angle.unit * temp.x * dir
  data
}

## interpolate segment or rectangle
## rectangle only use stepping as y

transformToRectInCircle <- function(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                      direction = c("clockwise", "anticlockwise"),
                      n = 5, mul = 0.05){

  direction <- match.arg(direction)
  data.back <- data
  inter.fun <- function(x, y) approx(x, y, n = n)
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
      values(data)$stepping <- disjointBins(data)
  }
  df <- as.data.frame(data)
  lst <- lapply(1:nrow(df), function(i){
    if(df[i,"width"] > 1){
      res.x <- as.integer(inter.fun(c(df[i,"start"], df[i,"end"]), c(0, 0))$x)
    }else{
      res.x <- c(df[i, "start"], df[i, "start"])
    }
    res <- data.frame(.biovizBase.new.x = c(res.x, rev(res.x)))
    N <- nrow(res)
    res$.biovizBase.group <- i
    res$.biovizBase.level <- c(rep(df[i, "stepping"] - 0.4, N/2),
                               rep(df[i, "stepping"] + 0.4, N/2))
    df.extra <- do.call("rbind", lapply(1:N, function(k) df[i,]))
    res <- cbind(res, df.extra)
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
  res <- transformToGenome(res.gr, space.skip)
  res <- transformToCircle(res, y = .y, radius = radius, trackWidth = trackWidth,
                   direction = direction, mul = mul)
  res
}

transformToBarInCircle <- function(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                     direction = c("clockwise", "anticlockwise")){

  direction <- match.arg(direction)
  if(length(y)){
    temp.y <- values(data)[,as.character(y)]
    if(!all(check.integer(temp.y)))
      stop("geom(bar) require specified 'y', must be integer which indicates the levels")
    values(data)$stepping <- temp.y
  }else{
      values(data)$stepping <- disjointBins(data)
  }
  df <- as.data.frame(data)
  lst <- lapply(1:nrow(df), function(i){
    res.x <- rep(df[i, "start"], 2)
    res <- data.frame(.biovizBase.new.x = res.x)
    N <- nrow(res)
    res$.biovizBase.group <- i
    res$.biovizBase.level <- c(df[i, "stepping"] - 0.4,
                               df[i, "stepping"] + 0.4)
    df.extra <- rbind(df[i, ], df[i, ])
    res <- cbind(res, df.extra)
    res
  })
  res <- do.call("rbind", lst)
  .y <- ".biovizBase.level"
  res.gr <- GRanges(res$seqnames, IRanges(start = res$.biovizBase.new.x,
                                          width = 1),
                    strand = res$strand)
  values(res.gr) <- subset(res, select = -c(start, end, width, strand,seqnames))
  seqlengths(res.gr) <- seqlengths(data)
  res <- transformToGenome(res.gr, space.skip)
  res <- transformToCircle(res, y = .y, radius = radius, trackWidth = trackWidth,
                   direction = direction)
  res
}

## ok, segment allow user to use a flexible y
transformToSegInCircle <- function(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                      direction = c("clockwise", "anticlockwise"),
                      n = 5){
  direction <- match.arg(direction)
  ## do the linear interpolatoin first
  inter.fun <- function(x, y) approx(x, y, n = n)  
  if(!length(y)){
    values(data)$stepping <- disjointBins(ranges(data))
  }
  df <- as.data.frame(data)
  lst <- lapply(1:nrow(df), function(i){
    res.x <- as.integer(inter.fun(c(df[i,"start"], df[i,"end"]), c(0, 0))$x)
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
  res <- transformToGenome(res.gr, space.skip)
  res <- transformToCircle(res, y = .y, radius = radius, trackWidth =trackWidth,
                   direction = direction)
  res
}

## for a special GRanges
transformToLinkInCircle <- function(data, linked.to, space.skip = 0.1, trackWidth = 10,
                                    radius = 10, link.fun = function(x, y, n = 100)
                                    bezier(x, y, evaluation = n),
                                    direction = c("clockwise", "anticlockwise")){
  direction <- match.arg(direction)
  if(missing(linked.to))
    stop("linked.to must be provided and be a GRanges")
  ## trace id
  N <- length(data)
  values(data)$.biovizBase.idx <- 1:N
  values(values(data)[, linked.to])$.biovizBase.idx <- 1:N
  obj <- transformToGenome(data, space.skip)
  obj <- transformToCircle(obj, y = 0, radius = radius, trackWidth =trackWidth,
                   direction = direction)
  obj <- obj[order(values(obj)$.biovizBase.idx)]  
  df <- as.data.frame(obj)
  linktodata <- values(data)[,linked.to]
  ## values(linktodata)$.biovizBase.idx <- 1:length(linktodata)
  ## missing y
  linktodata <- transformToGenome(linktodata, space.skip)
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
##   idx <- elementLengths(grl) > 0
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
                           width = NULL, strand = NULL){


  ## this losing seqinfo....
  colnms <- colnames(data)
  if(is.null(seqnames)){
    if("seqnames" %in% colnms)
      seqnames <- "seqnames"
    else
      stop("Please sepicify which column represent the seqnames")
  }

  if(is.null(start)){
    if("start" %in% colnms)
      start <- "start"
    else
      stop("Please sepicify which column represent the start")
  }

  if(is.null(end)){
    if("end" %in% colnms)
      end <- "end"
    else if(is.null(width)){
      if("width" %in% colnms)
        width <- "width"
      else
        stop("Must provide end or width column names")
    }}
  
  if(is.null(width)){
    if("width" %in% colnms)
        width <- "width"
  }

  if(!is.null(end))
    gr <- GRanges(data[,seqnames], IRanges(start = data[,start], end = data[,end]))
  else
    gr <- GRanges(data[,seqnames], IRanges(start = data[,start], width = data[,width]))

  if(is.null(strand)){
    if("strand" %in% colnms){
      strand <- "strand"
      strand(gr) <- data[,strand]
    }
  }
  values(gr) <- data[,!colnames(data) %in% c(start, end, width, seqnames, strand)]
  gr
}

check.integer <- function(x){
  sapply(x, function(x)
    !length(grep("[^[:digit:]]", as.character(x))))
}
