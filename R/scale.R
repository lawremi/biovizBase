## utils for creating proper scale breaks and labels for certain coordinate
setGeneric("getXScale", function(obj, ...) standardGeneric("getXScale"))
setMethod("getXScale", "GRanges", function(obj, type = c("default",
                                                  "all", "left", "right")){
  type <- match.arg(type)
  if(is_coord_truncate_gaps(obj)){
    message("using coord:truncate_gaps to parse x scale")
    res <- switch(type,
                  left = {
                    list(breaks = start(obj),
                         labels = start(values(obj)$.ori))
                  },
                  right = {
                    list(breaks = end(obj),
                         labels = end(values(obj)$.ori))
                  },
                  all = {
                    list(breaks = c(start(obj), end(obj)),
                         labels = c(start(values(obj)$.ori),
                           end(values(obj)$.ori)))
                  },
                  default = {
                    ## breaks = pretty(c(min(start(obj)), max(end(obj))), n = 5)
                    breaks = c(min(start(obj)), max(end(obj)))
                    labels = c(min(start(values(obj)$.ori)),
                      max(end(values(obj)$.ori)))
                    ## labels = c(min(start(values(obj)$.ori)),
                    ##   rep("~", length(breaks) -2),
                    ##   max(end(values(obj)$.ori)))
                    list(breaks = breaks, labels = labels)
                  })
  }else if(is_coord_genome(obj)){
    message("using coord:genome to parse x scale")
    res <- switch(type,
                  default = {
                      grl <- split(obj, seqnames(values(obj)$.ori))
                      res <- lapply(grl, function(gr){
                        mid <- (min(start(gr)) + max(end(gr)))/2
                        seqs <- unique(as.character(seqnames(values(gr)$.ori)))
                        df <- data.frame(breaks = mid, labels = seqs)
                      })
                      res <- do.call(rbind, res)
                      list(breaks = res$breaks,
                           labels = res$labels)
                  })
  }else{
    stop("coord other than truncate_gaps are not suported for getXScale yet")
  }
  res
})
