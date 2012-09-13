## going to be moved to biovizBase
setGeneric("getYLab", function(obj,...) standardGeneric("getYLab"))
setMethod("getYLab", "TranscriptDb", function(obj){
  md <- metadata(obj)
  ds <- md[md[,1] == "Data source",2]
  lb <- character()
  if(ds == "UCSC"){
    lb <- md[md[,1] == "UCSC Table",2]
  }
  if(ds == "BioMart"){
    lb <- md[md[,1] == "BioMart Database",2]
  }
  lb
})

## going to be moved to biovizBase
setGeneric("getXLab", function(obj,...) standardGeneric("getXLab"))
setMethod("getXLab", "GRanges", function(obj){
  chrs <- unique(as.character(seqnames(obj)))
  gms <- genome(obj)
  gm <- gms[chrs]
  if(length(chrs) > 1){
    xlab <- "Position"
  }else{
    xlab <- paste(ifelse(is.na(gm), chrs, paste0(gm, "::", chrs)))
    xlab <- paste("Position on", xlab)
  }
  xlab
})

setMethod("getXLab", "GRangesList", function(obj){
  chrs <- unique(as.character(seqnames(unlist(obj))))
  gms <- genome(obj)
  gm <- gms[chrs]
  if(length(chrs) > 1){
    xlab <- "Position"
  }else{
    xlab <- paste(ifelse(is.na(gm), chrs, paste0(gm, "::", chrs)))
    xlab <- paste("Position on", xlab)
  }
  xlab  
})

setMethod("getXLab", "GappedAlignments", function(obj){
  chrs <- unique(as.character(seqnames(obj)))
  gms <- genome(obj)
  gm <- gms[chrs]
  if(length(chrs) > 1){
    xlab <- "Position"
  }else{
    xlab <- paste(ifelse(is.na(gm), chrs, paste0(gm, "::", chrs)))
    xlab <- paste("Position on", xlab)
  }
  xlab
})

