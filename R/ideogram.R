##----------------------------------------------------------------##
##                         ideogram
##----------------------------------------------------------------##
getIdeogram <- function(genome,subchr = NULL,cytobands=TRUE){
  .gnm <- genome
  lst <- lapply(.gnm, function(genome){
    print(genome)
  ## to remove the "heavy dependency" we put require here.
  require(rtracklayer)
  if(!(exists("session")&&extends(class(session),"BrowserSession")))
    session <- rtracklayer::browserSession()
  if(missing(genome)){
    choices <- rtracklayer::ucscGenomes()[,1]
    res <- menu(choices,title="Please specify genome")
    genome <- as.character(choices[res])
  }
  if(cytobands){
    message("Loading ideogram...")
    tryres <- try(query <-
                  rtracklayer::ucscTableQuery(session,"cytoBand",genome))
    if(!inherits(tryres, "try-error")){
      rtracklayer::tableName(query) <- "cytoBand"
      df <- rtracklayer::getTable(query)
      gr <- GRanges(seqnames=df$chrom,
                    IRanges(start=df$chromStart,end=df$chromEnd))
      values(gr) <- df[,c("name","gieStain")]
      message("Loading ranges...")
   
      gr.r <- rtracklayer::GRangesForUCSCGenome(genome)
      suppressWarnings(seqlengths(gr) <- seqlengths(gr.r)[names(seqlengths(gr))])
      gr <- trim(gr)
      
    }else{
      message("cytoBand informatin is not available, only get ranges.")
      message("Loading ranges...")
      gr <- rtracklayer::GRangesForUCSCGenome(genome)
      message("Done")      
    }
  }else{
    message("Loading...")
    gr <- rtracklayer::GRangesForUCSCGenome(genome)
    message("Done")
  }
  if(length(subchr))
    gr <- gr[seqnames(gr) == subchr]

  gr <- sort(gr)
  gr
  })
  names(lst) <- .gnm
  if(length(lst) == 1)
    res <- lst[[1]]
  else 
    res <- lst
  res
}

isIdeogram <- function(obj){
  is(obj, "GRanges") && (c("gieStain", "name") %in% colnames(values(obj)))
}

isSimpleIdeogram <- function(obj){
  is(obj, "GRanges") && all(elementNROWS(split(obj, as.character(seqnames(obj)))) == 1)
}

