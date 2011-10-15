##----------------------------------------------------------------##
##                         ideogram
##----------------------------------------------------------------##
getIdeogram <- function(genome,subchr,cytobands=TRUE){
  ## to remove the "heavy dependency" we put require here.
  require(rtracklayer)
  if(!(exists("session")&&extends(class(session),"BrowserSession")))
    session <- browserSession()
  if(missing(genome)){
    choices <- ucscGenomes()[,1]
    res <- menu(choices,title="Please specify genome")
    genome <- as.character(choices[res])
  }
  if(cytobands){
    message("Loading...")
    query <- ucscTableQuery(session,"cytoBand",genome)
    tableName(query) <- "cytoBand"
    df <- getTable(query)
    gr <- GRanges(seqnames=df$chrom,
                  IRanges(start=df$chromStart,end=df$chromEnd))
    values(gr) <- df[,c("name","gieStain")]
    message("Done")
  }else{
    message("Loading...")
    gr <- GRangesForUCSCGenome(genome)
    message("Done")
  }
  if(!missing(subchr))
    gr <- gr[seqnames(gr) == subchr]
  ## gr <- sortChr(gr)
  ## sortChr is removed
  ##
  gr <- sort(gr)
  gr
}

isIdeogram <- function(obj){
  is(obj, "GRanges") && (c("gieStain", "name") %in% colnames(values(obj)))
}

isSimpleIdeogram <- function(obj){
  is(obj, "GRanges") && all(elementLengths(split(obj, as.character(seqnames(obj)))) == 1)
}

