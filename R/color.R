.DNABasesNColor <- function(){
  c(A = "#ABD9E9",
    T = "#2C7BB6",
    G = "#D7191C",
    C = "#FDAE61",
    N = "#FFFFBF")
}

.DNABasesColor <- function(){
  c(A = "#ABD9E9",
    T = "#2C7BB6",
    G = "#D7191C",
    C = "#FDAE61")
}

.DNAAlphabetColor <- function(pal = c("BluetoDarkOrange.18",
                                "DarkRedtoBlue.18")){

  pal <- match.arg(pal)
  cols <- dscale(DNA_ALPHABET, dichromat_pal(pal))
  names(cols) <- DNA_ALPHABET
  cols
}

.RNABasesNColor <- function(){
  c(A = "#ABD9E9",
    U = "#2C7BB6",
    G = "#D7191C",
    C = "#FDAE61",
    N = "#FFFFBF")
}

.RNABasesColor <- function(){
  c(A = "#ABD9E9",
    U = "#2C7BB6",
    G = "#D7191C",
    C = "#FDAE61")
}

## blind.pal.info
.RNAAlphabetColor <- function(){
  ## pal <- match.arg(pal)
  ## cols <- dscale(RNA_ALPHABET, dichromat_pal(pal))
  N <- length(RNA_ALPHABET)
  cols <- colorBlindSafePal("Categorical.12")(N, repeatable = TRUE)
  names(cols) <- RNA_ALPHABET
  cols
}


.IUPACCodeMapColor <- function(){
  N <- length(IUPAC_CODE_MAP)
  ## cols <- dscale(IUPAC_CODE_MAP, dichromat_pal(pal))
  cols <- colorBlindSafePal("Categorical.12")(N, repeatable = TRUE)
  names(cols) <- names(IUPAC_CODE_MAP)
  cols
}

.strandColor <- function(){
  c("+" = "#1B9E77",
    "-" = "#7570B3",
    "*" = "#D95F02")
}

## preset gpos1:100 and other widely used color
.cytobandColor <- function(){
  .getCytobandColor(c("gneg", "gvar", "stalk", "acen", paste0("gpos", 1:100)))
}

## got these codes and hints from Florian Hahne, from Gviz
.getCytobandColor <- function(type){
    type <- as.character(type)
    ocols <-  c(gneg = "grey100",  gpos100 = "grey0", gvar = "grey100", 
                stalk = "brown3", acen = "brown4")
    cols <- c(ocols[c("gneg", "stalk", "acen")],
              gpos=unname(ocols["gpos100"]), gvar=unname(ocols["gpos100"]))
    gpcols <- unique(grep("gpos", type, value=TRUE))
    crmp <- colorRampPalette(c(cols["gneg"], cols["gpos"]))(100)
    posCols <- setNames(crmp[as.integer(gsub("gpos", "", gpcols))], gpcols)
    return(c(cols, posCols))
}


.AminoAcidCodeColor <- function(){
  N <- length(AMINO_ACID_CODE)
  cols <- colorBlindSafePal("Categorical.12")(N, repeatable = TRUE)
  names(cols) <- names(AMINO_ACID_CODE)
  cols
}


.AAAlphabetColor <- function(){
  N <- length(AA_ALPHABET)
  cols <- colorBlindSafePal("Categorical.12")(N, repeatable = TRUE)
  names(cols) <- AA_ALPHABET
  cols
}

getBioColor <- function(type = c("DNA_BASES_N",
                          "DNA_BASES",
                          "DNA_ALPHABET",
                          "RNA_BASES_N",
                          "RNA_BASES",
                          "RNA_ALPHABET",
                          "IUPAC_CODE_MAP",
                          "AMINO_ACID_CODE",
                          "AA_ALPHABET",
                          "STRAND",
                          "CYTOBAND"),
                        source = c("option",
                          "default")){

  type <- match.arg(type)
  source <- match.arg(source)

  if(source == "default"){
    cols <- switch(type,
                   DNA_BASES_N = .DNABasesNColor(),
                   DNA_BASES = .DNABasesColor(),
                   DNA_ALPHABET = .DNAAlphabetColor(),
                   RNA_BASES_N = .RNABasesNColor(),
                   RNA_BASES = .RNABasesColor(),
                   RNA_ALPHABET = .RNAAlphabetColor(),
                   IUPAC_CODE_MAP = .IUPACCodeMapColor(),
                   AMINO_ACID_CODE = .AminoAcidCodeColor(),
                   AA_ALPHABET = .AAAlphabetColor(),
                   STRAND = .strandColor(),
                   CYTOBAND = .cytobandColor())
  }

  if(source == "option"){
    colname <- switch(type,
                      DNA_BASES_N = "DNABasesNColor",
                      DNA_BASES = "DNABasesColor",
                      DNA_ALPHABET = "DNAAlphabetColor",
                      RNA_BASES_N = "RNABasesNColor",
                      RNA_BASES = "RNABasesColor",
                      RNA_ALPHABET = "RNAAlphabetColor",
                      IUPAC_CODE_MAP = "IUPACCodeMapColor",
                      AMINO_ACID_CODE = "AminoAcidCodeColor",                     
                      AA_ALPHABET = "AAAlphabetColor",
                      STRAND = "strandColor",
                      CYTOBAND = "cytobandColor")
    cols <- getOption("biovizBase")[[colname]]
    ## res <- cols[base]
  }
  cols
}

plotColorLegend <- function(colors, labels, title){
  if(missing(labels))
    labels <- names(colors)
  n <- length(colors)
  plot(c(0, 1), c(1, -n), type = "n", xlab="", ylab="", axes = FALSE)
  if(!missing(title))
    text(0.5, 1, title)
  rect(0.5, -(1:n), 1, -(1:n)+0.5+0.25, col = colors)
  text(0, -(1:n) + 0.5, labels)
}




## this should include as many as I cann
colorBlindSafePal <- function(palette){
  namelist <- rownames(blind.pal.info)
  if(missing(palette)){
    return(blind.pal.info)
  }
  if(is.numeric(palette) && (palette < nrow(blind.pal.info)))
    name <- rownames(blind.pal.info)[palette]
  else{
    if(is.character(palette) && (palette %in% namelist))
      name <- palette
    else
      stop("Invalid palette name, please check colorBlindSafePal().")
  }
  pkg <- blind.pal.info[name,"pkg"]
  maxcolors <- as.numeric(blind.pal.info[name,"maxcolors"])
  getPalColor <- function(n, pal, pkg){
    cols <- switch(pkg,
                   RColorBrewer = {
                     brewer.pal(n, pal)
                   },
                   dichromat = {
                     colorschemes[[pal]][seq(n)]
                   })
  }
  function(n, repeatable = FALSE){
    if (n < 3) {
      warning("minimal value for n is 3,
    returning requested palette with 3 different levels\n")
      cols <- getPalColor(3, name, pkg)
    }

    if (n > maxcolors) {
      if(!repeatable){
        warning(paste("n too large, allowed maximum for palette", 
                      name, "is", maxcolors, 
                      "\nReturning the palette you asked for with that many colors\n"))
        cols <- getPalColor(maxcolors, name, pkg)
      }else{
        ## repeatable colors
        rpt <- n %/% maxcolors
        res <- n %% maxcolors
        full.cols <- getPalColor(maxcolors, name, pkg)
        cols <- c(rep(full.cols, rpt), full.cols[-maxcolors][1:res])      
      }
    }else{
      cols <- getPalColor(n, name, pkg)
    }
    cols
  }
}

genBrewerBlindPalInfo <- function(){
  res <- RColorBrewer::brewer.pal.info
  idx.seq <- which(rownames(res)
                   %in% c("Blues","BuPu","GnBu","Greens",
                          "Greys","Oranges","OrRd","PuBu","PuBuGn",
                          "PuRd","Purples","Reds","YlGn","YlGnBu",
                          "YlOrBr","YlOrRd"))
  idx.div <- which(rownames(res)
                   %in% c("PiYG","PRGn","PuOr","RdBu","RdYlBu"))
  idx.qual <- which(rownames(res) %in% c("Dark2","Paired","Set2"))
  res <- res[c(idx.seq, idx.div, idx.qual),]
  res[rownames(res) %in% c("Dark2", "Set2"), 1] <- 3
  res[rownames(res) %in% c("Paired"), 1] <- 4
  res[, c("maxcolors", "category"), drop=FALSE]
}


## return current pallete which is color-blind safe, if missing n,
## return list of supported pal, if not, use index to return a set
## color. should have qual/div/seq
.colorschemes.num <-function(){
  lst <- gregexpr("[0-9]",names(colorschemes))
  unlist(lapply(1:length(lst), function(i){
  pos <- lst[[i]]
  if(length(pos) == 1)
    pos <- c(pos, pos)
  as.numeric(substr(names(colorschemes)[i], pos[1], pos[2]))
}))}

## make a similar color pal.info like brewer.pal.info
genDichromatPalInfo <- function(){
  dichromat.pal.info <- data.frame(maxcolors = .colorschemes.num(),
                                   category = "div", stringsAsFactors = FALSE)
  rownames(dichromat.pal.info) <- names(colorschemes)
  idx.seq <-  rownames(dichromat.pal.info) %in%
  c("LightBluetoDarkBlue.10", "LightBluetoDarkBlue.7", "SteppedSequential.5")   
  idx.qual <- rownames(dichromat.pal.info) %in% "Categorical.12"
  dichromat.pal.info[idx.seq, 2] <- "seq"
  dichromat.pal.info[idx.qual, 2] <- "qual"
  dichromat.pal.info <- dichromat.pal.info[order(dichromat.pal.info[,2]),]
  dichromat.pal.info                   
}



genBlindPalInfo <- function(){
  blind.pal.info <- rbind(brewer.pal.blind.info, dichromat.pal.blind.info)
  id <- seq(nrow(blind.pal.info))
  pkg <- c(rep("RColorBrewer", nrow(brewer.pal.blind.info)),
           rep("dichromat", nrow(dichromat.pal.blind.info)))
  blind.pal.info$pkg <- pkg
  idx <- order(blind.pal.info$category, blind.pal.info$maxcolors)
  blind.pal.info <- blind.pal.info[idx,]
  blind.pal.info$pal.id <- id
  blind.pal.info
}

brewer.pal.blind.info <- genBrewerBlindPalInfo()
dichromat.pal.blind.info <- genDichromatPalInfo()
blind.pal.info <- genBlindPalInfo()

## show color
show_col2 <- function(colours) {
  n <- length(colours)
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  
  colors <- c(colours, rep(NA, nrow * ncol - length(colours)))
  colors <- matrix(colors, ncol = ncol, byrow = TRUE)

  colours <- c(colours, rep(NA, nrow * ncol - length(colours)))
  nms <- matrix(names(colours), ncol = ncol, byrow = TRUE)
  
  old <- par(pty = "s", mar = c(0, 0, 0, 0))
  on.exit(par(old))
  size <- max(dim(colors))
  plot(c(0, size), c(0, -size), type = "n", xlab="", ylab="", axes = FALSE)
  rect(col(colors) - 1, -row(colors) + 1, col(colors), -row(colors),
    col = colors)
  text(col(colors) - 0.5, -row(colors) + 0.5, nms)
}

showColor <- function(colors, label = c("color", "name")){
  label <- match.arg(label)
  colfun <- switch(label,
                color = scales::show_col,
                name = show_col2)
  colfun(colors)
}

