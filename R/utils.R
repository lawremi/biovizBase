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
    if(!length(obj))
        return(GRanges())
    if(!length(facets))
        facets <- as.formula(~seqnames)
    allvars <- all.vars(as.formula(facets))
    allvars.extra <- allvars[!allvars %in% c(".", "seqnames")]
    if(!is.null(group.name)){
        if(!group.name %in% colnames(values(obj)))
            stop(group.name, " is not in obj")
        grl <- splitByFacets(obj, facets)
        grl <- endoapply(grl, function(dt){
            res <- split(dt, values(dt)[,group.name])
            gps <- gaps(ranges(res))
            idx <- elementLengths(gps) > 0
            res.sub <- res[idx,]
            stds <- unlist(lapply(res.sub, function(x) as.character(strand(x))[1]))
            gps.sub <- gps[idx,]
            dfs <- values(unlist(res.sub))[cumsum(elementLengths(res.sub)),c(allvars.extra, "stepping"),
                                           drop = FALSE]

            ir <- unlist(gps.sub)
            if(length(ir)){
                gr <- GRanges(unique(seqnames(dt)), ir)
                values(gr) <- dfs[togroup(gps.sub),,drop = FALSE]
                strand(gr) <- stds[togroup(gps.sub)]
            }else{
                gr <- GRanges()
            }
            gr
        })
        res <- unlist(grl)
    }else{

    }
    if(length(res)){
        values(res)$type <- "gaps"
        res <- resize(res, width = width(res) + 2, fix = "center")
    }else{
        GRanges()
    }
    res
})



getIdeoGR <- function(gr){
    if(!is(gr, "GenomicRanges"))
        stop("require GenomicRanges")
    if(all(is.na(seqlengths(gr)))){
        warning("geom(ideogram) need valid seqlengths information for accurate mapping,
                 now use reduced information as ideogram... ")
        res <- range(reduce(gr, ignore = TRUE))
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


## test heterogeneity
is_homo <- function(grl){
    all(unlist(lapply(grl, function(gr){
        length(unique(as.character(seqnames(gr)))) == 1
    })))
}


genGenesymbolTable <- function(db, keys, keytype = "SYMBOL",
                               columns = c("CHR", "CHRLOC", "CHRLOCEND", "SYMBOL"), unique = TRUE){

    res <- select(db, keys = keys, keytype = keytype, columns = columns)
    chr <- res[, "CHR"]
    s <- substr(res[, "CHRLOC"], 1, 1)
    s[s != "-"] <- "+"

    st <- as.numeric(substr(res[, "CHRLOC"], 2, nchar(res[, "CHRLOC"])))
    ed <- as.numeric(substr(res[, "CHRLOCEND"], 2, nchar(res[, "CHRLOCEND"])))
    sym <- res[, "SYMBOL"]

    res.gr <- GRanges(chr, IRanges(start = st, end = ed), strand = s, symbol = sym)
    if(unique){
        res.gr.l <- split(res.gr, sym)
        syms <- sapply(res.gr.l, function(x) unique(x$symbol))
        res.gr.l <- range(res.gr.l)
        res.gr <- unlist(res.gr.l)
        names(res.gr) <- syms
    }
    res.gr
}


## setGeneric("fetch", function(obj, ...) standardGeneric("fetch"))
fetch <- function(obj, which, ..., gene.id,
                  truncate.gaps = FALSE,
                  truncate.fun = NULL,
                  ratio = 0.0025,
                  resize.extra = 10,
                  include.level = TRUE,
                  use.names = TRUE,
                  columns = c("tx_id", "tx_name","gene_id"),
                  type = c("all", "single", "exons.reduce", "exons.all",
                      "exons.geneid",  "gapped.pair", "raw", "all")){
    .txdb.type <- c("all", "single", "exons.reduce",
                    "exons.all", "exons.geneid")
    .bamfile.type <- c("gapped.pair", "raw", "all")

    ## type <- match.arg(type)
    if(length(type) >1){
        if(is(obj,"TxDb"))
            type <- "all"
        if(is(obj,"BamFile"))
            type <- "gapped.pair"
    }
    if(is(obj, "TxDb")){
        if(!(type %in% .txdb.type))
            stop("type for TxDb must be ", .txdb.type)
        if(is.list(which)){
            message("Parsing exons based on which(list) arguments")
            temp <- exons(obj, vals = which, columns = columns)
            which <- range(temp)
        }
        ## 1st set all the sequences to be inactive:
        isActiveSeq(obj)[seqlevels(obj)] <- FALSE
        ## Then set only "chr9" to be active:
        seqnms <- as.character(unique(seqnames(which)))
        isActiveSeq(obj)[seqnms] <- TRUE
        ## You can see which are active like this:
        if(type == "exons.geneid"){
            message("return exons within one gene...")
            ## this require gene id, only exons
            exons.gene <- exonsBy(obj, by = "gene")
            if(missing(gene.id))
                stop("Please specify the gene id")
            res <- exons.gene[[as.character(gene.id)]]
        }
        if(type == "exons.reduce"){
            message("Parsing exons...")
            exons <- exonsByOverlaps(obj, which, ...)
            exons <- reduce(exons)
            values(exons)$type <- "exon"
            res <- exons
        }
        if(type == "exons.all"){
            message("Parsing exons...")
            exons <- exonsByOverlaps(obj, which, ...)
            values(exons)$type <- "exon"
            res <- exons
        }
        if(type %in% c("all","single")){
            message("Parsing exons...")
            exons <- exonsBy(obj, "tx")
            message("Parsing cds...")
            cdss <- cdsBy(obj, "tx")
            message("Parsing transcripts...")
            tx <- transcripts(obj, columns = columns)
            tx <- subsetByOverlaps(tx, which)
            message("Aggregating...")
            txids <- values(tx)$tx_id
            lst <- lapply(txids, function(id){
                id <- as.character(id)
                tx_nm <- values(tx)[values(tx)$tx_id == id, "tx_name"]
                gene_id <- values(tx)[values(tx)$tx_id == id, "gene_id"]
                gene_id <- paste(unlist(gene_id), sep = ",")
                if(!length(gene_id))
                    gene_id <- NA
                exons.cur <- exons[[id]]
                values(exons.cur) <- data.frame(tx_id = id,
                                                tx_name = tx_nm,
                                                gene_id = gene_id,
                                                type = "exon")

                cds.cur <- cdss[[id]]
                seqs <- unique(as.character(seqnames(exons.cur)))
                exon_union <- reduce(exons.cur)
                ir.g <- IRanges(gaps(ranges(exons.cur)))
                if(length(ir.g)){
                    introns.cur <- GRanges(seqs,ir.g)
                    values(introns.cur) <- data.frame(tx_id = id,
                                                      tx_name = tx_nm,
                                                      gene_id = gene_id,
                                                      type = "intron")

                }else{
                    introns.cur <- GRanges()
                }
                if(!is.null(cds.cur)){
                    values(cds.cur) <- data.frame(tx_id = id,
                                                  tx_name = tx_nm,
                                                  gene_id = gene_id,
                                                  type = "cds")

                    cds_union <- reduce(cds.cur)
                    utrs <- setdiff(exon_union, cds_union)
                    if(length(utrs))
                        values(utrs) <- data.frame(tx_id = id,
                                                   tx_name = tx_nm,
                                                   gene_id = gene_id,
                                                   type = "utr")
                    else
                        utrs <- GRanges()
                    ir.g2 <- gaps(reduce(c(ranges(cds_union), ranges(utrs))))
                    if(length(ir.g2)){
                        gaps.cur <- GRanges(seqs,ir.g2)
                        values(gaps.cur) <- data.frame(tx_id = id,
                                                       tx_name = tx_nm,
                                                       gene_id = gene_id,
                                                       type = "gap")
                    }else{
                        gaps.cur <- GRanges()
                    }
                    gr <- c(exons.cur, introns.cur, cds.cur, gaps.cur, utrs)
                }else{
                    utrs <- exons.cur
                    values(utrs)$type <- factor("utr")
                    gaps.cur <- introns.cur
                    if(length(gaps.cur)){
                        values(gaps.cur)$type <- factor("gap")
                    }
                    gr <- c(exons.cur, introns.cur, utrs, gaps.cur)
                }
                gr
            })
            res <- do.call("c", lst)
            isActiveSeq(obj)[seqlevels(obj)] <- TRUE
            if(length(res))
                res <- keepSeqlevels(res, unique(as.character(seqnames(res))))
            else
                res <- GRanges()
        }
        if(type == "single"){
            if(length(res)){
                cds.s <- reduce(res[values(res)$type == "cds"])
                values(cds.s)$type <- factor("cds")
                exon.s <- reduce(res[values(res)$type == "exon"])
                values(exon.s)$type <- factor("exon")
                utr.s <- setdiff(exon.s, cds.s)
                values(utr.s)$type <- factor("utr")
                gap.s <- gaps(cds.s, start = min(start(cds.s)),
                              end = max(end(cds.s)))
                values(gap.s)$type <- factor("gap")
                res <- c(cds.s, utr.s, gap.s)
            }
        }
        message("Done")
    }
    if(is(obj, "GappedAlignments")){
        ## require(Rsamtools)
        ## require(Rsamtools)
        if(!missing(which))
            obj <- subsetByOverlaps(obj, which)
        if(!length(obj))
            stop("No entry in the GappedAlignments data")
        if(is.null(names(obj)))
            stop("Missing qname, please make use.names = TRUE, when reading GappedAlignments")
        grg <- grglist(obj)
        grg.u <- stack(grg, ".grl.name")
        message("extracting information...")
        values(grg.u)$junction <- rep(ifelse(elementLengths(grg)>1, TRUE, FALSE),
                                      times = elementLengths(grg))
        names(grg.u) <- NULL
        res <- grg.u
    }
    if(is(obj, "BamFile")){
        if(!(type %in% .bamfile.type))
            stop("type for TxDb must be ", .bamfile.type)
        ## require(Rsamtools)
        if(type == "gapped.pair"){
            message("Read GappedAlignments from BamFile...")
            ga <- readGAlignments(obj, param = ScanBamParam(which = which),
                                       use.names = use.names, ...)
            res <- fetch(ga)
        }

        if(type == "raw"){
            message("Read Raw from BamFile by calling scanBam...")
            res <- scanBamGRanges(obj, which, ...)
        }

        if(type == "all"){
            message("Read Raw from BamFile by calling scanBam...")
            res.mb <- scanBamGRanges(obj, which,
                                     flag = scanBamFlag(isFirstMateRead = TRUE))
            message("Read GappedAlignments from BamFile...")
            ga <- readGAlignments(obj, param = ScanBamParam(which = which),
                                       use.names = use.names, ...)
            res.gp <- fetch(ga)
            message("Combine...")
            nms <- values(res.mb)$qname
            ## fow now just, using isize
            isize <- values(res.mb)$issize
            names(isize) <- nms
            values(res.gp)$isize <- isize[values(res.gp)$qname]
            res <- res.gp
        }
    }
    if(truncate.gaps){
        if(is.null(truncate.fun)){
            if("gap" %in% unique(values(res)$type))
                idx <- values(res)$type %in% c("utr", "cds")
            res.s <- reduce(res[idx], ignore.strand = TRUE)
            truncate.fun <- shrinkageFun(gaps(res.s, min(start(res.s)), max(end(res.s))),
                                         maxGap(gaps(res.s, min(start(res.s)), max(end(res.s))),
                                                ratio = ratio))
        }
        res <- truncate.fun(res)
    }
    res
}

rectifySeqlevelsStyle <- function(x, y) {
    seqlevelsStyle(x) <- seqlevelsStyle(y)
    x
}

.transformSeqinfo <- function(obj){
    ss <- seqlengths(obj)
    res <- GRanges(seqnames(obj), IRanges(start = 1, end = ss))
    res <- keepSeqlevels(res, names(ss))
    seqlengths(res) <- ss
    res
}

