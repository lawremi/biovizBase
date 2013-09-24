setGeneric("crunch", function(obj, ...) standardGeneric("crunch"))
setMethod("crunch", "TranscriptDb", function(obj, which,
                                             columns = c("tx_id", "tx_name","gene_id"),
                                            type = c("all", "reduce")){


    if(is.list(which)){
        message("Parsing exons based on which(list) arguments")
        temp <- exons(obj, vals = which, columns = columns)
        which <- range(temp)
    }

    seqnms <- as.character(unique(seqnames(which)))

    seqlevels(obj, force = TRUE) <- seqnms
    on.exit(seqlevels0(obj))


    message("Parsing exons...")
    ## system.time(exons <- exonsBy(obj, "tx")) #takes 3.5s
    exons <- exonsByOverlaps(obj, which, columns = c("exon_id", "tx_id")) # 1.6
    txids <- unlist(exons$tx_id)
    idx <- togroup(exons$tx_id)
    exons <- exons[idx]
    exons <- exons[,1]
    exons$tx_id <- txids
    exons <- split(exons, exons$tx_id)

    message("Parsing cds...")
    cdss <- cdsByOverlaps(obj, which, columns = c("exon_id", "tx_id"))
    txids <- unlist(cdss$tx_id)
    idx <- togroup(cdss$tx_id)
    cdss <- cdss[idx]
    cdss <- cdss[,1]
    cdss$tx_id <- txids
    cdss <- split(cdss, cdss$tx_id)

    ## system.time(cdss <- cdsBy(obj, "tx")) #2.552s

    message("Parsing transcripts...")
    ## tx <- transcripts(obj, columns = columns)
    ## subsetByOverlaps(tx, which)
    
    tx <- transcriptsByOverlaps(obj, which, columns = columns)

    message("Parsing utrs and aggregating...")      
    txids <- as.character(values(tx)$tx_id)

    ## new method operate on GRangesList levels, and it is much faster to aggregate
    gids <- values(tx)$gene_id
    gids <- unlist(lapply(gids, function(x) do.call(paste0, list(as.list(x), collapse = ","))))

    tx_nm <- values(tx)$tx_name

    ## match table
    mt <- data.frame(txids, gids, tx_nm)
    rownames(mt) <- txids
    colnames(mt) <- c("tx_id", "gene_id", "tx_name")
    

    ## exons
    gr.exons <- unlist(exons)
    .nms <- as.character(gr.exons$tx_id)
    .gid.nms <- mt[.nms, "gene_id"]
    .tx.nms <- mt[.nms, "tx_name"]
    values(gr.exons) <- NULL
    values(gr.exons) <- data.frame(tx_id = .nms, tx_name = .tx.nms,
                                   gene_id = .gid.nms, type = "exon")
    values(gr.exons)$type <-     as.character(values(gr.exons)$type)
    names(gr.exons) <- NULL
    
    ## cds
    gr.cdss <- unlist(cdss)
    .nms <- as.character(gr.cdss$tx_id)
    .gid.nms <- mt[.nms, "gene_id"]
    .tx.nms <- mt[.nms, "tx_name"]
    values(gr.cdss) <- NULL
    values(gr.cdss) <- data.frame(tx_id = .nms, tx_name = .tx.nms,
                                   gene_id = .gid.nms, type = "cds")
    values(gr.cdss)$type <-     as.character(values(gr.cdss)$type)
    names(gr.cdss) <- NULL    
    ## intron
    irl.introns <- gaps(ranges(exons))
    ir.introns <- unlist(irl.introns)
    .nms <- names(ir.introns)
    .gid.nms <- mt[.nms, "gene_id"]
    .tx.nms <- mt[.nms, "tx_name"]
    gr.introns <- GRanges(seqnms, ir.introns, tx_id = .nms,
                          tx_name = .tx.nms, gene_id = .gid.nms,
                          type = "gap")
    names(gr.introns) <- NULL
    values(gr.introns)$type <-     as.character(values(gr.introns)$type)
    
    ## utrs
    irl.utrs <- setdiff(ranges(exons), ranges(cdss))
    ir.utrs <- unlist(irl.utrs)
    .nms <- names(ir.utrs)
    .gid.nms <- mt[.nms, "gene_id"]
    .tx.nms <- mt[.nms, "tx_name"]
    gr.utrs <- GRanges(seqnms, ir.utrs, tx_id = .nms,
                       tx_name = .tx.nms, gene_id = .gid.nms,
                       type = "utr")
    names(gr.utrs) <- NULL
    values(gr.utrs)$type <-     as.character(values(gr.utrs)$type)

    ## combine
    res <- c(gr.exons, gr.cdss, gr.introns, gr.utrs)
    if(!length(res))
        res <- GRanges()
    res$type <- factor(res$type)
    
    if(type == "reduce"){
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
            ## change it to *
            strand(res) <- "*"
        }
    }
    message("Done")
    res
})

setMethod("crunch", "GAlignments", function(obj, which){
    if(!missing(which))
        obj <- subsetByOverlaps(obj, which)
    if(!length(obj))
        stop("No entry in the GAlignments data")
    if(is.null(names(obj)))
        stop("Missing qname, please make use.name = TRUE, when reading GAlignments")
    grg <- grglist(obj)
    grg.u <- stack(grg, ".grl.name")
    message("extracting information...")
    values(grg.u)$junction <- rep(ifelse(elementLengths(grg)>1, TRUE, FALSE),
                                  times = elementLengths(grg))
    names(grg.u) <- NULL
    res <- grg.u
})

setMethod("crunch", "BamFile", function(obj, which, ...){
    if(!(type %in% .bamfile.type))
        stop("type for TranscriptDb must be ", .bamfile.type)
    ## require(Rsamtools)
    if(type == "gapped.pair"){
        message("Read GAlignments from BamFile...")
        ga <- readGAlignmentsFromBam(obj,
                                     param = ScanBamParam(which = which),
                                     use.names = use.name, ...)
        res <- crunch(ga)
    }

    if(type == "raw"){
        message("Read Raw from BamFile by calling scanBam...")
        res <- scanBamGRanges(obj, which, ...)
    }

    if(type == "all"){
        message("Read Raw from BamFile by calling scanBam...")
        res.mb <- scanBamGRanges(obj, which,
                                 flag = scanBamFlag(isFirstMateRead = TRUE))
        message("Read GAlignments from BamFile...")
        ga <- readGAlignmentsFromBam(obj,
                                     param = ScanBamParam(which = which),
                                     use.names = use.name, ...)
        res.gp <- crunch(ga)
        message("Combine...")
        nms <- values(res.mb)$qname
        ## fow now just, using isize
        isize <- values(res.mb)$issize
        names(isize) <- nms
        values(res.gp)$isize <- isize[values(res.gp)$qname]
        res <- res.gp
    }
    res
})

scanBamGRanges <- function(file, which, what = c("rname", "strand", "pos", "qwidth"),
                           ...){
    ##  check
    .names <- c("rname", "strand", "pos", "qwidth")
    if(!all(.names %in% what)){
        what <- unique(c(what, .names))
    }
    dfnm <- setdiff(what, .names)
    param <- ScanBamParam(which = which, what = what,...)
    message("scanBam...")
    bam <- scanBam(file, param = param)
    message("Coerce list to GRanges")
    .names <- c("rname", "strand", "pos", "qwidth")
    res <- lapply(1:length(bam), function(i){
        bm <- bam[[i]]
        df <- as(bm, "DataFrame")[,dfnm]
        if(length(bm$rname))
            GRanges(bm$rname, IRanges(start = bm$pos, width = bm$qwidth),
                    strand = bm$strand, df)
        else
            GRanges()
    })
    bam.gr <- do.call("c", res)
}


