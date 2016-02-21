setGeneric("crunch", function(obj, ...) standardGeneric("crunch"))
setMethod("crunch", "TxDb", function(obj, which,
                                     columns = c("tx_id", "tx_name","gene_id"),
                                     type = c("all", "reduce"),
                                     truncate.gaps = FALSE,
                                     truncate.fun = NULL,
                                     ratio = 0.0025){

    type <- match.arg(type)
    if(is.list(which)){
        message("Parsing exons based on which(list) arguments")
        temp <- exons(obj, vals = which, columns = columns)
        which <- range(temp)
    }

    seqnms <- as.character(unique(seqnames(which)))
    if(seqnms %in% seqlevels(obj)){
        seqlevels(obj, force = TRUE) <- seqnms
    }else{
        stop(seqnms, " is not matched with seqlevels of your object, please rename your 'which' arguments ")
    }
    ## seqlevels(obj, force = TRUE) <- seqnms
    on.exit(restoreSeqlevels(obj))  # needed only because TxDb are reference
                                    # objects (unlike R objects in general that
                                    # have a copy-on-change semantics)

    ## system.time(cdss <- cdsBy(obj, "tx")) #2.552s

    message("Parsing transcripts...")
    ## tx <- transcripts(obj, columns = columns)
    ## subsetByOverlaps(tx, which)

    tx <- transcriptsByOverlaps(obj, which, columns = columns)

    if(!length(tx)){
        message("No transcripts found at this region.")
        return(GRanges())
    }else{

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


        message("Parsing utrs...")
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
        message("------exons...")
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
        message("------cdss...")
        gr.cdss <- unlist(cdss)
        .nms <- as.character(gr.cdss$tx_id)
        .gid.nms <- mt[.nms, "gene_id"]
        .tx.nms <- mt[.nms, "tx_name"]
        if(length(gr.cdss)){
            values(gr.cdss) <- NULL
            values(gr.cdss) <- data.frame(tx_id = .nms, tx_name = .tx.nms,
                                          gene_id = .gid.nms, type = "cds")
            values(gr.cdss)$type <-     as.character(values(gr.cdss)$type)
            names(gr.cdss) <- NULL
        }else{
            gr.cdss <- GRanges()
        }

        ## intron
        message("------introns...")
        irl.introns <- gaps(ranges(exons))
        ir.introns <- unlist(irl.introns)
        if(length(ir.introns)){
            .nms <- names(ir.introns)
            .gid.nms <- mt[.nms, "gene_id"]
            .tx.nms <- mt[.nms, "tx_name"]
            gr.introns <- GRanges(seqnms, ir.introns, tx_id = .nms,
                                  tx_name = .tx.nms, gene_id = .gid.nms,
                                  type = "gap")
            names(gr.introns) <- NULL
            values(gr.introns)$type <-  as.character(values(gr.introns)$type)
        }else{
            gr.introns <- GRanges()
        }
        ## utrs
        message("------utr...")
        if(length(exons) && length(cdss)){
            txnms <- intersect(names(exons), names(cdss))
            irl.utrs <- setdiff(ranges(exons[txnms]), ranges(cdss[txnms]))
            ir.utrs <- unlist(irl.utrs)
            if(length(ir.utrs)){
                .nms <- names(ir.utrs)
                .gid.nms <- mt[.nms, "gene_id"]
                .tx.nms <- mt[.nms, "tx_name"]
                gr.utrs <- GRanges(seqnms, ir.utrs, tx_id = .nms,
                                   tx_name = .tx.nms, gene_id = .gid.nms,
                                   type = "utr")
                names(gr.utrs) <- NULL
                values(gr.utrs)$type <-     as.character(values(gr.utrs)$type)
            }else{
                gr.utrs <- GRanges()
            }
        }else{
            gr.utrs <- GRanges()
        }
        ## combine
        res <- c(gr.exons, gr.cdss, gr.introns, gr.utrs)

        return(reduceNtruncate(res, type=type, truncate.gaps=truncate.gaps,
                               truncate.fun=truncate.fun, ratio=ratio))
    }
})

setMethod("crunch", "GAlignments", function(obj, which,
                                            truncate.gaps = FALSE,
                                            truncate.fun = NULL,
                                            ratio = 0.0025){
    if(!missing(which))
        obj <- subsetByOverlaps(obj, which)
    if(!length(obj))
        stop("No entry in the GAlignments data")
    if(is.null(names(obj)))
        stop("Missing qname, please make use.names = TRUE, when reading GAlignments")
    grg <- grglist(obj)
    grg.u <- stack(grg, ".grl.name")
    message("extracting information...")
    values(grg.u)$junction <- rep(ifelse(elementNROWS(grg)>1, TRUE, FALSE),
                                  times = elementNROWS(grg))
    names(grg.u) <- NULL
    res <- grg.u
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
})

setMethod("crunch", "BamFile", function(obj, which, ...,
                                        type = c("gapped.pair", "raw", "all"),
                                        truncate.gaps = FALSE,
                                        truncate.fun = NULL,
                                        ratio = 0.0025){

    type <- match.arg(type)
    ## require(Rsamtools)
    if(type == "gapped.pair"){
        message("Read GAlignments from BamFile...")
        ga <- readGAlignments(obj, param = ScanBamParam(which = which),
                                   use.names = TRUE, ...)

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
        ga <- readGAlignments(obj, param = ScanBamParam(which = which),
                                   use.names = TRUE, ...)

        res.gp <- crunch(ga)
        message("Combine...")
        nms <- values(res.mb)$qname
        ## fow now just, using isize
        isize <- values(res.mb)$issize
        names(isize) <- nms
        values(res.gp)$isize <- isize[values(res.gp)$qname]
        res <- res.gp
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


####============================================================
##  crunch method from biovizBase
##
##  which can be a GRanges object or an object extending BasicFilter, or a
##  list of such filter objects.
####------------------------------------------------------------
setMethod("crunch", "EnsDb", function(obj, which,
                                      columns = c("tx_id",
                                                  "gene_name",
                                                  "gene_id"),
                                      type = c("all", "reduce"),
                                      truncate.gaps = FALSE,
                                      truncate.fun = NULL,
                                      ratio = 0.0025){
    type <- match.arg(type)
    ## which can be a single (!) GRanges an object extending BasicFilter or
    ## a list of BasicFilter objects.
    if(is(which, "list")){
        ## Has to be a list of objects extending BasicFilter.
        if(!all(vapply(which, is, logical(1L), "BasicFilter"))){
            stop("Which should be either a GRanges object, an object extending ",
                 "BasicFilter or a list of objects extending BasicFilter!")
        }
        exFilter <- which
    }else{
        if(!is(which, "GenomicRanges") & !is(which, "BasicFilter"))
                stop("Which should be either a GRanges object, an object extending ",
                     "BasicFilter or a list of objects extending BasicFilter!")
    }
    if(is(which, "GenomicRanges")){
        if(length(which) == 0){
            message("No transcripts found at this region.")
            return(GRanges())
        }
        if(length(which) > 1)
            stop("'which' has to be a single GRanges object.")
        if(!is.na(genome(which))){
            if(unname(genome(which)) != unique(unname(genome(obj))))
                stop("Genome versions do not fit! Argument 'which' has ",
                     unname(genome(which)), " argument 'obj' ",
                     unname(unique(genome(which))), "!")
        }
        ## Check if we've got the seqnames.
        if(!(seqlevels(which) %in% seqlevels(obj)))
            stop(seqlevels(which), " does not match any seqlevel in argument 'obj'!")
        exFilter <- GRangesFilter(which, condition="overlapping")
    }
    if(is(which, "BasicFilter"))
        exFilter <- which
    ## Check input argument 'columns':
    notAvailable <- !(columns %in% listColumns(obj))
    if(any(notAvailable)){
        if(all(notAvailable))
            stop("None of the columns specified by arguments 'columns' are available!",
                 " Allowed values are:", paste(listColumns(obj), collapse=", "), ".")
        ## Reducing to those which are allowed...
        warning("Columns ", paste(columns[notAvailable], collapse=", "), " are not",
                " available in the database and have been removed.")
        columns <- columns[!notAvailable]
    }
    ## Approach:
    ## 1) Get all exons in that region and retrieve also the tx_coding_seq_start.
    ## We're fetching the data just once and calculating everything we need from that,
    ## i.e. cds, utr and introns.
    message("Fetching data...", appendLF=FALSE)
    requiredCols <- c("tx_cds_seq_start", "tx_cds_seq_end", "exon_id", "tx_id")
    ## Forcing tx_id on the columns:
    columns <- unique(c(columns, "tx_id"))
    ## In order to solve also the "overlapping" condition I have first to fetch the transcripts
    ## in the region and their exons. That way we get, for a GRangesFilter or GRanges all
    ## transcripts that have an exon or an intron at the specified region.
    txInRegion <- transcripts(obj, filter=exFilter)
    if(length(txInRegion) == 0){
        message("No transcripts found at this region.")
        return(GRanges())
    }
    regExons <- exons(obj, filter=TxidFilter(unique(txInRegion$tx_id)),
                      columns=unique(c(requiredCols, columns)))
    ## Simple sanitizing: check if what we got is on the same chromosome!
    if(length(unique(as.character(seqnames(regExons)))) > 1)
        stop("Got features from different chromosomes! Please adjust argument",
             " 'which' in order to fetch only features from a single chromosome.")
    message("OK")
    if(length(regExons) > 0){
        message("Parsing exons...", appendLF=FALSE)
        ## Get an DataFrame that we can use as mcols for the GRanges:
        mDf <- unique(mcols(regExons)[, columns])
        ## Actually, with this I have all I need.
        ## 2) Define introns.
        regExonsList <- split(regExons, regExons$tx_id)
        message("OK\nDefining introns...", appendLF=FALSE)
        ints <- gaps(ranges(regExonsList))
        ints <- unlist(ints)
        if(length(ints) > 0){
            regIntrons <- GRanges(seqnames=unique(seqnames(regExons)), ranges=ints,
                                  strand="*",
                                  mDf[match(names(ints), mDf$tx_id), ])
            regIntrons$type <- "gap"
        }else{
            regIntrons <- GRanges()
        }
        message("OK\nDefining UTRs...", appendLF=FALSE)
        ## 3) Define UTRs.
        codingTx <- regExons[!is.na(regExons$tx_cds_seq_end)]
        if(length(codingTx) > 0){
            ## Define the whole CDS region per tx
            codingReg <- GRanges(seqnames=seqnames(codingTx),
                                 ranges=IRanges(codingTx$tx_cds_seq_start,
                                                codingTx$tx_cds_seq_end),
                                 strand=strand(codingTx),
                                 tx_id=codingTx$tx_id)
            codingTx <- split(codingTx, codingTx$tx_id)
            ## Subset to one CDS region per tx and split
            codingReg <- codingReg[match(names(codingTx), codingReg$tx_id)]
            codingReg <- split(codingReg, codingReg$tx_id)
            regUTRs <- unlist(setdiff(codingTx, codingReg))
            mcols(regUTRs) <- mDf[match(names(regUTRs), mDf$tx_id), ]
            regUTRs$type <- "utr"
            message("OK\nDefining CDS...", appendLF=FALSE)
            ## 4) Define CDS
            regCDSs <- unlist(intersect(codingTx, codingReg))
            mcols(regCDSs) <- mDf[match(names(regCDSs), mDf$tx_id), ]
            regCDSs$type <- "cds"
        }else{
            regUTRs <- GRanges()
            regCDSs <- GRanges()
        }
        message("OK\n", appendLF=FALSE)
        regExons <- regExons[, columns]
        regExons$type <- "exon"
        ensRes <- c(regExons, regIntrons, regUTRs, regCDSs)
    }else{
        warning("Did not find any transcript at the specified region!")
        return(GRanges())
    }
    return(reduceNtruncate(ensRes, type=type, truncate.gaps=truncate.gaps,
                           truncate.fun=truncate.fun, ratio=ratio))

})

## Just some simple helper function to avoid repeating code for TxDb and EnsDb.
reduceNtruncate <- function(res, type=c("all", "reduce"),
                            truncate.gaps=FALSE, truncate.fun=NULL,
                            ratio = 0.0025){
        message("aggregating...")
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
        if(truncate.gaps){
            message("truncating ...")
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
        message("Done")
        return(res)
}

