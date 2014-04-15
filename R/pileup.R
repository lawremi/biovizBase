##' This function summarize reads from bam files for nucleotides on
##' single base unit in a given region, this allows the downstream
##' mismatch summary analysis.
##'
##' It's a wrapper around \code{applyPileup} function in Rsamtools
##' package, more detailed control could be found under manual of
##' PileupParam function in Rsamtools. \code{pileupAsGRanges} function
##' return a GRanges object which including summary of nucleotides,
##' depth, bam file path. This object could be read directly into
##' \code{pileupGRangesAsVariantTable} function for mismatch
##' summary.
##' @title Summarize reads for certain region
##' @param bams A character which specify the bam file path.
##' @param regions A GRanges object specifying the region to be
##' summarized. This passed to \code{which} arguments in
##' \code{PileupParam}.
##' @param DNABases Nucleotide type you want to summarize in the result
##' and in specified order. It must be one or more of A,C,G,T,N.
##' @param ... Extra parameters passed to \code{PileupParam}.
##' @return A GRanges object, each row is one single base unit. and
##' elementMetadata contains summary about this position about all
##' nucleotides specified by DNABases. and \code{depth} for total
##' reads, \code{bam} for file path.
##' @author  Michael Lawrence, Tengfei Yin
pileupAsGRanges <- function(bams, regions,
                            DNABases = c("A", "C", "G", "T", "N"),
                            ...) {
  pileupFiles <- PileupFiles(bams)
  pileupParams <- PileupParam(which = regions, ...)
  ## what if it's multiple bam files
  bamNames <- names(bams)
  if (is.null(bamNames))
    bamNames <- bams
  pileupFun <- function(x) {
    grl <- lapply(seq_len(length(bamNames)), function(i) {
      seq <- x$seq[DNABases,i,]
      sn <- unique(names(x$seqnames))
      GRanges(seqnames = Rle(names(x$seqnames), x$seqnames),
              ranges = IRanges(x$pos, width = rep(1L, length(x$pos))),
              strand = Rle("+", length(x$pos)),
              t(seq),
              depth = colSums(seq), bam = Rle(bamNames[i], length(x$pos)))
    })
    do.call(c, grl[!is.null(grl)])
  }
  gr <- do.call(c, applyPileups(pileupFiles, pileupFun, param = pileupParams))
  seqinfo(gr) <- seqinfo(regions)
  gr
}

##' Compare to reference genome and compute mismatch summary for
##' certain region of reads.
##'
##' User need to make sure to pass the right reference genome to this
##' function to get the right summary. This function drop the position
##' has no reads and only keep the region with coverage in the
##' summary. The result could be used to show stacked barchart for
##' mismatch summary.
##' @title Mismatch summary
##' @param gr A GRanges object, with nucleotides summary, each base
##' take one column in elementMetadata or user can simply passed the
##' returned result from \code{pileupAsGRanges} function to this
##' function.
##' @param genome BSgenome object, need to be the reference genome.
##' @param DNABases Nucleotide types contained in passed GRanges
##' object. Default is A/C/G/T/N, it tries to match the column names
##' in elementMetadata to those default nucleotides. And treat the
##' matched column as base names.
##' @return A GRanges object. Containing the following elementMetadata
##' \describe{
##'  \itemize{
##'  \item {[ref] Nucleotide in reference genome.}
##'  \item {[read] Nucleotide contained in the reads at particular position,
##'  if multiple nucleotide, either matched or unmatched are found, they will be
##' summarized in different rows.}
##'  \item {[count] Count for read column.}
##'  \item{[match] Logical value, whether matched to reference genome or not}
##'  \item{[bam] Character indicate bam file path.}
##' }
##' }
##' @author Michael Lawrence, Tengfei Yin
pileupGRangesAsVariantTable <- function(gr, genome,
                                        DNABases = c("A", "C", "G", "T", "N")) {
  require(BSgenome)
  colnms <- colnames(values(gr))
  DNAbases <- colnms[colnms %in% DNABases]
  gr <- rectifySeqlevelsStyle(gr, genome)
  refBases <- getSeq(genome, gr, as.character = TRUE)
  counts <- as.matrix(as.data.frame(values(gr)[,DNABases]))
  refCounts <- counts[cbind(seq(nrow(counts)), match(refBases, DNABases))]
  variantsForBase <- function(base) {
    baseCounts <- counts[,base, drop=TRUE]
    keep <- baseCounts > 0
    if(sum(keep)){
      count <- baseCounts[keep]
      gr <- GRanges(seqnames(gr)[keep], ranges(gr)[keep], strand(gr)[keep],
                    ref = refBases[keep], read = base,
                    count = baseCounts[keep],
                    values(gr)[keep, setdiff(colnames(values(gr)), DNABases)])
      values(gr)$match <- values(gr)$ref == values(gr)$read
      gr
    }else{
      NULL
    }
  }
  lst <- lapply(colnames(counts), variantsForBase)
  idx <- !sapply(lst, is.null)
  res <- do.call("c", lst[idx])  ## why this is slow
  if(length(res))
    return(res[order(as.numeric(start(res)))])
  else
    GRanges()
}
