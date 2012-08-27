### =========================================================================
### Coverage utilities
### -------------------------------------------------------------------------

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Estimation of coverage
###

setGeneric("estimateCoverage",
           function(x, ...) standardGeneric("estimateCoverage"))

setMethod("estimateCoverage", "BamFile", function(x, maxBinSize = 2^14) {
  bai_path <- paste(index(x), ".bai", sep = "")
  bai <- file(bai_path, "rb")
  on.exit(close(bai))
  bytes <- readBin(bai, "raw", n = file.info(bai_path)[,"size"])
  bam_si <- seqinfo(x)

  bin_offsets <- .Call(scan_bam_bin_offsets, bytes)
  names(bin_offsets) <- seqlevels(bam_si)

  ## CONSTANTS
  UBSIZE <- 2^16 # hopefully is not going to change
  CBSIZE <- 19100 # determined empirically
  BINSIZES <- 2 ^ seq(29, 14, -3)
  TOTALSIZE <- 2 ^ 29
  NBIN <- sum(TOTALSIZE / BINSIZES)
  METABIN <- NBIN + 1
  
  makeBinRanges <- function() {
    unlist(seqapply(BINSIZES[BINSIZES <= maxBinSize], function(x) {
      chunks <- breakInChunks(TOTALSIZE, x)
      IRanges(start(chunks), end(chunks))
    }))
  }

  bin_ranges <- makeBinRanges()
  
  bin_offsets <- lapply(bin_offsets, function(x) {
    x[,x[1,] != METABIN] # drop the samtools metadata bin
  })
  
  off_scores <- lapply(seqlevels(bam_si), function(sn) {
    off <- t(bin_offsets[[sn]])
    colnames(off) <- c("bin", "start.coffset", "end.coffset",
                       "start.uoffset", "end.uoffset")
    coffdiff <- off[,"end.coffset"] - off[,"start.coffset"]
    off_diff <- ifelse(coffdiff > 0,
                       UBSIZE - off[,"start.uoffset"] +
                       UBSIZE * pmax(0, round(coffdiff / CBSIZE - 1)) +
                       off[,"end.uoffset"],
                       off[,"end.uoffset"] - off[,"start.uoffset"])
    ## sum bins with multiple chunks
    off_diff <- rowsum(as.numeric(off_diff), off[,"bin"])
    score <- numeric(NBIN)
    score[as.integer(rownames(off_diff)) + 1L] <- as.vector(off_diff)
    score[rep(BINSIZES, TOTALSIZE / BINSIZES) <= maxBinSize]
  })

  ## for performance, important to construct GRanges outside of loop
  sn_bin <- rep(seqlevels(bam_si), each = length(bin_ranges))
  bin_ranges <- restrict(rep(bin_ranges, length(bam_si)),
                         end = seqlengths(bam_si)[sn_bin],
                         keep.all.ranges = TRUE)
  gr <- GRanges(sn_bin, bin_ranges, score = unlist(off_scores),
                seqlengths = seqlengths(bam_si))
  gr[width(gr) > 0]
})
