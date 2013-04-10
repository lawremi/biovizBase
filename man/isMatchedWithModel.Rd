\name{isMatchedWithModel}
\alias{isMatchedWithModel}
\alias{isJunctionRead}
\title{Utils for Splicing Summary}
\description{
  Utilities used for summarizing isoforms
}
\usage{
isJunctionRead(cigar)
isMatchedWithModel(model, gr)
}
\arguments{
  \item{cigar}{
    A CIGAR string vector.
  }
  \item{model}{
    A GRanges object.
  }
  \item{gr}{
    A GRanges object.
  }
}
\value{
  Logical vectors.
}
\details{
  \code{isJunctionRead} simply parsing the CIGAR string to see if there
  is "N" in between and return a logical vector of the same length as
  cigar parameters, indicate it's junction read or not.

  \code{isMatchedWithModel} mapping \code{gr} to \code{model}, and counting overlapped
  cases for each row of model, If \code{gr} contains all the read, this
  will return a logical vector of the same length as \code{gr}, and
  indicate if each read is the support for this model. NOTICE: we only
  assume it's a full model, so each \code{model} here is simply one
  isoform. So we only treat the gaped reads which only overlapped with two
  consecutive exons in \code{model} as one support for it.
}
\examples{
library(Rsamtools)
bamfile <- system.file("extdata", "SRR027894subRBM17.bam",
  package="biovizBase")
## get index of junction read
which(isJunctionRead(cigar(readGAlignmentsFromBam(bamfile))))
##
model <- GRanges("chr1", IRanges(c(10, 20, 30, 40), width  = 5))
gr <- GRanges("chr1", IRanges(c(10, 10, 12, 22, 33), c(31, 40, 22, 32,
  44)))
isMatchedWithModel(model, gr)
}
\author{Tengfei Yin}

