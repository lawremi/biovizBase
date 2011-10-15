\name{GCcontent}
\alias{GCcontent}
\title{GC content computation for BSgenome}
\usage{
GCcontent(obj, ..., view.width, as.prob = TRUE)
}
\description{Compute GC content in a certain region of a BSgenome object}
\details{GC content is an interesting variable may be related to various
  biological questions. So we need a way to compute GC content in a
  certain region of a reference genome.  GCcontent function is a
  wrapper around getSeq function in BSgenome package and
  \code{letterFrequency}, \code{letterFrequencyInSlidingView}
  in Biostrings package.

  if the \code{view.width} is specified, the GC content will be computed
  in the sliding view}
\value{Numeric value indicate count or percentage}
\author{Tengfei Yin}
\arguments{\item{obj}{BSgenome object}
  \item{...}{Arguments passed to getSeq method for BSgenome package.}
  \item{view.width}{ Passed to \code{letterFrequencyInSlidingView}, the
    constant (e.g. 35, 48, 1000) size of the "window" to slide along
    \code{obj}.  The specified \code{letters} are tabulated in each window
    of length \code{view.width}. The rows of the result (see value)
    correspond to the various windows.}
  \item{as.prob}{
  If TRUE return percentage of GC content, otherwise
  return counts.}
}
\examples{
library(BSgenome.Hsapiens.UCSC.hg19)
GCcontent(Hsapiens, GRanges("chr1", IRanges(1e6, 1e6 + 1000)))
}

