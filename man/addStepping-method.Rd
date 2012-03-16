\name{addStepping-method}
\alias{addStepping}
\alias{addStepping,GenomicRanges-method}
\usage{
\S4method{addStepping}{GenomicRanges}(obj, group.name, extend.size = 0,
                       group.selfish = TRUE)
}
\title{Adding disjoint levels to a GenomicRanges object}
\description{Adding disjoint levels to a GenomicRanges object}
\details{This is a tricky question, for example, pair-end RNA-seq data
could be represented as two set of GenomicRanges object, one
indicates the read, one indicates the junction.  At the same time,
we need to make sure pair-ended read are shown on the same level,
and nothing falls in between. For better visualization of  the
data, we may hope to add invisible extended buffer to the reads, so
closely neighbored reads will be on the different levels.}
\arguments{\item{obj}{A GenomicRanges object}
\item{group.name}{Column name in the elementMetadata which specify
the grouping information of all the entries. If provided, this
will make sure all intervals belong to the same group will try to
be on the same level and nothing falls in between.}
\item{extend.size}{Adding invisible buffered region to the
GenomicRanges object, if it's 10, then adding 5 at both end. This
make the close neighbors assigned to the different levels and make
your eyes easy to identify.}
 \item{group.selfish}{
    Passed to \code{addStepping}, control whether to show each group as
    unique level or not. If set to \code{FALSE}, if two groups are not
    overlapped with each other, they will probably be layout in the same
    level to save space.
 }  
}
\value{A modified GenmicRanges object with \code{stepping} as one column.}
\docType{methods}
\author{Tengfei Yin}
\examples{
library(GenomicRanges)
set.seed(1)
N <- 500
## sample GRanges
gr <- GRanges(seqnames =
              sample(c("chr1", "chr2", "chr3", "chrX", "chrY"),
                     size = N, replace = TRUE),
              IRanges(
                      start = sample(1:300, size = N, replace = TRUE),
                      width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N,
                replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              group = sample(c("Normal", "Tumor"),
                size = N, replace = TRUE),
              pair = sample(letters, size = N,
                replace = TRUE))

## grouping and extending
head(addStepping(gr))
head(addStepping(gr, group.name = "pair"))
gr.close <- GRanges(c("chr1", "chr1"), IRanges(c(10, 20), width = 9))
addStepping(gr.close)
addStepping(gr.close, extend.size = 5)
}


