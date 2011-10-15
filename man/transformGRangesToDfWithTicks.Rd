\name{transformGRangesToDfWithTicks}
\alias{transformGRangesToDfWithTicks}
\title{transformGRangesToDfWithTicks}
\description{
  Transform GRanges object to a data.frame with ticks position stored in
  a list. This is different from simply coerce GRanges to a data frame,
  it is used for Grand Linear Compact View, which need to recompute the
  genomic coordinate and forcing them to be continued values.
}
\usage{
transformGRangesToDfWithTicks(gr, fixed.length = NULL)
}
\arguments{
  \item{gr}{
    A GRanges object.
  }
  \item{fixed.length}{
    Fixed length for arranging space. If NULL, use data range, which
    means it's simply concatenate the data based on seqnames. If it's
    provided, it's required that it's a numeric vector and names of the
    vector corresponding to the seqnames.
  }
}
\value{
  A list contains a data frame which with adjusted coordinates and a
  vector indicates the ticks position, which will be used for Grand
  Linear Compact View labeling in some graphics.
}
\details{
  For compact grand linear view, we need to recompute the coordinates
  and concatenate them in a linear manner, remove the x axis labels and relabel them with
  chromosome names. So this function serve the basis for doing such transformation.
}
\examples{
library(GenomicRanges)
gr <- GRanges("chr1", IRanges(seq(1,100, length.out = 10), width = 5))
transformGRangesToDfWithTicks(gr)
}
\author{Tengfei Yin}

