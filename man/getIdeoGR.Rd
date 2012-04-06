\name{getIdeoGR}
\alias{getIdeoGR}
\title{Get ideogram information}
\description{
  This function tries to parse ideogram information from
  \code{seqlengths} of a \code{GRanges} object. This information is
  usually used to plot chromosome background for kaytogram or esitmate
  proper lengths of chromosomes from data space for showing overview. 
}
\usage{
   getIdeoGR(gr)
}
\arguments{
  \item{gr}{
    A \code{GRanges} object with or without lengths information.
  }
}
\value{
  A ideogram \code{GRanges} object, each row indicate one single
  chromosome, and start with 1 and end with real chromosome length or
  estimated lengths.
}
\examples{
   library(GenomicRanges)
   data("hg19IdeogramCyto", package = "biovizBase")
   hg19IdeogramCyto
   ## without seqlengths, simply reduce
   getIdeoGR(hg19IdeogramCyto)
   ## with seqlengths
   gr <- GRanges("chr1",  IRanges(1,3))
   seqlevels(gr) <- c("chr1", "chr2", "chr3")
   nms <- c(100, 200, 300)
   names(nms) <- c("chr1", "chr2", "chr3")
   seqlengths(gr) <- nms
   gr
   getIdeoGR(gr)
}
\author{Tengfei Yin}

