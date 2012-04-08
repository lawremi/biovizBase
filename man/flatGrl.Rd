\name{flatGrl}
\alias{flatGrl}
\title{Transform GRangesList to GRanges}
\description{
  Transform GRangesList to GRanges.
}
\usage{
flatGrl(object, indName = "grl_name")
}
\arguments{
  \item{object}{
    a \code{GRangesList} object.
  }
  \item{indName}{
    column named by 'indName' that groups the records by their
    original element in 'object'.
  }
}
\value{
  A \code{GRanges} object.
}
\details{
  This method is different from default \code{stack}, it integrate
  elementMetadata of \code{GRangesList} to the final coerced \code{GRanges}.
}
\examples{
library(GenomicRanges)
gr1 <- GRanges("chr1", IRanges(1:10, width = 5))
gr2 <- GRanges("chr2", IRanges(1:10, width = 5))
obj <- GRangesList(gr1, gr2)
values(obj) <- data.frame(a = 1:2, b = letters[1:2])
stack(obj)
flatGrl(obj)
}
\author{Tengfei Yin}

