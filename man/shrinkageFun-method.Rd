\name{shrinkageFun-method}
\alias{shrinkageFun}
\alias{is_coord_truncate_gaps}
\usage{
## For IRanges
\S4method{shrinkageFun}{IRanges}(obj, max.gap = 1L)

## For GenomicRanges
\S4method{shrinkageFun}{GenomicRanges}(obj, max.gap = 1L)

is_coord_truncate_gaps(obj)
}
\title{Shrinkage function}
\description{Create a shrinkage function based on specified gaps and shrinkage rate.}
\details{
  \code{shrinkageFun} function will read in a GenomicRanges object
which represent the gaps, and return a function which works for
another GenomicRanges object, to shrink that object based on
previously specified gaps shrinking information.  You could use
this function to treat multiple tracks(e.g. GRanges) to make sure
they shrunk based on the common gaps and the same ratio.

\code{is_coord_truncate_gaps} is used to check if a \code{GRanges}
object is in "truncate_gaps" coordiantes or not.

}
\value{
A shrinkage function which could shrink a GenomicRanges
object, this function will add coord "truncate_gaps" and max.gap to metadata, ".ori" for oringal
data as extra column
}
\alias{shrinkageFun,GenomicRanges-method}
\alias{shrinkageFun,IRanges-method}
\author{Michael Lawrence, Tengfei Yin}
\arguments{\item{obj}{GenomicRanges object which represent gaps}
\item{max.gap}{Gaps to be kept, it's a fixed value, if this value
is bigger than certain gap interval, then that gap is not going to
be shrunk.}
}
\examples{
library(GenomicRanges)
gr1 <- GRanges("chr1", IRanges(start = c(100, 300, 600),
                               end = c(200, 400, 800)))

shrink.fun1 <- shrinkageFun(gaps(gr1), max.gap = maxGap(gaps(gr1), 0.1))
shrink.fun2 <- shrinkageFun(gaps(gr1, start(gr1), end(gr1)), max.gap = maxGap(gaps(gr1), 0.1))
shrink.fun3 <- shrinkageFun(gaps(gr1), max.gap = 0)
s1 <- shrink.fun1(gr1)
s2 <- shrink.fun2(gr1)
s3 <- shrink.fun3(gr1)
metadata(s1)$coord
values(s1)$.ori
is_coord_truncate_gaps(s1)
}
