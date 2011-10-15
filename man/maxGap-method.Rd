\name{maxGap-method}
\alias{maxGap}
\title{Estimated max gaps}
\usage{
\S4method{maxGap}{GenomicRanges}(obj, ratio = 0.0025)
}
\description{Compute an estimated max gap information, which could
be used as max.gap argument in \code{shringkageFun}.}
\details{This function tries to estimate an appropriate max gap
to be used for creating a shrinkage function.}
\value{A numeric value}
\author{Tengfei Yin}
\docType{methods}
\alias{maxGap,GenomicRanges-method}
\seealso{\code{\link{shrinkageFun}}}
\arguments{\item{obj}{GenomicRanges object}
\item{ratio}{Multiple by the range of the provided gaps as the max gap.}
}
\examples{require(GenomicRanges)
gr1 <- GRanges("chr1", IRanges(start = c(100, 300, 600),
end = c(200, 400, 800)))
gr2 <- GRanges("chr1", IRanges(start = c(100, 350, 550),
end = c(220, 500, 900)))
gaps.gr <- intersect(gaps(gr1, start = min(start(gr1))),
gaps(gr2, start = min(start(gr2))))
shrink.fun <- shrinkageFun(gaps.gr, max.gap = maxGap(gaps.gr))
shrink.fun(gr1)
shrink.fun(gr2)}

