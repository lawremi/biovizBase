\name{getXScale}
\alias{getXScale}
\alias{getXScale,GRanges-method}
\title{get x scale breaks and labels}
\description{
  get x scale breaks and labels for GRanges with different
  coordintes(currently only "truncate_gaps" and "genome" supported).
}
\usage{
 \S4method{getXScale}{GRanges}(obj, type = c("default", "all", "left", "right"))
}
\arguments{
  \item{obj}{
    a \code{GRanges} object. "coord" in metadata shows proper
    coordinates transformation for this object.
  }
  \item{type}{
    types of labels for transformed data.
  }
}
\value{
  list of breaks and labels.
}
\examples{
library(GenomicRanges)
gr1 <- GRanges("chr1", IRanges(start = c(100, 300, 600),
                               end = c(200, 400, 800)))
shrink.fun1 <- shrinkageFun(gaps(gr1), max.gap = maxGap(gaps(gr1), 0.15))
shrink.fun2 <- shrinkageFun(gaps(gr1), max.gap = 0)
s1 <- shrink.fun1(gr1)
getXScale(s1)
#  coord:genome
set.seed(1)
gr1 <- GRanges("chr1", IRanges(start = as.integer(runif(20, 1, 100)),
width = 5))
gr2 <- GRanges("chr2", IRanges(start = as.integer(runif(20, 1, 100)),
width = 5))
gr <- c(gr1, gr2)
gr.t <- transformToGenome(gr, space.skip = 1)
getXScale(gr.t)
}
\author{Tengfei Yin}

