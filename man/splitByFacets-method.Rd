\name{splitByFacets}
\alias{splitByFacets}
\alias{splitByFacets,GRanges,formula-method}
\alias{splitByFacets,GRanges,GRanges-method}
\alias{splitByFacets,GRanges,NULL-method}
\alias{splitByFacets,GRanges,missing-method}
\title{split a GRanges by formula}
\description{
  Split a GRanges by formula into GRangesList. Parse variables in
  formula and form a interaction factors then split the \code{GRanges}
  by the factos.
}
\usage{
\S4method{splitByFacets}{GRanges,formula}(object, facets)
\S4method{splitByFacets}{GRanges,GRanges}(object, facets)
\S4method{splitByFacets}{GRanges,NULL}(object, facets)
\S4method{splitByFacets}{GRanges,missing}(object, facets)
}
\arguments{
  \item{object}{
    \code{GRanges} object.
  }
  \item{facets}{
    formula object, such as . ~ seqnames. Or \code{GRanges} object, or \code{NULL}.
  }
}
\value{
  A \code{GRangesList}.
}
\details{
  if facets is formula, factors are crreated based on interaction of
  variables in formula, then split it with this factor. If facets is
  \code{GRanges} object, it first subset the original data by facets
  \code{GRanges}, then split by each region in the facets.  If facets is
  \code{NULL}, split just by seqnames. This function is used to
  perform computation in conjunction with facets argments in higher
  level graphic function.
}
\examples{
library(GenomicRanges)
N <- 1000
gr <- GRanges(seqnames = 
              sample(c("chr1", "chr2", "chr3"),
                     size = N, replace = TRUE),
              IRanges(
                      start = sample(1:300, size = N, replace = TRUE),
                      width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N, 
                replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"), 
                size = N, replace = TRUE),
              pair = sample(letters, size = N, 
                replace = TRUE))
facets <- sample ~ seqnames
splitByFacets(gr, facets)
splitByFacets(gr)
gr.sub <- GRanges("chr1", IRanges(c(1, 200, 250), width = c(50, 10,
30)))
splitByFacets(gr, gr.sub)
}
\author{Tengfei Yin}

