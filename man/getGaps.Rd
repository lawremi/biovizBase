\name{getGaps}
\alias{getGaps}
\alias{getGaps,GRanges-method}
\title{get gaps for a faceted and grouped GRanges object}
\description{
  Get gaps for a faceted and grouped GRanges object, for visualization
  purpose. Please see details below for motivation.
}
\usage{
\S4method{getGaps}{GRanges}(obj, group.name = NULL, facets = NULL)
}
\arguments{
  \item{obj}{
    A \code{GRanges} object.
  }
  \item{group.name}{
    group name, such as transcript ID, this is the group method within
    each panel of facets and gaps will be computed for each group of intervals.
  }
  \item{facets}{
    formula used for creating graphics, all variables must be present in
    the data.
  }
}
\value{
  A \code{GRanges} object representing gaps and for each row, the
  "stepping" column help later visualization and make sure gaps and
  intervals they are generated from are showed on the expected place. 
}
\details{
  Since faceting is a subset and group process in visualization
  stage, some statistical computation need to be taken place after
  that. This leaves some computation like computing gaps hard based on
  solely \code{GRanges} object. Extra information like facets formula
  and group method would help to generate gaps which make sure they are
  aligned on the same level and within the same panel for another
  \code{GRanges}. facets variables will be added to gaps \code{GRanges}
  along with group.name.
}
\examples{
set.seed(1)
N <- 100
library(GenomicRanges)
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

grl <- splitByFacets(gr, sample ~ seqnames)
gr <- unlist(endoapply(grl, addStepping))
gr.gaps <- getGaps(gr, group.name = "stepping", facets = sample ~ seqnames)
}
\author{Tengfei Yin}

