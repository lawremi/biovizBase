\name{estimateCoverage}
\alias{estimateCoverage}
\alias{estimateCoverage,BamFile-method}
\title{Estimation of Coverage}
\description{
  Estimation of Coverage
}
\usage{
\S4method{estimateCoverage}{BamFile}(x, maxBinSize = 2^14)
}
\arguments{
  \item{x}{
    A \code{BamFile} object.
  }
  \item{maxBinSize}{
    Max bin size.
  }
}
\value{
  A \code{GRanges} object.
}
\author{Michael Lawrence}

