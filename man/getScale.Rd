\name{getScale}
\alias{getScale}
\title{Get scale information from a GRanges}
\description{
  Trying to get scale information from a \code{GRanges} object, used for
  circular view for geom "scale".
}
\usage{
getScale(gr, unit = NULL, n = 100, type = c("M", "B", "sci"))
}
\arguments{
  \item{gr}{
    a \code{GRanges} object.
  }
  \item{unit}{
    A numeric value for scale unit. Default \code{NULL} use argument
    \code{n} to estimate the unit.
  }
  \item{n}{
    Integer value to indicate how many scale ticks to make.
  }
  \item{
    unit types to shown.
  }
}
\value{
 A \code{GRanges} object, with extra column: "type" indicate it's longer
 major ticks or shorter minor ticks. "scale.y"  indicates y height for
 major and minor ticks. Default ratio is 3:1.
}
\author{Tengfei Yin}

