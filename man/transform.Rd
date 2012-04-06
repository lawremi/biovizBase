\name{transformtoGenome}
\alias{transformtoGenome}
\alias{transformToCircle}
\alias{transformToRectInCircle}
\alias{transformToBarInCircle}
\alias{transformToSegInCircle}
\alias{transformToLinkInCircle}
\alias{transformDfToGr}
\alias{transformToDf-method}
\alias{transformToDf,GRanges-method}
\title{Transform GRanges to genome coordinates and circular layout}
\description{
  
}
\usage{
transformToGenome(data, space.skip = 0.1)

transformToCircle(data, x = NULL, y = NULL,
                  radius = 10, trackWidth =10,
                  direction = c("clockwise", "anticlockwise"),
                   mul = 0.05)

transformToRectInCircle(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                      direction = c("clockwise", "anticlockwise"),
                      n = 5, mul = 0.05)

transformToBarInCircle(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                     direction = c("clockwise", "anticlockwise"))

transformToSegInCircle(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                      direction = c("clockwise", "anticlockwise"), n =
                    5)

transformToLinkInCircle(data, linked.to, space.skip = 0.1, trackWidth = 10, radius = 10,
                      link.fun = function(x, y, n = 100) bezier(x, y, evaluation = n),
                      direction = c("clockwise", "anticlockwise"))

transformDfToGr(data, seqnames = NULL, start = NULL, end = NULL,
                           width = NULL, strand = NULL)

\S4method{transformToDf}{GRanges}(data)

}
\arguments{
  \item{data}{
    a \code{GRanges} object.

    for function \code{transformDfToGr} it's data.frame.
  }
  \item{x}{
    character for variable as x axis used for transformation.
  }
  \item{y}{
    character for variable as y axis used for transformation.
  }
  \item{space.skip}{
    numeric values indicates skipped ratio of whole space, not skipped
    space is identical between each space.
  }
  \item{radius}{
    numeric value, indicates radius when transform to a circle.
  }
  \item{trackWidth}{
    numeric value, for track width.
  }
  \item{direction}{
    "clockwise" or "counterclockwise", for layout or transform direction to circle.
  }
  \item{mul}{
    numeric value, passed to \code{expand_range} function, to control
    margin of y in the track.
  }
  \item{n}{
    integer value, control interpolated points numbers.
  }
  \item{linked.to}{
    a \code{GRanges} object, indicate the linked line's end point. 
  }
  \item{link.fun}{
    function used to generate linking lines.
  }
  \item{seqnames}{
    character indicate variable mapped to seqnames, default \code{NULL}
    use "seqnames".
  }
  \item{start}{
    character indicate variable mapped to start, default \code{NULL}
    use "start".
  }
  \item{end}{
    character indicate variable mapped to end, default \code{NULL}
  use "end".
    
  }
  \item{width}{
    character indicate variable mapped to width, default \code{NULL}
  use "width".
    
  }
  \item{strand}{
    character indicate variable mapped to strand, default \code{NULL}
    use "strand".
  }
}
\value{
  A \code{GRanges} object, with calculated new variables, including
  ".circle.x" for transformed x position, ".circle.y" for transformed y
  position, ".circle.angle" for transformed angle.
}
\examples{
library(GenomicRanges)
set.seed(1)
gr1 <- GRanges("chr1", IRanges(start = as.integer(runif(20, 1, 100)),
width = 5))
gr2 <- GRanges("chr2", IRanges(start = as.integer(runif(20, 1, 100)),
width = 5))
gr <- c(gr1, gr2)
gr.t <- transformToGenome(gr, space.skip = 1)
transformToCircle(gr)
transformToRectInCircle(gr)
transformToBarInCircle(gr)
transformToSegInCircle(gr)
values(gr1)$to.gr <- gr2
transformToLinkInCircle(gr1, linked.to = "to.gr")
}
\author{Tengfei Yin}

