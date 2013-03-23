\name{transformToGenome}
\alias{transformToGenome}
\alias{transformToGenome,GRanges-method}
\alias{transformToGenome,GRangesList-method}
\alias{transformToArch}
\alias{transformToArch,GRanges-method}
\alias{transformToCircle}
\alias{transformToRectInCircle}
\alias{transformToBarInCircle}
\alias{transformToSegInCircle}
\alias{transformToLinkInCircle}
\alias{transformDfToGr}
\alias{transformToDf}
\alias{transformToDf-method}
\alias{transformToDf,GRanges-method}
\alias{is_coord_genome}
\title{Transform GRanges to different coordinates and layout}
\description{
  Used for coordiante genome transformation, other transformation in
  circular view.
}
\usage{
\S4method{transformToGenome}{GRanges}(data, space.skip = 0.1, chr.weight
= NULL)
\S4method{transformToGenome}{GRangesList}(data, space.skip = 0.1,
chr.weight = NULL)

\S4method{transformToArch}{GRanges}(data, width = 1)
transformToCircle(data, x = NULL, y = NULL, ylim = NULL, 
                  radius = 10, trackWidth =10,
                  direction = c("clockwise", "anticlockwise"),
                   mul = 0.05)


transformToRectInCircle(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                      direction = c("clockwise", "anticlockwise"),
                      n = 100, mul = 0.05, chr.weight = NULL)

transformToBarInCircle(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                     direction = c("clockwise", "anticlockwise"),
                      n = 100, mul = 0.05, chr.weight = NULL)

transformToSegInCircle(data, y = NULL, space.skip = 0.1, trackWidth = 10, radius = 10,
                      direction = c("clockwise", "anticlockwise"), n =
                      100, chr.weight = NULL)


transformToLinkInCircle(data, linked.to, space.skip = 0.1, trackWidth = 10, radius = 10,
                      link.fun = function(x, y, n = 100) bezier(x, y, evaluation = n),
                      direction = c("clockwise", "anticlockwise"), chr.weight = NULL)

transformDfToGr(data, seqnames = NULL, start = NULL, end = NULL,
                            width = NULL, strand = NULL,
                            to.seqnames = NULL, to.start = NULL, to.end = NULL,
                            to.width = NULL, to.strand = NULL, linked.to
                            = to.gr)



\S4method{transformToDf}{GRanges}(data)

is_coord_genome(data)
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
  \item{ylim}{
    numeric range to control the ylimits on tracks when 'y' information
  is involved. 
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
    a column name of \code{GRanges} object, indicate the linked line's
  end point which represented as a \code{GRanges} too.. 
  }
  \item{link.fun}{
    function used to generate linking lines.
  }
  \item{seqnames}{
    character or integer values for column name or index indicate variable mapped to seqnames, default \code{NULL}
    use "seqnames".
  }
  \item{start}{
    character or integer values for column name or index indicate variable mapped to start, default \code{NULL}
    use "start".
  }
  \item{end}{
    character or integer values for column name or index indicate variable mapped to end, default \code{NULL}
  use "end".
  }
  \item{width}{
    character or integer values for column name or index indicate variable mapped to width, default \code{NULL}
  use "width".
    
  }
  \item{strand}{
    character or integer values for column name or index indicate variable mapped to strand, default \code{NULL}
    use "strand".
  }
  \item{to.seqnames}{
    character or integer values for column name or index indicate variable mapped to linked seqnames, default
  \code{NULL}, create GRanges without new GRanges attached as column. If
  this varialbe is not \code{NULL}, this mean you try to parse linked
  \code{GRanges} object.
  }
  \item{to.start}{
    character or integer values for column name or index indicate variable mapped to start of linked GRanges, default \code{NULL}
    use "to.start".
  }
  \item{to.end}{
    character or integer values for column name or index indicate variable mapped to end of linked GRanges, default \code{NULL}
    use "to.end".
  }
  \item{to.width}{
    character or integer values for column name or index indicate variable mapped to width of linked GRanges, default \code{NULL}
    use "to.width".
  }
  \item{to.strand}{
    character or integer values for column name or index indicate variable mapped to strand, default \code{NULL}
    use "to.strand" or just use *.
  }
  \item{chr.weight}{
    numeric vectors which sum to <1, the names of vectors has to be
  matched with seqnames in seqinfo, and you can only specify part of the
  seqnames, other lengths of chromosomes will be assined proportionally
  to their seqlengths, for example, you could specify chr1 to be 0.5, so
  the chr1 will take half of the space and other chromosomes squeezed to
  take left of the space. 
  }  
}
\value{
  A \code{GRanges} object, with calculated new variables, including
  ".circle.x" for transformed x position, ".circle.y" for transformed y
  position, ".circle.angle" for transformed angle.
}
\examples{
library(biovizBase)
library(GenomicRanges)
set.seed(1)
gr1 <- GRanges("chr1", IRanges(start = as.integer(runif(20, 1, 2e9)),
width = 5))
gr2 <- GRanges("chr2", IRanges(start = as.integer(runif(20, 1, 2e9)),
width = 5))
gr <- c(gr1, gr2)
gr.t <- transformToGenome(gr, space.skip = 0.1)
is_coord_genome(gr.t)
transformToCircle(gr.t)
transformToRectInCircle(gr)
transformToSegInCircle(gr)
values(gr1)$to.gr <- gr2
transformToLinkInCircle(gr1, linked.to = "to.gr")
}
\author{Tengfei Yin}

