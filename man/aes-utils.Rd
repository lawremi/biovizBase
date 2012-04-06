\name{parseArgsForAes}
\alias{parseArgsForAes}
\alias{parseArgsForNonAes}
\title{Utils for parsing (un)evaluated arguments}
\description{
  
}
\usage{
parseArgsForAes(args)
parseArgsForNonAes(args)
}
\arguments{
  \item{args}{
    arguments list.
  }
}
\value{
  For \code{parseArgsForAes} return a unevaluated arguments.

  For \code{parseArgsNonForAes} return a evaluated/quoted arguments.
}
\examples{
args <- alist(a = color, b = "b")
# parseArgsForAes(args)
}
\author{Tengfei Yin}

