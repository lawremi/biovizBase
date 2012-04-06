\name{getFormalNames}
\alias{getFormalNames}
\title{Get formals from functions}
\description{
  Get formals from functions, used for dispatching arguments inside.
}
\usage{
getFormalNames(..., remove.dots = TRUE)
}
\arguments{
  \item{...}{
    functions.
  }
  \item{remove.dots}{
    logical value, indicate remove dots in formals or not, default is \code{TRUE}.
  }
}
\value{
  A character vector for formal arguments.
}
\examples{
getFormalNames(plot, sapply)
}
\author{Tengfei Yin}

