\name{subsetArgsByFormals}
\alias{subsetArgsByFormals}
\title{Subset list of arguments by functions}
\description{
  find arguments matched by formals of passed functions, 
}
\usage{
subsetArgsByFormals(args, ..., remove.dots = TRUE)
}
\arguments{
  \item{args}{
    list of arguments with names indicate the formals.
  }
  \item{...}{
    functions used to parse formals.
  }
  \item{remove.dots}{
    logical value indicate whether to include dots in formals or not.
  }
}
\value{
  argumnets that matched with passed function.
}
\examples{
args <- list(x = 1:3, simplify = TRUE, b = "b")
subsetArgsByFormals(args, plot, sapply)
}
\author{Tengfei Yin}
