\name{strip_formula_dots}
\alias{strip_formula_dots}
\title{strip dots around a formula variables}
\description{
  strip dots around variables in a formula.
}
\usage{
  strip_formula_dots(formula)
}
\arguments{
  \item{formula}{
    A formula. such as coverage ~ seqnames, or ..coverage.. ~ seqnames.
  }
}
\value{
  A formula wihout ".." around for all variables.
}
\examples{
obj <- ..coverage.. ~ seqnames
strip_formula_dots(obj)
}
\author{Tengfei Yin}

