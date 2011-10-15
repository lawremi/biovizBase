\name{containLetters}
\alias{containLetters}
\title{Checking if string contains letters or not}
\usage{containLetters(obj, all=FALSE)
}
\description{Test if a string contain any letters}
\details{Useful when processing/sorting seqnames}
\value{Logical value}
\author{tengfei}
\arguments{\item{obj}{String}
\item{all}{If set to FALSE, return TRUE when any letters appears;
if all is set to TRUE, return TRUE only when the string is composed
of just letters.}
}
\examples{containLetters("XYZ123")
containLetters("XYZ123", TRUE)}

