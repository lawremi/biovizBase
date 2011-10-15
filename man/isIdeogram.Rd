\name{isIdeogram}
\alias{isIdeogram}
\title{Ideogram checking}
\usage{isIdeogram(obj)
}
\description{Check if an object is ideogram or not}
\details{Simply test if it's the result coming
from \code{getIdeogram} function or not, make sure it's a
\code{GRanges} and with extra column}
\value{A logical value to indicate it's a ideogram or not.}
\author{Tengfei Yin}
\arguments{\item{obj}{object}
}
\examples{data(hg19IdeogramCyto)
data(hg19Ideogram)
isIdeogram(hg19IdeogramCyto)
isIdeogram(hg19Ideogram)}

