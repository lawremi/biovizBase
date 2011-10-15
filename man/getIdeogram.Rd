\name{getIdeogram}
\alias{getIdeogram}
\title{Get ideogram.}
\usage{getIdeogram(genome, subchr, cytobands=TRUE)
}
\description{Get ideogram w/o cytoband for certain genome}
\details{This function require a network connection, it will parse the data
on the fly, function is a wrapper of some functionality from
rtracklayer package to get certain table like cytoBand, a full
table schema could be found
\url{http://genome.ucsc.edu/cgi-bin/hgTables} in UCSC genome
browser.This is useful for visualization of the whole genome or
single chromosome, you can see some examples in ggbio package.}
\value{A GRanges object.}
\author{Tengfei Yin}
\arguments{\item{genome}{Single specie names, which must be one of the result from
\code{ucscGenomes()$db}. If missing, will invoke a menu for users to choose
from.}
\item{subchr}{A character vector used to subset the result.}
\item{cytobands}{If TRUE, return ideogram with gieStain and name
column. If FALSE, simply return the genome information for each
chromosome.}
}
\examples{
\dontrun{hg19IdeogramCyto <- getIdeogram("hg19", cytoband = TRUE)}
}

