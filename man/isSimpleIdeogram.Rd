\name{isSimpleIdeogram}
\alias{isSimpleIdeogram}
\title{Simple ideogram checking}
\usage{isSimpleIdeogram(obj)
}
\description{Check if an object is a simple ideogram or not}
\details{test if it's GRanges or not, doesn't require cytoband
information. But it double check to see if there is only one entry
per chromosome}
\value{A logical value to indicate it's a simple ideogram or not.}
\author{tengfei}
\arguments{\item{obj}{object}
}
\examples{data(hg19IdeogramCyto)
data(hg19Ideogram)
isSimpleIdeogram(hg19IdeogramCyto)
isSimpleIdeogram(hg19Ideogram)}

