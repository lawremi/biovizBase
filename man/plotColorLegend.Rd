\name{plotColorLegend}
\alias{plotColorLegend}
\title{Show colors}
\usage{plotColorLegend(colors, labels, title)
}
\description{Plot color legend, simple way to check your default color scheme}
\details{Show color sheme as a legend style, labels}
\value{A graphic device showing color legend.}
\author{Tengfei Yin}
\arguments{\item{colors}{A character vector of colors}
\item{labels}{Labels to put aside colors, if \code{missing}, use names of
the colors character vector.}
\item{title}{Title for the color legend.}
}
\examples{cols <- getOption("biovizBase")$baseColor
plotColorLegend(cols, title = "strand legend")}

