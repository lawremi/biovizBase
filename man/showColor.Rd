\name{showColor}
\alias{showColor}
\title{Show colors}
\description{
  Show colors with color string or names of color vectors.
}
\usage{
showColor(colors, label = c("color", "name"))
}
\arguments{
  \item{colors}{
    A color character vector.
  }
  \item{label}{
    "color" label color with simple color string, and "name" label color
    with names of the vectors.
  }
}
\value{
  A plot.
}
\examples{
\dontrun{
showColor(getBioColor("CYTOBAND"))
}
}
\author{Tengfei Yin}

