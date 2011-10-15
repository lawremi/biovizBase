\name{colorBlindSafePal}
\alias{colorBlindSafePal}
\alias{genBrewerBlindPalInfo}
\alias{genDichromatPalInfo}
\alias{genBlindPalInfo}
\alias{blind.pal.info}
\alias{brewer.pal.blind.info}
\alias{dichromat.pal.blind.info}
\title{Color blind safe palette generator}
\description{
  This function help users create a function based on specified color
  blind safe palette. And the returned function could be used for color
  generation. 
}
\usage{
genDichromatPalInfo()
genBrewerBlindPalInfo()
genBlindPalInfo()
colorBlindSafePal(palette)
blind.pal.info
brewer.pal.blind.info
dichromat.pal.blind.info
}
\arguments{
  \item{palette}{
    A index numeric value or character. Please see
    \code{blind.pal.info}, the palette could be "pal_id" or names the
    row in which users could specify the palette you want to use.
  }
}
\value{
  A color generating function with arguments \code{n} and
  \code{repeatable}. \code{n} specifying how many different discrete colors
you want to get from them palette, and if repeatable turned on and set
to TRUE, you can specify n even larger than maximum color. The color
will be repeated following the same order.
}
\details{
  We get those color-blind safe palette based on
  \url{http://colorbrewer2.org/} and
  \url{http://geography.uoregon.edu/datagraphics/color_scales.htm}
  those color are implemented in two packages, RColorBrewer and
  dichromat. But RColorBrewer doesn't provide subset of color-blind safe
  palette info. And dichromat doesn't group color based on "quality",
  "sequential" and "divergent", so we pick those color manually and
  create a combined palette, blind.pal.info.

  \code{colorBlindSafePal} will return a function, this function will
  take two arguments, \code{n} and \code{repeatable}. if n is smaller than
  3(n >= 3 is required by RColorBrewer), we use 3 instead and return a
  color vector. If n is over the \code{maxcolors} column in blind.pal.info, we
  use \code{maxcolors} instead when \code{repeatable} set to \code{FALSE}, if
  \code{repeatable} set to \code{TRUE} we repeat the color of all the
  allowed colors(length equals \code{maxcolors}) in the same order. This has
  special case in certain graphics which is always displayed side by
  side and don't worry about the repeated colors being neighbors.

  \code{genBrewerBlindPalInfo} return \code{brewer.pal.blind.info} data
  frame containing all color-blind safe palettes from
  \code{brewer.pal.info} defined in RColorBrewer, but it's not only just
  subset of it, it also changes some \code{maxcolors} info.
  
  \code{genDichromatPalInfo} return \code{dichromat.pal.blind.info} data
  frame.
  
  \code{genBlindPalInfo}  return \code{blind.pal.info} data frame.
}
\examples{
\dontrun{
library(scales)
## brewer subse of only color blind palette
brewer.pal.blind.info
genBrewerBlindPalInfo()
## dichromat info
dichromat.pal.blind.info
genDichromatPalInfo()
## all color blind palette, adding id/pkg.
blind.pal.info
## with no parameters, just return blind.pal.info
colorBlindSafePal()
mypal <- colorBlindSafePal(20)
## or pass character name
mypal <- colorBlindSafePal("Set2")
mypal12 <- colorBlindSafePal(22)
show_col(mypal(12, repeatable = FALSE)) # warning
show_col(mypal(11, repeatable = TRUE))  # no warning, and repeat
show_col(mypal12(12))
}
}
\author{Tengfei Yin}

