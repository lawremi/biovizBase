##' biovizBase is a package which provides utilities and color scheme
##' for higher level graphic package which aim to visualize biological
##' data especially genetic data.
##'
##' This package provides default color scheme for nucleotide, strand,
##' amino acid, try to pass colorblind checking as possible as we can.
##' And also provide giemsa stain result color scheme used to show
##' cytoband. This package also provides utilites to manipulate and
##' summarize raw data to get them ready to be visualized.
##' @name biovizBase-package
##' @aliases biovizBase
##' @docType package
NULL

##' Gene symbols with position
##'
##' This data set provides genen symbols in human with position
##' and starnd information, stored as GRanges object.
##' @name genesymbol
##' @docType data
##' @usage data(genesymbol)
##' @format GRanges: length of 16914 
##' @keywords datasets
##' @examples
##' data(genesymbol)
##' head(genesymbol)
NULL

##' Hg19 ideogram without cytoband information
##'
##' This data set provides hg19 genome information wihout cytoband
##' information.
##' @name hg19Ideogram
##' @docType data
##' @usage data(hg19Ideogram)
##' @format GRanges: length of 24
##' @keywords datasets
##' @examples
##' data(hg19Ideogram)
##' hg19Ideogram
NULL

##' Hg19 ideogram with cytoband information
##'
##' This data set provides hg19 genome information with cytoband
##' information.
##' @name hg19IdeogramCyto
##' @docType data
##' @usage data(hg19IdeogramCyto)
##' @format GRanges: length of 862
##' @keywords datasets
##' @examples
##' data(hg19IdeogramCyto)
##' hg19IdeogramCyto
NULL

##' Subset of RNA editing sites in hg19
##'
##' This data set provides a subset(500 sites only) of hg19 RNA editing sites, and
##' originally from DARNED \url{http://darned.ucc.ie/} for the hg19
##' assembly. 
##' @name darned_hg19_subset500
##' @docType data
##' @usage data(darned_hg19_subset500)
##' @format GRanges: length of 500
##' @keywords datasets
##' @examples
##' data(darned_hg19_subset500)
##' darned_hg19_subset500
NULL


