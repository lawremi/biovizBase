\name{crunch}
\alias{crunch}
\alias{crunch,TxDb-method}
\alias{crunch,EnsDb-method}
\alias{crunch,GAlignments-method}
\alias{crunch,BamFile-method}
\title{Fetching GRanges from various data source}
\description{
  Fetching Granges from various data source, currently supported objects
  are TxDb, EnsDb, GAlignments and BamFile.
}
\usage{
\S4method{crunch}{TxDb}(obj, which, columns = c("tx_id", "tx_name","gene_id"),
       type = c("all", "reduce"), truncate.gaps = FALSE,
       truncate.fun = NULL, ratio = 0.0025)
\S4method{crunch}{EnsDb}(obj, which, columns = c("tx_id", "gene_name","gene_id"),
       type = c("all", "reduce"), truncate.gaps = FALSE,
       truncate.fun = NULL, ratio = 0.0025)
\S4method{crunch}{GAlignments}(obj, which, truncate.gaps = FALSE,
       truncate.fun = NULL, ratio = 0.0025)
\S4method{crunch}{BamFile}(obj, which, ..., type = c("gapped.pair", "raw", "all"),
       truncate.gaps = FALSE, truncate.fun = NULL, ratio = 0.0025)
}
\arguments{
  \item{obj}{
    supported objects are \code{TxDb}, \code{GAlignments} and \code{BamFile}.
  }
  \item{which}{
    GRanges object. For TxDb object, could aslo be a list.
    For \code{EnsDb}, it can also be a single filter object
    extending \code{\link[ensembldb]{BasicFilter-class}} or a
    \code{list} of such objects.
  }
  \item{columns}{
    columns to include in the output.
  }
  \item{type}{
    default 'all' is to show the full model, 'reduce' is to show a
    single model.
  }
  \item{truncate.gaps}{
    logical value, default \code{FALSE}. Whether to truncate gaps or
    not.
  }
  \item{truncate.fun}{
    shrinkage function.
  }
  \item{ratio}{
    numeric value, shrinking ratio.
  }
  \item{...}{
    arguments passed to function \code{readGAlignments}.
  }
}
\value{
  GRanges object.
}
\examples{
library(biovizBase)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
obj <- txdb
temp <- crunch(txdb, which = genesymbol["BRCA1"], type = "all")
temp <- crunch(txdb, which = genesymbol["BRCA1"], type = "reduce")

## Fetching data from a EnsDb object.
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
gr <- genesymbol["BRCA1"]
seqlevels(gr) <- sub(seqlevels(gr), pattern="chr", replacement="")
temp <- crunch(edb, which=gr)

## Alternatively, use the GenenameFilter from the ensembldb package:
temp <- crunch(edb, which=GenenameFilter("BRCA1"))

}
\author{Tengfei Yin}

