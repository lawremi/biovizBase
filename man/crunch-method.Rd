\name{crunch}
\alias{crunch}
\alias{crunch,TranscriptDb-method}
\alias{crunch,GAlignments-method}
\alias{crunch,BamFile-method}
\title{Fetching GRanges from various data source}
\description{
  Fetching Granges from various data source, currently supported objects
  are TranscriptDb, GAligment and BamFile. 
}
\usage{
\S4method{crunch}{TranscriptDb}(obj, which,
                 columns = c("tx_id", "tx_name","gene_id"),
                 type = c("all", "reduce"))
\S4method{crunch}{GAlignments}(obj, which)
\S4method{crunch}{BamFile}(obj, which, ...)
}
\arguments{
  \item{obj}{
    supported objects are \code{TranscriptDb}, \code{GAligment} and \code{BamFile}. 
  }
  \item{which}{
    GRanges object. For TranscriptDb object, could aslo be a list. 
  }
  \item{columns}{
    columns to include in the output.
  }
  \item{type}{
    default 'all' is to show the full model, 'reduce' is to show a
    single model. 
  }
  \item{...}{
    arguments passed to function \code{readGAlignmentsFromBam}.
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
temp <- crunch(txdb, which = genesymbol["BRCA1"], type = "single")
}
\author{Tengfei Yin}

