\name{pileupAsGRanges}
\alias{pileupAsGRanges}
\title{Summarize reads for certain region}
\usage{pileupAsGRanges(bams, regions, DNABases=c("A", "C", "G", "T", "N"), ...)
}
\description{This function summarize reads from bam files for nucleotides on
single base unit in a given region, this allows the downstream
mismatch summary analysis.}
\details{It's a wrapper around \code{applyPileup} function in Rsamtools
package, more detailed control could be found under manual of
PileupParam function in Rsamtools. \code{pileupAsGRanges} function
return a GRanges object which including summary of nucleotides,
depth, bam file path. This object could be read directly into
\code{pileupGRangesAsVariantTable} function for mismatch
summary.}
\value{A GRanges object, each row is one single base unit. and
elementMetadata contains summary about this position about all
nucleotides specified by DNABases. and \code{depth} for total
reads, \code{bam} for file path.}
\author{Michael Lawrence, Tengfei Yin}
\arguments{\item{bams}{A character which specify the bam file path.}
\item{regions}{A GRanges object specifying the region to be
summarized. This passed to \code{which} arguments in
\code{PileupParam}.}
\item{DNABases}{Nucleotide type you want to summarize in the result
and in specified order. It must be one or more of A,C,G,T,N.}
\item{...}{Extra parameters passed to \code{PileupParam}.}
}
\examples{
\dontrun{
library(Rsamtools)
data(genesymbol)
library(BSgenome.Hsapiens.UCSC.hg19)    
bamfile <- system.file("extdata", "SRR027894subRBM17.bam", package="biovizBase")
test <- pileupAsGRanges(bamfile, region = genesymbol["RBM17"])
test.match <- pileupGRangesAsVariantTable(test, Hsapiens)
head(test[,-7])
head(test.match[,-5])
}
}

