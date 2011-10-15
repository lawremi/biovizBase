\name{pileupGRangesAsVariantTable}
\alias{pileupGRangesAsVariantTable}
\title{Mismatch summary}
\usage{pileupGRangesAsVariantTable(gr, genome, DNABases=c("A", "C", "G", "T", "N"))
}
\description{Compare to reference genome and compute mismatch summary for
certain region of reads.}
\details{User need to make sure to pass the right reference genome to this
function to get the right summary. This function drop the position
has no reads and only keep the region with coverage in the
summary. The result could be used to show stacked barchart for
mismatch summary.}
\value{A GRanges object. Containing the following elementMetadata
\itemize{
\item {ref}{Nucleotide in reference genome.}
\item {read}{Nucleotide contained in the reads at particular position,
if multiple nucleotide, either matched or unmatched are found, they will be
summarized in different rows.}
\item {count}{Count for read column.}
\item{match}{Logical value, whether matched to reference genome or not}
\item{bam}{Character indicate bam file path.}
}}
\author{Michael Lawrence, Tengfei Yin}
\arguments{\item{gr}{A GRanges object, with nucleotides summary, each base
take one column in elementMetadata or user can simply passed the
returned result from \code{pileupAsGRanges} function to this
function.}
\item{genome}{BSgenome object, need to be the reference genome.}
\item{DNABases}{Nucleotide types contained in passed GRanges
object. Default is A/C/G/T/N, it tries to match the column names
in elementMetadata to those default nucleotides. And treat the
matched column as base names.}
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

