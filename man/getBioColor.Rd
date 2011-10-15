\name{getBioColor}
\alias{getBioColor}
\title{Color scheme getter for biological data}
\description{
  This function tries to get default color scheme either from fixed
  palette or options for specified data set, usually just biological
  data.
}
\usage{
getBioColor(type = c("DNA_BASES_N", "DNA_BASES", "DNA_ALPHABET",
                     "RNA_BASES_N", "RNA_BASES", "RNA_ALPHABET",
                     "IUPAC_CODE_MAP", "AMINO_ACID_CODE", "AA_ALPHABET",
                     "STRAND", "CYTOBAND"),
            source = c("option", "default"))
}
\arguments{
  \item{type}{
    Color set based on which you want to get.
  }
  \item{source}{
    "option" tries to get color scheme from \code{Options}. This allow
    user to edit the color globally.
    "default" gets fixed color scheme.
  }
}
\value{
  A character of vector contains color(rgb), and the names of the vector
  is originally from the name of different data set. e.g. for DNA_BASES,
  it's just \code{A,C,T,G}. This allow users to get color for a vector
  of specified names. Please see the examples below.
}
\details{
  Most data set specified in the type argument are defined in Biostrings
  package, such as "DNA_BASES", "DNA_ALPHABET", "RNA_BASES",
  "RNA_ALPHABET", "IUPAC_CODE_MAP", "AMINO_ACID_CODE", "AA_ALPHABET", please check the
  manual for more details.

  "DNA_BASES_N" is just "DNA_BASES" with extra "N" used in certain
  cases, like the result from \code{applyPileup} in Rsamtools package.
  We start with the five most used nucleotides, \code{A,T,C,G,N}, In
  genetics, \code{GC-content} usually has special biological
  significance because GC pair is bound by three hydrogen bonds instead
  of two like AT pairs. So it has higher thermostability which could
  result in different significance, like higher annealing temperature in
  PCR. So we hope to choose warm color for \code{G,C} and cold color for
  \code{A,T}, and a color in between to represent \code{N}.  They are
  chosen from diverging color set of \code{color brewer}. So we should
  be able to easily tell the GC enriched region. This set of color also
  passed vischeck for colorblind people.

  In \code{GRanges} object, we have \code{strand} which contains three
  levels, +, -, *. We are using a qualitative color set from \code{Color
    Brewer}. This color set pass the colorblind test. It should be a safe
  color set to use to color strand.

  For most default color scheme if they are under 18, we are trying to
  use package dichromat to set color for color blind people. But for
  data set that contains more than 18 objects, it's not possible to
  assign colorblind-safe color to them anymore. So we need to repeat
  some color. It should not matter too much, because even normal people
  cannot tell the difference anymore.

  Here are the definition for the data sets.

  \describe{
    \item{DNA_BASES}{
      Contains \code{A,C,T,G}.
    }
    \item{DNA_ALPHABET}{
     This alphabet contains all letters from the IUPAC Extended Genetic
     Alphabet (see "?IUPAC_CODE_MAP") + the gap ("-") and the hard
     masking ("+") letters.      
    }
    \item{DNA_BASES_N}{
      Contains \code{A,C,T,G,N}
    }
    \item{RNA_BASES_N}{
      Contains \code{A,C,U,G,N}
    }
    \item{RNA_BASES}{
      Contains \code{A,C,T,G}
    }
    \item{RNA_ALPHABET}{
     This alphabet contains all letters from the IUPAC Extended Genetic
     Alphabet (see ?IUPAC_CODE_MAP) where "T" is replaced by "U"
     + the gap ("-") and the hard masking ("+") letters. 
    }
    \item{IUPAC_CODE_MAP}{
      The "IUPAC_CODE_MAP" named character vector contains the mapping
     from the IUPAC nucleotide ambiguity codes to their meaning.
   }
    \item{AMINO_ACID_CODE}{
      Single-Letter Amino Acid Code (see "?AMINO_ACID_CODE"). 
    }
    \item{AA_ALPHABET}{
      This alphabet contains all letters from the Single-Letter Amino
     Acid Code (see "?AMINO_ACID_CODE") + the stop ("*"), the gap
     ("-") and the hard masking ("+") letters      
    }
    \item{STRAND}{
      Contains \code{"+", "-", "*"}
    }
    \item{CYTOBAND}{
      Contains giemsa stain results:\code{gneg, gpos25, gpos50, gpos75,
	gpos100, gvar, stalk, acen}. Color defined in package geneplotter.
    }
  }
}
\examples{
opts <- getOption("biovizBase")
opts$DNABasesNColor[1] <- "red"
options(biovizBase = opts)
## get from option(default)
getBioColor("DNA_BASES_N")
## get default fixed color
getBioColor("DNA_BASES_N", source = "default")
seqs <- c("A", "C", "T", "G", "G", "G", "C")
## get colors for a sequence.
getBioColor("DNA_BASES_N")[seqs]
}
\author{Tengfei Yin}

