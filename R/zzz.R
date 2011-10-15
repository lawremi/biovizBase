.onLoad <- function(libname, pkgname){
  options(biovizBase = list(
            DNABasesNColor = getBioColor("DNA_BASES_N", source = "default"),
            DNABasesColor = getBioColor("DNA_BASES", source = "default"),
            DNAAlphabetColor = getBioColor("DNA_ALPHABET", source = "default"),
            RNABasesNColor = getBioColor("RNA_BASES_N", source = "default"),
            RNABasesColor = getBioColor("RNA_BASES", source = "default"),
            RNAAlphabetColor = getBioColor("RNA_ALPHABET", source = "default"),
            IUPACCodeMapColor = getBioColor("IUPAC_CODE_MAP", source = "default"),
            AminoAcidCodeColor = getBioColor("AMINO_ACID_CODE", source = "default"),                  AAAlphabetColor = getBioColor("AA_ALPHABET", source = "default"),
            strandColor = getBioColor("STRAND", source = "default"),
            cytobandColor = getBioColor("CYTOBAND", source = "default")
            ))
}
