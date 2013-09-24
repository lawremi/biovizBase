library(biovizBase)
## ------------------------------------------------------------
## getBioColor
## ------------------------------------------------------------
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
getBioColor("DNA_BASES")
getBioColor("DNA_ALPHABET")
getBioColor("RNA_BASES_N")
getBioColor("RNA_BASES")
getBioColor("RNA_ALPHABET")
getBioColor("IUPAC_CODE_MAP")
getBioColor("AMINO_ACID_CODE")
getBioColor("AA_ALPHABET")
getBioColor("STRAND")
getBioColor("CYTOBAND")

## done
## ------------------------------------------------------------
## ColorBlindPal
## ------------------------------------------------------------
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
showColor(mypal(12, repeatable = FALSE)) # warning
show_col(mypal(11, repeatable = TRUE))  # no warning, and repeat
show_col(mypal12(12))
showColor(getBioColor("CYTOBAND"), "name")
showColor(getBioColor("DNA_BASES"), "name")
plotColorLegend(getBioColor("DNA_BASES"))

mypalFun <- colorBlindSafePal(21)   
par(mfrow = c(1, 3))
showColor(mypalFun(4))
library(dichromat)
showColor(dichromat(mypalFun(4), "deutan"))
showColor(dichromat(mypalFun(4), "protan"))

getOption("biovizBase")$cytobandColor

par(mfrow = c(1, 3))
cols <- getBioColor("DNA_BASES_N", "default")
showColor(cols, "name")
cols.deu <- dichromat(cols, "deutan")
names(cols.deu) <- names(cols)
cols.pro <- dichromat(cols, "protan")
names(cols.pro) <- names(cols)
showColor(cols.deu, "name")
showColor(cols.pro, "name")


library(GenomicRanges)
set.seed(1)
N <- 500
gr <- GRanges(seqnames =
              sample(c("chr1", "chr2", "chr3", "chrX", "chrY"),
                     size = N, replace = TRUE),
              IRanges(
                      start = sample(1:300, size = N, replace = TRUE),
                      width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N,
                replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              group = sample(c("Normal", "Tumor"),
                size = N, replace = TRUE),
              pair = sample(letters, size = N,
                replace = TRUE))



head(addSteppings(gr))
head(addSteppings(gr, group.name = "pair"))
gr.close <- GRanges(c("chr1", "chr1"), IRanges(c(10, 20), width = 9))
addSteppings(gr.close)
addSteppings(gr.close, extend.size = 5)


gr.temp <- GRanges("chr1", IRanges(start = c(100, 250),
                                 end = c(200, 300)))
maxGap(gaps(gr.temp, start = min(start(gr.temp))))
maxGap(gaps(gr.temp, start = min(start(gr.temp))), ratio = 0.5)



gr1 <- GRanges("chr1", IRanges(start = c(100, 300, 600),
                               end = c(200, 400, 800)))
shrink.fun1 <- shrinkageFun(gaps(gr1), max.gap = maxGap(gaps(gr1), 0.15))
shrink.fun2 <- shrinkageFun(gaps(gr1), max.gap = 0)
library(ggbio)
p1 <- qplot(gr1)
p2 <- qplot(shrink.fun1(gr1))
p3 <- qplot(shrink.fun2(gr1))
tracks(p1, p2, p3)

gr2 <- GRanges("chr1", IRanges(start = c(100, 350, 550),
                               end = c(220, 500, 900)))
gaps.gr <- intersect(gaps(gr1, start = min(start(gr1))),
                     gaps(gr2, start = min(start(gr2))))
shrink.fun <- shrinkageFun(gaps.gr, max.gap = maxGap(gaps.gr))
p1 <- qplot(gr1)
p2 <- qplot(gr2)
p3 <- qplot(shrink.fun(gr1))
p4 <- qplot(shrink.fun(gr2))
tracks(p1, p2, p3, p4)

library(Rsamtools)
data(genesymbol)
library(BSgenome.Hsapiens.UCSC.hg19)    
bamfile <- system.file("extdata", "SRR027894subRBM17.bam", package="biovizBase")
test <- pileupAsGRanges(bamfile, region = genesymbol["RBM17"])
test.match <- pileupGRangesAsVariantTable(test, Hsapiens)

## library(rtracklayer)
library(biovizBase)
hg19IdeogramCyto <- getIdeogram("hg19", cytoband = TRUE)
hg19Ideogram <- getIdeogram("hg19", cytoband = FALSE)
unknowIdeogram <- getIdeogram()

head(ucscGenomes()$db)

data(hg19IdeogramCyto)
head(hg19IdeogramCyto)
data(hg19Ideogram)
head(hg19Ideogram)

## to be removed
library(biovizBase)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb@organism
wh <- genesymbol["BRCA1"]
genome(txdb)
str(txdb)
obj <- txdb
## decrease from 18s to 13s, then to 4.7s
system.time(temp <- crunch(txdb, which = genesymbol["BRCA1"], type = "all"))
system.time(temp <- crunch(txdb, which = genesymbol["BRCA1"], type = "single"))

library("org.Hs.eg.db")
hdb <- org.Hs.eg.db
k <- keys(hdb, keytype = "SYMBOL")
res <- genGenesymbolTable(hdb, keys = c("BRCA1", "BRCA2"))
res <- genGenesymbolTable(hdb, keys = c("BRCA1", "BRCA2"), unique  = FALSE)

