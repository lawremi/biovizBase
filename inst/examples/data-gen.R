library(biovizBase)
library(GenomicRanges)
## with no 
hg19IdeogramCyto <- getIdeogram("hg19", cytobands = TRUE)
seqinfo(hg19Ideogram)
hg19Ideogram <- getIdeogram("hg19", cytobands = FALSE)
## then for genesymbol
library("org.Hs.eg.db")
ls("package:org.Hs.eg.db")
available.dbschemas()

symbol_id <- mappedRkeys(org.Hs.egSYMBOL)
head(symbol_id)
length(symbol_id)


gene_id <- mget(symbol_id, revmap(org.Hs.egSYMBOL), ifnotfound=NA)
gene_id <- unlist2(gene_id)
gene_id <- gene_id[!is.na(gene_id)]
chr <- toTable(org.Hs.egCHR[gene_id])
chrloc <- toTable(org.Hs.egCHRLOC[gene_id])
chrlocend <- toTable(org.Hs.egCHRLOCEND[gene_id])
identical(chrloc$gene_id, chrlocend$gene_id)
chrs <- chrloc
head(chrlocend)
chrs$end_location <- chrlocend$end_location
head(chrs)
res <- chrs
sum(duplicated(gene_id))# no duplicate
idx <- match(res$gene_id, gene_id)
res$symbol <- names(gene_id[idx])
sum(with(res, start_location > 0 & end_location < 0))# good
sum(with(res, start_location < 0 & end_location > 0))# good
## merge(merge(data.frame(gene_id=gene_id, symbol=names(gene_id)), chrloc),
## chrlocend)
## df <- merge(merge(merge(data.frame(gene_id = gene_id,
##                                    symbol = names(gene_id)), chr), chrloc), chrlocend)
## do I need to subset it? No.
## add extra column to indicate sense or antisense.
##  FOXD4L2,
res$sense <- with(res, start_location>=0&end_location>=0)
head(res)
res.sense <- res[res$sense,]
sum(with(res.sense, start_location > end_location))
res.antisense <- res[!res$sense,]
sum(with(res.antisense, start_location < end_location))
head(res)
## need to make sure you get this right
res$strand <- strand(!res$sense)
head(res)
sum(with(res, start_location == end_location))
## df2 <- subset(df, start_location>=0&end_location>=0)
## strand <- ifelse(df2$start_location>df2$end_location, "-", "+")
## strand[df2$start_location == df2$end_location] <- "*"
## id <- strand == "-"
## tmp <- df2[id,"start_location"]
## df2[id,"start_location"] <- df2[id,"end_location"]
## df2[id,"end_location"] <- tmp
## df2$strand <- strand
## df2 <- df2[!duplicated(df2),]
## for ensembl
ensembl_id <- unlist(mget(gene_id, org.Hs.egENSEMBL))
head(ensembl_id)
idx <- match(res$gene_id, names(ensembl_id))
res$ensembl_id <- ensembl_id[idx]
head(res)
library(GenomicRanges)
genesymbol <- with(res,GRanges(seqnames = paste("chr",Chromosome,sep = ""),
                              IRanges(abs(start_location), abs(end_location)),
                               strand = strand,
                              symbol = symbol,
                               ensembl_id = ensembl_id))
names(genesymbol) <- values(genesymbol)$symbol
save(genesymbol, file = "~/Datas/rdas/genesymbol.rda")
save(genesymbol, file = "~/Codes/gitrepos/biovizBase/data/genesymbol.rda")
save(hg19IdeogramCyto, file = "~/Codes/gitrepos/biovizBase/data/hg19IdeogramCyto.rda")
save(hg19Ideogram, file = "~/Codes/gitrepos/biovizBase/data/hg19Ideogram.rda")

id <- unlist2(mget("FOXD4L2", revmap(org.Hs.egSYMBOL)))
id <- revmap(org.Hs.egSYMBOL)[["FOXD4L2"]]
org.Hs.egCHRLOC[[id]]
org.Hs.egCHRLOCEND[[id]]

## new data generation method, with seqlengths
ideogramCyto <- getIdeogram(c("hg19", "hg18", "mm10", "mm9"))
ideogram <- getIdeogram(c("hg19", "hg18", "mm10", "mm9"), cytoband = FALSE)
ideo <- ideogram
ideoCyto <- ideogramCyto
save(ideo,  file = "../../../../bioc-devel/biovizBase/data/ideo.rda")
save(ideoCyto, file = "/Users/tengfei/Code/svnrepos/bioc-devel/biovizBase/data/ideoCyto.rda")
hg19 <- ideoCyto$hg19
ideoCyto$hg18 <- keepSeqlevels(ideoCyto$hg18, paste0("chr", c(1:22, "X", "Y")))
ideoCyto$hg19 <- keepSeqlevels(ideoCyto$hg19, paste0("chr", c(1:22, "X", "Y")))
ideoCyto$mm10 <- keepSeqlevels(ideoCyto$mm10, paste0("chr", c(1:19, "X", "Y")))
ideoCyto$mm9 <- keepSeqlevels(ideoCyto$mm9, paste0("chr", c(1:19, "X", "Y")))

## for circular view
crc1 <- system.file("extdata", "crc1-missense.csv", package = "biovizBase")
crc1 <- read.csv(crc1)
library(GenomicRanges)
mut.gr <- with(crc1,GRanges(Chromosome, IRanges(Start_position, End_position),
                            strand = Strand))
values(mut.gr) <- subset(crc1, select = -c(Start_position, End_position, Chromosome))
data("hg19Ideogram", package = "biovizBase")
seqs <- seqlengths(hg19Ideogram)
## subset_chr
chr.sub <- paste("chr", 1:22, sep = "")
## levels tweak
seqlevels(mut.gr) <- c(chr.sub, "chrX")
mut.gr <- keepSeqlevels(mut.gr, chr.sub)
seqs.sub <- seqs[chr.sub]
## remove wrong position
bidx <- end(mut.gr) <= seqs.sub[match(as.character(seqnames(mut.gr)),
                                      names(seqs.sub))]
mut.gr <- mut.gr[which(bidx)]
## assign_seqlengths
seqlengths(mut.gr) <- seqs.sub
## reanme to shorter names
new.names <- as.character(1:22)
names(new.names) <- paste("chr", new.names, sep = "")
new.names
mut.gr.new <- renameSeqlevels(mut.gr, new.names)
head(mut.gr.new)

hg19Ideo <- hg19Ideogram
hg19Ideo <- keepSeqlevels(hg19Ideogram, chr.sub)
hg19Ideo <- rena

rearr  <- read.csv(system.file("extdata", "crc-rearrangment.csv", package = "biovizBase"))
## start position
gr1 <- with(rearr, GRanges(chr1, IRanges(pos1, width = 1)))
## end position
gr2 <- with(rearr, GRanges(chr2, IRanges(pos2, width = 1)))
## add extra column
nms <- colnames(rearr)
.extra.nms <- setdiff(nms, c("chr1", "chr2", "pos1", "pos2"))
values(gr1) <- rearr[,.extra.nms]
## remove out-of-limits data
seqs <- as.character(seqnames(gr1))
.mx <- seqlengths(hg19Ideo)[seqs]
idx1 <- start(gr1) > .mx
seqs <- as.character(seqnames(gr2))
.mx <- seqlengths(hg19Ideo)[seqs]
idx2 <- start(gr2) > .mx
idx <- !idx1 & !idx2
gr1 <- gr1[idx]
seqlengths(gr1) <- seqlengths(hg19Ideo)
gr2 <- gr2[idx]
seqlengths(gr2) <- seqlengths(hg19Ideo)

values(gr1)$to.gr <- gr2
## rename to gr
gr <- gr1
values(gr)$rearrangements <- ifelse(as.character(seqnames(gr))
                                    == as.character(seqnames((values(gr)$to.gr))),
                                    "intrachromosomal", "interchromosomal")
crc.gr <- gr
crc.gr
mut.gr <- mut.gr.new
hg19sub <- hg19Ideo
save(crc.gr, mut.gr, hg19sub, file = "/Users/tengfei/Code/svnrepos/bioc-devel/biovizBase/data/CRC.rda")


snp <- read.table("/Users/tengfei/Downloads/plink.assoc.txt", header = TRUE)
N <- 2000
set.seed(123)
snp <- snp[sample(1:nrow(snp), size = N, replace = FALSE),]
rownames(snp) <- NULL
head(snp)
write.table(snp, 
            file = "/Users/tengfei/Code/svnrepos/bioc-devel/biovizBase/inst/extdata/plink.assoc.sub.txt", row.names = FALSE)
## https://github.com/stephenturner/qqman
