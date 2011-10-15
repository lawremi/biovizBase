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
