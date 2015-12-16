## Simple test function to test the crunch method on EnsDb
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

test_crunch_ensdb <- function(){
    zbtb <- genes(edb, filter=GenenameFilter("ZBTB16"))

    Zens <- crunch(edb, zbtb)
    ## crunch without which but using a filter.
    Test <- crunch(edb, GenenameFilter("ZBTB16"))

    ## Testing exceptions:
    checkException(crunch(edb, which="bla"))
    checkException(crunch(edb, which=list(GenenameFilter("ZBTB16"), "bla")))

    Test <- crunch(edb, zbtb)
    ## Expect a single gene.
    checkEquals(length(unique(Test$gene_id)), 1)
    strand(zbtb) <- "*"
    ## Expect two genes.
    Test <- crunch(edb, zbtb)
    checkEquals(length(unique(Test$gene_id)), 3)

    ## Testing columns
    wantCols <- c("gene_biotype", "tx_biotype")
    Test <- crunch(edb, zbtb, columns=wantCols)
    checkTrue(all(wantCols %in% colnames(mcols(Test))))

    ## Now let's try a tricky one...
    Test <- crunch(edb, which=GenenameFilter("ALDOA"))
    gr <- GRanges(seqnames=16, IRanges(30768000, 30770000), strand="*")
    Test <- crunch(edb, which=gr)
}

