newDataAfterFacetByGr <- function(gr, which, id.name){
  gr <- sort(gr)
  which <- sort(which)
  ## suppose which is a GRanges
  res <- lapply(seq_len(length(which)),function(i){
    res <- subsetByOverlaps(gr, which[i])
    values(res)$.bioviz.facetid <- i
    res
  })
  do.call(c, res)
}

strip_formula_dots <- function(formula){
  allvars <- all.vars(as.formula(formula))
  idx <- ggplot2_is_calculated_aes(allvars)
  if(sum(idx))
    formula[[which(idx)+1]] <- as.name(unlist(ggplot2_strip_dots(allvars[idx])))
  formula
}

setGeneric("splitByFacets", function(object, facets, ...) standardGeneric("splitByFacets"))

setMethod("splitByFacets", c("GRanges", "formula"), function(object, facets){
  .checkFacetsRestrict(facets, object)
  facets <- strip_formula_dots(facets)
  allvars <- all.vars(as.formula(facets))
  if(length(allvars) > 1){
  if(allvars[1] == "."){
    res <- split(object, seqnames(object))
  }else{
    if(allvars[1] != "strand")
      fts <- interaction(as.character(values(object)[,allvars[1]]),as.character(seqnames(object)))
    else
      fts <- interaction(as.character(strand(object)),as.character(seqnames(object)))
    res <- split(object, fts)    
  }}
  if(length(allvars) == 1)
    res <- split(object, seqnames(object))
  res <- res[elementNROWS(res)>0]
  res
})

setMethod("splitByFacets", c("GRanges", "GRanges"), function(object, facets){
  if(length(facets)){
      object <- newDataAfterFacetByGr(object, facets)
      grl <- split(object, values(object)$.bioviz.facetid)
    }else{
      stop("facets is an empty GenomicRanges object")
    }
  grl <- grl[elementNROWS(grl)>0]
  grl
})


setMethod("splitByFacets", c("GRanges", "missing"), function(object, facets){
  res <- split(object, seqnames(object))
  res <- res[elementNROWS(res)>0]
  res
})

setMethod("splitByFacets", c("GRanges", "NULL"), function(object, facets){
  res <- split(object, seqnames(object))
  res <- res[elementNROWS(res)>0]  
  res
})


isFacetByOnlySeq <- function(facets){
  allvars <- all.vars(as.formula(facets))
  allvars <- allvars[!allvars %in% c(".", "seqnames")]
  if(!length(allvars))
    TRUE
  else
    FALSE
}





## that's for linear or default!
.checkFacetsRestrict <- function(facets, object){
  facets <- strip_formula_dots(facets)
  allvars <- all.vars(as.formula(facets))
  if(length(allvars) == 1){
    if(allvars[1] != "seqnames")
      stop("Column of facets formula can only be seqnames, such as . ~ seqnames, in default restrict mode, you can only change row varaibles")
  }
  if(length(allvars) > 1){
    if(allvars[2] != "seqnames" & allvars[2] != "."){
      stop("Column of facets formula can only be seqnames or '.', such as . ~ seqnames, in default restrict mode, you can only change row varaibles")
    }
    if(allvars[1] != "."){
      if(!allvars[1] %in% c(colnames(values(object)), "strand"))
        stop(allvars[1]," doesn't exists in data columns")
    }
  }
}

