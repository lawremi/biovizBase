(TeX-add-style-hook "intro"
 (lambda ()
    (LaTeX-add-labels
     "fig:dichromat"
     "fig:cytoband"
     "fig:cytoband-show-col"
     "fig:strand"
     "fig:ne"
     "fig:base-legend"
     "fig:shrink-single"
     "fig:shrink-two")
    (TeX-add-symbols
     '("software" 1)
     '("Rcode" 1)
     '("Rclass" 1)
     '("Rfunarg" 1)
     '("Rmethod" 1)
     '("Rpackage" 1)
     '("Robject" 1)
     '("Rfunction" 1)
     "R"
     "IRanges"
     "biovizBase"
     "ggbio"
     "visnab")
    (TeX-run-style-hooks
     "verbatim"
     "hyperref"
     "latex2e"
     "art10"
     "article"
     "10pt")))

