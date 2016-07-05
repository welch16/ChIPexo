(TeX-add-style-hook "chipexo_QC_vs2"
 (lambda ()
    (LaTeX-add-bibliographies
     "chip_exo_paper")
    (LaTeX-add-labels
     "sec:intro"
     "sec:results"
     "tab:qc"
     "tab:qc_sig"
     "sec:conc"
     "sec:methods"
     "mod"
     "fig:chip_diagram"
     "fig:comp"
     "fig:scc_exo"
     "fig:qcdiagram"
     "fig:enrich"
     "fig:strand"
     "fig:eval"
     "fig:methods_comp"
     "fig:design"
     "fig:reso_all")
    (TeX-add-symbols
     '("RW" 1)
     '("SK" 1)
     '("pname" 1)
     "beginsupplement"
     "sig")
    (TeX-run-style-hooks
     "inputenc"
     "utf8"
     "soul"
     "listings"
     "gensymb"
     "xr"
     "framed"
     "multirow"
     "tabu"
     "verbatim"
     "float"
     "tikz"
     "graphicx"
     "amssymb"
     "amsthm"
     "amsmath"
     "url"
     "xcolor"
     "latex2e"
     "bmcart10"
     "bmcart")))

