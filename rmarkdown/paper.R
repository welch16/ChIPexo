## ----table1,include = TRUE,echo = FALSE,eval = TRUE----------------------
  load("../data/paper_tables/Landick_rif_summary.RData")
  load("../data/paper_tables/Landick_aero_summary.RData")
  load("../data/paper_tables/pugh_ctcf.RData")
  load("../data/paper_tables/Carroll_mm9_FoxA1.RData")
  load("../data/paper_tables/Meijsing_hg19_GR.RData")
  load("../data/paper_tables/Carroll_hg19_ER.RData")
  do.call(rbind,list(landick_aero,
                     landick_rif,
                     pugh_ctcf,
                     carroll_mm9_FoxA1,
                     meijing_hg19_GR,
                     carroll_hg19_ER))

## ----comp_param , include = FALSE,echo = FALSE,eval = TRUE---------------

  bin_size <- 150
  frag_len <- 150


## ----comp_param,include = FALSE,echo = FALSE, eval = TRUE----------------

  id_dist <- 20
  ext <- 15


## ----table_chipseq,include = TRUE,echo = FALSE,eval = TRUE---------------
  load("../data/paper_tables/Landick_rif_chipseq_summary.RData")
  do.call(rbind,chipseq)

