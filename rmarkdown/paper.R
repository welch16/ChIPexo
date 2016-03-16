## ----table1,include = FALSE,echo = FALSE,eval = TRUE---------------------
  library(data.table)
  load("../data/paper_tables/Landick_rif_summary.RData")
  load("../data/paper_tables/Landick_aero_summary.RData")
  load("../data/paper_tables/pugh_ctcf.RData")
  load("../data/paper_tables/Carroll_mm9_FoxA1.RData")
  load("../data/paper_tables/Meijsing_hg19_GR.RData")
  load("../data/paper_tables/Carroll_hg19_ER.RData")
  landick_aero[,Organism := "E.Coli"]
  landick_aero[,Replicate := c(1,2,1,2)]
  landick_aero[,TF := "Sig70"]
  landick_aero[ ,cond := c(rep("Aerobic",2),rep("Anaerobic", 2))]
  landick_rif[,Organism := "E.Coli"]
  landick_rif[,Replicate := c(1,1,2,2)]
  landick_rif[,TF := "Sig70"]
  landick_rif[ ,cond := rep(c("Rif-0min","Rif-20min"),each = 2)]
  pugh_ctcf[,Organism := "Human"]
  pugh_ctcf[,Replicate := 1]
  pugh_ctcf[,TF := "CTCF"]
  pugh_ctcf[ ,cond := "HeLA"]
  carroll_mm9_FoxA1[,Organism := "Mouse"]
  carroll_mm9_FoxA1[,Replicate := c(3,1,2)]
  carroll_mm9_FoxA1[,TF := "Fox A1"]
  carroll_mm9_FoxA1[ ,cond := "Mouse Liver"]
  meijing_hg19_GR[,Organism := "Human"]
  meijing_hg19_GR[,Replicate := 1]
  meijing_hg19_GR[, TF := "GR"]
  meijing_hg19_GR[,cond := c("IMR90","K562","U20S")]
  carroll_hg19_ER[,Organism := "Human"]  
  carroll_hg19_ER[,Replicate := c(1,3,2)]
  carroll_hg19_ER[,TF := "ER"]
  carroll_hg19_ER[ ,cond := "MCF-7"]
  A <- do.call(rbind,list(landick_aero,
                     landick_rif[order(Replicate)],
                     pugh_ctcf[order(Replicate)],
                     carroll_mm9_FoxA1[order(Replicate)],
                     meijing_hg19_GR[order(Replicate)],
                     carroll_hg19_ER[order(Replicate)]))
  A[,nsc := unlist(nsc)]
  setcolorder(A, c("Organism","TF","cond","Replicate","nreads","pbc","nsc","files"))
  A[,files := NULL]
  A[,nreads := prettyNum(nreads,big.mark = ",")]
  A[,pbc := round(pbc,4)]
  A[,nsc := round(nsc,4)]
  setnames(A,names(A),c("Organism","IP/TF","Condition/Cell","Rep.","Depth","PBC","NSC"))

## ----table2,include = TRUE,echo = FALSE, eval = TRUE---------------------
  A

## ----comp_param , include = FALSE,echo = FALSE,eval = TRUE---------------

  bin_size <- 150
  frag_len <- 150


## ----comp_param,include = FALSE,echo = FALSE, eval = TRUE----------------

  id_dist <- 20
  ext <- 15


## ----methods_param , include = FALSE,echo = FALSE,eval = TRUE------------

  fdr <- 0.01
  topM <- 300


## ----table_chipseq1,include = TRUE,echo = FALSE,eval = TRUE--------------
  load("../data/paper_tables/Landick_rif_chipseq_summary.RData")
  do.call(rbind,chipseq)

