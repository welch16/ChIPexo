## ----table1,include = FALSE,echo = FALSE,eval = TRUE---------------------
  library(data.table)
  load("../data/paper_tables/Landick_rif_summary.RData")
  load("../data/paper_tables/Landick_aero_summary.RData")
  load("../data/paper_tables/pugh_ctcf.RData")
  load("../data/paper_tables/Carroll_mm9_FoxA1.RData")
  load("../data/paper_tables/Meijsing_hg19_GR.RData")
  load("../data/paper_tables/Carroll_hg19_ER.RData")
  load("../data/paper_tables/Venters_hg19_TBP.RData")
  load("../data/paper_tables/ChIPnexus_dm3.RData")
  load("../data/paper_tables/ChIPnexus_hg19.RData")
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
  venters_hg19_TBP[,Organism := "Human"]
  venters_hg19_TBP[,cond := "K562"]
  venters_hg19_TBP[,Replicate := 1:3]
  venters_hg19_TBP[,TF := "TBP"]
  chipnexus_dm3[, Organism := "D. Melanogaster"]
  chipnexus_dm3[, cond := c(rep("Embryo",4),rep("S2",4))]
  chipnexus_dm3[, Replicate := rep(1:2,4)]
  chipnexus_dm3[, TF := c(rep("Dorsal",2),rep("Twist",2),rep("Max",2),rep("MyC",2))]
  chipnexus_K562[ , Organism := "Human"]
  chipnexus_K562[ , cond  := "K562"]
  chipnexus_K562[ , Replicate := 1:2]
  chipnexus_K562[ , TF := "TBP"]
  ## this creates the table
  A <- do.call(rbind,list(landick_aero,
                     landick_rif[order(Replicate)],
                     pugh_ctcf[order(Replicate)],
                     carroll_mm9_FoxA1[order(Replicate)],
                     meijing_hg19_GR[order(Replicate)],
                     carroll_hg19_ER[order(Replicate)],
                     venters_hg19_TBP[order(Replicate)]))
  B <- do.call(rbind,list(chipnexus_dm3,
                     chipnexus_K562))
  format_tab <- function(A){
    A[,nsc := unlist(nsc)]
    setcolorder(A, c("Organism","TF","cond","Replicate","nreads","pbc","nsc","files"))
    A[,files := NULL]
    A[,nreads := prettyNum(nreads,big.mark = ",")]
    A[,pbc := round(pbc,4)]
    A[,nsc := round(nsc,4)]
    setnames(A,names(A),c("Organism","IP/TF","Condition/Cell","Rep.","Depth","PBC","NSC"))
    return(A)}

  A <- format_tab(A)
  B <- format_tab(B)

## ----table2,include = TRUE,echo = FALSE, eval = TRUE---------------------
  A

## ----comp_param , include = FALSE,echo = FALSE,eval = TRUE---------------

  bin_size <- 150
  frag_len <- 150


## ----table3,include = TRUE,echo = FALSE, eval = TRUE---------------------
  B

## ----all_param,include = FALSE,echo = FALSE,eval = TRUE------------------

  Ntimes <- prettyNum(1e4,big.mark = ",")
  Msamp <- prettyNum(1e3,big.mark = ",")
  thresh1 <- 10


## ----comp_param,include = FALSE,echo = FALSE, eval = TRUE----------------

  id_dist <- 20
  ext <- 15


## ----methods_param , include = FALSE,echo = FALSE,eval = TRUE------------

  fdr <- 0.01
  topM <- 300


## ----table_chipseq1,include = TRUE,echo = FALSE,eval = TRUE--------------
  load("../data/paper_tables/Landick_rif_chipseq_summary.RData")
  chipseq$set <- chipseq$set[,Organism := "E.Coli"]
  chipseq$set <- chipseq$set[,Replicate := c(1,1,2,2)]
  chipseq$set <- chipseq$set[,TF := "Sig70"]
  chipseq$set <- chipseq$set[ ,cond := rep(c("Rif-0min","Rif-20min"),each = 2)]
  chipseq$pet <- chipseq$pet[,Organism := "E.Coli"]
  chipseq$pet <- chipseq$pet[,Replicate := c(1,1,2,2)]
  chipseq$pet <- chipseq$pet[,TF := "Sig70"]
  chipseq$pet <- chipseq$pet[ ,cond := rep(c("Rif-0min","Rif-20min"),each = 2)]
  chipseq <- lapply(chipseq,format_tab)
  chipseq <- do.call(rbind,chipseq)
  chipseq <- chipseq[, Protocol := rep(c("PE","SE"),each = 4)]
  chipseq

