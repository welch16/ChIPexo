## ----<style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------------------
  BiocStyle::latex()

## ----init,include=FALSE,echo=FALSE,eval=TRUE,results='hide'-----------------------------

  library(knitr)
  library(data.table)
  library(devtools)

  load_all("/p/keles/ChIPexo/volume4/ChIPexoQual")



## ----depth,include=FALSE,echo=FALSE,eval=TRUE,results='hide'----------------------------

  datadir = "/p/keles/ChIPexo/volume3/Analysis"
  rensamples = c("H3k27ac","H3k4me1-rep1","H3k4me1-rep2")
  rendir = file.path(datadir,"Ren",rensamples,"data")
  carrollsamples = c(paste0("ER-rep",1:3),paste0("FoxA1-rep",1:3))
  carrolldir = file.path(datadir,"Carroll",rep(c("human","mouse"),each=3),
    carrollsamples,"data")
  rifsamples = paste0(rep(c("BetaPrimeFlag","Beta","Sig70"),each=4),"-rif",
    rep(c(0,20),each=2),"-rep",1:2)
  statsamples = c("Input",paste0(rep(c("Sig70","SigmaS"),each=4),rep(c("-exp-","-stat-"),each=2),
    "rep",rep(1:2,4)))
  landicksamples = c(rifsamples,statsamples)
  landickdir = file.path(datadir,"Landick",c(rep("rif-treatment",12),rep("stat-vs-exp",9)),
    landicksamples,"data")

  load_depth <- function(direc,samp){
    load(file.path(direc,paste0(samp,"_depth.RData")))
    return(depth)}

  rendepths = mapply(load_depth,rendir,rensamples)
  names(rendepths) = rensamples

  carrolldepths = mapply(load_depth,carrolldir,carrollsamples)
  names(carrolldepths) = carrollsamples

  landickdepths = mapply(load_depth,landickdir,landicksamples)
  names(landickdepths) = landicksamples


## ----pbc,include=FALSE,echo=FALSE,eval=TRUE,results='hide'------------------------------

  load_PBC <- function(direc,samp){    
    load(file.path(direc,paste0(samp,"_PBC.RData")))
  return(pcr_coeff)}

  renpbc= mapply(load_PBC,rendir,rensamples)
  names(renpbc) = rensamples

  carrollpbc = mapply(load_PBC,carrolldir,carrollsamples)
  names(carrollpbc) = carrollsamples

  landickpbc = mapply(load_PBC,landickdir,landicksamples)
  names(landickpbc) = landicksamples



## ----cross_corr,include=FALSE,echo=FALSE,eval=TRUE,results='hide'-----------------------
  
  load_cc <- function(direc,samp){
    load(file.path(direc,paste0(samp,"_cross.corr.RData")))
  return(cc)}

  rencc= mapply(load_cc,rendir,rensamples,SIMPLIFY=FALSE)
  names(rencc) = rensamples

  carrollcc = mapply(load_cc,carrolldir,carrollsamples,SIMPLIFY=FALSE)
  names(carrollcc) = carrollsamples

  landickcc = mapply(load_cc,landickdir,landicksamples,SIMPLIFY=FALSE)
  names(landickcc) = landicksamples

  nsc <- function(cc)cc[,max(cross.corr)/min(cross.corr)]

  rennsc = sapply(rencc,nsc)
  carrollnsc = sapply(carrollcc,nsc)
  landicknsc = sapply(landickcc,nsc)

  create_table <- function(samples,depth,pbc,nsc)
  {
    return(data.table(dataset =samples,depth = as.character(depth) , pbc = pbc,
                      nsc = nsc))
  }
 
  landicksamples = gsub("PrimeFlag","Pf",landicksamples)
  landick = create_table(landicksamples,landickdepths,landickpbc,landicknsc)
  ren = create_table(rensamples,rendepths,renpbc,rennsc)
  carroll = create_table(carrollsamples,carrolldepths,carrollpbc,carrollnsc)

  kkable <- function(dt)
  {
    out = dt[,(c("depth","pbc","nsc")):=
      list(prettyNum(depth,big.mark=","),
           round(pbc,3) , round(nsc,3)),
      by = dataset]
    return(kable(as.data.frame(out)))
  }


## ----include=TRUE,echo=FALSE,eval=TRUE,results='asis'-----------------------------------
  kkable(landick[grep("rif",dataset,invert=TRUE)][-1])

## ----include=TRUE,echo=FALSE,eval=TRUE,results='asis'-----------------------------------
  kkable(landick[grep("rif",dataset)])

## ----include=TRUE,echo=FALSE,eval=TRUE,results='asis'-----------------------------------
  kkable(carroll[grep("Fox",dataset)])

## ----include=TRUE,echo=FALSE,eval=TRUE,results ='asis'----------------------------------
  kkable(carroll[grep("Fox",dataset,invert=TRUE)])

## ----include=TRUE,echo=FALSE,eval=TRUE,results='asis'-----------------------------------
  kkable(ren)

## ----echo=FALSE,out.width='5cm',fig.show='hold'-----------------------------------------
ggplot(rencc[[1]],aes(shift,cross.corr))+geom_line()+geom_abline(slope=0,intercept=0,linetype=2)+scale_y_continuous(limits = c(0,1))+geom_vline(xintercept=0,linetype=2)
ggplot(carrollcc[[1]],aes(shift,cross.corr))+geom_line()+geom_abline(slope=0,intercept=0,linetype=2)+scale_y_continuous(limits = c(0,1))+geom_vline(xintercept=0,linetype=2)
ggplot(carrollcc[[6]],aes(shift,cross.corr))+geom_line()+geom_abline(slope=0,intercept=0,linetype=2)+scale_y_continuous(limits = c(0,1))+geom_vline(xintercept=0,linetype=2)

