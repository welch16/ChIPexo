## ----<style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------------------
  BiocStyle::latex()

## ----init,include=FALSE,echo=FALSE,eval=TRUE,results='hide'-----------------------------

  library(knitr)
  library(data.table)
  library(devtools)

  load_all("/p/keles/ChIPexo/volume4/ChIPexoQual")



## ----depth,include=FALSE,echo=FALSE,eval=TRUE,results='hide'----------------------------

  datadir = "/p/keles/ChIPexo/volume3/Analysis"
  pughsamples = "CTCF"
  pughdir = file.path(datadir,"Pugh",pughsamples,"data")
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

  pughdepths = mapply(load_depth,pughdir,pughsamples)
  names(pughdepths) = pughsamples

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

  pughpbc = mapply(load_PBC,pughdir,pughsamples)
  names(pughpbc) =pughsamples

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

  pughcc = mapply(load_cc,pughdir,pughsamples,SIMPLIFY=FALSE)
  names(pughcc) = pughsamples

  rencc= mapply(load_cc,rendir,rensamples,SIMPLIFY=FALSE)
  names(rencc) = rensamples

  carrollcc = mapply(load_cc,carrolldir,carrollsamples,SIMPLIFY=FALSE)
  names(carrollcc) = carrollsamples

  landickcc = mapply(load_cc,landickdir,landicksamples,SIMPLIFY=FALSE)
  names(landickcc) = landicksamples

  nsc <- function(cc)cc[between(shift,1,300),max(cross.corr)/min(cross.corr)]

  pughnsc = sapply(pughcc,nsc)
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
  pugh = create_table(pughsamples,pughdepths,pughpbc,pughnsc)

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

## ----include = TRUE,echo=FALSE,eval =TRUE,results ='asis'-------------------------------
  kkable(pugh)

## ----include=TRUE,echo=FALSE,eval=TRUE,results='asis'-----------------------------------
  kkable(carroll[grep("Fox",dataset)])

## ----include=TRUE,echo=FALSE,eval=TRUE,results ='asis'----------------------------------
  kkable(carroll[grep("Fox",dataset,invert=TRUE)])

## ----include=TRUE,echo=FALSE,eval=TRUE,results='asis'-----------------------------------
  kkable(ren)

## ----include=FALSE,echo=FALSE,eval=TRUE,results='hide'----------------------------------

stack_cc = function(cc_list)
{
  nms = names(cc_list)
  if(is.null(nms)){
    nms = 1:length(cc_list)
  }
  cc_list= mapply(function(cc,nm){
    cc[,name:=nm]
    return(cc)
  },cc_list,nms,SIMPLIFY = FALSE)
  return(do.call(rbind,cc_list))  
}

allrencc = stack_cc(rencc)
allhumantf = stack_cc(c(carrollcc[1:3],pughcc))
allmouse = stack_cc(carrollcc[4:6])
landick1 = stack_cc(landickcc[14:17])
landick2 = stack_cc(landickcc[5:8])
landick3 = stack_cc(landickcc[1:4])


## ----echo=FALSE,out.width='5cm',fig.show='hold',warning=FALSE---------------------------
ggplot(allrencc,aes(shift,cross.corr,colour = name))+geom_line()+
  geom_abline(slope=0,intercept=0,linetype=2)+
  scale_y_continuous(limits = c(0,.01))+geom_vline(xintercept=0,linetype=2)+
  theme(legend.position = "top") +
  scale_x_continuous(limits = c(0,300))

ggplot(allhumantf,aes(shift,cross.corr,colour = name))+geom_line()+
  geom_abline(slope=0,intercept=0,linetype=2)+
  scale_y_continuous(limits = c(0,.5))+geom_vline(xintercept=0,linetype=2)+
  theme(legend.position = "top")+
  scale_x_continuous(limits = c(0,300))

ggplot(allmouse,aes(shift,cross.corr,colour = name))+geom_line()+
  geom_abline(slope=0,intercept=0,linetype=2)+
  scale_y_continuous(limits = c(0,.25))+geom_vline(xintercept=0,linetype=2)+
  theme(legend.position = "top") +
  scale_x_continuous(limits = c(0,300))
     

## ----echo=FALSE,out.width='5cm',fig.show='hold',warning=FALSE---------------------------
  ggplot(landick1,aes(shift,cross.corr,colour = name))+geom_line()+
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_vline(xintercept=0,linetype=2)+
  theme(legend.position = "top") +
  scale_x_continuous(limits = c(0,300))+scale_y_continuous(limits = c(0,.35))


ggplot(landick2,aes(shift,cross.corr,colour = name))+geom_line()+
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_vline(xintercept=0,linetype=2)+
  theme(legend.position = "top") +
  scale_x_continuous(limits = c(0,300))+scale_y_continuous(limits = c(0,.2))


ggplot(landick3,aes(shift,cross.corr,colour = name))+geom_line()+
  geom_abline(slope=0,intercept=0,linetype=2)+
  geom_vline(xintercept=0,linetype=2)+
  theme(legend.position = "top") +
  scale_x_continuous(limits = c(0,300))+scale_y_continuous(limits = c(0,.15))               

