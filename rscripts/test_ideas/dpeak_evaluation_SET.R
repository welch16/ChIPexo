rm(list = ls())

library(ggplot2)
library(GenomicAlignments)
library(data.table)
library(grid)
library(gridExtra)

dr <- "/p/keles/ChIPexo/volume7/Landick/K12"
files <- list.files(dr,recursive = TRUE)
files <- files[grep("rif_treatment",files)]
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]
## files <- files[grep("PET",files,invert = TRUE)] ## only keep SET
files <- files[grep("Input",files,invert = TRUE)] ## only chip
files <- files[c(4,8,12)]


annot_dir <- "/p/keles/ChIPexo/volume6/K12/annotations"
afiles <- list.files(annot_dir)
afiles <- afiles[grep("bed",afiles)]

annots <- lapply(file.path(annot_dir,afiles),read.table)
annots[[1]] <- IRanges(start = annots[[1]][,2],width = 1)
annots[[2]] <- IRanges(end = annots[[2]][,3],width = 1)
annots <- reduce(c(annots[[1]],annots[[2]]))
width(annots) <- 1
annots <- sort(annots)

size <- 250
anchors <- annots[which(diff(mid(annots)) > size)]
center <- mid(anchors)
start(anchors) <- center -2* size 
end(anchors) <- center + 2*size 
start(anchors[1]) <- 1
anchors <- anchors[countOverlaps(anchors,annots) > 1]
anchors <- reduce(anchors)

peaks <- ChIPUtils::gr2dt(GRanges(seqnames = "U00096",anchors))

temp_peakfile <- tempfile(pattern = "peak",fileext = ".txt")

write.table(peaks,file = temp_peakfile,row.names = FALSE,col.names = FALSE,quote = FALSE)

## reads chunk
load_reads <- function(file,PET)
{
  if(PET){
    out <- readGAlignmentPairs(file,param = NULL)
    out1 <- left(out)
    out2 <- right(out)
    sqnms <- seqnames(out1)
    st <- start(out1)
    en <- end(out2)
    out <- GRanges(seqnames =sqnms,
            ranges =  IRanges(start = st , end = en),strand = strand(out1))
  }else{
    out <- readGAlignments(file,param = NULL)
    out <- as(out,"GRanges")
  }
  return(out)
}

reads <- mcmapply(load_reads,file.path(dr,files),c(FALSE,TRUE,FALSE),mc.cores = 3,SIMPLIFY = FALSE)
names(reads) <- sapply(strsplit(basename(files),"_"),function(x)x[1])

reads <- mapply(function(x,PET){
  if(!PET){
    x <- resize(x,1)
  }
  return(x)},reads,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)

##

library(dpeak)

get_sites <- function(file,PET,peakfile){
  dpeak <- dpeakRead(peakfile = peakfile, readfile = file,fileFormat = "bam",parallel = TRUE,
                     nCore = 24,PET = PET)
  fit <- dpeakFit(dpeak,maxComp = 5,nCore = 24)
  return(fit)
}

## fits <- mapply(get_sites,file.path(dr,files),c(FALSE,TRUE,FALSE),
##                MoreArgs = list(temp_peakfile),SIMPLIFY = FALSE)
## save(fits, file = "data/dpeak_experiment_fits.RData")

load("data/dpeak_experiment_fits.RData")

sapply(reads,length)


fit_to_sites <- function(fit){
  temp_bed <- tempfile(pattern = "binding",fileext = ".bed")
  export(fit,type = "bed",filename = temp_bed)
  out <- data.table(read.table(temp_bed,skip = 1))
  return(out)
}

sites <- lapply(fits,fit_to_sites)

predictions <- lapply(sites,function(x)
  IRanges(start = mid(x[,IRanges(start = V2,end = V3)]) + 1,width = 1))
names(predictions) <- names(reads)


## Resolution

nreads <- lapply(reads,function(x)countOverlaps(anchors,ranges(x)))
nreads <- data.table(do.call(cbind,nreads))
if(ncol(peaks) == 4)peaks <- cbind(peaks,nreads)

peaks[, nevents := countOverlaps(anchors,annots)]
peaks[,table(nevents)]

resolution <- function(pred,anchors,annot)
{
  anchors <- split(anchors,1:length(anchors))
  out <- lapply(anchors,function(x,pred,annot){
    annot <- subsetByOverlaps(annot,x)
    pred <- subsetByOverlaps(pred,x)
    out <- sapply(mid(annot),function(y)min(abs(y - mid(pred))))
    if(length(pred) < length(annot)){
      out <- sort(out)[1:length(pred)]
    }
    return(out)},pred,annot)
}

res <- mclapply(predictions,resolution,anchors,annots,mc.cores = 8 )

res_boxplot <- function(res,rif,repl,minreads,maxDist,peaks)
{
  idx <- which(peaks[[1]] > minreads)
  res1 <- data.table(seq = "ChIPexo",reso =  do.call(c,res[[1]][idx]))
  res2 <- data.table(seq = "ChIPseq_SET",reso  = do.call(c,res[[3]][idx]))
  res3 <- data.table(seq = "ChIPseq_PET",reso  = do.call(c,res[[2]][idx]))
  res1 <- res1[reso <= maxDist]
  dt <- rbind(res1,res2,res3)
  dt[,rif := rif]
  dt[,rep := repl]
  return(dt)
}

resoDT <- res_boxplot(res,"20min","rep2",500,80,peaks)

ggplot(resoDT,aes(seq , reso,colour = seq))+
  geom_boxplot()+coord_cartesian(ylim = c(0,150))+facet_grid(rif ~ rep)+theme_bw()+
  scale_color_brewer(palette = "Set1")+theme(legend.position = "none")
dev.off()


## resoDT <- list()
## resoDT[[1]] <- res_boxplot("edsn1311","edsn1396","0min","rep-1",200,80,peaks,res)
## resoDT[[2]] <- res_boxplot("edsn1314","edsn1398","20min","rep-1",500,80,peaks,res)
## resoDT[[3]] <- res_boxplot("edsn1317","edsn1400","0min","rep-2",500,80,peaks,res)
## resoDT[[4]] <- res_boxplot("edsn1320","edsn1402","20min","rep-2",500,80,peaks,res)
## resoDT <- do.call(rbind,resoDT)


## pdf(file = "figs/resolution_experiment_fixed_regions.pdf")
## ggplot(resoDT,aes(seq , reso))+geom_boxplot()+coord_cartesian(ylim = c(0,1e2))+facet_grid(rif ~ rep)+theme_bw()
## dev.off()

## ## sensitivity
sensitivity <- function(pred,anchors,annot,ext)
{
  anchors <- split(anchors,1:length(anchors))
  out <- lapply(anchors,function(x,pred,annot,ext){
    ## pred <- subsetByOverlaps(pred,x)
    annot <- subsetByOverlaps(annot,x)
    output <- list()
    output[["npred"]] <- length(subsetByOverlaps(pred,x))
    output[["nannot"]] <- length(annot)
    output[["eventDist"]] <- mean(diff(sort(mid(annot))))
    output[["nidenDC"]] <- sum(sapply(mid(annot),function(x)
    as.numeric(length(which(abs(x - mid(pred)) <= ext)) >= 1)))
    output[["niden"]] <- sum( sapply(mid(annot),function(x)min(abs(x - mid(pred)))) <= ext)
    output[["sens"]] <- output[["niden"]] / output[["nannot"]]
    return(output)
  },pred,annot,ext)
  return(out)
}


library(MASS)

sens_plot <- function(set1,set2,rif,repl,minreads,peaks,sens)
{
  idx <- which(peaks[[set1]] > minreads)
  sens1 <- sens[[set1]][idx]
  sens2 <- sens[[set2]][idx]
  nms <- names(sens1[[1]])
  sens1 <- data.table(do.call(rbind,lapply(sens1,unlist)))
  sens2 <- data.table(do.call(rbind,lapply(sens2,unlist)))
  sens <- rbind(sens1[,seq := "ChIPexo"],sens2[,seq:="ChIPseq_SET"])
  sens[,rif := rif]
  sens[,repl := repl]

  ggplot(sens[!(sens == 0 & eventDist >200)],aes(eventDist,sens,colour = seq))+
    geom_point(aes(size = factor(nannot)))+
    geom_smooth(se = FALSE,method = "rlm")+
    scale_color_brewer(palette = "Set1")+theme_bw()+
    theme(legend.position = "top")+coord_cartesian(xlim = c(0,3e2))+
    xlab("Average distance between annotated events")+ylab("Sensitivity")
}

ext <- 30
minreads <- 6e3
sens <- mclapply(predictions,sensitivity,anchors,annots,ext,mc.cores = 8)
pdf(file = "figs/profiles/resolution_sensitivity_experiment.pdf")
ggplot(resoDT,aes(seq , reso,colour = seq))+geom_boxplot()+coord_cartesian(ylim = c(0,1e2))+
  facet_grid(rif ~ rep)+theme_bw()+scale_color_brewer(palette = "Set1")+theme(legend.position = "none")
sens_plot("edsn1320","edsn1402","20min","rep-2",minreads,peaks,sens)
dev.off()

sensDT <- list()
sensDT[[1]] <- sens_plot("edsn1311","edsn1396","0min","rep-1",200,peaks,sens)

sensDT[[2]] <- sens_plot("edsn1314","edsn1398","20min","rep-1",500,peaks,sens)
sensDT[[3]] <- sens_plot("edsn1317","edsn1400","0min","rep-2",500,peaks,sens)
sensDT[[4]] <- sens_plot("edsn1320","edsn1402","20min","rep-2",700,peaks,sens)
sensDT <- do.call(rbind,sensDT)




## plot chunk

## plot_regions <- function(reads,predictions,name,anchors,annots)
## {
##    depth <- length(reads)
##    reads <- split(reads,strand(reads))

##    covF <- coverage(ranges(reads[["+"]]))
##    covR <- coverage(ranges(reads[["-"]]))
##    anchors <- split(anchors,1:length(anchors))

##    covplot <- function(anchor,covF,covR,annots,pred,name){
##      covF <- covF[anchor]
##      covR <- covR[anchor]
##      if(length(covF) > 0){
##        dt1 <- data.table(x = start(anchor):end(anchor),y = 1e9 * as.vector(covF) / depth,strand = "F")
##      }else{
##        dt1 <- data.table(x = start(anchor):end(anchor),y = 0,strand = "F")
##      }
##      if(length(covR) > 0){
##        dt2 <- data.table(x = start(anchor):end(anchor),y = -1e9 * as.vector(covR)/depth,strand = "R")
##      }else{
##        dt2 <- data.table(x = start(anchor):end(anchor),y = 0 ,strand = "R")
##      }
##      dt <- rbind(dt1,dt2)
##      dt[,name := name]
##      ss <- mid(subsetByOverlaps(annots,anchor))
##      pp <- mid(subsetByOverlaps(pred,anchor))

##      p <- ggplot(dt,aes(x,y,colour = strand))+geom_step()+scale_color_brewer(palette = "Set1")+
##        theme_bw()+theme(legend.position = "none")+
##        xlab("genomic coordinates")+ylab("normalized counts")

##      if(length(ss) > 0){
##        p <- p + geom_vline(xintercept = ss,linetype = 2,colour = "black")
##      }

##      if(length(pp) > 0){
##        p <- p + geom_vline(xintercept =  pp,linetype = 2, colour = "red")
##      }

##      p <- p +facet_grid(name ~. )
     
##      return(p)

##    }

##    plots <- mclapply(anchors,covplot,covF,covR,annots,predictions,name,mc.cores = detectCores())

##    return(plots)

## }

## plots <- mapply(plot_regions,reads,predictions,names(reads),MoreArgs = list(anchors,annots),SIMPLIFY = FALSE)


## pdf(file = "figs/profiles/EColi_methods_comp.pdf",width = 9,height = 12)
## for(k in 1:length(plots[[1]])){
##   grid.arrange(plots[[1]][[k]],
##                plots[[2]][[k]],
##                plots[[3]][[k]],nrow = 3)
## }
## dev.off()

## pdf(file = "figs/profiles/EColi_ChIPexo_dpeak_experiment.pdf",width = 9,height = 18)
## for(k in 1:length(plots[[1]])){
##   grid.arrange(plots[[1]][[k]],
##                plots[[2]][[k]],
##                plots[[3]][[k]],
##                plots[[4]][[k]],nrow = 4)
## }
## dev.off()



## pdf(file = "figs/profiles/EColi_ChIPseq_SET_dpeak_experiment.pdf",width = 9,height = 18)
## for(k in 1:length(plots[[1]])){
##   grid.arrange(plots[[5]][[k]],
##                plots[[6]][[k]],
##                plots[[7]][[k]],
##                plots[[8]][[k]],nrow = 4)
## }
## dev.off()







