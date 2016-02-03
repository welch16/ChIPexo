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
files <- files[grep("PET",files,invert = TRUE)] ## only keep SET
files <- files[grep("Input",files,invert = TRUE)] ## only chip

annot_dir <- "/p/keles/ChIPexo/volume6/K12/annotations"
afiles <- list.files(annot_dir)
afiles <- afiles[grep("bed",afiles)]

annots <- lapply(file.path(annot_dir,afiles),read.table)
annots[[1]] <- IRanges(start = annots[[1]][,2],width = 1)
annots[[2]] <- IRanges(end = annots[[2]][,3],width = 1)
annots <- reduce(c(annots[[1]],annots[[2]]))
width(annots) <- 1
annots <- sort(annots)

size <- 300
anchors <- annots[which(diff(mid(annots)) > size)]
center <- mid(anchors)
start(anchors) <- center - size * 2
end(anchors) <- center + size * 2
start(anchors[1]) <- 1
anchors <- anchors[countOverlaps(anchors,annots) > 1]

dt <- ChIPUtils::gr2dt(GRanges(seqnames = "U00096",anchors))

temp_peakfile <- tempfile(pattern = "peak",fileext = ".txt")

write.table(dt,file = temp_peakfile,row.names = FALSE,col.names = FALSE,quote = FALSE)

## reads <- lapply(reads,function(x){
##   x <- resize(x,1)
##   return(x)})

library(dpeak)

get_sites <- function(file,peakfile){
  dpeak <- dpeakRead(peakfile = peakfile, readfile = file,fileFormat = "bam",parallel = TRUE,
                     nCore = 24)
  fit <- dpeakFit(dpeak,maxComp = 5,nCore = 24)
  temp_bed <- tempfile(pattern = "binding",fileext = ".bed")
  export(fit,type = "bed",filename = temp_bed)
  out <- data.table(read.table(temp_bed,skip = 1))
  return(out)
}

sites <- lapply(file.path(dr,files),get_sites,temp_peakfile)

predictions <- lapply(sites,function(x)
  IRanges(start = mid(x[,IRanges(start = V2,end = V3)]),width = 1))


reads <- mclapply(file.path(dr,files),readGAlignments,param = NULL,mc.cores = 4)
reads <- lapply(reads,as,"GRanges")
names(reads) <- sapply(strsplit(basename(files),"_"),function(x)x[1])


reads <- lapply(reads,function(x){
  x <- resize(x,1)
  return(x)})


plot_regions <- function(reads,predictions,name,anchors,annots)
{
   depth <- length(reads)
   reads <- split(reads,strand(reads))

   covF <- coverage(ranges(reads[["+"]]))
   covR <- coverage(ranges(reads[["-"]]))
   anchors <- split(anchors,1:length(anchors))

   covplot <- function(anchor,covF,covR,annots,pred,name){
     covF <- covF[anchor]
     covR <- covR[anchor]
     if(length(covF) > 0){
       dt1 <- data.table(x = start(anchor):end(anchor),y = 1e9 * as.vector(covF) / depth,strand = "F")
     }else{
       dt1 <- data.table(x = start(anchor):end(anchor),y = 0,strand = "F")
     }
     if(length(covR) > 0){
       dt2 <- data.table(x = start(anchor):end(anchor),y = -1e9 * as.vector(covR)/depth,strand = "R")
     }else{
       dt2 <- data.table(x = start(anchor):end(anchor),y = 0 ,strand = "R")
     }
     dt <- rbind(dt1,dt2)
     dt[,name := name]
     ss <- mid(subsetByOverlaps(annots,anchor))
     pp <- mid(subsetByOverlaps(pred,anchor))

     p <- ggplot(dt,aes(x,y,colour = strand))+geom_step()+scale_color_brewer(palette = "Set1")+
       theme_bw()+theme(legend.position = "none")+
       xlab("genomic coordinates")+ylab("normalized counts")

     if(length(ss) > 0){
       p <- p + geom_vline(xintercept = ss,linetype = 2,colour = "black")
     }

     if(length(pp) > 0){
       p <- p + geom_vline(xintercept =  pp,linetype = 2, colour = "red")
     }

     p <- p +facet_grid(name ~. )
     
     return(p)

   }

   plots <- mclapply(anchors,covplot,covF,covR,annots,predictions,name,mc.cores = detectCores())

   return(plots)

}

plots <- mapply(plot_regions,reads,predictions,names(reads),MoreArgs = list(anchors,annots),SIMPLIFY = FALSE)


pdf(file = "figs/profiles/EColi_ChIPexo_dpeak_experiment.pdf",width = 9,height = 18)
for(k in 1:length(plots[[1]])){
  grid.arrange(plots[[1]][[k]],
               plots[[2]][[k]],
               plots[[3]][[k]],
               plots[[4]][[k]],nrow = 4)
}
dev.off()



pdf(file = "figs/profiles/EColi_ChIPseq_SET_dpeak_experiment.pdf",width = 9,height = 18)
for(k in 1:length(plots[[1]])){
  grid.arrange(plots[[5]][[k]],
               plots[[6]][[k]],
               plots[[7]][[k]],
               plots[[8]][[k]],nrow = 4)
}
dev.off()

## pdf("figs/profiles/EColi_annots_rif_SE_ChIPseq.pdf",width = 9,height = 20)
## for(k in 1:length(plots[[1]])){
##   grid.arrange(plots[[1]][[k]],
##                plots[[2]][[k]],
##                plots[[3]][[k]],
##                plots[[4]][[k]],nrow = 4)
## }
## dev.off()


## pdf("figs/profiles/EColi_annots_rif_more1BS_SE_ChIPseq.pdf",width = 9,height = 20)
## for(k in which(countOverlaps(anchors,annots) > 1)){
##   grid.arrange(plots[[1]][[k]],
##                plots[[2]][[k]],
##                plots[[3]][[k]],
##                plots[[4]][[k]],nrow = 4)
## }
## dev.off()






