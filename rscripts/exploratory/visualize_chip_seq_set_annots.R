rm(list = ls())

library(ggplot2)
library(GenomicAlignments)
library(data.table)
library(grid)
library(gridExtra)

dr <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_SET"
files <- list.files(dr,recursive = TRUE)
files <- files[grep("rif_treatment",files)]
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)][-1]

reads <- mclapply(file.path(dr,files),readGAlignments,param = NULL,mc.cores = 4)
reads <- lapply(reads,as,"GRanges")
names(reads) <- sapply(strsplit(basename(files),"_"),function(x)x[1])


reads <- lapply(reads,function(x){
  x <- resize(x,1)
  return(x)})

annot_dir <- "/p/keles/ChIPexo/volume6/K12/annotations"
files <- list.files(annot_dir)
files <- files[grep("bed",files)]

annots <- lapply(file.path(annot_dir,files),read.table)
annots[[1]] <- IRanges(start = annots[[1]][,2],width = 1)
annots[[2]] <- IRanges(end = annots[[2]][,3],width = 1)
annots <- reduce(c(annots[[1]],annots[[2]]))
width(annots) <- 1
annots <- sort(annots)

size <- 300
anchors <- annots[which(diff(mid(annots)) > size)]
center <- mid(anchors)
start(anchors) <- center - size
end(anchors) <- center + size
start(anchors[1]) <- 1


plot_regions <- function(reads,anchors,annots)
{
   depth <- length(reads)
   reads <- split(reads,strand(reads))

   covF <- coverage(ranges(reads[["+"]]))
   covR <- coverage(ranges(reads[["-"]]))
   anchors <- split(anchors,1:length(anchors))

   covplot <- function(anchor,covF,covR,annots){
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
     yy <- mid(subsetByOverlaps(annots,anchor))

     p <- ggplot(dt,aes(x,y,colour = strand))+geom_step()+scale_color_brewer(palette = "Set1")+
       theme_bw()+theme(legend.position = "none")+geom_vline(xintercept = yy ,linetype = 2)+
       xlab("genomic coordinates")+ylab("normalized counts")

     return(p)

   }

   plots <- mclapply(anchors,covplot,covF,covR,annots,mc.cores = detectCores())

   return(plots)

}

plots <- lapply(reads,plot_regions,anchors,annots)

pdf("figs/profiles/EColi_annots_rif_SE_ChIPseq.pdf",width = 9,height = 20)
for(k in 1:length(plots[[1]])){
  grid.arrange(plots[[1]][[k]],
               plots[[2]][[k]],
               plots[[3]][[k]],
               plots[[4]][[k]],nrow = 4)
}
dev.off()


pdf("figs/profiles/EColi_annots_rif_more1BS_SE_ChIPseq.pdf",width = 9,height = 20)
for(k in which(countOverlaps(anchors,annots) > 1)){
  grid.arrange(plots[[1]][[k]],
               plots[[2]][[k]],
               plots[[3]][[k]],
               plots[[4]][[k]],nrow = 4)
}
dev.off()





