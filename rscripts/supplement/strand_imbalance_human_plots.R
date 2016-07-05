
rm(list = ls())

library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(scales)
library(GenomicRanges)

load("data/CTCF_strand_imbalance.RData")
load("data/TBP_strand_imbalance.RData")
load("data/ER_strand_imbalance.RData")


density_comp <- function(samp1,samp2,lab1,lab2,oneSite = FALSE)
{ 
  ov <- findOverlaps(samp1,samp2)
  samp1 <- samp1[queryHits(ov)]
  samp2 <- samp2[subjectHits(ov)]
  dt1 <- data.table(as.data.frame(mcols(samp1)))
  dt2 <- data.table(as.data.frame(mcols(samp2)))
  dt <- rbind(dt1[,samp := lab1],dt2[,samp := lab2])
  r <- brewer.pal(3,name = "Set1")[c(1,3)]
  dt[,fsr := fwd / (fwd + bwd)]
  if(oneSite){
    dt <- dt[nsites == 1]
  }else{
    dt <- dt[nsites >= 1]
  }
  out <- ggplot(dt,aes(fsr,colour = samp))+stat_density(geom = "line")+
    theme_bw()+theme(legend.position = "top")+
    scale_color_manual(name = "",values = r)+
    geom_vline(xintercept = .5,linetype = 2,colour = "black")+
    xlim(0,1)+ylab("Density")+xlab("Ratio of forward strand reads")
  return(out)
}


### CTCF commons peaks ChIP-exo against ChIP-seq - rep1
### and ChIP-seq - rep2 by separate

ctcf_dens <- list()
ctcf_dens[[1]] <- density_comp(ctcf_peaks[[1]],ctcf_peaks[[2]],
                           "ChIP-exo","ChIP-seq (SE) Rep-1",TRUE)
ctcf_dens[[2]] <- density_comp(ctcf_peaks[[1]],ctcf_peaks[[3]],
                           "ChIP-exo","ChIP-seq (SE) Rep-2",TRUE)

### TBP commons peaks

tbp_dens <- list()
tbp_dens[[1]] <- density_comp(tbp_peaks[[1]],tbp_peaks[[4]],
                           "ChIP-exo Rep-1","ChIP-seq (SE) Rep-1",TRUE)
tbp_dens[[2]] <- density_comp(tbp_peaks[[2]],tbp_peaks[[4]],
                           "ChIP-exo Rep-2","ChIP-seq (SE) Rep-1",TRUE)
tbp_dens[[3]] <- density_comp(tbp_peaks[[3]],tbp_peaks[[4]],
                           "ChIP-exo Rep-3","ChIP-seq (SE) Rep-1",TRUE)
tbp_dens[[4]] <- density_comp(tbp_peaks[[1]],tbp_peaks[[5]],
                           "ChIP-exo Rep-1","ChIP-seq (SE) Rep-2",TRUE)
tbp_dens[[5]] <- density_comp(tbp_peaks[[2]],tbp_peaks[[5]],
                           "ChIP-exo Rep-2","ChIP-seq (SE) Rep-2",TRUE)
tbp_dens[[6]] <- density_comp(tbp_peaks[[3]],tbp_peaks[[5]],
                           "ChIP-exo Rep-3","ChIP-seq (SE) Rep-2",TRUE)

### ER common peaks

er_dens <- list()
er_dens[[1]] <- density_comp(er_peaks[[1]],er_peaks[[4]],
                           "ChIP-exo Rep-1","ChIP-seq (SE) Rep-1",TRUE)
er_dens[[2]] <- density_comp(er_peaks[[2]],er_peaks[[5]],
                           "ChIP-exo Rep-2","ChIP-seq (SE) Rep-2",TRUE)
er_dens[[3]] <- density_comp(er_peaks[[3]],er_peaks[[6]],
                           "ChIP-exo Rep-3","ChIP-seq (SE) Rep-3",TRUE)

figs_dir <- "figs/supplement"
pdf(file.path(figs_dir,"CTCF_imbalance_1Site.pdf"))
u <- lapply(ctcf_dens,print)
dev.off()

pdf(file.path(figs_dir,"TBP_imbalance_1Site.pdf"))
u <- lapply(tbp_dens,print)
dev.off()

pdf(file.path(figs_dir,"ER_imbalance_1Site.pdf"))
u <- lapply(er_dens,print)
dev.off()


### CTCF commons peaks ChIP-exo against ChIP-seq - rep1
### and ChIP-seq - rep2 by separate

ctcf_dens <- list()
ctcf_dens[[1]] <- density_comp(ctcf_peaks[[1]],ctcf_peaks[[2]],
                           "ChIP-exo","ChIP-seq (SE) Rep-1",FALSE)
ctcf_dens[[2]] <- density_comp(ctcf_peaks[[1]],ctcf_peaks[[3]],
                           "ChIP-exo","ChIP-seq (SE) Rep-2",FALSE)

### TBP commons peaks

tbp_dens <- list()
tbp_dens[[1]] <- density_comp(tbp_peaks[[1]],tbp_peaks[[4]],
                           "ChIP-exo Rep-1","ChIP-seq (SE) Rep-1",FALSE)
tbp_dens[[2]] <- density_comp(tbp_peaks[[2]],tbp_peaks[[4]],
                           "ChIP-exo Rep-2","ChIP-seq (SE) Rep-1",FALSE)
tbp_dens[[3]] <- density_comp(tbp_peaks[[3]],tbp_peaks[[4]],
                           "ChIP-exo Rep-3","ChIP-seq (SE) Rep-1",FALSE)
tbp_dens[[4]] <- density_comp(tbp_peaks[[1]],tbp_peaks[[5]],
                           "ChIP-exo Rep-1","ChIP-seq (SE) Rep-2",FALSE)
tbp_dens[[5]] <- density_comp(tbp_peaks[[2]],tbp_peaks[[5]],
                           "ChIP-exo Rep-2","ChIP-seq (SE) Rep-2",FALSE)
tbp_dens[[6]] <- density_comp(tbp_peaks[[3]],tbp_peaks[[5]],
                           "ChIP-exo Rep-3","ChIP-seq (SE) Rep-2",FALSE)

### ER common peaks

er_dens <- list()
er_dens[[1]] <- density_comp(er_peaks[[1]],er_peaks[[4]],
                           "ChIP-exo Rep-1","ChIP-seq (SE) Rep-1",FALSE)
er_dens[[2]] <- density_comp(er_peaks[[2]],er_peaks[[5]],
                           "ChIP-exo Rep-2","ChIP-seq (SE) Rep-2",FALSE)
er_dens[[3]] <- density_comp(er_peaks[[3]],er_peaks[[6]],
                           "ChIP-exo Rep-3","ChIP-seq (SE) Rep-3",FALSE)

figs_dir <- "figs/supplement"
pdf(file.path(figs_dir,"CTCF_imbalance.pdf"))
u <- lapply(ctcf_dens,print)
dev.off()

pdf(file.path(figs_dir,"TBP_imbalance.pdf"))
u <- lapply(tbp_dens,print)
dev.off()

pdf(file.path(figs_dir,"ER_imbalance.pdf"))
u <- lapply(er_dens,print)
dev.off()
