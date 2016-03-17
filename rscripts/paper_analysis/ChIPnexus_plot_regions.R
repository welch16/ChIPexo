
rm(list = ls())
library(data.table)
library(devtools)
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(reshape2)
library(scales)
library(RColorBrewer)


bam_dir <- "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"
files <- list.files(bam_dir)
files <- files[grep("K562",files)]
files <- files[grep("bai",files,invert = TRUE)]


load_all("~/Desktop/Docs/Code/ChIPexoQual")
nexus <- lapply(file.path(bam_dir,files),create_exo_experiment,calc_summary = TRUE,parallel = TRUE,mc.cores = 20)
names(nexus) <- gsub(".sort.bam","",files)

stats <- lapply(nexus,summary_stats)
tags <- lapply(nexus,reads)

tags <- lapply(tags,function(x)resize(x,1))

sapply(stats,nrow)

stats <- lapply(stats,function(x)x[f > 0 & r > 0])
sapply(stats,nrow)

stats <- lapply(stats,function(x)x[npos > 100])
sapply(stats,nrow)

stat2gr <- function(x)x[,GRanges(seqnames = seqnames,
  ranges = IRanges(start = start, end = end))]                                 

regions <- lapply(stats,stat2gr)

common <- intersect(regions[[1]],regions[[2]])


plot_regions <- function(reads,name,anchors)
{
  
  depth <- length(reads)
  reads <- split(reads,strand(reads))

  covF <- coverage(reads[["+"]])
  covR <- coverage(reads[["-"]])
  anchors <- split(anchors,1:length(anchors))

  covplot <- function(anchor,covF,covR,name)
  {
    covF <- covF[anchor][[as.character(seqnames(anchor))]]
    covR <- covR[anchor][[as.character(seqnames(anchor))]]
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

    p <- ggplot(dt,aes(x,y,colour = strand))+geom_step()+scale_color_brewer(palette = "Set1")+
      theme_bw()+theme(legend.position = "none")+
      xlab("genomic coordinates")+ylab("normalized counts")

    p <- p +facet_grid(name ~. )+ggtitle( paste0(as.character(seqnames(anchor)),
                                                 ":",prettyNum(start(anchor),big.mark = ","),
                                                 "-",prettyNum(end(anchor),big.mark = ",")))
     
    return(p)
  }

  anchors <- as.list(anchors)
  plots <- mclapply(anchors,covplot,covF,covR,name,mc.cores = 22)

   return(plots)

}

plots <- mapply(plot_regions,tags,names(tags),MoreArgs = list(common),SIMPLIFY = FALSE)

library(gridExtra)

figsdir <- "figs/ChIPnexus"
pdf(file = file.path(figsdir,"ChIPnexus_K562_TBP_profiles.pdf"),height = 6,width = 8)
for(i in 1:length(plots[[1]]))grid.arrange(plots[[1]][[i]],plots[[2]][[i]],nrow = 2)
dev.off()
  
