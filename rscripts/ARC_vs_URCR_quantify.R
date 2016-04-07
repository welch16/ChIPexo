
rm(list = ls())

library(devtools)
library(viridis)
library(hexbin)
library(data.table)
library(GenomicAlignments)
library(parallel)
library(scales)

data_dir <- "/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles"
data_dir <- "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"

files <- list.files(data_dir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]
 
files <- files[grep("embryo_dors",files)]
load_all("~/Desktop/Docs/Code/ChIPexoQual")

exo <- lapply(file.path(data_dir,files),
  create_exo_experiment,parallel = TRUE,height = 1)

stats <- lapply(exo,summary_stats)

enrichment_plots <- function(stat)
{
  r <- viridis(100,option = "D")

  plots <- list()
  plots[[1]] <- ggplot(stat,aes(npos,depth))+stat_binhex(bins = 75)+
    scale_fill_gradientn(colours = r,trans = "log10",
      labels = trans_format("log10",math_format(10^.x)))+
    theme_bw()+theme(legend.position = "top")+
      xlim(0,500)+ylim(0,3000)
  plots[[2]] <- ggplot(stat,aes(ave_reads,cover_rate))+
    stat_binhex(bins = 75)+
    scale_fill_gradientn(colours = r,trans = "log10",
      labels = trans_format("log10",math_format(10^.x)))+
    theme_bw()+theme(legend.position = "top")+
      xlim(0,4)+ylim(0,1)
  plots[[3]] <- ggplot(stat,aes(npos,width))+stat_binhex(bins = 75)+
    scale_fill_gradientn(colours = r,trans = "log10",
      labels = trans_format("log10",math_format(10^.x)))+
    theme_bw()+theme(legend.position = "top")+
      xlim(0,500)+ylim(0,1200)
  plots[[4]] <- ggplot(stat,aes(depth,width))+stat_binhex(bins = 75)+
    scale_fill_gradientn(colours = r,trans = "log10",
      labels = trans_format("log10",math_format(10^.x)))+
    theme_bw()+theme(legend.position = "top")+
      xlim(0,3000)+ylim(0,1200)
  return(plots)
}

plots <- lapply(stats,enrichment_plots)

library(gridExtra)

figs_dir <- "figs/test_quantify"
labels <- c("Rep-3","Rep-1","Rep-2")
labels <- as.character(1:3)
labels <- c("1311","1314","1317","1320")
labels <- as.character(seq(931,937,by = 2))
labels <- c("Rep1","Rep2")

pdf(file.path(figs_dir,"ChIPnexus_embryo_dorsal.pdf"))
p <- plots[[1]]
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],nrow = 2,top = labels[1])
p <- plots[[2]]
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],nrow = 2,top = labels[2])
## p <- plots[[3]]
## grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],nrow = 2,top = labels[3] )
## p <- plots[[4]]
## grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],nrow = 2,top = labels[4] )
dev.off()

density_by_npos <- function(stat,what,npos_val)
{
  out <- lapply(npos_val,function(x)
    ggplot(stat[npos > x],aes_string(x = what))+
      stat_density(geom = "line",kernel = "rectangular")+xlim(0,10)
                )
  return(out)
}

pdf("aa.pdf")
density_by_npos(stats[[1]],"ave_reads",c(1,5,10,20,50,100,200))
dev.off()




