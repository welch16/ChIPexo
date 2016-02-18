
rm(list = ls())
library(viridis)
library(data.table)
library(GenomicAlignments)
library(viridis)
library(devtools)
library(parallel)
load_all("~/Desktop/Docs/Code/ChIPexoQual")

indir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(indir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert= TRUE)]
mc <- detectCores()

expts <- lapply(file.path(indir,files),create_exo_experiment,calc_summary = TRUE,height = 1,parallel = TRUE , mc = detectCores())
names(expts) <- plyr::mapvalues(gsub(".sort.bam","",files),from = c("ERR336935","ERR336942","ERR336956"),
                                to = c("Rep-3","Rep-1","Rep-2"))

stats <- lapply(expts,summary_stats)
stats <- mapply(function(x,y)x[,repl := y],stats,names(expts),SIMPLIFY = FALSE)
stats <- do.call(rbind,stats)


library(ggplot2)
library(scales)
library(data.table)
library(RColorBrewer)
library(hexbin)
library(viridis)


r <- viridis::viridis(1e3, option = "D")

figs_dir <- "figs/ChIPexoQuals"

pdf(file = file.path(figs_dir,"Carroll_FoxA1_mouse_enrichment.pdf"),width =9 ,height = 5)
ggplot(stats,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+facet_grid( ~ repl)+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,.85))+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = "log10")+
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")
dev.off()


