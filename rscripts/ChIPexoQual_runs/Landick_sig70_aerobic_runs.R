
rm(list = ls())
library(devtools)
##library(ChIPUtils)
load_all("~/Desktop/Docs/Code/ChIPexoQual")

indir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic"
files <- list.files(indir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]

expts <- mclapply(file.path(indir,files),create_exo_experiment,calc_summary = TRUE,height = 1,mc.cores =4)
lapply(expts,function(A)A@summary_stats[,sum(depth)]/depth(A))

stats <- lapply(expts,summary_stats)
stats <- mapply(function(x,y)x[,condition := y],stats,1:4,SIMPLIFY = FALSE)
stats <- do.call(rbind,stats)

stats[,repl := ifelse(condition %in% 1:2,"Rep-1","Rep-2")]
stats[,condition := ifelse(condition %in% c(1,3),"Aerobic","Anaerobic")]


library(ggplot2)
library(scales)
library(data.table)
library(RColorBrewer)
library(hexbin)
library(viridis)


r <- viridis::viridis(1e3, option = "D")

figs_dir <- "figs/ChIPexoQuals"

pdf(file = file.path(figs_dir,"Sig70_aero_enrichment.pdf"))
ggplot(stats,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+facet_grid(repl ~ condition)+
  scale_x_continuous(limits = c(0,10))+
  scale_y_continuous(limits = c(0,.5))+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = "log10")+
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")
dev.off()


