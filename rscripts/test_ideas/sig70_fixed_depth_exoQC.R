
rm(list = ls())

library(devtools)
library(parallel)
library(GenomicAlignments)
library(data.table)
library(viridis)
library(hexbin)

load_all("~/Desktop/Docs/Code/ChIPexoQual")

indir <- "/p/keles/ChIPexo/volume6/K12/saturation/ChIPexo/"

files <- list.files(indir,recursive = TRUE)
files <- files[grep("bam",files)]
files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("samples",files)]
files <- files[grep("12345",files)]

exo <- mclapply(file.path(indir,files),create_exo_experiment,calc_summary = TRUE,parallel = FALSE,mc.cores = 24)
names(exo) <- basename(files)

stats <- lapply(exo,summary_stats)

stats <- mapply(function(x,y)x[,file := y],stats,names(stats),SIMPLIFY = FALSE)

STAT <- do.call(rbind,stats)
STAT[,file := gsub(".bam","",file)]
STAT[,edsn := sapply(strsplit(file,"_",fixed = TRUE),function(x)x[1])]
setkey(STAT,edsn)

library(scales)
r <- viridis::viridis(100,option = "D")

pdf(file = "figs/sig70_rif_enrichment.pdf",height = 12,width = 12)
p <- ggplot(STAT["edsn1311"],aes(ave_reads,cover_rate))+stat_binhex(bins = 50)+
  facet_wrap( ~ file,ncol = 3)+scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+xlim(0,3)+ylim(0,1)+
  theme_bw()+theme(legend.position = "top")+
  xlab("Average read coverage")+
  ylab("Unique read coverage rate")
p
p %+% STAT["edsn1314"]
p %+% STAT["edsn1317"]
p %+% STAT["edsn1320"]
dev.off()
