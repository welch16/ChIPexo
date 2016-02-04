
rm(list = ls())

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ChIPUtils)
library(GenomicAlignments)

prom_dir <- "/p/keles/ChIPexo/volume6/gold_standard/Landick"
pfiles <- list.files(prom_dir)
pfiles <- pfiles[grep(".bed",pfiles)]
promoters <- lapply(file.path(prom_dir,pfiles),read.table)
promoters <- lapply(promoters,data.table)

promoters[[1]][,str := "R"]
promoters[[2]][,str := "F"]

sites <- do.call(rbind,promoters)
sites[, id := 1:nrow(sites)]
setnames(sites,names(sites),c("seqnames","start","end","strand","id"))

sites[,seqnames := "U00096"]
sites[,strand := ifelse(strand == "F","+","-")]

####################################################################

## load predictions

dpeak <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_R1/ChIP-exo_sigma70_exp_phase_dpeak_fit.bed"
dpeak <- data.table(read.table(dpeak))


apex <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_Apex/931_ecoli_Apex.pp.bed"
apex <- data.table(read.table(apex))


gem <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_GEM_one_sample_R1/GEM_sigma70_exp_phase_R1_GEM_events.txt"
gem <- read.delim(gem,sep = "\t")

####################################################################

## convert to GRanges all
dpeak[, V1 := "U00096"]
setnames(dpeak,names(dpeak),c("seqnames","start","end"))
dpeak <- dpeak[!is.na(start)]
dpeak <- dt2gr(dpeak)

apex[ , V1 := "U00096"]
apex_group <- apex[,(V4)]
apex[,V4 := NULL]
setnames(apex,names(apex),c("seqnames","start","end"))
apex <- dt2gr(apex)

gem_all <- gem
gem <- rownames(gem)

gem <- strsplit(gem,":",fixed = TRUE)

gem <- GRanges(seqnames = sapply(gem,function(x)x[1]),
               ranges = IRanges(
                 start = as.numeric(sapply(gem,function(x)x[2])),
                 width = 1))

resolution <- function(site,positions)
{
  out <- min(abs(site - positions))
  return(out)
}

dpeak_res <- sapply(sites[,(start)],resolution,start(dpeak))
apex_res <- sapply(sites[,(start)],resolution,start(apex))
gem_res <- sapply(sites[,(start)],resolution , start(gem))

res1 <- data.table(reso = dpeak_res,method = "dPeak")
res2 <- data.table(reso = apex_res, method = "Mace")
res3 <- data.table(reso = gem_res,method = "Gem")

reso <- rbind(res1,res2,res3)
reso[ , method := factor(method, levels = c("dPeak","Mace","Gem"))]

pdf(file = "figs/for_paper/algorithm_resolution.pdf",width = 5 , height = 5)
ggplot(reso[reso < 200], aes(method , reso,colour = method))+geom_boxplot()+
  xlab("Algorithm")+ylab("Resolution")+coord_cartesian(ylim = c(0,150))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0),legend.position = "none")+
  scale_color_brewer(palette = "Set2")
dev.off()


