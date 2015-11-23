
rm(list = ls())

library(ChIPUtils)
library(reshape2)
library(ggplot2)
library(scales)
library(RColorBrewer)

hg19 <- system.file("extdata","chrom.sizes","hg19.chrom.sizes",
  package = "ChIPUtils")
sizes <- data.table(read.table(hg19))
shift <- 1:200

exo_file <- "/p/keles/ChIPexo/volume3/PughData/CTCF.bam"
set_file <- "/p/keles/ChIPexo/volume3/PughData/encode-Uw-Helas3-Ctcf-rep1.bam"

exo <- create_reads(exo_file)
set <- create_reads(set_file)


scc_exo <- strand_cross_corr(exo,chrom.sizes = sizes,parallel = TRUE,shift = 1:300)
scc_set <- strand_cross_corr(set,chrom.sizes = sizes,parallel = TRUE,shift = 1:300) 

dt <- rbind(scc_exo[,seq := "exo"],scc_set[,seq := "set"])

pdf(file ="figs/for_paper/scc_ctcf.pdf",width = 8,height = 4)
ggplot(dt,aes(shift , cross.corr,colour = seq))+geom_line(size = 1)+
  theme_bw()+scale_color_manual(name = "",values = brewer.pal(3,"Set1")[c(1,3)])+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))
dev.off()
