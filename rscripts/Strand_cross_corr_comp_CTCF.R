
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

exo_file <- "/p/keles/ChIPexo/volume4/pugh_data/raw_data/CTCF_hg18.sort.bam"
set_file1 <- "/p/keles/ChIPexo/volume4/pugh_data/encode/encode-Uw-Helas3-Ctcf-rep1.bam"
set_file2 <- "/p/keles/ChIPexo/volume4/pugh_data/encode/encode-Uw-Helas3-Ctcf-rep2.bam"

exo <- create_reads(exo_file)
set1 <- create_reads(set_file1)
set2 <- create_reads(set_file2)


scc_exo <- strand_cross_corr(exo,chrom.sizes = sizes,parallel = TRUE,shift = 1:300)
scc_set1 <- strand_cross_corr(set1,chrom.sizes = sizes,parallel = TRUE,shift = 1:300) 
scc_set2 <- strand_cross_corr(set2,chrom.sizes = sizes,parallel = TRUE,shift = 1:300) 
 
scc_exo[,seq := "ChIP-exo"]
scc_set1[,seq := "ChIP-Seq (SE)"]
scc_set2[,seq := "ChIP-Seq (SE)"]


dt <- rbind(scc_exo[, repl := 1],
            scc_set1[,repl := 1],
            scc_set2[,repl := 2])

r1 <- brewer.pal(name = "Set1",5)
r2 <- brewer.pal(name = "Dark2",5)
r <- c(r2[c(1,3:5)],r1[2:3])

r <- brewer.pal(name = "Set1",3)

pdf(file ="figs/for_paper/scc_ctcf.pdf",width = 8,height = 5)
ggplot(dt,
  aes(shift , cross.corr,colour = as.factor(repl),linetype = seq))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Replicate",values = r[c(1,2,2)],
                     guide = guide_legend(nrow = 1))+
  scale_linetype_manual(name = "Protocol",values = c(1,2,2),
         guide = guide_legend(nrow = 1))+ylab("Strand Cross-Corr. (SCC)")+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))
dev.off()



pdf(file ="figs/for_paper/scc_ctcf2.pdf",width = 8,height = 4)
ggplot(dt,
  aes(shift , cross.corr,linetype = as.factor(repl),colour = seq))+
  geom_line(size = 1)+
  theme_bw()+
  scale_color_manual(name = "Replicate",values = r[c(1,3,3)],guide = FALSE)+
  scale_linetype_manual(name = "Replicate",values = c(1,2,2))+
  ylab("Strand Cross-Corr. (SCC)")+
  theme_bw()+scale_x_continuous(breaks = seq(10,300,by = 40))+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  facet_grid(  . ~ seq)+
  geom_vline(xintercept = 36,colour = "darkgrey",linetype = 3,size = 1)+
  geom_vline(xintercept = 90,colour = "black",linetype = 3,size = 1)
dev.off()
