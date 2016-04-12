
rm(list = ls())


library(devtools)
library(data.table)
library(GenomicAlignments)

load_all("~/Desktop/Docs/Code/ChIPUtils")

indir <- "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"
files <- list.files(indir)
files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("K562",files,invert = FALSE)]

mc <- 2
reads <- lapply(file.path(indir,files),create_reads)#,mc.cores = mc)
names(reads) <- gsub(".sort.bam","",files)
nreads <- sapply(reads,nreads)
pbc <- do.call(c,mclapply(reads,PBC,mc.cores = mc))
 
fsr <- sapply(reads,function(x){
  tab <- summary(x)
  fwd <- sum(tab[,as.numeric(readsF)])
  bwd <- sum(tab[,as.numeric(readsR)])
  fwd / (fwd + bwd)
})


human.size <-system.file("extdata","chrom.sizes","hg19.chrom.sizes",package = "ChIPUtils")
human.size <- data.table(read.table(human.size))
scc <- lapply(reads,strand_cross_corr,shift = 1:300,chrom.size = human.size, parallel = TRUE)

nsc <- sapply(scc,function(x)x[,max(cross.corr) / min(cross.corr)])

chipnexus_K562 <- data.table(files = names(reads),nreads,pbc,nsc)

save(chipnexus_K562,file = "data/paper_tables/ChIPnexus_hg19.RData")

scc <- mapply(function(x,y)x[,expt := y],scc,names(reads),SIMPLIFY = FALSE)
scc <- do.call(rbind,scc)

library(scales)

pdf(file = "figs/ChIPnexus_scc_K562.pdf")
ggplot(scc,aes(shift,cross.corr,colour = expt))+geom_line()+
    theme_bw()+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Dark2",name = "")+
  guides(colour = guide_legend(nrow = 4))
dev.off()
