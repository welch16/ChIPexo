
rm(list = ls())
library(ChIPUtils)
library(data.table)
library(GenomicAlignments)

indir <- "/p/keles/ChIPexo/volume4/pugh_data"
files <- "CTCF.bam"

reads <- lapply(file.path(indir,files),create_reads)
names(reads) <- gsub(".bam","",files)
nreads <- sapply(reads,nreads)
pbc <- sapply(reads,PBC)

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

pugh_ctcf <- data.table(files = names(reads),nreads,pbc,nsc)

save(pugh_ctcf,file = "data/paper_tables/pugh_ctcf.RData")
