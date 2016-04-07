rm(list = ls())
library(ChIPUtils)
library(data.table)
library(GenomicAlignments)

indir <- "/p/keles/ChIPexo/volume4/venters_data/sortbam"
files <- list.files(indir)
files <- files[grep(".sort.bam",files)]
files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("TBP",files)]

reads <- lapply(file.path(indir,files),create_reads)#,mc.cores = 3)
names(reads) <- gsub(".sort.bam","",files)
nreads <- sapply(reads,nreads)
pbc <- do.call(c,mclapply(reads,PBC,mc.cores = 3))

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

venters_hg19_TBP <- data.table(files = names(reads),nreads,pbc,nsc)

save(venters_hg19_TBP,file = "data/paper_tables/Venters_hg19_TBP.RData")
