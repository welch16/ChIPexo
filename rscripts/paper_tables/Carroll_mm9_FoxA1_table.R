
rm(list = ls())
library(ChIPUtils)
library(data.table)
library(GenomicAlignments)

indir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(indir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]

reads <- mclapply(file.path(indir,files),create_reads,mc.cores = 3)
names(reads) <- gsub(".sort.bam","",files)
nreads <- sapply(reads,nreads)
pbc <- do.call(c,mclapply(reads,PBC,mc.cores = 3))

fsr <- sapply(reads,function(x){
  tab <- summary(x)
  fwd <- sum(tab[,as.numeric(readsF)])
  bwd <- sum(tab[,as.numeric(readsR)])
  fwd / (fwd + bwd)
})

mouse.size <-system.file("extdata","chrom.sizes","mm9.chrom.sizes",package = "ChIPUtils")
mouse.size <- data.table(read.table(mouse.size))
scc <- lapply(reads,strand_cross_corr,shift = 1:300,chrom.size = mouse.size, parallel = TRUE)

nsc <- sapply(scc,function(x)x[,max(cross.corr) / min(cross.corr)])

carroll_mm9_FoxA1 <- data.table(files = names(reads),nreads,pbc,nsc)

save(carroll_mm9_FoxA1,file = "data/paper_tables/Carroll_mm9_FoxA1.RData")
