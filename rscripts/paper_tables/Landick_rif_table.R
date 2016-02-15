
rm(list = ls())
library(ChIPUtils)
library(data.table)
library(GenomicAlignments)

indir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment"
files <- list.files(indir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]

reads <- mclapply(file.path(indir,files),create_reads,mc.cores = 4)
names(reads) <- gsub(".sort.bam","",files)

nreads <- sapply(reads,nreads)
pbc <- sapply(reads,PBC)

fsr <- sapply(reads,function(x){
  tab <- summary(x)
  fwd <- tab[,as.numeric(readsF)]
  bwd <- tab[,as.numeric(readsR)]
  fwd / (fwd + bwd)
})

ecoli.size <- data.table(read.table("/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/peak_inputs/ecoli.size"))
setnames(ecoli.size,names(ecoli.size),c("chr","size"))
scc <- mclapply(reads,strand_cross_corr,shift = 1:300,chrom.size = ecoli.size, parallel = FALSE,mc.cores = 4)

nsc <- lapply(scc,function(x)x[,max(cross.corr) / min(cross.corr)])

landick_rif <- data.table(files = names(reads),nreads,pbc,nsc)

save(landick_rif,file = "data/paper_tables/Landick_rif_summary.RData")
