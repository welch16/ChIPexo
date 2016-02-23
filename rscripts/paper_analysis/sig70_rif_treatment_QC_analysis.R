
rm(list = ls())

library(GenomicAlignments)
library(ChIPUtils)
library(data.table)

indir <- "/p/keles/ChIPexo/volume7/Landick/K12"
figsdir <- "figs/ChIPseq_QC"
outdir <- "data/ChIPseq_QC"

files <- list.files(indir,recursive = TRUE)

files <- files[grep("seq",files)]
files <- files[grep("rif",files)]
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("Input",files,invert = TRUE)]

pet <- files[grep("PET",files)]
set <- files[-grep("PET",files)]

set <- mclapply(file.path(indir,set),create_reads,mc.cores = 4)
pet <- mclapply(file.path(indir,pet),create_reads,is_PET = TRUE,mc.cores = 4)

pet_nreads <- sapply(pet,nreads)*2
set_nreads <- sapply(set,nreads)

## fix nreads for pet

pet_pbc <- sapply(pet,PBC)
set_pbc <- sapply(set,PBC)

pet_fsr <- sapply(pet,function(x){
  tab <- summary(x)
  fwd <- tab[,as.numeric(readsF)]
  bwd <- tab[,as.numeric(readsR)]
  fwd / (fwd + bwd)
})

set_fsr <- sapply(set,function(x){
  tab <- summary(x)
  fwd <- tab[,as.numeric(readsF)]
  bwd <- tab[,as.numeric(readsR)]
  fwd / (fwd + bwd)
})


ecoli.size <- data.table(read.table("/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/peak_inputs/ecoli.size"))
setnames(ecoli.size,names(ecoli.size),c("chr","size"))
pet_scc <- mclapply(pet,strand_cross_corr,shift = 1:300,chrom.size = ecoli.size, parallel = FALSE,mc.cores = 4)
set_scc <- mclapply(set,strand_cross_corr,shift = 1:300,chrom.size = ecoli.size, parallel = FALSE,mc.cores = 4)

pet_nsc <- lapply(pet_scc,function(x)x[,max(cross.corr) / min(cross.corr)])
set_nsc <- lapply(set_scc,function(x)x[,max(cross.corr) / min(cross.corr)])


## pet_scc <- mapply(function(x,y)x[,file := y],pet_nsc,files[grep("PET",files)],SIMPLIFY = FALSE)
## pet_scc <- do.call(rbind,pet_scc)
## set_scc <- mapply(function(x,y)x[,file := y],set_nsc,files[-grep("PET",files)],SIMPLIFY = FALSE)
## set_scc <- do.call(rbind,set_scc)

## scc <- rbind(pet_scc,set_scc)


pet <- data.table(files = basename(files[grep("PET",files)]),nreads = pet_nreads,pbc = pet_pbc,nsc = pet_nsc)
set <- data.table(files = basename(files[grep("SET",files)]),nreads = set_nreads,pbc = set_pbc,nsc = set_nsc)

chipseq = list(pet,set)
names(chipseq) <- c("pet","set")

save(chipseq,file = "data/paper_tables/Landick_rif_chipseq_summary.RData")
