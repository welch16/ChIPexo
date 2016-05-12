
rm(list = ls())

library(parallel)
library(data.table)
library(GenomicAlignments)
library(devtools)

load_all("~/Desktop/Docs/Code/ChIPUtils")

mc <- detectCores()

read_dir <- "/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles"

files <- list.files(read_dir,pattern = "chipseq.sort",full.names = TRUE)
files <- files[grep("bai",files,invert = TRUE)]

reads <- mclapply(files,create_reads,FALSE,mc.cores = 3)

reads <- create_reads(bamfile,isPET)
sizefile <- "hg19"

if(tolower(sizefile) %in% c("hg19","mm9","mm10","dm3")){
  sizedir <- system.file("extdata","chrom.sizes", package = "ChIPUtils")
  sizefiles <- list.files(sizedir)
  sizefile <- sizefiles[grep(sizefile,sizefiles)]
  sizefile <- file.path(sizedir,sizefile)
  rm(sizedir,sizefiles)
}

sizes <- data.table(read.table(sizefile,header = FALSE))

scc <- lapply(reads,strand_cross_corr,shift = 1:300,
   chrom.sizes = sizes,parallel = TRUE)

mapply(write.table,scc,file.path("data/SCC_curves",
  gsub(".sort.bam","_SCC.txt",basename(files))),
       MoreArgs = list(quote = FALSE,
   sep = "\t",row.names = FALSE,col.names = TRUE)) 

