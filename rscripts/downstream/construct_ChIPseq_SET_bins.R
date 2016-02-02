

rm(list = ls())

library(mosaics)
library(GenomicAlignments)
library(parallel)

frag_len <- 150
bin_size <- 150

out_dir <- "/p/keles/ChIPexo/volume6/K12/downstream"
in_dir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_SET/rif_treatment"


out_dir <- file.path(out_dir,"ChIPseq_SET")
out_dir <- file.path(out_dir,"bins")

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)


files <- c("edsn1396_Sig70.sort.bam","edsn1398_Sig70.sort.bam",
           "edsn1400_Sig70.sort.bam","edsn1402_Sig70.sort.bam")

A <- mclapply(file.path(in_dir,files),
              constructBins,
              fileFormat = "bam",
              outfileLoc = out_dir,PET = FALSE,
              fragLen = frag_len,binSize = bin_size,mc.cores = 4 )

