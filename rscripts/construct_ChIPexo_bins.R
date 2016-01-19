
rm(list = ls())

library(mosaics)
library(GenomicAlignments)
library(parallel)

frag_len <- 150
bin_size <- 150

out_dir <- "/p/keles/ChIPexo/volume6/resolution"
in_dir <- "/p/keles/ChIPexo/volume7/Landick/ChIPexo/rif_treatment"

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

out_dir <- file.path(out_dir,"ChIPexo")
check_create(out_dir)

out_dir <- file.path(out_dir,"bins")
check_create(out_dir)

files <- c("edsn1311_Sig70.sort.bam","edsn1314_Sig70.sort.bam",
           "edsn1317_Sig70.sort.bam","edsn1320_Sig70.sort.bam")

A <- mclapply(file.path(in_dir,files),
              constructBins,
              fileFormat = "bam",
              outfileLoc = out_dir,PET = FALSE,
              fragLen = frag_len,binSize = bin_size,mc.cores = 4 )

