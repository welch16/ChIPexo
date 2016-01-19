

rm(list = ls())

library(mosaics)
library(GenomicAlignments)

frag_len <- 150
bin_size <- 150

out_dir <- "/p/keles/ChIPexo/volume6/resolution"
in_dir <- "/p/keles/ChIPexo/volume7/Landick/ChIPseq_SET/rif_treatment"

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)

out_dir <- file.path(out_dir,"ChIPseq_SET")
check_create(out_dir)

out_dir <- file.path(out_dir,"bins")
check_create(out_dir)

files <- c("edsn1396_Sig70.sort.bam","edsn1398_Sig70.sort.bam",
           "edsn1400_Sig70.sort.bam","edsn1402_Sig70.sort.bam")

A <- mclapply(file.path(in_dir,files),
              constructBins,
              fileFormat = "bam",
              outfileLoc = out_dir,PET = FALSE,
              fragLen = frag_len,binSize = bin_size,mc.cores = 4 )

