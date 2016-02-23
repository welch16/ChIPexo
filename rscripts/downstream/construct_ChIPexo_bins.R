
rm(list = ls())

library(mosaics)
library(GenomicAlignments)
library(parallel)

frag_len <- 150
bin_size <- 150

out_dir <- "/p/keles/ChIPexo/volume6/K12/downstream"
in_dir1 <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment"
in_dir2 <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic"

out_dir <- file.path(out_dir,"ChIPexo")
out_dir <- file.path(out_dir,"bins")

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

files1 <- c("edsn1311_Sig70.sort.bam","edsn1314_Sig70.sort.bam",
           "edsn1317_Sig70.sort.bam","edsn1320_Sig70.sort.bam")
files2 <- c("edsn931_Sig70.sort.bam","edsn933_Sig70.sort.bam")


A <- mclapply(c(file.path(in_dir1,files1),file.path(in_dir2,files2)),
              constructBins,
              fileFormat = "bam",
              outfileLoc = out_dir,PET = FALSE,
              fragLen = frag_len,binSize = bin_size,mc.cores = 6 )

