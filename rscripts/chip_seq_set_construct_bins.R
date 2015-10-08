
rm(list = ls())

library(mosaics)

bin_dir <- "/p/keles/ChIPexo/volume6/results/mosaics_peaks/Landick/ChIPseq_SET/bins"
bam_dir <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET"

files <- list.files(bam_dir)

files <- files[grep("edsn",files)]
files <- files[grep("filter",files)]
files <- files[grep("bai",files,invert = TRUE)]

binSize <- 150
fragLen <- 150

out <- lapply(file.path(bam_dir,files),
       constructBins,
       fileFormat = "bam",
       outfileLoc = bin_dir,
       PET = FALSE,
       fragLen = 150,
       binSize = 150)
