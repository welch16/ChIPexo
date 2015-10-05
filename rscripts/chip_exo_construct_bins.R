
rm(list = ls())

library(mosaics)

bin_dir <- "/p/keles/ChIPexo/volume6/results/mosaics_peaks/Landick/ChIPexo/bins"
bam_dir <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"

files <- list.files(bam_dir)

files <- files[grep("edsn",files)]
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
