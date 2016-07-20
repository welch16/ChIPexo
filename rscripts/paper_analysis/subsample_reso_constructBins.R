
rm(list = ls())

library(mosaics)
library(parallel)

seed <- "12345"

frag_len <- 150
bin_size <- 150

indir <- "/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity"

files <- list.files(indir,full.names = TRUE,recursive = TRUE)
files <- files[grep("subsample_1",basename(files))]

files <- files[grep(seed,files)]
files <- files[grep("fragL",files,invert = TRUE)]

is_pet <- vapply(files,function(x)grepl("PET",x),TRUE)
names(is_pet) <- NULL

file_format <- vapply(files,function(x)grepl(".eland",x),TRUE)
file_format <- ifelse(file_format,"eland_result","sam")
names(file_format) <- NULL

out <- mcmapply(constructBins,
                infile = files,
                fileFormat = file_format,
                outfileLoc = dirname(files),
                PET = is_pet,
                fragLen = frag_len,
                binSize = bin_size,SIMPLIFY = FALSE,mc.cores = 4)
                  
