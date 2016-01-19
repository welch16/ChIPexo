
rm(list = ls())

library(mosaics)
library(GenomicAlignments)

frag_len <- 150
bin_size <- 150

out_dir <- "/p/keles/ChIPexo/volume6/resolution/inputs"
in_dir <- "/p/keles/ChIPexo/volume7/Landick"

file <- "edsn1369_Input.sort.bam"

constructBins(infile = file.path(in_dir,"ChIPseq_PET","rif_treatment",file),
              fileFormat = "bam",
              outfileLoc = file.path(out_dir,"ChIPseq_PET"),PET = TRUE,
              fragLen = frag_len,binSize = bin_size)

constructBins(infile = file.path(in_dir,"ChIPseq_SET","rif_treatment",file),
              fileFormat = "bam",
              outfileLoc = file.path(out_dir,"ChIPseq_SET"),PET = FALSE,
              fragLen = frag_len,binSize = bin_size)

