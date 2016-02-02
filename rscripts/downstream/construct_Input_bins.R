
rm(list = ls())

library(mosaics)
library(GenomicAlignments)

frag_len <- 150
bin_size <- 150
in_dir <- "/p/keles/ChIPexo/volume7/Landick/K12"
out_dir <- "/p/keles/ChIPexo/volume6/K12/downstream"
file <- "edsn1369_Input.sort.bam"

constructBins(infile = file.path(in_dir,"ChIPseq_PET","rif_treatment",file),
              fileFormat = "bam",
              outfileLoc = file.path(out_dir,"ChIPseq_PET"),PET = TRUE,
              fragLen = frag_len,binSize = bin_size)

constructBins(infile = file.path(in_dir,"ChIPseq_SET","rif_treatment",file),
              fileFormat = "bam",
              outfileLoc = file.path(out_dir,"ChIPseq_SET"),PET = FALSE,
              fragLen = frag_len,binSize = bin_size)

