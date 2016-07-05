
rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(GenomicRanges)
library(parallel)
library(BSgenome.Ecoli.NCBI.20080805)

peakdir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPseq_PET/peaks"
fdr <- "FDR5"
genome <- "NC_010473"


peakdir <- file.path(peakdir,fdr)
files <- list.files(peakdir,include.dirs = TRUE,full.names = TRUE)
peaks <- lapply(files,fread)

peaks <- lapply(peaks,function(x)x[ V8 > 500])

peaks <- lapply(peaks,function(x)x[,GRanges(seqnames = genome,
     ranges = IRanges(start = V2,end = V3))])

peaks <- reduce(do.call(c,peaks))
peaks <- sort(shift(peaks,1))

save(peaks, file = "data/sig70_ChIPSeq_PE_unified_peaks.RData")

sequences <- getSeq(BSgenome.Ecoli.NCBI.20080805,peaks)

nms <- paste0(as.character(seqnames(peaks)),":",
              start(peaks),"-",end(peaks))

fasta_formats <- paste0(">",nms,"\n",sequences)

out_dir <- "/p/keles/ChIPexo/volume6/K12/meme_sig70"

write.table(fasta_formats,
  file.path(out_dir,"ChIPseq_PE_sequences.s"),
  quote = FALSE,row.names = FALSE,col.names = FALSE)



