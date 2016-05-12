
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
peaks <- lapply(files,read.table)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[ V8 > 1])

peaks <- lapply(peaks,function(x)x[,GRanges(seqnames = genome,
    ranges = IRanges(start = V2,end = V3))])

sequences <- mclapply(peaks,function(x)
    getSeq(BSgenome.Ecoli.NCBI.20080805,x),mc.cores =4)

nms <- lapply(peaks,function(x){
  paste0(as.character(seqnames(x)),":",start(x),"-",end(x))})

fasta_formats <- mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)

out_dir <- "/p/keles/ChIPexo/volume6/K12/meme"
 
mapply(write.table,fasta_formats,
       file.path(out_dir,
                 gsub("peaks.txt","sequences.s",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))



