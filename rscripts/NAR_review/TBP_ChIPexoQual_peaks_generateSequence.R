
rm(list = ls())

library(ChIPexoQual)
library(magrittr)

dr = "/p/keles/ChIPexo/volume4"

files = list.files(dr,recursive = TRUE,full.names = TRUE)
files = files[grep("TBP",files)]
files = files[grep("bam",files)]
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("txt",files,invert = TRUE)]
files = files[grep("chipseq",files,invert = TRUE)]

library(parallel)

options("mc.cores" = 22)

reads = files %>% lapply(readGAlignments,param = NULL)
reads = reads %>% mclapply(as,"GRanges")

names(reads) = c(paste0("Exo",seq_len(3)),paste0("Nexus",seq_len(2)))

exo = lapply(reads,function(x)ExoData(reads =  x))

peakfiles = list.files(dr,full.names = TRUE,recursive = TRUE,pattern = "peaks")
peakfiles = peakfiles[grep("venters",peakfiles)]
peakfiles = peakfiles[grep("encode",peakfiles,invert = TRUE)]

library(readr)
library(dplyr)
library(data.table)

peaks = mclapply(peakfiles,read_delim,delim = " ",col_names = FALSE)
peaks = lapply(peaks,function(x){
    x = x %>% select(X1,X2,X3)
    setnames(x,names(x),c("seqnames","start","end"))
    x = as.data.table(x)
  return(ChIPUtils::dt2gr(x))})

## peak columns:
## chrID peakStart peakStop peakSize logAveP logMinP aveLogP aveChipCount maxChipCount map GC 

## join all peak regions together
all_peaks = Reduce(c,peaks) %>% reduce

exo_peaks = mclapply(exo,subsetByOverlaps,all_peaks)
readlength = reads %>% sapply(function(x)x %>% width %>% median)

## common sense filter

## remove chrM
exo_peaks = exo_peaks %>% mclapply(function(x)x[as.character(seqnames(x)) != "chrM"])

## width analysis

exo_peaks %>% lapply(function(x)width(x) %>% summary)

rl = median(readlength)

K = 3
exo_peaks = exo_peaks %>% lapply(function(x)x[width(x) >= K * rl])

library(BSgenome.Hsapiens.UCSC.hg19)

sequences = mclapply(exo_peaks,function(x)
   getSeq(Hsapiens,x),mc.cores = 3)

sequences = mclapply(sequences,function(x)as.character(x),mc.cores =3)

nms = lapply(exo_peaks,function(x){
  paste0(as.character(seqnames(x)),":",start(x),"-",end(x))})

fasta_formats = mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)

exo_dir = "/p/keles/ChIPexo/volume4/tbp_analysis/sequences"
mapply(write.table,fasta_formats,
       file.path(exo_dir,
                 gsub(".sort.bam","_exo_peak_sequences.fna",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))
