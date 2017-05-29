
rm(list = ls())

library(ChIPexoQual)
library(magrittr)
library(readr)

dr = "data/figures/fig5"

files = list.files(dr,recursive = TRUE,full.names = TRUE)
files = files[grep("carroll",files)]
files = files[grep("Fox",files)]
files = files[grep("stats",files)]

library(parallel)
options("mc.cores" = 22)


stats = files %>% mclapply(read_tsv)

## files = files[grep("sort",files)]
## files = files[grep("bai",files,invert = TRUE)]
## files = files[grep("txt",files,invert = TRUE)]


## reads = files %>% mclapply(readGAlignments,param = NULL)
## reads = reads %>% mclapply(as,"GRanges")
## names(reads) = paste0("Rep",c(3,1,2))

## exo = lapply(reads,function(x)ExoData(reads =  x))

exo = stats %>% lapply(function(x)
    GRanges(seqnames = x$seqnames,
            ranges = IRanges(
                start = x$start,end = x$end)))
names(exo) = paste0("Rep",seq_len(3))


peakdr = "/p/keles/ChIPexo/volume4"
peakfiles = list.files(peakdr,full.names = TRUE,recursive = TRUE,pattern = "peaks")
peakfiles = peakfiles[grep("mouse",peakfiles)]


library(dplyr)
library(data.table)

peaks = mclapply(peakfiles,read_delim,delim = " ",col_names = FALSE)
peaks = lapply(peaks,function(x){
    x = x %>% select(X1,X2,X3)
    setnames(x,names(x),c("seqnames","start","end"))
    x = as.data.table(x)
  return(ChIPUtils::dt2gr(x))})
names(peaks) = paste0("Rep",c(3,1,2))

## peak columns:
## chrID peakStart peakStop peakSize logAveP logMinP aveLogP aveChipCount maxChipCount map GC 

## join all peak regions together
all_peaks = Reduce(c,peaks) %>% reduce

exo_peaks = mclapply(exo,subsetByOverlaps,all_peaks)

readlength = exo %>% sapply(function(x)min(width(x)))


#    reads %>% sapply(function(x)x %>% width %>% median)

## common sense filter

## remove chrM
exo_peaks = exo_peaks %>% mclapply(function(x)x[as.character(seqnames(x)) != "chrM"])

## width analysis

exo_peaks %>% lapply(function(x)width(x) %>% summary)

## $Rep3
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   34.00   36.00   53.00   80.97   98.00 3399.00 

## $Rep1
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    35.0    60.0   233.0   254.8   400.0  3390.0 

## $Rep2
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    35.0    39.0    86.0   147.8   234.0  3388.0 

rl = median(readlength)

K = 3
exo_peaks = exo_peaks %>% lapply(function(x)x[width(x) >= K * rl])

library(BSgenome.Mmusculus.UCSC.mm9)

sequences = mclapply(exo_peaks,function(x)
   getSeq(Mmusculus,x),mc.cores = 3)

sequences = mclapply(sequences,function(x)as.character(x),mc.cores =3)

nms = lapply(exo_peaks,function(x){
  paste0(as.character(seqnames(x)),":",start(x),"-",end(x))})

fasta_formats = mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)

exo_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo"

mapply(write.table,fasta_formats,
       file.path(exo_dir,"sequences",
                 gsub("_stats.tsv","_exo_peak_sequences2.fna",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))

mapply(write.table,fasta_formats,
       file.path(exo_dir,"sequences",
                 gsub(".sort.bam","_exo_peak_sequences.fna",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))
