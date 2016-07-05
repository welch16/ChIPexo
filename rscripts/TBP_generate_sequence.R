
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)

load_all("~/Desktop/Docs/Code/ChIPexoQual")

data_dir <- "data/ChIPexo_QC_runs"
files <- list.files(data_dir,full.names = TRUE)

files <- files[grep("K562",files)]
files <- files[grep("TBP",files)]
files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("TBP",files)]
files <- files[-7]

load_file <- function(x){
  load(x)
  return(ext_stats[["stats"]])
}

files <- files[grep("samp",files,invert = TRUE)]

stats <- mclapply(files,load_file,mc.cores = 22)
names(stats) <- gsub(".RData","",basename(files))

rl <- 50

## 1 - wider regions
stats <- lapply(stats,function(x)x[between(width,rl,2e3)])

## 2 - removed single - stranded regions
stats <- lapply(stats,function(x)x[f > 0 & r > 0])

## 3 - deeper regions
stats <- lapply(stats,function(x){
  mm <- x[,median(depth)]
  return(x[depth > mm])})


regions <- lapply(stats,function(x)
                  ChIPUtils::dt2gr(x[,2:4,with = FALSE]))
regions <- lapply(regions,sort)

library(BSgenome.Hsapiens.UCSC.hg19)

sequences <- mclapply(regions,function(x)
   getSeq(Hsapiens,x),mc.cores = 5)

sequences <- mclapply(sequences,function(x)as.character(x),mc.cores = 5)

nms <- lapply(regions,function(x){
  paste0(as.character(seqnames(x)),":",start(x),"-",end(x))})

fasta_formats <- mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)

out_dir <- "/p/keles/ChIPexo/volume4/tbp_analysis/sequences"

mapply(write.table,fasta_formats,
       file.path(out_dir,
                 gsub(".RData","sequences.s",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))

