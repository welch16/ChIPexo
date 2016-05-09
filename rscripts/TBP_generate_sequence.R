
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)

load_all("~/Desktop/Docs/Code/ChIPexoQual")

data_dir1 <- "/p/keles/ChIPexo/volume4/venters_data/sortbam"
files <- list.files(data_dir1,
       pattern = "sort.bam",full.names = TRUE,include.dirs = TRUE)

files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("TBP",files)]

exo <- lapply(files,create_exo_experiment,parallel = TRUE)

exo_reads <- lapply(exo,reads)
rl_summary <- lapply(exo_reads,function(x)table(width(x)))
rl <- 45


stats <- lapply(exo,summary_stats)
### filter regions

## 1 - wider regions
stats <- lapply(stats,function(x)x[between(width,3 * rl,1e3)])

## most of the regions 
## > stats <- lapply(stats,function(x)x[between(width,2 * rl,1e3)])
## > sapply(stats,nrow)
## [1] 32326  7447  5500
## > stats <- lapply(stats,function(x)x[between(width,3 * rl,1e3)])
## > sapply(stats,nrow)
## [1] 2146 2578 1832

regions <- lapply(stats,function(x)
                  ChIPUtils::dt2gr(x[,2:4,with = FALSE]))
regions <- lapply(regions,sort)

library(BSgenome.Hsapiens.UCSC.hg19)

sequences <- mclapply(regions,function(x)
   getSeq(Hsapiens,x),mc.cores = 3)

sequences <- mclapply(sequences,function(x)as.character(x),mc.cores =3)

nms <- lapply(regions,function(x){
  paste0(as.character(seqnames(x)),":",start(x),"-",end(x))})

fasta_formats <- mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)


out_dir <- "/p/keles/ChIPexo/volume4/tbp_analysis/sequences/chipexo"

mapply(write.table,fasta_formats,
       file.path(out_dir,
                 gsub("sort.bam","sequences.fna",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))

#### for chip nexus

data_dir2 <- "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"
files <- list.files(data_dir2,
       pattern = "sort.bam",full.names = TRUE,include.dirs = TRUE)

files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("TBP",files)]

nexus <- lapply(files,create_exo_experiment,parallel = TRUE)

nexus_reads <- lapply(nexus,reads)

rl_summary <- lapply(exo_reads,function(x)table(width(x)))
rl <- 45

stats <- lapply(nexus,summary_stats)
### filter regions

## 1 - wider regions
stats <- lapply(stats,function(x)x[between(width,3 * rl,1e3)])

## most of the regions 
## > stats <- lapply(stats,function(x)x[between(width,2 * rl,1e3)])
## > sapply(stats,nrow)
## [1] 32326  7447  5500
## > stats <- lapply(stats,function(x)x[between(width,3 * rl,1e3)])
## > sapply(stats,nrow)
## [1] 2146 2578 1832

stats <- lapply(stats,function(x)x[depth > 100])
stats <- lapply(stats,function(x)x[npos > 20])


regions <- lapply(stats,function(x)
                  ChIPUtils::dt2gr(x[,2:4,with = FALSE]))
regions <- lapply(regions,sort)

sequences <- mclapply(regions,function(x)
   getSeq(Hsapiens,x),mc.cores = 3)

sequences <- mclapply(sequences,function(x)as.character(x),mc.cores =3)

nms <- lapply(regions,function(x){
  paste0(as.character(seqnames(x)),":",start(x),"-",end(x))})

fasta_formats <- mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)


out_dir <- "/p/keles/ChIPexo/volume4/tbp_analysis/sequences/chipnexus"

mapply(write.table,fasta_formats,
       file.path(out_dir,
                 gsub("sort.bam","sequences.fna",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))
