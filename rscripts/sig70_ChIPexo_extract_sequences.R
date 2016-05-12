
rm(list = ls())

library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(BSgenome.Ecoli.NCBI.20080805)
genome <- "NC_010473"

devtools::load_all("~/Desktop/Docs/Code/ChIPexoQual")


read_dir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo"
read_files <- list.files(read_dir,recursive = TRUE,
                         full.names = TRUE,pattern = "sort.bam")
read_files <- read_files[grep("bai",read_files,invert = TRUE)]

exo <- mclapply(read_files,create_exo_experiment,parallel = FALSE,mc.cores = 8)
stats <- lapply(exo,summary_stats)

## filter a high quality set of regions

rl <- 50
stats <- lapply(stats,function(x)x[between(width , rl,2e3)])
stats <- lapply(stats,function(x)x[f > 0  & r > 0])

minNpos <- 15
stats <- lapply(stats,function(x)x[npos > minNpos])

minDepth <- 500
stats <- lapply(stats,function(x)x[depth > minDepth])

gr <- lapply(stats,function(x)
   GRanges(seqnames = genome, ranges =              
  ranges(ChIPUtils::dt2gr(x[,2:4,with = FALSE]))))

nms <- lapply(stats,function(x)x[,(match)])

sequences <- mclapply(gr,function(x)
    getSeq(BSgenome.Ecoli.NCBI.20080805,x),mc.cores = 8 )

fasta_formats <- mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)

out_dir <- "/p/keles/ChIPexo/volume6/K12/meme"
 
mapply(write.table,fasta_formats,
       file.path(out_dir,
                 gsub(".sort.bam","_sequences.s",basename(read_files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))

