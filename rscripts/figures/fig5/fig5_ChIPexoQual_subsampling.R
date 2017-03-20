#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--readsfile",action = "store_true",type = "character",
                help = "Bam file with the aligned reads from the ChIP-exo/nexus experiment"),
    make_option("--depths",action = "store_true",type = "character",
                help = "Number of reads (in millions) used to subsample, separated by comma"),
    make_option("--nregions",action = "store_true",type = "numeric",default = 1e3,
                help = "Number of regions used to calculate QC scores"),
    make_option("--ntimes",action = "store_true",type = "numeric",default = 1e3,
                help = "Number of times the bootstrap procedure is repeared
                        to calculate QC scores"),
    make_option("--statsfile",action = "store_true",type = "character",
                default = tempfile(pattern = "stats",fileext = ".tsv"),
                help = "Name of the file where the summary stats are stored."),   
    make_option("--scoresfile",action = "store_true",type = "character",
                default = tempfile(pattern = "scores",fileext =".tsv"),
               help = "Name of the file where the qc scores are saved"),
    make_option("--cores",action = "store_true",type = "numeric",default = 22,
                help = "Number of cores used for parallel")
)


opt = parse_args(OptionParser(option_list = optList))

library(base,quietly = TRUE)
library(magrittr,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(GenomicAlignments,quietly = TRUE)
library(GenomicRanges,quietly = TRUE)
library(parallel,quietly = TRUE)
library(tidyverse,quietly = TRUE)
library(ChIPexoQual,quietly = TRUE)
library(readr,quietly = TRUE)


options(mc.cores = opt$cores)

reads = readGAlignments(opt$readsfile,param = NULL)
reads = reads %>% as("GRanges")

opt$depths = opt$depths %>% strsplit(",") %>% unlist %>% as.numeric

exo = ExoDataSubsampling(reads = reads,sample.depth = opt$depths * 1e6)

depths = opt$depths %>% as.character %>% paste0("M")

##

stats = exo %>% lapply(function(x)as.data.frame(x) %>% as.tbl)
qcscores = exo %>% lapply(function(x)paramDist(x) %>% as.data.frame %>% as.tbl)

message("Writing files")

statfiles = depths %>% lapply(function(x)gsub(".tsv",paste0("_",x,".tsv"),opt$statsfile))

scoresfiles = depths %>% lapply(function(x)gsub(".tsv",paste0("_",x,".tsv"),opt$scoresfile))

aa = mapply(write_tsv,stats,statfiles)
bb = mapply(write_tsv,qcscores,scoresfiles)

