
rm(list = ls())

library(ChIPUtils)
library(dplyr)
library(readr)
library(GenomicAlignments)
library(dplyr)
library(GenomicRanges)


chrom.sizes = "/p/keles/SOFTWARE/hg19.chrom.sizes"

dr = "/p/keles/ChIPexo/volume4/pugh_data"
files = list.files(dr,full.names = TRUE,pattern = "bam",recursive = TRUE)
files = files[grep("bai",files,invert =TRUE)]

exofiles = files[grep("mace",files)]
setfiles = files[grep("encod",files)]
setfiles = setfiles[grep("input",setfiles,invert = TRUE)]

sizes = read_delim(chrom.sizes,delim = "\t",col_names = FALSE)

exo = lapply(exofiles,create_reads)
set = lapply(setfiles,create_reads)

shift = 300

library(parallel)

sizes = as.data.table(sizes)

scc = lapply(c(exo,set), strand_cross_corr,shift = seq_len(shift),
   chrom.sizes = sizes,parallel = TRUE)

scc = mapply(function(x,y){as.tbl(x) %>% mutate(sample = y)},
             scc,c(paste0("ChIP-exo",seq_along(exofiles)),
                   paste0("ChIP-seq",seq_along(setfiles))),
             SIMPLIFY = FALSE) %>%
    bind_rows

write_tsv(scc,"data/figures/fig1/fig1F_scc.tsv")
