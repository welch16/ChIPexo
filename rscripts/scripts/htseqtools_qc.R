#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

## defines option parse and help caption for script
optList = list(
  make_option("--chipfile" , action = "store_true",type = "character",
              help = "File corresponding to the ChIP-seq/exo/nexus file used"),
  make_option("--outfile", action = "store_true",type = "character",default = tempfile(),
              help = "File where the SCC is going to be stored"),
  make_option("--ncluster",action = "store_true",type = "numeric",default = 12,
              help = "Number of clusters used for parallel processing.")
 )

opt = parse_args(OptionParser(option_list = optList))

## opt$chipfile = "/p/keles/ChIPexo/volume4/carroll_data/mouse/ERR336942.sort.bam"

stopifnot(file.exists(opt$chipfile))

library(htSeqTools)
library(GenomicAlignments)
library(parallel)
library(BiocParallel)
library(magrittr)

reads = readGAlignments(opt$chipfile,param = NULL) %>% as("GRanges")

SSD = reads %>% ssdCoverage %>% round(4)
gini = reads %>% giniCoverage(mc.cores = 22) %>% round(4)

library(dplyr)

write.csv(tibble(SSD,GINI = gini[1],ADJ = gini[2]),
            file = opt$outfile,quote =FALSE,row.names =FALSE)


