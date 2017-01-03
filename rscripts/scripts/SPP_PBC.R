#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

## defines option parse and help caption for script
optList = list(
  make_option("--chipfile" , action = "store_true",type = "character",
              help = "File corresponding to the ChIP-seq/exo/nexus file used"),
  make_option("--sccfile", action = "store_true",type = "character",
              help = "File corresponding to the SCC curve"),
  make_option("--outfile", action = "store_true",type = "character",default = tempfile(),
              help = "File where the SCC is going to be stored"),
  make_option("--ncluster",action = "store_true",type = "numeric",default = 12,
              help = "Number of clusters used for parallel processing.")

 )


opt = parse_args(OptionParser(option_list = optList))

stopifnot(file.exists(opt$chipfile))

library(spp)
library(snow)
library(parallel)
library(readr)

tags = read.bam.tags(opt$chipfile)

table.chip.data = mclapply(tags$tags,table,mc.cores = opt$ncluster)

nUniq = sum(sapply(table.chip.data,function(x) sum(x==1)))
nTotal = sum(sapply(table.chip.data,length))

pbc = nUniq / nTotal

library(dplyr)
library(GenomicAlignments)

ga = readGAlignments(opt$chipfile,param = NULL)
w = width(ga)
rl = floor(mean(w))

scc = read_delim(opt$sccfile,delim = "\t")
scc = scc %>% mutate(cross.corr = ifelse(cross.corr < 0,0,cross.corr))

scc = scc %>% filter(cross.corr > 0)

nsc = scc %>% summarize(max(cross.corr ) /min(cross.corr))
nsc = nsc[[1]]

c1 = scc %>% summarise(max(cross.corr))
c2 = scc %>% summarise(min(cross.corr))
c3 = scc %>% filter(shift == rl) %>% select(cross.corr)

rsc = (c1[[1]] - c2[[1]]) / (c3[[1]] - c2[[1]])


write_delim(tibble(pbc,read_length = rl , nsc , rsc),opt$outfile , delim = "\t")
